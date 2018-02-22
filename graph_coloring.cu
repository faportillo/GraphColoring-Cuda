#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include <thrust/count.h>
#include <thrust/logical.h>
#include <thrust/functional.h>

#define MAXBLOCKS 1<<30

using namespace std;

//find min color
__device__
int min_color(int v,int n,int* colorMask){
	int i = 1;
	while((colorMask[i] == v) && (i < n)){
		i++;
	}

	return i;
}

__global__
void colorTopoKernel(int n, int* NNZ, int* preSum, int* colIndex, int* colors, bool* changed, bool* colored){

    int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x*gridDim.x;
    int c;
    int colorMask[256];
    for (int i = 0; i < 256; i++){
    	colorMask[i] = 512;
    }

	for(int i = index; i < n; i+= stride){

        if (!colored[i]){
        	for (int k = preSum[i]; k < preSum[i + 1]; k++){
			    int j = colIndex[k];
			    colorMask[colors[j]] = i;
			}
			c = min_color(i, n, colorMask);
            colors[i] = c;
            colored[i] = true;
            changed[i] = true;
        }
	}
}

//check for collisions due to coloring adjacent nodes same color
__global__
void checkCollisions(int n, int* NNZ, int* preSum, int* colIndex, int* colors, bool* colored){
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x*gridDim.x;

	for(int i = index; i < n; i+=stride){
		for (int k = preSum[i]; k < preSum[i + 1]; k++){
		    int j = colIndex[k];
		    if ((colors[i] == colors[j]) && ( i < j)){
		    	colored[i] = false;
		    }
		}
	}
}

//topology driven parallel graph coloring algorithm
void colorTopo(int n, int* NNZ, int* preSum, int* colIndex, int* colors){
	bool* changed;
    bool* colored;

    cudaMallocManaged(&colored, sizeof(bool)*n);
    cudaMallocManaged(&changed, sizeof(bool)*n);

    thrust::fill(colored, colored+n, false);
    thrust::fill(colors, colors+n,0);

    do{
    	thrust::fill(changed, changed+n, false);
    	int nt = 256;
    	int nb = min((n + nt -1)/nt, MAXBLOCKS);
    	colorTopoKernel<<<nb,nt>>>(n, NNZ, preSum, colIndex, colors, changed, colored);
    	cudaDeviceSynchronize();
    	checkCollisions<<<nb,nt>>>(n, NNZ,preSum, colIndex,  colors,colored);
    	cudaDeviceSynchronize();
    }while(thrust::any_of(changed,changed+n,thrust::identity<bool>()));
    
    cudaFree(colored);
    cudaFree(changed);
}

//jones plassmann luby algorithm for graphcoloring
__global__ 
void colorJPLKernel(int n, int c, int* NNZ, int* preSum, int* colIndex,int* randoms, int* colors){

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x*gridDim.x;
    
    if (index < n){
    	for (int i = index; i < n; i += stride){
		bool f = true;

		if ((colors[i] != -1 )){
			continue;
			
		}

		int ir = randoms[i];
        //instead of looping through 
		for (int k = preSum[i]; k < preSum[i + 1]; k++){
			int j = colIndex[k];
			int jc = colors[j];
			if (((jc != -1) && (jc != c)) || (i == j)){
				continue;
			}
			int jr = randoms[j];
			if (ir <= jr){
				f = false;
			}
		}
		if (f){
			colors[i] = c;
		}
	}
    }
	
}

//jones plassmann luby algorithm for graphcoloring
void colorJPL(int n, int* NNZ, int* preSum, int* colIndex, int* colors){
    int* randoms;
    cudaMallocManaged(&randoms, sizeof(int)* n);
    for (int i = 0; i  < n; i++){
    	randoms[i] = rand(); 
    }

    thrust::fill(colors, colors+n, -1);

    for(int c = 1; c< n+1; c++){
    	int nt = 256;
    	int nb = min((n + nt -1)/nt, MAXBLOCKS);

    	colorJPLKernel<<<nb,nt>>>(n,c,NNZ,preSum,colIndex, randoms,colors);
        cudaDeviceSynchronize();
    	int left = (int)thrust::count(colors, colors+n, -1);
    	if (left ==0){
    		break;
    	}
    }

    cudaFree(randoms);
}

// Counts the number of unique colors in a solution
int CountColors(int V, int* color){
   int num_colors = 0;
   set<int> seen_colors;

   for (int i = 0; i < V; i++) {
      if (seen_colors.find(color[i]) == seen_colors.end()) {
         seen_colors.insert(color[i]);
         num_colors++;
      }  
   }

   return num_colors;
}

// Returns true if the color assignment is valid for the graph
bool IsValidColoring(bool* graph, int V, int* color) {
   for (int i = 0; i < V; i++) {
      for (int j = 0; j < V; j++) {
         if (graph[i * V + j]) {
            if (i != j && color[i] == color[j]) {
               printf("Vertex %d and Vertex %d are connected and have the same color %d\n", i, j, color[i]);
               return false;
            }
            if (color[i] < 1) {
               printf("Vertex %d has invalid color %d\n", i, color[i]);
               return false;
            }
         }
      }
   }

   return true;
}

// Read DIMACS graphs
// Assumes input nodes are numbered starting from 1
void ReadColFile(const char filename[], bool** graph, int* V) {
   string line;
   ifstream infile(filename);
   if (infile.fail()) {
      printf("Failed to open %s\n", filename);
      return;
   }

   int num_rows, num_edges;

   while (getline(infile, line)) {
      istringstream iss(line);
      string s;
      int node1, node2;
      iss >> s;
      if (s == "p") {
         iss >> s; // read string "edge"
         iss >> num_rows;
         iss >> num_edges;
         *V = num_rows;
         *graph = new bool[num_rows * num_rows];
         memset(*graph, 0, num_rows * num_rows * sizeof(bool));
         continue;
      } else if (s != "e")
         continue;
      
      iss >> node1 >> node2;

      // Assume node numbering starts at 1
      (*graph)[(node1 - 1) * num_rows + (node2 - 1)] = true;
      (*graph)[(node2 - 1) * num_rows + (node1 - 1)] = true;
   }
   infile.close();
}

// Read MatrixMarket graphs
// Assumes input nodes are numbered starting from 1
void ReadMMFile(const char filename[], bool** graph, int* V) {
   string line;
   ifstream infile(filename);
   if (infile.fail()) {
      printf("Failed to open %s\n", filename);
      return;
   }

   // Reading comments
   while (getline(infile, line)) {          
      istringstream iss(line);
      if (line.find("%") == string::npos)
         break;
   }

   // Reading metadata
   istringstream iss(line);
   int num_rows, num_cols, num_edges;
   iss >> num_rows >> num_cols >> num_edges;

   *graph = new bool[num_rows * num_rows];
   memset(*graph, 0, num_rows * num_rows * sizeof(bool));
   *V = num_rows;

   // Reading nodes
   while (getline(infile, line)) {          
      istringstream iss(line);
      int node1, node2, weight;
      iss >> node1 >> node2 >> weight;
      
      // Assume node numbering starts at 1
      (*graph)[(node1 - 1) * num_rows + (node2 - 1)] = true;
      (*graph)[(node2 - 1) * num_rows + (node1 - 1)] = true;
   }
   infile.close();
}

//store sparse graph in compressed sparse row format
void CSRConvert(bool** graph, int rows, int** NNZ, int* preSum, int** colIndex, int* counter){

    //assume square matrix
    int cols = rows;
    *counter = 0;
    int rowElem[rows];

    for (int i = 0; i < rows; i++){
        rowElem[i] = 0;
    }
    //initialize preSum
    preSum[0] = 0;

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            if ((*graph)[i*rows + j] == false){
                continue;
            }
            else{
                //reallocate size of NNZ and colIndex
               *NNZ = (int*)realloc(*NNZ, sizeof(int)*(*counter + 1));
                 
                (*NNZ)[*counter] = 1;

                *colIndex = (int*)realloc(*colIndex, sizeof(int) * (*counter +1));
                (*colIndex)[*counter] = j;

                //preSum[counter + 1] = preSum[counter] + prevRowCount;
                rowElem[i]++;
                *counter += 1;
            }
        }
    }

    for (int i = 0; i < rows +1; i++){
        preSum[i+1] = preSum[i] + rowElem[i];
    }
}

//Assignment API
void GraphColoringGPU(const char filename[], int** color){
    bool* graph;
    int V;

    if (string(filename).find(".col") != string::npos)
      ReadColFile(filename, &graph, &V);
    else{
      ReadMMFile(filename, &graph, &V);
    }

    //convert the sparse array into compact sparse row format
    int *NNZ = (int*)malloc(sizeof(int)); 
    int *preSum = (int*)malloc(sizeof(int) * (V + 1));
    int *colIndex= (int*)malloc(sizeof(int));
    int counter = 0;
    CSRConvert(&graph, V, &NNZ, preSum, &colIndex, &counter);
    
    //migrate values to GPU
    int* Ao;
    int* Av;
    int* Ac;
    int* colors;
    cudaMallocManaged(&Ao, sizeof(int)*(V+1));
    cudaMallocManaged(&Av, sizeof(int)*counter);
    cudaMallocManaged(&Ac, sizeof(int)*counter);
    cudaMallocManaged(&colors, sizeof(int)* V);
    
    
    for(int i = 0; i < counter; i++){
    	Av[i] = NNZ[i];
    	Ac[i] = colIndex[i];
    }
    
    
    //printf("offset values : ");
    for (int i = 0; i < V + 1; i++){
    	Ao[i] = preSum[i];

    }
    
    // colorJPL(V, Av, Ao, Ac,colors);
    // printf("JPL coloring found solution with %d colors\n", CountColors(V, colors));
    // printf("Valid coloring: %d\n", IsValidColoring(graph, V, colors));

    colorTopo(V,Av,Ao,Ac,colors);

    printf("Topo coloring found solution with %d colors\n", CountColors(V, colors));
    printf("Valid coloring: %d\n", IsValidColoring(graph, V, colors));
    
    free(NNZ);
    free(preSum);
    free(colIndex);
    cudaFree(Ao);
    cudaFree(Av);
    cudaFree(Ac);

}
