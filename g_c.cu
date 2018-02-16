/*
  Felix A. Portillo
  EEC289Q HW3
*/

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <set>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>

//Can i just do <thrust> ????
#include <thrust/count.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

using namespace std;

void GraphColoringGPU(const char filename[], int** color);
// Read MatrixMarket graphs
// Assumes input nodes are numbered starting from 1
void ReadMMFile(const char filename[], bool** graph, int* V) 
{
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


// Read DIMACS graphs
// Assumes input nodes are numbered starting from 1
void ReadColFile(const char filename[], bool** graph, int* V) 
{
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
   //Print input matrix
   for (int i=0; i<num_rows; i++){
      for(int j=0; j<num_rows; j++){
        cout << (*graph)[i+num_rows+j]<<"   ";
      }
      cout << "\n\n";
   }
   infile.close();
}  

/*
  API deliverable. Put code here instead of in main!!
*/
void GraphColoringGPU(const char filename[], int** color){
  bool* graph;
  int V; //num rows
  int* color;
  
  //Experiment with thrust vectors
  thrust::device_vector<bool> dev_graph(graph,graph+V*V);
  int deg = thrust::count(thrust::device, dev_graph.begin(),dev_graph.end(),1);
  
   
   if (string(filename[1]).find(".col") != string::npos)
      ReadColFile(filename[1], &graph, &V);
   else if (string(filename[1]).find(".mm") != string::npos) 
      ReadMMFile(filename[1], &graph, &V);
   else
      

}

/*
  Generate random list based on probabilities
*/
__global__
void Graph_Coloring_GPU(bool* graph,int **color, int V){
  
}

int main(void){
    int **color;
   GraphColoringGPU(argv[1],color);
   /*if (string(argv[1]).find(".col") != string::npos)
      ReadColFile(argv[1], &graph, &V);
   else if (string(argv[1]).find(".mm") != string::npos) 
      ReadMMFile(argv[1], &graph, &V);
   else
      return -1;*/

  return 0;
}