#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>

using namespace std;


// Read MatrixMarket graphs
// Assumes input nodes are numbered starting from 1
void ReadMMFile(const char filename[], bool** graph, int* V, int *N) 
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

// Counts the number of unique colors in a solution
int CountColors(int V, int* color)
{
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
bool IsValidColoring(bool* graph, int V, int* color) 
{
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

/* A utility function to print solution */
void PrintSolution(int* color, int V) 
{
   printf("Solution Exists:"
         " Following are the assigned colors \n");
   for (int i = 0; i < V; i++)
      printf(" %d ", color[i]);
   printf("\n");
}

/* A utility function to check if the current color assignment
   is safe for vertex v */
bool IsSafe (int v, bool* graph, int V, int* color, int c)
{
   for (int i = 0; i < V; i++)
      if (graph[v * V + i] && c == color[i])
         return false;

   return true;
}

/* A recursive utility function to solve m coloring problem */
bool FindColoring(bool* graph, int V, int m, int* color, int v)
{
   /* base case: If all vertices are assigned a color then
      return true */
   if (v == V)
      return true;

   /* Consider this vertex v and try different colors */
   for (int c = 1; c <= m; c++)
   {
      /* Check if assignment of color c to v is fine*/
      if (IsSafe(v, graph, V, color, c))
      {
         color[v] = c;

         /* recur to assign colors to rest of the vertices */
         if (FindColoring(graph, V, m, color, v+1) == true)
            return true;

         /* If assigning color c doesn't lead to a solution
            then remove it */
         color[v] = 0;
      }
   }

   /* If no color can be assigned to this vertex then return false */
   return false;
}

// Brute-force method
void GraphColoring(bool* graph, int V, int** color)
{
   *color = new int[V];

   // Initialize all color values as 0. This initialization is needed
   // correct functioning of IsSafe()
   for (int i = 0; i < V; i++)
      (*color)[i] = 0;

   for (int m = 1; m <= V; m++) { 
      if (FindColoring(graph, V, m, *color, 0)) {
         break;
      }
   }
}

// Assigns colors (starting from 0) to all vertices and prints
// the assignment of colors
void GreedyColoring(bool* graph, int V, int** color)
{
   *color = new int[V];
   int* solution = *color;

   // Assign the first color to first vertex
   solution[0]  = 1;

   // Initialize remaining V-1 vertices as unassigned
   for (int u = 1; u < V; u++)
      solution[u] = 0;  // no color is assigned to u

   // Assign colors to remaining V-1 vertices
   for (int u = 1; u < V; u++)
   {
      set<int> used_colors; // colors already assigned to adjacent vertices
      // Process all adjacent vertices and flag their colors
      // as unavailable
      for (int i = 0; i < V; i++)
         if (graph[u * V + i])
            if (solution[i] != 0)
               used_colors.insert(solution[i]);

      
      // Find the first available color
      int cr;
      for (cr = 1; cr < V; cr++) 
         if (used_colors.find(cr) == used_colors.end())
            break; // color not in used set

      solution[u] = cr; // Assign the found color
   }
} 

int main(int argc, char* argv[])
{
   bool* graph;
   int V; //Num rows
   int* color;
   
   if (string(argv[1]).find(".col") != string::npos)
      ReadColFile(argv[1], &graph, &V);
   else if (string(argv[1]).find(".mm") != string::npos) 
      ReadMMFile(argv[1], &graph, &V);
   else
      return -1;

 //  GraphColoring(graph, V, &color);
 //  printf("Brute-foce coloring found solution with %d colors\n", CountColors(V, color));
 //  printf("Valid coloring: %d\n", IsValidColoring(graph, V, color));

   GreedyColoring(graph, V, &color);
   printf("Greedy coloring found solution with %d colors\n", CountColors(V, color));
   printf("Valid coloring: %d\n", IsValidColoring(graph, V, color));
   
   return 0;
}
