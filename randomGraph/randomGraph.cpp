#include "Snap.h"
#include <iostream>

//sing namespace TSnap;

PUNGraph Graph = TUNGraph::New();


void generateRandomGraph(int n, double p)
{


  for(int i = 0; i < n; i++)
  {
    Graph->AddNode(i);
    //printf("Node number is: %d\n", i);
  }

  printf("Simulation not done\n");
  printf("Number of edges is: %d\n", Graph->GetEdges());

  for(int i = 0; i < n; i++)
  {
    for(int j = i; j < n - i; j++)
    {
      //srand(time(NULL));
      int rndNum = rand() % 100 + 1;
      printf("this is the random number: %d", rndNum);
      if(p * 100 > rndNum)
      {

        printf("rand was: %d and so Edge added\n", rndNum);
        Graph->AddEdge(i , j);
      }
    }
  }

  printf("Simulation done\n");
}

int main()
{
  generateRandomGraph(10, 0.8);
  //Graph = TSnap::GenRndGnm<PNGraph>(1000000, 2000000, false);
  //Graph->DelEdge(253,254);
  printf("Nodes: %d\n", Graph->GetNodes());
  printf("Edges: %d\n", Graph->GetEdges());
  //printf("\nThe length is: %d\n", TSnap::GetShortPath(Graph, 253, 127));
  return 0;
}
