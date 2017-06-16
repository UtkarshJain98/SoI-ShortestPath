#include "Snap.h"
#include <iostream>

//sing namespace TSnap;

PNGraph Graph = TNGraph::New();


/*void generateRandomGraph(int n, int p, int m)
{
  printf("%d", n);
  for(int i = 0; i < n; i++)
  {
    //Graph->AddNode(i);
    printf("%d", i);
  }
  for(int i = 0; i < m; i++)
  {
    //Graph->AddEdge(rand() % 1000 + 1);
  }
} */

int main()
{
  //generateRandomGraph(10, 5, 10);
  Graph = TSnap::GenRndGnm<PNGraph>(1000000, 2000000, false);
  //Graph->DelEdge(253,254);
  printf( "%d", Graph->GetEdges() ); //<< endl;
  printf("\nThe length is: %d\n", TSnap::GetShortPath(Graph, 253, 127));
  return 0;
}
