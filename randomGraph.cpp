#include "Snap.h"
#include <iostream>
#include <random>

/*Generating a random undirected graph called Graph  */
PUNGraph Graph = TUNGraph::New();


/* Function to generate a random graph with n number of
nodes. Edges are randomly created between all nodes with
a probability of p.*/
void generateRandomGraph(int n, double p)
{
  for(int i = 0; i < n; i++)
  {
    Graph->AddNode(i);
    //printf("Node number is: %d\n", i);
  }

  //printf("Simulation not done\n");
  //printf("Number of edges is: %d\n", Graph->GetEdges());

  for(int i = 0; i < n; i++)
  {
    for(int j = i +1 ; j <= n -1 ; j++)
    {
      //srand(time(NULL));
      //int rndNum = rand() % 100 + 1;
      //printf("this is the random number: %d", rndNum);
      std::random_device rd;
      std::mt19937 gen(rd());
      std::bernoulli_distribution rv_b(p);

      if( rv_b(gen) )
      {
        printf("and so Edge added\n");
        Graph->AddEdge(i , j);
      }
    }
    //printf("%d runs done",i);
  }
  printf("Simulation done\n");

}

int main()
{
  generateRandomGraph(10, 0.5);   //Generating a random undirected graph. (number of nodes, probability of edges)
  //Graph = TSnap::GenRndGnm<PNGraph>(1000000, 2000000, false);
  //Graph->DelEdge(253,254);
  printf("Nodes: %d\n", Graph->GetNodes());
  printf("Edges: %d\n", Graph->GetEdges());
  //printf("\nThe length is: %d\n", TSnap::GetShortPath(Graph, 253, 127));
  return 0;
}


/*void floydWarshall(int[][] pairs, int n)
{
  int numberOfPairs = sizeof(pairs);
  for(int a = 0; a < numberOfPairs; a++)
  {
    for(int k = 0; k < n; k++)
    {
      for(int i =0; i < n; i++)
      {
        for(int j = 0; j < n; j++)
        {

        }
      }
    }
  }
}

void dijkstra (int a, int b, int n)
{
  int[n] NodeValues;
  for(int i = a; i < n; i++)
  {
    for(int j = i; j < n; j++)
    {
      if(IsEdge(i , j))
      {
        if(NodeValues[j] > NodeValues[i] + 1)
          {
            NodeValues[j] = NodeValues[i] + 1;

          }
      }
    }
  }
} */

//Finds the shortest path between the source
//node and all the possibles nodes whose
//number is greater than the source node
//saves all the paths in a set/vector and
// returns number of edges / number of edges
int shortestPath(int source)
{
  for(int i = source; i < Graph->GetNodes() -1; i++)
  {
    //finding shortest path b/w source and source + i
    for(int j = i; j < source + i; j++)
    {
      
    }
  }
}
