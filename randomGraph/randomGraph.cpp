#include "Snap.h"
#include <iostream>
#include <random>
#include <vector>

using namespace std;

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
      random_device rd;
      mt19937 gen(rd());
      bernoulli_distribution rv_b(p);

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
// returns number of edges divide by
// number of edges
int shortestPath(int source)
{
  vector<int> MyVector; //should be <int, int> ??
  double dist[Graph->GetEdges()];
  double prev[Graph->GetEdges()];
  int n = Graph->GetEdges();
  int v = source;
  for(int v = source; v < n; v++)
  {
    dist[v] = 0;
    prev[v] = 0;
    if(find(MyVector.begin(), MyVector.end(), v) == MyVector.end())
    {
      MyVector.push_back(v);
    }
  }

  while(MyVector.size() != 0)
  {
    //select a u to remove
    int i = source;
    for(i = source; i < n; i++)
    {
      if(IsEdge(source, i))
        break;
    }
    MyVector.erase(remove(MyVector.begin()), MyVector.end(), i), MyVector.end());

    int noOfEdges = 0;
    for(int a = 0; a < n; a++)
    {
      if(IsEdge(source, a))
      {
        noOfEdges ++;
        int alt = dist[i] + 1;
        if(alt < dist[v])
        {
          dist[v] = alt;
          prev[v] = alt;
        }
      }
    }

  }


}
