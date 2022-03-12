#include <Rcpp.h>
#include <ctime>
#include <algorithm>
#include <queue>
#include <utility>
#include <vector>
#include <cmath>
#include <stack>
#include <map>
#include <string>




using namespace Rcpp;

typedef std::pair<double, int> iPair;




/*// [[Rcpp::export]]
NumericMatrix randomize(NumericMatrix & m)
{
  int rows = m.nrow();
  int cols = m.ncol();

  NumericMatrix rand_mat(rows,cols);

  for(int i=0; i<cols; i++)
  {
    NumericVector x = m(_,i);
    std::random_shuffle(x.begin(), x.end());
    rand_mat(_,i)=x;
  }

  return rand_mat;
} */


/*Calculates and returns pairwise distances for given a dataset */

// [[Rcpp::export]]
NumericMatrix calcPairwiseDist (NumericMatrix & centers){
  int outrows = centers.nrow();
  double d;
  NumericMatrix out(outrows,outrows);

  for (int i = 0 ; i < outrows - 1; i++)
  {

    NumericVector v1= centers.row(i);
    for (int j = i + 1  ; j < outrows ; j ++)
    {
      d = sum(pow(v1-centers.row(j),2.0));
      out(j,i)=d;
      out(i,j)=d;
    }
  }
  return (out) ;
}



/*Creates an adjacent matrix representation of a tree from a given parent vector */
std::vector<std::vector<int>> get_tree(std::vector<int>& parent)
{
  int nodes_size = parent.size();
  std::vector<std::vector<int>> tree(nodes_size);

  //NumericMatrix tree(nodes_size, nodes_size);

  for(int i=1; i<nodes_size; i++)
  {
    //tree(i,parent[i]) = 1;
    //tree(parent[i],i) = 1;

    tree[i].push_back(parent[i]);
    tree[parent[i]].push_back(i);

  }

  return tree;
}


/*Computes the Minimum Spanning Tree given the pairwise distances in dist_mat */
// [[Rcpp::export]]
std::vector<std::vector<int>> calculate_mst(NumericMatrix & dist_mat)
{
  int nodes_size = dist_mat.nrow();
  std::priority_queue< iPair, std::vector <iPair> , std::greater<iPair> > pq;


  int included = 0;
  int src = 0; // Taking vertex 0 as source

  // Create a vector for keys and initialize all
  // keys as infinite (INF)
  std::vector<double> key(nodes_size, DBL_MAX);

  // To store parent array which in turn store MST
  std::vector<int> parent(nodes_size, -1);

  // To keep track of vertices included in MST
  std::vector<bool> inMST(nodes_size, false);

  // Insert source itself in priority queue and initialize
  // its key as 0.
  pq.push(std::make_pair(0, src));
  key[src] = 0;


  /* Looping till all nodes have been inclusted in MST */
  while (included < nodes_size)//
  {

    // The first vertex in pair is the minimum key
    // vertex, extract it from priority queue.
    // vertex label is stored in second of pair (it
    // has to be done this way to keep the vertices
    // sorted key (key must be first item

    // in pair)

    int u = pq.top().second;
    pq.pop();

    if(inMST[u])
    {
      continue;
    }


    inMST[u] = true;  // Include vertex in MST
    included++;


    // 'i' is used to get all adjacent vertices of a vertex

    for (int i = 0; i<nodes_size; i++)
    {


      //  If i is not in MST and weight of (u,i) is smaller
      // than current key of i
      if (inMST[i] == false && key[i] > dist_mat(u,i))
      {

        key[i] = dist_mat(u,i);
        pq.push(std::make_pair(key[i], i));
        parent[i] = u;
      }

    }

  }

  return get_tree(parent);

}



/* Performs breadth first search on a given tree and returns a pair containing the last node
* in the search and its distance from the start node
*/
std::pair<int,int>  bfs(std::vector<std::vector<int>> & tree, int start)
{
  int nodes_size = tree.size();
  std::vector<int> distance(nodes_size, -1);
  std::queue<int> node_queue;


  // Mark the current node as visited and enqueue it
  distance[start] = 0;
  node_queue.push(start);



  while(!node_queue.empty())
  {
    // Dequeue a vertex from queue and print it
    int s = node_queue.front();
    node_queue.pop();

    // Get all adjacent vertices of the dequeued
    // vertex s. If an adjacent vertex has not been visited,
    // then mark it visited and enqueue it
    int degree = tree[s].size();
    for (int i =0; i<degree; i++)
    {
      if (distance[tree[s][i]]==-1)
      {
        distance[tree[s][i]] = distance[s] + 1;
        node_queue.push(tree[s][i]);
      }
    }
  }

  int end_node, max_dist=-1;

  /*Get the node that is furthest from start node. */
  for(int i=0; i<nodes_size;i++)
  {
    if(max_dist < distance[i])
    {
      max_dist = distance[i];
      end_node = i;
    }
  }

  return std::make_pair(end_node, max_dist);

}





/*Computes the length of longest path in a given tree */

// [[Rcpp::export]]
int get_longest_path_statistic(std::vector<std::vector<int>>& tree)
{
  std::pair<int, int> end = bfs(tree, 1);
  end = bfs(tree, end.first);
  return end.second;
}


/* Returns ceiling of num divided by 2 */
int ceiling2(int num, int div)
{
  return (num/div) + (num % div);
}


/* Returns ceiling of num divided by 2
int ceiling(int num)
{
  return (num/2) + (num % 2);
}
*/


NumericVector get_degrees(std::vector<std::vector<int>> & tree)
{
  int nodes = tree.size();
  NumericVector degrees(nodes);

  for(int i=0; i<nodes; i++)
  {
    degrees[i] = tree[i].size();
  }
  return degrees;
}


/* Computes the tree dimension of a given tree */
// [[Rcpp::export]]
List tree_dimension(std::vector<std::vector<int>> & tree)
{
  int total_nodes = tree.size();

  //NumericVector degrees = rowSums(tree);
  NumericVector degrees = get_degrees(tree);
  LogicalVector total_d1 = (degrees == 1);
  int terminals = sum(total_d1) - 2;

  int lp = get_longest_path_statistic(tree);
  double tree_dimension = 1.0 + ((total_nodes-lp-1) + ceiling2(terminals,2))/(double)(lp+1);

  List output = List::create(_["leafs"]=terminals+2, _["diameter"]=lp, _["td"]=tree_dimension);
  return output;

}

/*
//[[Rcpp::export]]
NumericMatrix DFtoNM(DataFrame x) {
  NumericMatrix y = internal::convert_using_rfunction(x, "as.matrix");
  return y;
}
*/


/*Converts edge list output from fast 'emst' to an adjacency list */
// [[Rcpp::export]]
std::vector<std::vector<int>> to_adj_mat(NumericMatrix edges)
{
  //Number of edges in MST
  int size = edges.nrow();

  //Adjacency list
  std::vector<std::vector<int>> tree(size+1);


  int from=0;  // 'from' column in edge list matrix
  int to = 1;  // 'to' column in edge list matrix

  for(int i=0; i<size; i++)
  {

    tree[edges(i,from)].push_back(edges(i,to));
    tree[edges(i,to)].push_back(edges(i,from));
  }
  return tree;

}


/* Function returns list of edges as a vector; for igraph MST plot*/
// [[Rcpp::export]]
List get_edges(NumericMatrix& mat, String MST)
{
  std::vector<std::vector<int>> tree;
  if(MST=="exact")
  {
    NumericMatrix dist_mat = calcPairwiseDist(mat);
    tree = calculate_mst(dist_mat);

  }


  else if(MST=="boruvka")
  {
    //Obtain environment containing function
    Rcpp::Environment package_env("package:mlpack");
    //Rcpp::Environment::namespace_env(pkg)

    // Make function callable from C++
    Rcpp::Function boruvka_mst = package_env["emst"];

    List emst_output = boruvka_mst(mat);
    NumericMatrix edges = emst_output["output"];

    //Transform MST format to adjacency matrix
    tree = to_adj_mat(edges);
  }

  NumericVector edges;
  int total_nodes = tree.size();

  for(int i=0; i<total_nodes; i++)
  {
    int list_size = tree[i].size();

    for(int j=0; j<list_size; j++)
    {
      if(i<tree[i][j]){
        edges.push_back(i+1);
        edges.push_back(tree[i][j]+1);
      }
    }
  }

  List output = List::create(_["tree"] = tree, _["edges"]=edges);
  return output;
}




/*Computes all tdt statistics */

// [[Rcpp::export]]
List getStatistics(NumericMatrix& mat, int sample_size, String MST, bool returnTree)
{
    // Adjacency list to hold minimum spanning tree
   std::vector<std::vector<int>> tree;

   if(MST=="exact")
   {
      NumericMatrix dist_mat = calcPairwiseDist(mat);
      tree = calculate_mst(dist_mat);
   }


   else if(MST=="boruvka")
   {
      //Obtain environment containing function
      Rcpp::Environment package_env("package:mlpack");
      //Rcpp::Environment::namespace_env(pkg)

      // Make function callable from C++
      Rcpp::Function boruvka_mst = package_env["emst"];

      //Call the fast emst function and obtain output as list. Obtain the edge list matrix
      List emst_output = boruvka_mst(mat);
      NumericMatrix edges = emst_output["output"];

      //Transform MST format to adjacency matrix
      tree = to_adj_mat(edges);
   }


   else
   {
     stop("MST not recognized.");
   }

  List td_stats = tree_dimension(tree);
  double td = td_stats["td"];
  double max_td = 1.0 + ((sample_size-3) + ceiling2(sample_size-3,2))/(double)(3);
  double min_td =1.0;

  double effect =(log(max_td) - log(td)) / (log(max_td) - log(min_td));
  double s = effect * sample_size;

  //Statistics to be returned as a list
  List output;
  //Return tree (as edges) if the returnTree argument is true. This is to avoid overhead of unnecessarily returning tree when
  //computing distribution parameters
  if(returnTree)
  {
     NumericVector tree_edges= get_edges(mat,MST)["edges"];
     output = List::create(_["td"]=td, _["stat"]=s, _["effect"]=effect, _["leafs"]=td_stats["leafs"], _["diameter"]=td_stats["diameter"],_["mst"]=tree_edges);
  }

  //Return output without tree (edges)
  else
  {
     output = List::create(_["td"]=td, _["stat"]=s, _["effect"]=effect, _["leafs"]=td_stats["leafs"], _["diameter"]=td_stats["diameter"]);

  }
  return output;
}




/* Computes null distributions for s statistic*/
// [[Rcpp::export]]
NumericVector computeDists(NumericMatrix& data, int perm, int sample_size, Function g, String MST)
{

  NumericVector td_dist(perm);
  NumericMatrix randomized_data;
  List stats;


  for(unsigned int i=0; i<perm; i++)
  {
    randomized_data = g(data.nrow(), data.ncol(),i,perm);
    td_dist[i] = getStatistics(randomized_data,sample_size, MST,false)["stat"];

  }

  return td_dist;
}





// NumericMatrix convert_to_tree(NumericMatrix& mat)
// {
//   int total_nodes = mat.nrow();
//   std::vector<std::vector<int>> tree;
//   NumericMatrix mtree(total_nodes,total_nodes);
//   NumericMatrix dist_mat = calcPairwiseDist(mat);
//   tree = calculate_mst(dist_mat);
//
//   for(int i=0; i<total_nodes; i++)
//   {
//     int list_size = tree[i].size();
//     for(int j=0; j<list_size; j++)
//     {
//       mtree(i,tree[i][j])=1;
//       mtree(tree[i][j],i)=1;
//     }
//   }
//
//
//
//   return mtree;
// }

//Custom dfs for subtree cover
void subtree_dfs(int vertex, List& tree, NumericVector& visited, NumericVector& inSubtree, int & pure_edges, StringVector& labels)
{

  if(!visited[vertex])
  {
    visited[vertex]=true;
  }

  NumericVector adjacent_nodes = tree[vertex];
  for(int i=0; i<adjacent_nodes.size(); i++)
  {
    if(!visited[adjacent_nodes[i]])
    {
      subtree_dfs(adjacent_nodes[i],tree, visited, inSubtree, pure_edges, labels);
      if(inSubtree[adjacent_nodes[i]])
      {
        inSubtree[vertex] = 1;
        if(labels[vertex]==labels[adjacent_nodes[i]])
          pure_edges ++;

      }
    }
  }
}



/*Computes minimum subtree cover for a set of nodes s on a tree */
// [[Rcpp::export]]
List minSubtreeCover (List& tree, NumericVector s, StringVector labels)
{
  int s_size = s.size();
  int tree_size = tree.size();
  int pure_edges = 0;
  NumericVector visited(tree_size);
  NumericVector is_in_subtree(tree_size);

  for(int i=0; i<tree_size;i++)
  {
    visited[i]=0;
    is_in_subtree[i]=0;
  }

  for(int i=0; i<s_size; i++)
  {
    is_in_subtree[s[i]] = 1;

  }

  subtree_dfs(s[0],tree, visited, is_in_subtree, pure_edges, labels);

  List output = List::create(_["subtree_nodes"]= (int)sum(is_in_subtree==1), _["pure_edges"]=pure_edges);
  return output;
  //return sum(is_in_subtree==1);

}


