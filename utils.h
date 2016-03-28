#ifndef _PPLL_UTILS_H_
#define _PPLL_UTILS_H_

#include <sys/time.h>
#include <ctime>
#include <bits/stdc++.h>
#include <cilk/cilk.h>
using namespace std;

namespace {
inline int CountAdjacentNodes(int node, vector<int> &offset,
			 	              vector<int> &adjacency) {
	return offset[node+1] - offset[node];
}
}  // namespace

namespace utils {

double GetCurrentTimeSec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

// If (a,b) and (b,a) are present in edges, one of them is removed.
void RemoveDoubleEdges(vector<pair<int,int> > &edges) {
	cilk_for(int i = 0; i < edges.size(); ++i)
		if (edges[i].first > edges[i].second) swap(edges[i].first, edges[i].second);
	sort(edges.begin(), edges.end());
	edges.erase(unique(edges.begin(), edges.end()), edges.end());
}

// For each edge (a,b), an edge (b,a) is added to edges.
// In case (a,b) and (b,a) have been already present, nothing is added.
void AddBackEdges(vector<pair<int,int> > &edges) {
	RemoveDoubleEdges(edges);
	int sz = edges.size();
	for (int i = 0; i < sz; ++i)
		edges.push_back(make_pair(edges[i].second, edges[i].first));
}

// Takes a list of edges and forms an adjacency representation.
// Representation stored within offset and adjacency.
// edge_array should cointain edges of type (a,b) - UNIQUE
void EdgeArrayToAdjacencyGraph(vector<pair<int,int> > &edge_array,
							   vector<int> &offset, vector<int> &adjacency) {
	int nodes = 0;
	for(int i = 0; i < edge_array.size(); ++i) {
		if (nodes < edge_array[i].first) nodes = edge_array[i].first;
		if (nodes < edge_array[i].second) nodes = edge_array[i].second;
	}
	nodes++; // we allow 0 to be used as a label for a node

	vector< pair<int,int> > symm = edge_array;

	sort(symm.begin(), symm.end());
	for (int i = 0; i < symm.size(); ++i)
		adjacency.push_back(symm[i].second);

	int curr_dx = 0;
	for (int i = 0; i < nodes; ++i) {
		int cnt = 0;
		while (curr_dx+cnt < symm.size() && symm[curr_dx+cnt].first == i)
			cnt++;
		offset.push_back(curr_dx);
		curr_dx += cnt;
	}

	offset.push_back(adjacency.size());
}

// Rearranges the graph by sorting the vertex-order by their
// degree property. Stores in inv the inverse mapping of the nodes.
void RearrangeGraph(vector<int> &offset, vector<int> &adjacency, vector<int> &inv) {
	for (int i = 0; i <= offset.size(); ++i)
		inv[i] = i;
	int nodes = offset.size() - 1;
	vector<int> res_adj, res_off;
	vector<pair<int,int> > data;
	vector<int> rinv(nodes, 0);

	for (int i = 0; i < nodes; ++i) 
		data.push_back(make_pair(CountAdjacentNodes(i, offset, adjacency), i));
	sort(data.rbegin(), data.rend());

	for (int i = 0; i < nodes; ++i)
		rinv[data[i].second] = i;

	int curr = 0;
	for (int i = 0; i < data.size(); ++i) {
		inv[i] = data[i].second;
		res_off.push_back(curr);
		for (int k = 0; k < data[i].first; ++k)
			res_adj.push_back( rinv[adjacency[offset[data[i].second] + k]] );
		curr += data[i].first;
	}

	res_off.push_back(curr);
	offset = res_off;
	adjacency = res_adj;
}
}  // namespace utils

#endif  // _PPLL_UTILS_H_