#ifndef _PPLL_LABELING_H_
#define _PPLL_LABELING_H_

#include <xmmintrin.h>
#include <bits/stdc++.h>
#include <atomic>
#include <assert.h>
#include <vector>
#include <cilk/cilk.h>
//#include "ligra2/utils.h"
#include "utils.h"

using namespace std;

//60
const int kNumBitParallelRoots = 256;
//const int kNumBitParallelRoots_ARRSZ = 60; // UPPER BOUND
const int MAX_DIST = 250; // should be bigger than the diameter


typedef unsigned char distance_type;

struct label{
	int f;
	distance_type s;
	label(int p1, int p2):f(p1),s(p2){};
};

struct index_t {
	    uint8_t bpspt_d[kNumBitParallelRoots];
	    uint64_t bpspt_s[kNumBitParallelRoots][2];  // [0]: S^{-1}, [1]: S^{0}
	    uint32_t *spt_v;
	    uint8_t *spt_d;
} __attribute__((aligned(64)));  // Aligned for cache lines

class Labeling {
public:
	Labeling() { 
	}	

	// Create a new labeling object over the underlying graph, given by
	// node_offset and adjacency_list. make_own_labels is used to specify
	// if the object is taking it's labels from another object or if its
	// going to create its own labels.
	explicit Labeling(vector<int> &node_offset,
					  vector<int> adjacencyList,
					  bool make_own_labels = 1) {
		labels = NULL;
		index_ = NULL;
		nodes = node_offset.size() - 1;
		edges = adjacencyList.size();
		offset = node_offset;
		adjacency = adjacencyList;
		this->make_own_labels = make_own_labels;
		kNumRoots = nodes;
	}

	virtual void AllocateMemory() {
		if (make_own_labels) {
			labels = new vector< label >[nodes+5];
			used = new bool[nodes];
			cilk_for(int i = 0; i < nodes; ++i)
				used[i] = false;
		}
		work   = new distance_type[nodes+5];
		work_c = new distance_type[nodes+5];
		inv    = new int[nodes+5];
		for (int i = 0; i <= nodes; ++i) {
			work[i] = MAX_DIST;
			work_c[i] = MAX_DIST;
			inv[i] = -1;
		}
	}

	virtual void BuildLabels() { }

	virtual void clearTempMemSpec(){}

	void takeLabelData(Labeling *L) {
		labels = L->labels;
		used = L->used;
		index_ = L->index_;
	}

	double GetAvgLabelSize() {
		if (labels == NULL) return 0;
		int s = 0;
		for (int i = 0; i < nodes; ++i) {
			s += labels[i].size();
		}
		return s / (double)nodes;
	}


	int GetDistanceBitParallel(int v, int w) {
		  const index_t &idx_v = index_[v];
		  const index_t &idx_w = index_[w];
		  int d = MAX_DIST;

		  for (int i = 0; i < kNumBitParallelRoots; ++i) {
		    int td = idx_v.bpspt_d[i] + idx_w.bpspt_d[i];
		    if (td - 2 <= d) {
		      td +=
		          (idx_v.bpspt_s[i][0] & idx_w.bpspt_s[i][0]) ? -2 :
		          ((idx_v.bpspt_s[i][0] & idx_w.bpspt_s[i][1]) | (idx_v.bpspt_s[i][1] & idx_w.bpspt_s[i][0]))
		          ? -1 : 0;

		      if (td < d) d = td;
		    }
		  }

		  return d;
	}

	int GetDistance(int n1, int n2) {
		int ret = INT_MAX;
		if (index_ != NULL)
			ret = GetDistanceBitParallel(n1, n2);
		if (labels == NULL)
			return ret;

		for (int i = 0; i < labels[n1].size(); ++i) work_c[labels[n1][i].f] = min(work_c[labels[n1][i].f], labels[n1][i].s);

		for (int i = 0; i < labels[n2].size(); ++i) {
			int node = labels[n2][i].f;
			int d = labels[n2][i].s;
			if (ret > work_c[node] + d)
				ret = work_c[node] + d;
		}

		for (int i = 0; i < labels[n1].size(); ++i) work_c[labels[n1][i].f] = MAX_DIST;	
		if (ret >= MAX_DIST) return INT_MAX;
		return ret;
	}

	bool PruneByDistance(int n1, int n2, int d) {
		for (int i = 0; i < labels[n1].size(); ++i) work_c[labels[n1][i].f] = min(work_c[labels[n1][i].f], labels[n1][i].s);

		int ret = INT_MAX;
		for (int i = 0; i < labels[n2].size(); ++i) {
			if (work_c[labels[n2][i].f] + labels[n2][i].s <= d) {
				for (int i = 0; i < labels[n1].size(); ++i) work_c[labels[n1][i].f] = MAX_DIST;	
				return 1;
			}
		}

		for (int i = 0; i < labels[n1].size(); ++i) work_c[labels[n1][i].f] = MAX_DIST;	
		return false;
	}

	int QueryDistance(int n1, int n2) {
		return GetDistance(n1, n2);
	}

	vector< label > getLabels(int node) {
		return labels[node];
	}

	void clearTempMem() {
		delete [] work;
		delete [] work_c;
		delete [] inv;
	}

	void clear() {
		clearTempMem();
		delete [] labels;
	}	

	void fixVertexOrder() {
		cilk_for (int i = 0; i < nodes; ++i) {
			for (int j = 0; j < labels[i].size(); ++j) {
				int dx = j;
				while (dx > 0 && labels[i][dx].f < labels[i][dx-1].f) {
					swap(labels[i][dx], labels[i][dx-1]);
					dx--;
				}
			}
		}
	}

	
	int QueryDistanceCacheEfficient(int v, int w) {
	  const index_t &idx_v = index_[v];
	  const index_t &idx_w = index_[w];
	  int d = MAX_DIST;

	  for (int i = 0; i < kNumBitParallelRoots; ++i) {
	    int td = idx_v.bpspt_d[i] + idx_w.bpspt_d[i];
	    if (td - 2 <= d) {
	      td +=
	          (idx_v.bpspt_s[i][0] & idx_w.bpspt_s[i][0]) ? -2 :
	          ((idx_v.bpspt_s[i][0] & idx_w.bpspt_s[i][1]) | (idx_v.bpspt_s[i][1] & idx_w.bpspt_s[i][0]))
	          ? -1 : 0;

	      if (td < d) d = td;
	    }
	  }

	  int it1 = 0, it2 = 0;
	  int sz1 = labels[v].size();
	  int sz2 = labels[w].size();
	  label *arr_v = labels[v].data();
	  label *arr_w = labels[w].data();
	  for ( ; ; ) {
		  	if (it1 == sz1 || it2 == sz2) break;
		    int v1 = arr_v[it1].f, v2 = arr_w[it2].f;
		    if (v1 == v2) {
		      int td = arr_v[it1].s + arr_w[it2].s;
		      if (td < d) d = td;
		      it1++;;
		      it2++;
		    } else {
		      if (v1 < v2) it1++;
		      if (v1 > v2) it2++;
		    }
		}
		if (d == MAX_DIST) return INT_MAX;
		return d;
	}

protected:
	int nodes;
	int edges;
	bool make_own_labels;
	bool *used;
	vector<int>offset, adjacency;
	vector< label > *labels;
	index_t *index_;
	distance_type *work,*work_c;
	int *inv;
	int kNumRoots;
	double initTime;
};

#endif  // _PPLL_LABELING_H_


































































































































































































































































































































































































