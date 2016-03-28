#ifndef _PPLL_LABELING_PARALLEL_PRUNED_LABELING_H_
#define _PPLL_LABELING_PARALLEL_PRUNED_LABELING_H_

#include <condition_variable>
#include <cilk/cilk.h>

#define ANALYSIS_MODE

const int MAX_THREADS = 82;

class ParallelPrunedLabeling: public Labeling{
public:
	explicit ParallelPrunedLabeling(vector<int> &node_offset, vector<int> adjacencyList, bool ownLabel = 1)
		: Labeling(node_offset, adjacencyList, ownLabel) {

		}

	void BuildLabels() {
		AllocateMemory();

		mutex m, memlock, storelock;
		int currNode = 0;
		int rootsLeft = kNumRoots;

		bool m_init[MAX_THREADS];
		int *m_q[MAX_THREADS];
		bool *m_visi[MAX_THREADS];
		int *m_good[MAX_THREADS];
		int *m_work[MAX_THREADS];

		stack<int> mem_id;
		for (int i = MAX_THREADS - 1; i >= 0; --i) {
			mem_id.push(i);
			m_init[i] = false;
		}

		int *marking = new int[nodes+5];
		cilk_for(int i = 0; i < nodes; ++i)
			marking[i] = false;
		
		bool overload_mode = false;
		vector< pair<int, pair<int,int> > > store;
		int limit = nodes+1;
		rootsLeft++;
			
		cilk_for(int prun = 0; prun < limit; ++prun){

			m.lock();
			int root = currNode;
			int skip = used[root] || (rootsLeft <= 0);
			currNode++; 
			rootsLeft--;
			m.unlock();

			if (skip || root >= nodes) continue;

			memlock.lock();
			int thread_mem_id = mem_id.top();
			mem_id.pop();
			memlock.unlock();
			if (m_init[thread_mem_id] == false) {
				m_init[thread_mem_id] = true;
				m_q[thread_mem_id] = new int[nodes+5];
				m_visi[thread_mem_id] = new bool[nodes+5];
				m_good[thread_mem_id] = new int[nodes+5];
				m_work[thread_mem_id] = new int[nodes+5];

				cilk_for(int i = 0; i < nodes; ++i) {
					m_work[thread_mem_id][i] = MAX_DIST;
					m_visi[thread_mem_id][i] = false;
				}
			}

			int *q = m_q[thread_mem_id];
			bool *visi = m_visi[thread_mem_id];
			int *good = m_good[thread_mem_id];
			int *work = m_work[thread_mem_id];

			int qlo = 0;
			int qsz = 0;
			q[qsz++] = root;
			visi[root] = 1;

			while (!CAS(marking+root, 0, 1));
			for (int i = 0; i < labels[root].size(); ++i)
				work[labels[root][i].f] = labels[root][i].s;
			marking[root] = 0;
			
			for (int d = 0; qlo < qsz; ++d) {
				int mem = qsz;
				for (int i = 0; i < mem-qlo; ++i) good[i] = 1;
				
				cilk_for(int i = qlo; i < mem; ++i) {
					if (used[q[i]]) {good[i-qlo] = 0; continue; }

					if (index_ != NULL) {
						index_t &idx_r = index_[root];
						index_t &idx_v = index_[q[i]];

						for (int k = 0; k < kNumBitParallelRoots; ++k) {
				            int td = idx_r.bpspt_d[k] + idx_v.bpspt_d[k];
				            if (td - 2 <= d) {
				              td +=
				                  (idx_r.bpspt_s[k][0] & idx_v.bpspt_s[k][0]) ? -2 :
				                  ((idx_r.bpspt_s[k][0] & idx_v.bpspt_s[k][1]) |
				                   (idx_r.bpspt_s[k][1] & idx_v.bpspt_s[k][0]))
				                  ? -1 : 0;
				              if (td <= d) {good[i-qlo] = 0; break;}
				            }
			         	}
					}

					if (good[i-qlo]==0) continue;

					while (!CAS(marking+q[i], 0, 1));
					for (int j = 0; j < labels[q[i]].size() && good[i-qlo]; j++) {
						if (d >= work[labels[q[i]][j].f] + labels[q[i]][j].s){ good[i-qlo] = 0; break; }
					}
					if (good[i-qlo])
						labels[q[i]].push_back(label(root, d));
					marking[q[i]] = 0;
				}
				
				
				
				for (int i = qlo; i < mem; ++i) {
					int u = q[i];

					if (!good[i-qlo]) continue;
					

					int cnte = offset[u+1] - offset[u];
					for (int i = 0; i < cnte; ++i) {
						int nxt = adjacency[offset[u] + i];
						if (visi[nxt]) continue;
						visi[nxt] = 1;
						q[qsz++] = nxt;
					}
					
				}
				qlo = mem;
				
			}
			used[root] = 1;



			while (!CAS(marking+root, 0, 1));
			for (int i = 0; i < labels[root].size(); ++i)
				work[labels[root][i].f] = MAX_DIST;
			marking[root] = 0;
			for (int i = 0; i < qsz; ++i)
				visi[q[i]] = false;

			memlock.lock();	
			mem_id.push(thread_mem_id);
			memlock.unlock();
		}
	}
};

#endif  // _PPLL_LABELING_PARALLEL_PRUNED_LABELING_H_