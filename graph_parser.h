#ifndef _PPLL_GRAPH_PARSER_H_
#define _PPLL_GRAPH_PARSER_H_

#include <bits/stdc++.h>
#include <fstream>
using namespace std;

namespace {
void ReadBasicGraph(int &nodes, vector< pair<int,int> > &edges, FILE *in) {
	int m;
	char buff[100];
	fscanf(in,"%d",&m);
	cout << m << endl;
	nodes = 0;
	for (int i = 0; i < m; ++i) {
		int a,b; fscanf(in,"%d%d",&a,&b);
		edges.push_back(make_pair(a,b));
		nodes = max(nodes, a + 1);
		nodes = max(nodes, b + 1);
	}
}
}  // namespace

void ReadGraph(int &nodes, vector< pair<int,int> > &edges, int &num_roots) {
	string form;
	cin >> form;

	num_roots = -1;

	if (form[0] >= '0' && form[0] <= '9') {
		// The number of bit-parallel roots was provided on the first line.
		num_roots = 0;
		for (int i = 0; i < form.size(); ++i)
			num_roots = num_roots * 10 + form[i] - '0';
		cout << num_roots<<endl;
		string inputFile;
		cin >> inputFile;
		FILE *in = fopen(inputFile.c_str(),"r");
		char buff[20];
		fscanf(in,"%s",buff);
		ReadBasicGraph(nodes, edges, in);
	} else if (form == "BASIC") {
		ReadBasicGraph(nodes, edges, stdin);
		cout << nodes << " " << edges.size() << endl;
	} else {
		cout << "Unsupported graph format was given within the input-graph." << endl;
		exit(1);
	}
}

#endif  // _PPLL_GRAPH_PARSER_H_