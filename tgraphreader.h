#ifndef TGRAPHREADER_H_
#define TGRAPHREADER_H_

//#include <vector>
//#include <map>
#include <algorithm>
#include <sys/types.h>

#include "btree_set.h"
#include "btree_map.h"
using namespace btree;
using namespace std;

class TGraphReader {
public:
        uint nodes;
        uint edges;
        uint changes;
        uint maxtime;
        
       	//TGraphReaderEventList* tgraph;
		btree_map< uint, btree_map< uint, btree_set<uint> > > tgraph;
        //TGraphReaderReverseList *revgraph;
		btree_map< uint, btree_set <uint> > revgraph;

        TGraphReader(uint n, uint e, uint c, uint t) {
                nodes = n;
                edges = e;
                changes = c;
                maxtime = t;
        }
        
        void addChange(uint u, uint v, uint t) {
            pair<btree_set<uint>::iterator, bool> info;

            info = tgraph[u][v].insert(t);

            // por si el datasets tiene intervalos de la forma I1 = [a,b) I2 = [b, c)
            // Guardamos solo el intervalo [a,c) ... (solo pasa si est‡ corrupto el dataset)
            if (info.second == false) {
                tgraph[u][v].erase(info.first);
            }

        }
        
        void addReverseEdge(uint v, uint u) {
        	revgraph[v].insert(u);
        }
};



#endif
