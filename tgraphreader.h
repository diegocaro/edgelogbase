#ifndef TGRAPHREADER_H_
#define TGRAPHREADER_H_

#include <vector>
#include <map>
#include <algorithm>
#include <sys/types.h>

using namespace std;

class TGraphReaderEventList {
public:
        map< uint, vector<uint> > timepoints;
        
        void addEvent(uint v, uint t) {
                timepoints[v].push_back(t);
        }
        
        void sort() {
          map< uint, vector<uint> >::iterator it;
          for (it=timepoints.begin(); it!=timepoints.end(); ++it) {
            std::sort(it->second.begin(), it->second.end());
          }
          
        }
        
};


class TGraphReaderReverseList {
public:
	vector <uint> neighbors;
};


class TGraphReader {
public:
        uint nodes;
        uint edges;
        uint changes;
        uint maxtime;
        
        TGraphReaderEventList* tgraph;
        TGraphReaderReverseList *revgraph;

        TGraphReader(uint n, uint e, uint c, uint t) {
                nodes = n;
                edges = e;
                changes = c;
                maxtime = t;
                
                
                tgraph = new TGraphReaderEventList[nodes];
                revgraph = new TGraphReaderReverseList[nodes];
        }
        
        void addChange(uint u, uint v, uint t) {
                tgraph[u].addEvent(v,t);
        }
        
        void setCapacity(uint u, uint v, uint n) {
          tgraph[u].timepoints[v].reserve(n);
        }

        void addReverseEdge(uint v, uint u) {
        	revgraph[v].neighbors.push_back(u);
        }

        void setCapacityReverse(uint u, uint n) {
          revgraph[u].neighbors.reserve(n);
        }

};



#endif
