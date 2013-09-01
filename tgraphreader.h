#ifndef TGRAPHREADER_H_
#define TGRAPHREADER_H_

#include <vector>
#include <map>
#include <sys/types.h>

using namespace std;

class TGraphReaderEventList {
public:
        map< uint, vector<uint> > timepoints;
        
        void addEvent(uint v, uint t) {
                timepoints[v].push_back(t);
        }
};



class TGraphReader {
public:
        uint nodes;
        uint edges;
        uint changes;
        uint maxtime;
        
        TGraphReaderEventList* tgraph;


        TGraphReader(uint n, uint e, uint c, uint t) {
                nodes = n;
                edges = e;
                changes = c;
                maxtime = t;
                
                
                tgraph = new TGraphReaderEventList[nodes];
        }
        
        void addChange(uint u, uint v, uint t) {
                tgraph[u].addEvent(v,t);
        }

};



#endif