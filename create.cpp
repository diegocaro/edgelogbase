#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <assert.h>
#include <list>
#include <unistd.h>

#include "tgraph.h"
#include "tgraphreader.h"
#include "debug.h"

struct opts {
	enum CP_FORMAT c; //bit data structure
	char *outfile;
};

int readopts(int argc, char **argv, struct opts *opts) {
	int o;
	
	
	// Default options
	opts->c = PFOR;

	while ((o = getopt(argc, argv, "c:")) != -1) {
		switch (o) {
			case 'c':
			if(strcmp(optarg, "S9")==0) {
				INFO("Using S9");
				opts->c = S9;
			}
			else if(strcmp(optarg, "S16")==0) {
				INFO("Using S16");
				opts->c = S16;
			}
			else if(strcmp(optarg, "VBYTE")==0) {
				INFO("Using VBYTE");
				opts->c = VBYTE;
			}
			else if(strcmp(optarg, "RICE")==0) {
				INFO("Using RICE");
				opts->c = RICE;
			}
			else if(strcmp(optarg, "PFOR")==0) {
				INFO("Using PFOR");
				opts->c = PFOR;
			}
			break;
			default: /* '?' */
			break;
		}
	}
	
        if (optind >= argc || (argc-optind) < 1) {
		fprintf(stderr, "%s [-c S9,S16,VBYTE,RICE,PFOR] <outputfile> \n", argv[0]);
		fprintf(stderr, "Expected argument after options\n");
		exit(EXIT_FAILURE);
        }
	
	opts->outfile = argv[optind];
	
	return optind;

}

TGraphReader* readcontacts() {
	uint nodes, edges, lifetime, contacts;
	uint u,v,a,b;

	vector < map<uint, list<uint> > > btable;
	vector < set<uint> > revgraph; //edges are in the form v,u
	scanf("%u %u %u %u", &nodes, &edges, &lifetime, &contacts);

	for(uint i = 0; i < nodes; i++) {
		map<uint, list<uint> > t;
		btable.push_back(t);

		set<uint> s;
		revgraph.push_back(s);
	}

	uint c_read = 0;
	while( EOF != scanf("%u %u %u %u", &u, &v, &a, &b)) {
		c_read++;
		if(c_read%500000==0) fprintf(stderr, "Processing %.1f%%\r", (float)c_read/contacts*100);

		btable[u][a].push_back(v);

		//reverse node
		revgraph[v].insert(u);

		if (b == lifetime-1) continue;

		btable[u][b].push_back(v);
	}
	fprintf(stderr, "Processing %.1f%%\r", (float)c_read/contacts*100);
	assert(c_read == contacts);


	TGraphReader *tgraphreader = new TGraphReader(nodes,edges,2*contacts,lifetime);


	map<uint, list<uint> >::iterator it;
	// temporal graph
	for(uint i = 0; i < nodes; i++) {
		if(i%10000==0) fprintf(stderr, "Copying temporal %.1f%%\r", (float)i/nodes*100);

		for( it = btable[i].begin(); it != btable[i].end(); ++it) {
			for(list<uint>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2 ) {
				tgraphreader->addChange(i, *it2, it->first);
			}

      //vector<uint>().swap(it->second);
      (it->second).clear();
      
		}
    btable[i].clear(); 
	}
  
  vector < map<uint, list<uint> > >().swap(btable);
  btable.clear();
  
	//reverse neighbors
	set<uint>::iterator its;
	for(uint i = 0; i < nodes; i++) {
		if(i%10000==0) fprintf(stderr, "Copying reverse %.1f%%\r", (float)i/nodes*100);
    
		for( its = revgraph[i].begin(); its != revgraph[i].end(); ++its) {
			tgraphreader->addReverseEdge(i, *its);
		}
    revgraph[i].clear();
	}
  vector < set<uint> >().swap(revgraph);
  revgraph.clear();


	return tgraphreader;
}
int main(int argc, char *argv[]) {
        struct opts opts;
	int optind;
        TGraph tg;
	
	optind = readopts(argc, argv, &opts);
        /*
         *     uint nodes, edges, changes, maxtime;
        uint u,v,t,o;
        uint p;

        //scanf("%d %d %d %d", &nodes, &edges, &changes, &maxtime);
        scanf("%u %u %u", &nodes,&changes, &maxtime);
        
        TGraphReader tgraphreader(nodes,edges,changes,maxtime);
        
        p = 0;
        while ( EOF != scanf("%u %u %u %u", &u, &v, &t, &o) ) {
                if (p%10000==0) fprintf(stderr, "Loading: %0.2f%%\r", (float)p/changes*100);
                p++;
                
                tgraphreader.addChange(u,v,t);
        }
        
        fprintf(stderr, "\n");
        */
	TGraphReader *tgraphreader;

	tgraphreader = readcontacts();

        tg.set_policy(opts.c);
        tg.create(*tgraphreader);
        
        INFO("Saving structure...");
	ofstream f;
	f.open(opts.outfile, ios::binary);
	tg.save(f);
	f.close();
        
        
        
        return 0;
}
