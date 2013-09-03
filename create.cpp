#include <cstdio>
#include <fstream>
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

int main(int argc, char *argv[]) {
        uint nodes, edges, changes, maxtime;
        uint u,v,t,o;
        uint p;
        struct opts opts;
	int optind;
        TGraph tg;
	
	optind = readopts(argc, argv, &opts);
        
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
        tg.set_policy(opts.c);
        tg.create(tgraphreader);
        
        INFO("Saving structure...");
	ofstream f;
	f.open(opts.outfile, ios::binary);
	tg.save(f);
	f.close();
        
        
        
        return 0;
}
