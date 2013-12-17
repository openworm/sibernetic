#include "owWorldSimulation.h"
#include <stdio.h>

//extern bool load_to_file = false;
extern bool load_from_file = false;
int main(int argc, char **argv)
{
	if(argc == 1)
		run( argc, argv);
	else{
		bool graph = true;
		bool load_to = false;
		for(int i = 1; i<argc; i++){
			if(strncmp(argv[i], "-no_g", 5) == 0)//run without graphics
				graph = false;
			if(strncmp(argv[i], "-l_to", 5) == 0){//run load config to file mode
				graph = false;
				load_to = true;
			}
			if(strncmp(argv[i], "-l_from", 7) == 0){//run load config from file mode
				graph = true;
				load_from_file = true;
			}
		}
		run( argc, argv, graph, load_to);
	}
	return 0;
}
