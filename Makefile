BIN = bam2bedgraph
CXXFLAGS= $(if $(DEBUG),-g -O0,-O2) -march=native -pipe -std=c++1y -Wall -Wextra -pedantic -I/usr/local/include/bamtools -I/usr/include/bamtools
LDLIBS= -lbamtools -lboost_program_options -lboost_filesystem -lboost_system -lboost_regex
-include $(BIN).d

$(BIN) : $(BIN).cc
	$(LINK.cc) $^ $(LDLIBS) -o $@
