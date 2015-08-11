#!/usr/bin/env bash
set -eu
export SCRIPT=$(cd "$(dirname "$0")" && pwd -P)/$(basename "$0") 
export SRC=/${SCRIPT%.sh}.cc
export BINDIR=~/.cache/bin
export BIN=$BINDIR/$HOSTNAME/$SRC
cd "$(dirname "$0")" 

make SHELL='bash -eu' -s -f - "$BIN" -j <<'EOF'
-include $(BIN).d
CXXFLAGS= $(if $(DEBUG),-g -O0,-O2) -march=native -pipe -std=c++1y -Wall -Wextra -pedantic -I/usr/local/include/bamtools -I/usr/include/bamtools
LDLIBS= -lbamtools -lboost_program_options -lboost_filesystem -lboost_system -lboost_regex
$(BIN) : $(SRC) $(SCRIPT)
	mkdir -m 700 -p "$(BINDIR)" "$$(dirname "$@")" && \
	TMP=$$(TMPDIR=$$(dirname "$@") mktemp -t "tmpXXXXXX") && \
	$(LINK.cc) -MMD -MP -MF "$$TMP.d" -MQ "$@" "$<" $(LDLIBS) -o "$$TMP" && \
	perl -i -pe "s/<stdin>//" "$$TMP.d" && \
	mv "$$TMP.d" "$@.d" && \
	mv "$$TMP" "$@"
EOF

exec -a "$0" "$BIN" "$@"
exit $?
