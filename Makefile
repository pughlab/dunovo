CFLAGS = -Wall -shared -fPIC

all: local kalign
.PHONY: all

local:
	gcc $(CFLAGS) align.c -o libalign.so
	gcc $(CFLAGS) swalign.c -o libswalign.so -lm
	gcc $(CFLAGS) seqtools.c -o libseqtools.so
	gcc $(CFLAGS) consensus.c -o libconsensus.so
.PHONY: local

kalign:
	if [ -f kalign/Makefile ]; then make -C kalign; fi
.PHONY: kalign

clean: clean_local clean_kalign
.PHONY: clean

clean_kalign:
	if [ -f kalign/Makefile ]; then make -C kalign clean; fi
.PHONY: clean_kalign

clean_local:
	rm -f libalign.so libswalign.so libseqtools.so libconsensus.so
.PHONY: clean_local
