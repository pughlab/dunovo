CFLAGS = -Wall -shared -fPIC

all:
	gcc $(CFLAGS) align.c -o libalign.so
	gcc $(CFLAGS) swalign.c -o libswalign.so -lm
	gcc $(CFLAGS) seqtools.c -o libseqtools.so
	gcc $(CFLAGS) consensus.c -o libconsensus.so
	g++ -std=c++11 -Wall seqan_align.cpp -o seqan_align
