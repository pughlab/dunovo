CFLAGS = -Wall -shared -fPIC

all: c seqan

c:
	gcc $(CFLAGS) align.c -o libalign.so
	gcc $(CFLAGS) swalign.c -o libswalign.so -lm
	gcc $(CFLAGS) seqtools.c -o libseqtools.so
	gcc $(CFLAGS) consensus.c -o libconsensus.so

seqan:
	g++ -std=c++14 $(CFLAGS) seqan_align.cpp -o libseqan_align.so

clean:
	rm -f libalign.so libswalign.so libseqtools.so libconsensus.so seqan_align
