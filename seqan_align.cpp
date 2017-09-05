#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;

int align(int nseq, char *seqs[]) {

    Align<String <Dna5> > align;
    resize(rows(align), nseq);
    for (int i = 0; i < nseq; i++) {
        assignSource(row(align, i), seqs[i]);
    }

    //TODO: Check the scoring parameters are appropriate.
    Score<int> score(0, -1, -1, -2);
    globalMsaAlignment(align, score);
    std::cout << align << std::endl;

    return 0;
}

int main(int argc, char *argv[]) {
    for (int i = 1; i < argc; i++) {
        argv[i-1] = argv[i];
    }
    align(argc-1, argv);
    return 0;
}
