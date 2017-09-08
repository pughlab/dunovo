#include <iostream>
#include <stdlib.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/sequence.h>
#include <seqan/graph_msa.h>

using namespace seqan;

char **align(int nseq, char *seqs[]) {

    Align<String<Dna5>> align;
    resize(rows(align), nseq);
    for (int i = 0; i < nseq; i++) {
        assignSource(row(align, i), seqs[i]);
    }

    //TODO: Check that the scoring parameters are appropriate.
    globalMsaAlignment(align, EditDistanceScore());

    // Convert the Align rows to char *'s and store back in seqs.
    typedef typename Row<Align<String<Dna5>>>::Type TRow;
    for (int i = 0; i < nseq; i++) {
        // Each row is type TRow, but also functions as a Gaps. This is why isGap accepts it.
        TRow arow = row(align, i);
        int len = (int)length(arow);
        char *new_seq = (char *)malloc(sizeof(char) * len+1);
        int offset = 0;
        for (int j = 0; j < len; j++) {
            if (isGap(arow, j)) {
                new_seq[j] = '-';
                offset--;
            } else {
                new_seq[j] = seqs[i][j+offset];
            }
        }
        new_seq[len] = '\0';
        seqs[i] = new_seq;
    }

    return seqs;
}

int main(int argc, char *argv[]) {
    for (int i = 1; i < argc; i++) {
        argv[i-1] = argv[i];
    }
    char **aligned_seqs = align(argc-1, argv);
    for (int i = 0; i < argc-1; i++) {
        std::cout << aligned_seqs[i] << std::endl;
    }
    return 0;
}
