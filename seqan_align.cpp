#include <iostream>
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

    char activity = 's';
    if (activity == 's') {
        // Print unaligned sequences, translating Dna5 ints into chars.
        typedef typename StringSetType<Align<String<Dna5>>>::Type TStringSet;
        TStringSet strings = stringSet(align);
        for (int i = 0; i < nseq; i++) {
            char *s = (char *)toCString(strings[i]);
            char b;
            for (int j = 0; j < 40; j++) {
                switch ((int)s[j]) {
                    case 0:
                        b = 'A';
                        break;
                    case 1:
                        b = 'C';
                        break;
                    case 2:
                        b = 'G';
                        break;
                    case 3:
                        b = 'T';
                        break;
                    case 4:
                        b = 'N';
                        break;
                    default:
                        b = '?';
                }
                std::cout << b;
            }
            std::cout << std::endl;
        }
    } else if (activity == 'r') {
        // Access each row of the alignment.
        typedef typename Row<Align<String<Dna5>>>::Type TRow;
        for (int i = 0; i < nseq; i++) {
            TRow arow = row(align, i);
            // std::cout << getValue(arow, 1) << std::endl;
        }
    }
    // StringSet<String<Dna5> > strings = stringSet(align);
    //     seq[i] = toCString(row(align, i));
    //     seqan::Gaps<seqan::String> arow = row(align, i);
    //     String<char> seq(arow);
    //     std::cout << toCString(seq) << std::endl;

    return seqs;
}

int main(int argc, char *argv[]) {
    for (int i = 1; i < argc; i++) {
        argv[i-1] = argv[i];
    }
    char **aligned_seqs = align(argc-1, argv);
    // for (int i = 0; i < argc-1; i++) {
    //     std::cout << aligned_seqs[i] << std::endl;
    // }
    return 0;
}
