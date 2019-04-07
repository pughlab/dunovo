#!/usr/bin/env python
import os
import sys
import errno
import ctypes
import argparse
PY3 = sys.version_info.major >= 3

# Locate the library file.
LIBFILE = 'libconsensus.so'
script_dir = os.path.dirname(os.path.realpath(__file__))
library_path = os.path.join(script_dir, LIBFILE)
if not os.path.isfile(library_path):
  ioe = IOError('Library file "'+LIBFILE+'" not found.')
  ioe.errno = errno.ENOENT
  raise ioe

consensus = ctypes.cdll.LoadLibrary(library_path)
consensus.get_consensus.restype = ctypes.c_char_p
consensus.get_consensus_duplex.restype = ctypes.c_char_p
consensus.build_consensus_duplex_simple.restype = ctypes.c_char_p

ARG_DEFAULTS = {'alignment':sys.stdin}
DESCRIPTION = "Get the consensus of a set of aligned sequences."


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)
  parser.add_argument('alignment', type=argparse.FileType('r'),
    help='The aligned sequences, in FASTA format (but no multi-line sequences).')
  return parser


def main(argv):
  parser = make_argparser()
  args = parser.parse_args(argv[1:])
  sequences = []
  line_num = 0
  for line in args.alignment:
    line_num += 1
    if line_num % 2 == 0:
      sequences.append(line.rstrip('\r\n'))
  cons = get_consensus(sequences)
  print(cons)


# N.B.: The quality scores must be aligned with their accompanying sequences.
def get_consensus(align, quals=[], cons_thres=-1.0, min_reads=0, qual_thres=' ', gapped=False):
  cons_thres_c = ctypes.c_double(cons_thres)
  if PY3:
    qual_thres_val = ord(qual_thres)
  else:
    qual_thres_val = qual_thres
  qual_thres_c = ctypes.c_char(qual_thres_val)
  n_seqs = len(align)
  if gapped:
    gapped_c = 1
  else:
    gapped_c = 0
  assert not quals or len(quals) == n_seqs, 'Different number of sequences and quals.'
  seq_len = None
  for seq in (align + quals):
    if seq_len is None:
      seq_len = len(seq)
    else:
      if seq_len != len(seq):
        raise AssertionError('All sequences in the alignment must be the same length: {}bp != {}bp.'
                             '\nAlignment:\n{}'.format(seq_len, len(seq), '\n'.join(align)))
  align_c = str_pylist_to_str_carray(align, length=n_seqs)
  if quals:
    quals_c = str_pylist_to_str_carray(quals, length=n_seqs)
  else:
    quals_c = 0
  cons = consensus.get_consensus(align_c, quals_c, n_seqs, seq_len, cons_thres_c, min_reads,
                                 qual_thres_c, gapped_c)
  if PY3:
    return str(cons, 'utf8')
  else:
    return cons


# N.B.: The quality scores must be aligned with their accompanying sequences.
def get_consensus_duplex(align1, align2, quals1=[], quals2=[], cons_thres=-1.0, min_reads=0,
                         qual_thres=' ', method='iupac'):
  assert method in ('iupac', 'freq')
  if PY3:
    method_bytes = bytes(method, 'utf8')
    qual_thres_val = ord(qual_thres)
  else:
    method_bytes = method
    qual_thres_val = qual_thres
  qual_thres_c = ctypes.c_char(qual_thres_val)
  cons_thres_c = ctypes.c_double(cons_thres)
  n_seqs1 = len(align1)
  n_seqs2 = len(align2)
  assert (not quals1 and not quals2) or (quals1 and quals2)
  assert not quals1 or len(quals1) == n_seqs1
  assert not quals2 or len(quals2) == n_seqs2
  seq_len = None
  for seq in (align1 + align2 + quals1 + quals2):
    if seq_len is None:
      seq_len = len(seq)
    else:
      assert seq_len == len(seq), 'All sequences in the alignment must be the same length.'
  align1_c = str_pylist_to_str_carray(align1, length=n_seqs1)
  align2_c = str_pylist_to_str_carray(align2, length=n_seqs1)
  if quals1:
    quals1_c = str_pylist_to_str_carray(quals1, length=n_seqs1)
  else:
    quals1_c = 0
  if quals2:
    quals2_c = str_pylist_to_str_carray(quals2, length=n_seqs1)
  else:
    quals2_c = 0
  cons = consensus.get_consensus_duplex(align1_c, align2_c, quals1_c, quals2_c, n_seqs1, n_seqs2,
                                        seq_len, cons_thres_c, min_reads, qual_thres_c, method_bytes)
  if PY3:
    return str(cons, 'utf8')
  else:
    return cons


def build_consensus_duplex_simple(cons1_raw, cons2_raw, gapped=False):
  assert len(cons1_raw) == len(cons2_raw)
  if PY3:
    cons1_bytes = bytes(cons1_raw, 'utf8')
    cons2_bytes = bytes(cons2_raw, 'utf8')
  else:
    cons1_bytes = cons1_raw
    cons2_bytes = cons2_raw
  cons1_c = ctypes.c_char_p(cons1_bytes)
  cons2_c = ctypes.c_char_p(cons2_bytes)
  if gapped:
    gapped_c = 1
  else:
    gapped_c = 0
  cons = consensus.build_consensus_duplex_simple(cons1_c, cons2_c, gapped_c)
  if PY3:
    return str(cons, 'utf8')
  else:
    return cons


def str_pylist_to_str_carray(str_pylist, length=None, encoding='utf8'):
  if length is None:
    length = len(str_pylist)
  str_carray = (ctypes.c_char_p * length)()
  for i, str_raw in enumerate(align):
    if PY3:
      str_bytes = bytes(str_raw, encoding)
    else:
      str_bytes = str_raw
    str_carray[i] = ctypes.c_char_p(str_bytes)
  return str_carray


if __name__ == '__main__':
  sys.exit(main(sys.argv))
