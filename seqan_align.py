#!/usr/bin/env python
import os
import sys
import errno
import ctypes
import argparse

# Locate the library file.
LIBFILE = 'libseqan_align.so'
script_dir = os.path.dirname(os.path.realpath(__file__))
library_path = os.path.join(script_dir, LIBFILE)
if not os.path.isfile(library_path):
  ioe = IOError('Library file "'+LIBFILE+'" not found.')
  ioe.errno = errno.ENOENT
  raise ioe

seqan_align = ctypes.cdll.LoadLibrary(library_path)


def make_argparser():
  parser = argparse.ArgumentParser(description='Align a set of sequences.')
  parser.add_argument('input', type=argparse.FileType('r'), default=sys.stdin, nargs='?',
    help='Input sequences.')
  return parser


def main(argv):
  parser = make_argparser()
  args = parser.parse_args(argv[1:])
  seqs = []
  for line_raw in args.input:
    line = line_raw.rstrip('\r\n')
    if line.startswith('>'):
      continue
    else:
      seqs.append(line)
  alignment = align(seqs)
  print
  for i, seq in enumerate(alignment):
    print '>seq{}\n{}'.format(i, seq)


def align(seqs):
  num_seqs = len(seqs)
  seqs_c = strlist_to_c(seqs)
  seqan_align.align(seqs_c, num_seqs)
  return seqs_c


def strlist_to_c(strlist):
  c_strs = (ctypes.c_char_p * len(strlist))()
  for i, s in enumerate(strlist):
    c_strs[i] = ctypes.c_char_p(s)
  return c_strs


if __name__ == '__main__':
  sys.exit(main(sys.argv))
