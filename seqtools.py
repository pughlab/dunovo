import os
import sys
import errno
import ctypes
PY3 = sys.version_info.major >= 3

# Locate the library file.
LIBFILE = 'libseqtools.so'
script_dir = os.path.dirname(os.path.realpath(__file__))
library_path = os.path.join(script_dir, LIBFILE)
if not os.path.isfile(library_path):
  ioe = IOError('Library file "'+LIBFILE+'" not found.')
  ioe.errno = errno.ENOENT
  raise ioe

seqtools = ctypes.cdll.LoadLibrary(library_path)
seqtools.get_revcomp.restype = ctypes.c_char_p
seqtools.transfer_gaps.restype = ctypes.c_char_p


def get_revcomp(seq_raw):
  if PY3:
    seq_bytes = bytes(seq_raw, 'utf8')
  else:
    seq_bytes = bytes(seq_raw)
  return seqtools.get_revcomp(seq_bytes)


def get_diffs_frac_simple(consensus, family):
  consensus_c = pystr_to_cstr(consensus)
  family_c = str_pylist_to_str_carray(family)
  seqtools.get_diffs_frac_simple.restype = ctypes.POINTER(ctypes.c_double * len(family))
  diffs = seqtools.get_diffs_frac_simple(consensus_c, family_c, len(family))
  return tuple(diffs.contents)


def get_diffs_frac_binned(consensus, family, bins):
  seq_len = None
  consensus_c = pystr_to_cstr(consensus)
  family_c = (ctypes.c_char_p * len(family))()
  for i, seq in enumerate(family):
    if seq_len:
      if seq_len != len(seq):
        return None
    else:
      seq_len = len(seq)
    family_c[i] = pystr_to_cstr(seq)
  double_array_pointer = ctypes.POINTER(ctypes.c_double * bins)
  seqtools.get_diffs_frac_binned.restype = ctypes.POINTER(double_array_pointer * len(family))
  diffs_binned_c = seqtools.get_diffs_frac_binned(consensus_c, family_c, len(family), seq_len, bins)
  diffs_binned = []
  for diffs_c in diffs_binned_c.contents:
    diffs_binned.append(diffs_c.contents)
  return diffs_binned


def transfer_gaps(aligned, seq, gap_char_in='-', gap_char_out='-'):
  aligned_c = pystr_to_cstr(aligned)
  seq_c = pystr_to_cstr(seq)
  if PY3:
    gap_char_in_bytes = bytes(gap_char_in, 'utf8')
    gap_char_out_bytes = bytes(gap_char_out, 'utf8')
  else:
    gap_char_in_bytes = bytes(gap_char_in)
    gap_char_out_bytes = bytes(gap_char_out)
  gap_char_in_c = ctypes.c_char(gap_char_in_bytes)
  gap_char_out_c = ctypes.c_char(gap_char_out_bytes)
  seq_aligned = seqtools.transfer_gaps(aligned_c, seq_c, gap_char_in_c, gap_char_out_c)
  if PY3:
    return str(seq_aligned, 'utf8')
  else:
    return str(seq_aligned)


def transfer_gaps_multi(seqs, aligned, gap_char_in='-', gap_char_out='-'):
  if PY3:
    gap_char_in_bytes = bytes(gap_char_in, 'utf8')
    gap_char_out_bytes = bytes(gap_char_out, 'utf8')
  else:
    gap_char_in_bytes = bytes(gap_char_in)
    gap_char_out_bytes = bytes(gap_char_out)
  gap_char_in_c = ctypes.c_char(gap_char_in_bytes)
  gap_char_out_c = ctypes.c_char(gap_char_out_bytes)
  n_seqs = len(seqs)
  assert n_seqs == len(aligned), ('Unequal number of gapped and ungapped sequences ({} vs {} '
                                  'sequences, respectively)'.format(len(aligned), n_seqs))
  seqs_c = str_pylist_to_str_carray(seqs, length=n_seqs)
  aligned_c = str_pylist_to_str_carray(aligned, length=n_seqs)
  seqtools.transfer_gaps_multi.restype = ctypes.POINTER(ctypes.c_char_p * n_seqs)
  output_c = seqtools.transfer_gaps_multi(n_seqs, aligned_c, seqs_c, gap_char_in_c, gap_char_out_c)
  output = []
  for seq_raw in output_c.contents:
    if PY3:
      seq = str(seq_raw, 'utf8')
    else:
      seq = str(seq_raw)
    output.append(seq)
  return output


def str_pylist_to_str_carray(str_pylist, length=None, encoding='utf8'):
  if length is None:
    length = len(str_pylist)
  str_carray = (ctypes.c_char_p * length)()
  for i, pystr in enumerate(str_pylist):
    str_carray[i] = pystr_to_cstr(pystr, encoding=encoding)
  return str_carray


def pystr_to_cstr(pystr, encoding='utf8'):
  if PY3:
    pystr_bytes = bytes(pystr, encoding)
  else:
    pystr_bytes = bytes(pystr)
  return ctypes.c_char_p(pystr_bytes)
