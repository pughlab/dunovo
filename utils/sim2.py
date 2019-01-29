#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
import re
import os
import sys
import bisect
import random
import logging
import numbers
import tempfile
import argparse
import subprocess
script_path = os.path.realpath(__file__)
root_dir = os.path.dirname(os.path.dirname(script_path))
sys.path.insert(0, root_dir)
from bfx import getreads
from pcr import sim as pcr

PY3 = sys.version_info[0] >= 3

WGSIM_ID_REGEX = r'^(.+)_(\d+)_(\d+)_\d+:\d+:\d+_\d+:\d+:\d+_([0-9a-f]+)/[12]$'
USAGE = """%(prog)s [options] ref.fa [--frag-file frags.fq] -1 reads_1.fa -2 reads_2.fa
or     %(prog)s [options] ref.fa --stdout > reads.fa
or     %(prog)s [options] --frag-file frags.fq -1 reads_1.fa -2 reads_2.fa"""
DESCRIPTION = """Simulate a duplex sequencing experiment."""

RAW_DISTRIBUTION = (
  #  0     1     2     3     4     5     6     7     8     9
  # Low singletons, but then constant drop-off. From pML113 (see 2015-09-28 report).
  #  0,  100,   36,   31,   27,   22,   17,   12,    7,  4.3,
  #2.4,  1.2,  0.6,  0.3,  0.2, 0.15,  0.1, 0.07, 0.05, 0.03,
  # High singletons, but then a second peak around 10. From Christine plasmid (2015-10-06 report).
  #    0,  100, 5.24, 3.67, 3.50, 3.67, 3.85, 4.02, 4.11, 4.20,
  # 4.17, 4.10, 4.00, 3.85, 3.69, 3.55, 3.38, 3.15, 2.92, 2.62,
  # 2.27, 2.01, 1.74, 1.56, 1.38, 1.20, 1.02, 0.85,
  # Same as above, but low singletons, 2's, and 3's (rely on errors to fill out those).
     0,    1,    2,    3, 3.50, 3.67, 3.85, 4.02, 4.11, 4.20,
  4.17, 4.10, 4.00, 3.85, 3.69, 3.55, 3.38, 3.15, 2.92, 2.62,
  2.27, 2.01, 1.74, 1.56, 1.38, 1.20, 1.02, 0.85,
)


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION, usage=USAGE)
  io = parser.add_argument_group('I/O')
  io.add_argument('ref', metavar='ref.fa', nargs='?',
    help='Reference sequence. Omit if giving --frag-file.')
  io.add_argument('-1', '--reads1', type=argparse.FileType('w'),
    help='Write final mate 1 reads to this file.')
  io.add_argument('-2', '--reads2', type=argparse.FileType('w'),
    help='Write final mate 2 reads to this file.')
  io.add_argument('-o', '--out-format', choices=('fastq', 'fasta'), default='fasta',
    help='Default: %(default)s')
  io.add_argument('--stdout', action='store_true',
    help='Print interleaved output reads to stdout.')
  io.add_argument('-m', '--mutations', type=argparse.FileType('w'),
    help='Write a log of the PCR and sequencing errors introduced to this file. Will overwrite any '
         'existing file at this path.')
  io.add_argument('-b', '--barcodes', type=argparse.FileType('w'),
    help='Write a log of which barcodes were ligated to which fragments. Will overwrite any '
         'existing file at this path.')
  io.add_argument('--frag-file',
    help='The path of the FASTQ file of fragments. If --ref is given, these will be generated with '
         'wgsim and kept (normally a temporary file is used, then deleted). Note: the file will be '
         'overwritten! If --ref is not given, then this should be a file of already generated '
         'fragments, and they will be used instead of generating new ones.')
  io.add_argument('-Q', '--fastq-qual', default='I',
    help='The quality score to assign to all bases in FASTQ output. Give a character or PHRED '
         'score (integer). A PHRED score will be converted using the Sanger offset (33). Default: '
         '"%(default)s"')
  params = parser.add_argument_group('Simulation Parameters')
  params.add_argument('-n', '--n-frags', type=int, default=1000,
    help='The number of original fragment molecules to simulate. The final number of reads will be '
         'this multiplied by the average number of reads per family. If you provide fragments with '
         '--frag-file, the script will still only read in the number specified here. Default: '
         '%(default)s')
  params.add_argument('-r', '--read-len', type=int, default=100,
    help='Default: %(default)s')
  params.add_argument('-f', '--frag-len', type=int, default=400,
    help='Default: %(default)s')
  params.add_argument('-s', '--seq-error', type=float, default=0.001,
    help='Sequencing error rate per base (0-1 proportion, not percent). Default: %(default)s')
  params.add_argument('-p', '--pcr-error', type=float, default=0.001,
    help='PCR error rate per base (0-1 proportion, not percent). Default: %(default)s')
  params.add_argument('-c', '--cycles', type=int, default=25,
    help='Number of PCR cycles to simulate. Default: %(default)s')
  params.add_argument('-e', '--efficiency-decline', type=float, default=1.05,
    help='Rate at which the PCR replication efficiency declines.')
  params.add_argument('-i', '--indel-rate', type=float, default=0.15,
    help='Fraction of errors which are indels. Default: %(default)s')
  params.add_argument('-E', '--extension-rate', dest='ext_rate', type=float, default=0.3,
    help='Probability an indel is extended. Default: %(default)s')
  params.add_argument('-B', '--bar-len', type=int, default=12,
    help='Length of the barcodes to generate. Default: %(default)s')
  params.add_argument('-I', '--invariant', default='TGACT',
    help='The invariant linker sequence between the barcode and sample sequence in each read. '
         'Default: %(default)s')
  log = parser.add_argument_group('Logging')
  log.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  log.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  log.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  log.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  misc = parser.add_argument_group('Misc')
  misc.add_argument('-S', '--seed', type=int,
    help='Random number generator seed. By default, a random, 32-bit seed will be generated and '
         'logged to stdout.')
  return parser


def main(argv):
  # Parse and interpret arguments.
  parser = make_argparser()
  args = parser.parse_args(argv[1:])
  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()
  if not (args.ref or args.frag_file):
    parser.print_usage()
    fail('You must provide either a reference or fragments file.')
  if args.ref:
    if not os.path.isfile(args.ref):
      fail('Error: reference file {!r} not found.'.format(args.ref))
    if not os.path.getsize(args.ref):
      fail('Error: reference file {!r} empty (0 bytes).'.format(args.ref))
  else:
    if not (args.reads1 and args.reads2):
      fail('Error: must provide output --reads1 and --reads2 files.')
  if args.seed is None:
    seed = random.randint(0, 2**31-1)
    logging.info('seed: {}\n'.format(seed))
  else:
    seed = args.seed
  random.seed(seed)
  if args.stdout:
    reads1 = sys.stdout
    reads2 = sys.stdout
  else:
    reads1 = args.reads1
    reads2 = args.reads2
  if isinstance(args.fastq_qual, numbers.Integral):
    assert args.fastq_qual >= 0, '--fastq-qual cannot be negative.'
    fastq_qual = chr(args.fastq_qual + 33)
  elif isinstance(args.fastq_qual, basestring):
    assert len(args.fastq_qual) == 1, '--fastq-qual cannot be more than a single character.'
    fastq_qual = args.fastq_qual
  else:
    raise AssertionError('--fastq-qual must be a positive integer or single character.')
  qual_line = fastq_qual * args.read_len

  invariant_rc = pcr.get_revcomp(args.invariant)

  # Create a temporary directory to do our work in. Then work inside a try so we can finally remove
  # the directory no matter what exceptions are encountered.
  tmpfile = tempfile.NamedTemporaryFile(prefix='wgdsim.frags.', delete=False)
  tmpfile.close()
  try:
    # Step 1: Use wgsim to create fragments from the reference.
    if args.frag_file:
      frag_path = args.frag_file
    else:
      frag_path = tmpfile.name
    if args.ref:
      #TODO: Check exit status
      #TODO: Check for wgsim on the PATH.
      # Set error and mutation rates to 0 to just slice sequences out of the reference without
      # modification.
      run_command('wgsim', '-e', '0', '-r', '0', '-d', '0', '-R', args.indel_rate, '-S', seed,
                  '-N', args.n_frags, '-X', args.ext_rate, '-1', args.frag_len,
                  args.ref, frag_path, os.devnull)

    # NOTE: Coordinates here are 0-based (0 is the first base in the sequence).
    extended_dist = extend_dist(RAW_DISTRIBUTION)
    proportional_dist = compile_dist(extended_dist)
    n_frags = 0
    for raw_fragment in getreads.getparser(frag_path, filetype='fastq'):
      n_frags += 1
      if n_frags > args.n_frags:
        break
      chrom, id_num, start, stop = parse_read_id(raw_fragment.id)
      barcode1 = pcr.get_rand_seq(args.bar_len)
      barcode2 = pcr.get_rand_seq(args.bar_len)
      barcode2_rc = pcr.get_revcomp(barcode2)
      #TODO: Vary the size of the fragment.
      #      Could add ~100bp to frag_len arg to wgsim, then randomly select a subsequence here.
      raw_frag_full = barcode1 + args.invariant + raw_fragment.seq + invariant_rc + barcode2

      # Step 2: Determine how many reads to produce from each fragment.
      # - Use random.random() and divide the range 0-1 into segments of sizes proportional to
      #   the likelihood of each family size.
      # bisect.bisect() finds where an element belongs in a sorted list, returning the index.
      # proportional_dist is just such a sorted list, with values from 0 to 1.
      n_reads = bisect.bisect(proportional_dist, random.random())

      # Step 3: Introduce PCR errors.
      # - Determine the mutations and their frequencies.
      #   - Could get frequency from the cycle of PCR it occurs in.
      #     - Important to have PCR errors shared between reads.
      # - For each read, determine which mutations it contains.
      #   - Use random.random() < mut_freq.
      tree = pcr.build_good_pcr_tree(args.cycles, n_reads, args.efficiency_decline, 1000)
      # Add errors to all children of original fragment.
      subtree1 = tree.child1
      subtree2 = tree.child2
      #TODO: Only simulate errors on portions of fragment that will become reads.
      frag_len = len(raw_frag_full)
      pcr.add_pcr_errors(subtree1, '+', frag_len, args.pcr_error, args.indel_rate, args.ext_rate)
      pcr.add_pcr_errors(subtree2, '-', frag_len, args.pcr_error, args.indel_rate, args.ext_rate)
      pcr.apply_pcr_errors(tree, raw_frag_full)
      fragments = pcr.get_final_fragments(tree)
      pcr.add_mutation_lists(tree, fragments, [])

      # Step 4: Introduce sequencing errors.
      for fragment in fragments.values():
        for mutation in pcr.generate_mutations(args.read_len, args.seq_error, args.indel_rate,
                                               args.ext_rate):
          fragment['mutations'].append(mutation)
          fragment['seq'] = pcr.apply_mutation(mutation, fragment['seq'])

      # Print barcodes to log file.
      if args.barcodes:
        args.barcodes.write('{}-{}\t{}\t{}\n'.format(chrom, id_num, barcode1, barcode2_rc))
      # Print family.
      for frag_id in sorted(fragments.keys()):
        fragment = fragments[frag_id]
        read_id = '{}-{}-{}'.format(chrom, id_num, frag_id)
        # Print mutations to log file.
        if args.mutations:
          read1_muts = pcr.get_mutations_subset(fragment['mutations'], 0, args.read_len)
          read2_muts = pcr.get_mutations_subset(fragment['mutations'], 0, args.read_len,
                                                revcomp=True, seqlen=len(fragment['seq']))
          if fragment['strand'] == '-':
            read1_muts, read2_muts = read2_muts, read1_muts
          pcr.log_mutations(args.mutations, read1_muts, read_id+'/1', chrom, start, stop)
          pcr.log_mutations(args.mutations, read2_muts, read_id+'/2', chrom, start, stop)
        frag_seq = fragment['seq']
        read1_seq = frag_seq[:args.read_len]
        read2_seq = pcr.get_revcomp(frag_seq[len(frag_seq)-args.read_len:])
        if fragment['strand'] == '-':
          read1_seq, read2_seq = read2_seq, read1_seq
        if args.out_format == 'fasta':
          reads1.write('>{}\n{}\n'.format(read_id, read1_seq))
          reads2.write('>{}\n{}\n'.format(read_id, read2_seq))
        elif args.out_format == 'fastq':
          reads1.write('@{}\n{}\n+\n{}\n'.format(read_id, read1_seq, qual_line))
          reads2.write('@{}\n{}\n+\n{}\n'.format(read_id, read2_seq, qual_line))

  finally:
    try:
      os.remove(tmpfile.name)
    except OSError:
      pass


def run_command(*command, **kwargs):
  """Run a command and return the exit code.
  run_command('echo', 'hello')
  If "echo" keyword argument is set to True, this will print the command to stdout first."""
  command_strs = map(str, command)
  command_line = '$ '+' '.join(command_strs)+'\n'
  logging.info(command_line)
  if kwargs.get('echo'):
    print(command_line)
  devnull = open(os.devnull, 'w')
  try:
    exit_status = subprocess.call(map(str, command), stderr=devnull)
  except OSError:
    exit_status = None
  finally:
    devnull.close()
  return exit_status


def extend_dist(raw_dist, exponent=1.25, min_prob=0.00001, max_len_mult=2):
  """Add an exponentially decreasing tail to the distribution.
  It takes the final value in the distribution and keeps dividing it by
  "exponent", adding each new value to the end. It will not add probabilities
  smaller than "min_prob" or extend the length of the list by more than
  "max_len_mult" times."""
  extended_dist = list(raw_dist)
  final_sum = sum(raw_dist)
  value = raw_dist[-1]
  value /= exponent
  while value/final_sum >= min_prob and len(extended_dist) < len(raw_dist)*max_len_mult:
    extended_dist.append(value)
    final_sum += value
    value /= exponent
  return extended_dist


def compile_dist(raw_dist):
  """Turn the human-readable list of probabilities defined at the top into
  proportional probabilities.
  E.g. [10, 5, 5] -> [0.5, 0.75, 1.0]"""
  proportional_dist = []
  final_sum = sum(raw_dist)
  current_sum = 0
  for magnitude in raw_dist:
    current_sum += magnitude
    proportional_dist.append(current_sum/final_sum)
  return proportional_dist


def parse_read_id(read_id):
  match = re.search(WGSIM_ID_REGEX, read_id)
  if match:
    chrom = match.group(1)
    start = match.group(2)
    stop = match.group(3)
    id_num = match.group(4)
  else:
    chrom, id_num, start, stop = read_id, None, None, None
  return chrom, id_num, start, stop


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')

if __name__ == '__main__':
  sys.exit(main(sys.argv))
