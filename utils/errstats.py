#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import re
import os
import sys
import copy
import errno
import string
import random
import logging
import argparse
import collections
# sys.path hack to access lib package in root directory.
sys.path.insert(0, os.path.dirname(sys.path[0]))
# sys.path hack to allow overriding installed pyBamParser.
if os.environ.get('PYTHONPATH'):
  sys.path.insert(1, os.environ.get('PYTHONPATH'))
try:
  import pyBamParser.bam
except ImportError:
  pass
from utillib import simplewrap
from kalign import kalign
import consensus as consensuslib
import seqtools
import swalign
PY3 = sys.version_info.major >= 3


if PY3:
  trans_fxn = str.maketrans
else:
  trans_fxn = string.maketrans
IUPAC_TO_N_TABLE = trans_fxn('rymkbdhvRYMKBDHV', 'NNNNNNNNNNNNNNNN')
REVCOMP_TABLE = trans_fxn('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
GAP_WIN_LEN = 4

DESCRIPTION = """Tally statistics on errors in reads, compared to their (single-stranded) \
consensus sequences. Output is one tab-delimited line per single-read alignment (one mate within \
one strand (order) within one family (barcode)).
A "unique error" is a class of error defined by its reference coordinate and the erroneous base.
A single unique error may occur several times in the same family, if it happened on multiple reads.
The first four output columns are always:
  1. barcode
  2. order
  3. mate (0 or 1)
  4. number of reads (family size)
The following columns depend on the argument to --out-format:
--out-format reads (default):
  5. number of unique errors that were observed in more than one read
  6-end. number of errors in each read
--out-format errors1 (aka --all-repeats):
  5-end. count of how many reads each unique error was observed in
--out-format errors2:
  5. read ids (comma-delimited)
  6. consensus sequence GC content (0-1 proportion)
  7. number of unique errors
  8-end. count of how many reads each unique error was observed in

Format of --overlap-stats is tab-delimited statistics on each mate:
  1. barcode
  2. order
  3. mate
  4. "True"/"False": did we find this read's opposite mate in the alignment?
  5. length of the overlap region
  6. length of the non-overlap region (in this mate)
  7. number of unique errors in the overlap region
  8. number of unique errors outside the overlap, but aligned to the reference
  9. number of unique errors with no reference coordinate
  10. number of unique errors that appeared on both mates in the pair (duplicates)"""


def make_argparser():

  # Need to use argparse.RawDescriptionHelpFormatter to preserve formatting in the
  # description of columns in the tsv output. But to still accommodate different
  # terminal widths, dynamic wrapping with simplewrap will be necessary.
  wrap = simplewrap.Wrapper().wrap
  parser = argparse.ArgumentParser(description=wrap(DESCRIPTION),
                                   formatter_class=argparse.RawDescriptionHelpFormatter)

  parser.add_argument('input', metavar='families.msa.tsv', nargs='?', type=argparse.FileType('r'),
    default=sys.stdin,
    help='Aligned families (output of align_families.py). Omit to read from stdin.')
  parser.add_argument('-f', '--out-format', default='reads', choices=('reads', 'errors1', 'errors2'),
    help='Default: %(default)s')
  parser.add_argument('-R', '--all-repeats', dest='out_format', action='store_const', const='errors1',
    help='Backward compatibility shorthand for "--out-format errors1".')
  parser.add_argument('-D', '--duplex', action='store_true',
    help='Use the full duplex family to determine the consensus sequence, not just the single-'
         'stranded family. In this case, the second output column can be ignored (it will always '
         'be "ab"). The third column then represents the mate, but the duplex mate, not the input '
         'mate. The reads for both strands will be aligned to the corresponding duplex consensus '
         'at once, and errors calculated for the entire alignment.')
  parser.add_argument('-a', '--alignment', dest='human', action='store_true',
    help='Print human-readable output, including a full alignment with consensus bases masked '
         '(to highlight errors).')
  parser.add_argument('-r', '--min-reads', type=int, default=1,
    help='Minimum number of reads to form a consensus (and thus get any statistics). '
         'Default: %(default)s')
  parser.add_argument('-1', '--mate1', dest='mate_offset', action='store_const', default=0, const=1,
    help='Use 1-based indexing for mate numbering (1 and 2 instead of 0 and 1).')
  #TODO:
  # parser.add_argument('-c', '--cons-thres', type=float, default=0.5)
  parser.add_argument('-q', '--qual-thres', type=int, default=0,
    help='PHRED quality score threshold for consensus making. NOTE: This should be the same as was '
         'used for producing the reads in the bam file, if provided! Default: %(default)s')
  parser.add_argument('-Q', '--qual-errors', action='store_true',
    help='Don\'t count errors with quality scores below the --qual-thres in the error counts.')
  parser.add_argument('-I', '--no-indels', dest='indels', action='store_false', default=True,
    help='Don\'t count indels. Note: the indel counting isn\'t good right now. It counts every '
         'base of the indel as a separate error, so a 3bp deletion counts as 3 errors. '
         'Default is to count them, for backward compatibility. Some day I may remove this footgun.')
  parser.add_argument('-d', '--dedup', action='store_true',
    help='Figure out whether there is overlap between mates in read pairs and deduplicate errors '
         'that appear twice because of it. Requires --bam.')
  parser.add_argument('-b', '--bam',
    help='The final single-stranded consensus reads, aligned to a reference. Used to find overlaps.')
  parser.add_argument('-s', '--seed', type=int, default=0,
    help='The random seed. Used to choose which error to keep when deduplicating errors in '
         'overlaps. Default: %(default)s')
  parser.add_argument('-o', '--overlap-stats', type=argparse.FileType('w'),
    help='Write statistics on overlaps and errors in overlaps to this file. Warning: will '
         'overwrite any existing file.')
  parser.add_argument('-v', '--validate-kalign', action='store_true',
    help='When using --duplex, check that Kalign is returning aligned sequences in the same order '
         'they were given.')
  parser.add_argument('-L', '--dedup-log', type=argparse.FileType('w'),
    help='Log overlap error deduplication to this file. Warning: Will overwrite any existing file.')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-S', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.ERROR)
  parser.add_argument('-V', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('--debug', dest='volume', action='store_const', const=logging.DEBUG)

  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  if args.human and args.out_format == 'errors2':
    fail('Error: --alignment invalid with --out-format errors2.')

  if args.dedup:
    if not args.bam:
      fail('--dedup requires a --bam file to be supplied.')
    try:
      pyBamParser.bam
    except NameError:
      fail('Failed to import pyBamParser (required for --dedup).')

  if args.qual_errors:
    error_qual_thres = args.qual_thres
  else:
    error_qual_thres = 0

  logging.info('Calculating consensus sequences and counting errors..')
  single_strand_families = double_strand_families = 0
  family_stats = {}
  for family in parse_families(args.input):
    barcode = family['bar']
    if args.dedup:
      family_stats[barcode] = {'ab':[{}, {}], 'ba':[{}, {}]}
    add_consensi(family, args.qual_thres)
    if is_double_stranded(family):
      double_strand_families += 1
      if args.duplex:
        transform_to_duplex_family(family, validate=args.validate_kalign)
    else:
      single_strand_families += 1
      if args.duplex:
        continue
    for order in ('ab', 'ba'):
      for mate in (0, 1):
        if args.duplex and order == 'ba':
          # Since ab.1 == ba.2 and ab.2 == ba.1 in duplex families, only process ab.
          # This also means the mate value will still be accurate, meaning the duplex mate.
          continue
        seq_align = family[order][mate]['seqs']
        qual_align = family[order][mate]['quals']
        consensus = family[order][mate]['consensus']
        ids = family[order][mate]['ids']
        num_seqs = len(seq_align)
        # Replace non-"ACGTN-" bases with N, to make it easier to ignore them.
        consensus_mod = consensus.translate(IUPAC_TO_N_TABLE)
        error_types = get_family_errors(seq_align, qual_align, consensus_mod, error_qual_thres,
                                        count_indels=args.indels)
        overlap = collections.defaultdict(int)
        family_stat = {'num_seqs':num_seqs, 'consensus':consensus, 'errors':error_types,
                       'overlap':overlap, 'ids':ids}
        if args.dedup:
          family_stats[barcode][order][mate] = family_stat
        elif num_seqs >= args.min_reads:
          print_errors(barcode, order, mate+args.mate_offset, family_stat, args.out_format,
                       args.human, seq_align, qual_align)

  total = single_strand_families + double_strand_families
  logging.info('Processed {} families: {:0.2f}% single-stranded, {:0.2f}% double-stranded.'
               .format(total, 100*single_strand_families/total, 100*double_strand_families/total))

  if args.dedup:
    logging.info('Deduplicating errors in overlaps..')
    dedup_all_errors(args.bam, family_stats, args.dedup_log)
    for barcode in family_stats:
      for order in ('ab', 'ba'):
        for mate in (0, 1):
          if args.duplex and order == 'ba':
            continue
          family_stat = family_stats[barcode][order][mate]
          if family_stat['num_seqs'] < args.min_reads:
            continue
          print_errors(barcode, order, mate+args.mate_offset, family_stat, args.out_format)
          if args.overlap_stats:
            print_overlap_stats(barcode, order, mate, args.overlap_stats, family_stat['overlap'])


def parse_families(infile):
  """Parse a families.msa.tsv file.
  Yields a data structure for each family:
  family = {
    'bar': barcode,                          # family 'AAACCGACACAGGACTAGGGATCA'
    'ab': (                                    # order ab
            {                                    # mate 1
              'seqs':  [seq1,   seq2,   seq3],     # sequences
              'quals': [quals1, quals2, quals3],   # quality scores
              'ids':   [id1,    id2,    id3],      # read ids
            },
            {                                    # mate 2
              'seqs':  [seq1,   seq2],
              'quals': [quals1, quals2],
              'ids':   [id1,    id2]
            },
          ),
    'ba': (                                    # order ba
            {                                    # mate 1
              'seqs':  [seq1,   seq2],
              'quals': [quals1, quals2]
              'ids':   [id1,    id2]
            },
            {                                    # mate 2
              'seqs':  [],
              'quals': [],
              'ids':   [],
            },
          )
  }
  That is, each family is a dict with the 'bar' key giving the barcode sequence, and a key for both
  orders ('ab' and 'ba'). The value for each order is a tuple of 2 values, one for each mate. Each
  value in the tuple is itself a 2-tuple containing the aligned bases and quality scores.
  Examples:
  Getting the sequences for mate 1 of order "ab":
  seq_align  = family['ab'][0]['seqs']
  Getting the quality scores:
  qual_align = family['ab'][0]['quals']
  Getting the sequences for mate 2 of order "ba":
  seq_align = family['ba'][1]['seqs']
  """
  last_barcode = None
  empty_family = {
    'bar':None,
    'ab': (
      {'seqs':[], 'quals':[], 'ids':[]},
      {'seqs':[], 'quals':[], 'ids':[]}
    ),
    'ba': (
      {'seqs':[], 'quals':[], 'ids':[]},
      {'seqs':[], 'quals':[], 'ids':[]}
    )
  }
  family = copy.deepcopy(empty_family)
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    barcode, order, mate_str, name, seq, quals = fields
    mate = int(mate_str)-1
    fields = name.split()
    if fields and fields[0]:
      read_id = fields[0]
    else:
      read_id = '.'
    if barcode != last_barcode:
      if last_barcode is not None:
        family['bar'] = last_barcode
        yield family
      family = copy.deepcopy(empty_family)
      last_barcode = barcode
    family[order][mate]['seqs'].append(seq)
    family[order][mate]['quals'].append(quals)
    family[order][mate]['ids'].append(read_id)
  family['bar'] = barcode
  yield family


def add_consensi(family, qual_thres):
  for order in ('ab', 'ba'):
    for mate in (0, 1):
      subfamily = family[order][mate]
      subfamily['consensus'] = consensuslib.get_consensus(subfamily['seqs'], subfamily['quals'],
                                                          qual_thres=chr(qual_thres+32), gapped=True)


def get_duplex_consensi(family):
  consensi = []
  for (order1, mate1), (order2, mate2) in (('ab', 0), ('ba', 1)), (('ab', 1), ('ba', 0)):
    alignment = swalign.smith_waterman(family[order1][mate1]['consensus'].replace('-', ''),
                                       family[order2][mate2]['consensus'].replace('-', ''))
    consensi.append(consensuslib.build_consensus_duplex_simple(alignment.query, alignment.target))
  return consensi


def is_double_stranded(family):
  return (family['ab'][0]['seqs'] and family['ba'][1]['seqs'] and
          family['ab'][1]['seqs'] and family['ba'][0]['seqs'])


def transform_to_duplex_family(family, validate=False):
  # Get the duplex consensus sequences (mate 1 and 2), and replace the single-strand consensus
  # sequences with them.
  duplex_consensi = get_duplex_consensi(family)
  duplex_mates = ((('ab', 0), ('ba', 1)), (('ab', 1), ('ba', 0)))
  for mate_num, duplex_mate in enumerate(duplex_mates):
    for order, mate in duplex_mate:
      family[order][mate]['consensus'] = duplex_consensi[mate_num]
  # Realign the reads and quality scores to the duplex consensus sequences.
  for duplex_mate in duplex_mates:
    seq_align = []
    qual_align = []
    ids = []
    for order, mate in duplex_mate:
      seq_align  += family[order][mate]['seqs']
      qual_align += family[order][mate]['quals']
      ids += family[order][mate]['ids']
      consensus = family[order][mate]['consensus']
    consensus, seq_align, qual_align = realign_family_to_consensus(consensus, seq_align, qual_align,
                                                                   validate=validate)
    for order, mate in duplex_mate:
      family[order][mate]['seqs'] = seq_align
      family[order][mate]['quals'] = qual_align
      family[order][mate]['consensus'] = consensus
      family[order][mate]['ids'] = ids


def realign_family_to_consensus(consensus, family, quals, validate=False):
  unaligned_seqs = [consensus]
  unaligned_seqs.extend([seq.replace('-', '') for seq in family])
  aligned_seqs = kalign.align(unaligned_seqs)
  if validate:
    for input, output in zip(unaligned_seqs, aligned_seqs):
      assert input == output.replace('-', ''), (
        "Kalign may have returned sequences in a different order than they were given. Failed on "
        "sequence: {}".format(input)
      )
  aligned_consensus = aligned_seqs[0]
  aligned_family = aligned_seqs[1:]
  aligned_quals = seqtools.transfer_gaps_multi(quals, aligned_family, gap_char_out=' ')
  return aligned_consensus, aligned_family, aligned_quals


def get_gc_content(seq):
  gc = len(re.findall(r'[GC]', seq, re.I))
  return gc/len(seq)


def get_family_errors(seq_align, qual_align, consensus, qual_thres, count_indels=False):
  if not (seq_align and qual_align):
    return None
  errors = get_alignment_errors(consensus, seq_align, qual_align, qual_thres,
                                count_indels=count_indels)
  error_types = group_errors(errors)
  return list(error_types)


def print_errors(barcode, order, mate, family_stat, out_format, human=False,
                 seq_align=None, qual_align=None):
  errors_per_seq, repeated_errors, error_repeat_counts = tally_errors(family_stat['errors'],
                                                                      family_stat['num_seqs'])
  if human:
    masked_alignment = mask_alignment(seq_align, family_stat['errors'])
    for seq, seq_errors in zip(masked_alignment, errors_per_seq):
      print('{} errors: {}'.format(seq, seq_errors))
    if out_format == 'errors1':
      print('{} errors: {}, repeat errors: {}\n'.format(family_stat['consensus'],
                                                        sum(errors_per_seq),
                                                        ', '.join(map(str, error_repeat_counts))))
    else:
      print('{} errors: {}, repeat errors: {}\n'.format(family_stat['consensus'],
                                                        sum(errors_per_seq),
                                                        repeated_errors))
  else:
    fields = [barcode, order, mate, family_stat['num_seqs']]
    if out_format == 'reads':
      fields.append(repeated_errors)
      fields.extend(errors_per_seq)
    elif out_format == 'errors1':
      fields.extend(error_repeat_counts)
    elif out_format == 'errors2':
      fields.append(','.join(family_stat['ids']))
      fields.append('{:0.3f}'.format(get_gc_content(family_stat['consensus'])))
      fields.append(len(error_repeat_counts))
      fields.extend(error_repeat_counts)
    print(*fields, sep='\t')


def print_overlap_stats(barcode, order, mate, stats_fh, stats):
  columns = [barcode, order, mate+1]
  columns.append(stats['paired'])
  columns.append(stats['overlap_len'])
  columns.append(stats['non_overlap_len'])
  columns.append(stats['overlap_errors'])
  non_overlap_errors = stats['total_errors'] - stats['overlap_errors'] - stats['nonref_errors']
  columns.append(non_overlap_errors)
  columns.append(stats['nonref_errors'])
  columns.append(stats['duplicates'])
  stats_fh.write('\t'.join(map(str, columns))+'\n')


def get_alignment_errors(consensus, seq_align, qual_align, qual_thres, count_indels=False):
  qual_thres_char = chr(qual_thres+32)
  num_seqs = len(seq_align)
  errors = []
  last_bases = [None] * num_seqs
  running_indels = [None] * num_seqs
  for coord, (cons_base, bases, quals) in enumerate(zip(consensus, zip(*seq_align), zip(*qual_align))):
    for seq_num, (base, qual) in enumerate(zip(bases, quals)):
      #TODO: Figure out how to deal with quality scores for indels and consensus N's.
      if base != cons_base:
        # Mismatch between the read and consensus, and quality is above the threshold.
        if base == '-':
          # We're in a deletion.
          if not count_indels:
            continue
          running_indel = running_indels[seq_num]
          if running_indel and running_indel['type'] == 'del':
            # We were already tracking a deletion.
            running_indel['alt'] += 1
          else:
            # We weren't already tracking a deletion.
            # Either there was no indel being tracked, or we were tracking an insertion.
            running_indels[seq_num] = {'type':'del', 'seq':seq_num, 'coord':coord, 'alt':1}
            if running_indel and running_indel['type'] == 'ins':
              errors.append(running_indel)
        elif cons_base == '-':
          # We're in an insertion.
          if not count_indels:
            continue
          running_indel = running_indels[seq_num]
          if running_indel and running_indel['type'] == 'ins':
            # We were already tracking an insertion.
            running_indel['alt'] += base
          else:
            # We weren't already tracking an insertion.
            # Either there was no indel being tracked, or we were tracking a deletion.
            running_indels[seq_num] = {'type':'ins', 'seq':seq_num, 'coord':coord, 'alt':base}
            if running_indel and running_indel['type'] == 'del':
              errors.append(running_indel)
        else:
          # We're in an SNV.
          if running_indels[seq_num]:
            # Finish up any indels currently being tracked.
            errors.append(running_indels[seq_num])
            running_indels[seq_num] = None
          if cons_base != 'N' and qual > qual_thres_char:
            errors.append({'type':'SNV', 'seq':seq_num, 'coord':coord+1, 'alt':base})
  for indel in running_indels:
    if indel is not None:
      errors.append(indel)
  return errors


def get_gap_quality_score(quals, indel_coord):
  """Calculate what should be considered the quality score of the indel.
  Currently, this is just the simple average of the quality scores of the two bases bookending the
  indel.
  This should give the same results as the algorithm in `get_gap_qual()` of consensus.c.
  `quals` is the quality scores. Gaps should be indicated by a space character.
  `indel_coord` is any (1-based) coordinate inside the indel (any position where there's a '-')."""
  left_i = right_i = indel_coord-1
  while left_i >= 0 and quals[left_i] == ' ':
    left_i -= 1
  while right_i < len(quals) and quals[right_i] == ' ':
    right_i += 1
  score_sum = 0
  weight_sum = 0
  left_weight = right_weight = GAP_WIN_LEN
  scores = []
  weights = []
  while left_weight > 0 or right_weight > 0:
    if left_weight > 0 and left_i >= 0:
      if quals[left_i] != ' ':
        scores.insert(0, ord(quals[left_i]))
        weights.insert(0, left_weight)
        score_sum += left_weight * ord(quals[left_i])
        weight_sum += left_weight
        left_weight -= 1
      left_i -= 1
    if right_weight > 0 and right_i < len(quals):
      if quals[right_i] != ' ':
        scores.append(ord(quals[right_i]))
        weights.append(right_weight)
        score_sum += right_weight * ord(quals[right_i])
        weight_sum += right_weight
        right_weight -= 1
      right_i += 1
  return int(round(score_sum/weight_sum))


def group_errors(errors):
  """Group errors by coordinate and base."""
  last_error = None
  current_types = []
  for error in sorted(errors, key=lambda error: error['coord']):
    if (last_error is not None and last_error['coord'] == error['coord'] and
        last_error['alt'] == error['alt']):
      current_types.append(error)
    else:
      if current_types:
        yield tuple(current_types)
      current_types = [error]
    last_error = error
  if current_types:
    yield tuple(current_types)


def tally_errors(error_types, num_seqs):
  errors_per_seq = [0] * num_seqs
  repeated_errors = 0
  error_repeat_counts = []
  for error_type in error_types:
    error_repeat_counts.append(len(error_type))
    if len(error_type) > 1:
      repeated_errors += 1
    for error in error_type:
      errors_per_seq[error['seq']] += 1
  return errors_per_seq, repeated_errors, error_repeat_counts


def mask_alignment(seq_alignment, error_types):
  masked_alignment = [['.'] * len(seq) for seq in seq_alignment]
  for error_type in error_types:
    for error in error_type:
      if 'type' not in error or error['type'] == 'SNV':
        masked_alignment[error['seq']][error['coord']-1] = error['alt']
  return [''.join(seq) for seq in masked_alignment]


def dedup_all_errors(bam_path, family_stats, dedup_log):
  pair = [None, None]
  for read in pyBamParser.bam.Reader(bam_path):
    barcode, order, mate = get_read_identifiers(read)
    try:
      pair_stats = family_stats[barcode][order]
    except KeyError:
      fail('Read pair found in BAM but not in alignment:\nbar: {}, order: {}'
           .format(barcode, order))
    pair_stats[mate]['overlap']['found'] = True
    pair_stats[mate]['overlap']['paired'] = False
    # Skip if it's a secondary alignment or a supplementary alignment, or if it's not mapped in
    # the proper pair.
    flags = read.get_flag()
    if flags & (256+2048) or not flags & 2:
      continue
    if pair[mate]:
      # We already have this mate for this pair.
      # We must be on a new pair now, and the matching mate for the last one is missing.
      logging.debug('Failed to complete the pair for {}'.format(pair[mate].get_read_name()))
      pair = [None, None]
      pair[mate] = read
    else:
      other_mate = mate ^ 1
      if pair[other_mate]:
        barcode2, order2, mate2 = get_read_identifiers(pair[other_mate])
        if barcode2 == barcode and order2 == order:
          # It's a matching pair.
          pair_stats[0]['overlap']['paired'] = True
          pair_stats[1]['overlap']['paired'] = True
          pair[mate] = read
          dedup_pair(pair, pair_stats, dedup_log)
          pair = [None, None]
        else:
          # The reads do not match; they're from different pairs.
          # We must be on a new pair now, and the matching mate for the last one is missing.
          logging.debug('Failed to complete the pair for {}.{}'.format(barcode2, order2))
          pair = [None, None]
          pair[mate] = read
      else:
        # The pair is empty ([None, None]).
        # This happens on the first loop, and after a pair has been completed on the previous loop.
        logging.debug('Pair for {} empty ([None, None]).'.format(barcode))
        pair[mate] = read
  if (pair[0] and not pair[1]) or (pair[1] and not pair[0]):
    read = pair[0] or pair[1]
    logging.debug('Failed to complete the pair for {}'.format(read.get_read_name()))


def get_read_identifiers(read):
  name = read.get_read_name()
  barcode, order = name.split('.')
  flags = read.get_flag()
  if flags & 64:
    mate = 0
  elif flags & 128:
    mate = 1
  else:
    raise ValueError('Neither flag 64 nor 128 are set: {}'.format(read.get_flag()))
  return barcode, order, mate


def dedup_pair(pair, pair_stats, dedup_log=None):
  """We've gathered a pair of reads and the errors in them. Now correlate the data between them."""
  #TODO: Make this work for indels.
  dedup_log and dedup_log.write('{} ({} read pairs)\n'.format(pair[0].get_read_name(),
                                                              pair_stats[0]['num_seqs']))
  edges = get_edges(pair)
  overlap_len, non_overlap_lens = get_overlap_len(edges)
  errors_by_ref_coord, nonref_errors = convert_pair_errors(pair, pair_stats)
  count_errors_by_location(errors_by_ref_coord, nonref_errors, edges, pair_stats)
  new_errors_lists = null_duplicate_errors(errors_by_ref_coord, pair_stats, dedup_log)
  if dedup_log and (nonref_errors[0] or nonref_errors[1]):
    log_nonref_errors(nonref_errors, dedup_log)
  dedup_log and dedup_log.write('\n')
  for mate in (0, 1):
    pair_stats[mate]['overlap']['total_errors'] = len(pair_stats[mate]['errors'])
    new_errors_lists[mate].extend(nonref_errors[mate])
    pair_stats[mate]['errors'] = filter(lambda e: e is not None, new_errors_lists[mate])
    pair_stats[mate]['overlap']['non_overlap_len'] = non_overlap_lens[mate]
    pair_stats[mate]['overlap']['overlap_len'] = overlap_len


def convert_pair_errors(pair, pair_stats):
  """Convert the coordinate of each error type into reference coordinates.
  Returns two data structures:
  1. errors which have a reference coordinate
  - A list of two dicts, one per mate.
    - Each dict is indexed by a 2-tuple: (the reference coordinate, the erroneous base)
      - Each value of the dicts is an error_type, as created by group_errors().
  2. errors with no reference coordinate (failed conversion: read.to_ref_coord() returns None)
  - A list of two lists, one per mate.
    - Each value of each list is an error_type."""
  errors_by_ref_coord = [{}, {}]
  nonref_errors = [[], []]
  for mate in (0, 1):
    read = pair[mate]
    error_types = pair_stats[mate]['errors']
    for i, error_type in enumerate(error_types):
      error = error_type[0]
      if 'type' in error and error['type'] != 'SNV':
        #TODO: Make indel-compatible.
        continue
      read_coord = error['coord']
      if read.is_seq_reverse_complement():
        alt = get_revcomp(error['alt'])
      else:
        alt = error['alt']
      ref_coord = read.to_ref_coord(read_coord)
      if ref_coord is None:
        nonref_errors[mate].append(error_type)
      else:
        errors_by_ref_coord[mate][(ref_coord+1, alt)] = error_type
  return errors_by_ref_coord, nonref_errors


def count_errors_by_location(errors_by_ref_coord, nonref_errors, edges, pair_stats):
  for mate in (0, 1):
    for error_type in nonref_errors[mate]:
      pair_stats[mate]['overlap']['nonref_errors'] += 1
    for ref_coord, base in errors_by_ref_coord[mate].keys():
      if (edges[mate]['overlap_start'] and edges[mate]['overlap_end'] and
          edges[mate]['overlap_start'] <= ref_coord <= edges[mate]['overlap_end']):
        pair_stats[mate]['overlap']['overlap_errors'] += 1
      else:
        pair_stats[mate]['overlap']['non_overlap_errors'] += 1


def null_duplicate_errors(errors_by_ref_coord, pair_stats, dedup_log=None):
  """Correlate errors between the two mates of the pair & replace one of each duplicate pair with None.
  Whenever the same error (identified by reference coordinate and base) appears on both reads in the
  pair, randomly choose one and replace it with None.
  Returns a new data structure, a tuple containing two lists, one for each mate.
  Each list is simply a list of error_types (from group_errors()), except some elements which are
  None."""
  new_errors_lists = ([], [])
  errors1 = set(errors_by_ref_coord[0].keys())
  errors2 = set(errors_by_ref_coord[1].keys())
  all_errors = list(errors1.union(errors2))
  for ref_coord, alt in sorted(all_errors):
    these_error_types = [None, None]
    for mate in (0, 1):
      error_type = errors_by_ref_coord[mate].get((ref_coord, alt))
      these_error_types[mate] = error_type
    if these_error_types[0] and these_error_types[1]:
      # The same error occurred at the same position in both reads. Keep only one of them.
      mate = random.choice((0, 1))
      these_error_types[mate] = None
      pair_stats[0]['overlap']['duplicates'] += 1
      pair_stats[1]['overlap']['duplicates'] += 1
      dedup_log and dedup_log.write('omitting error {} {} from mate {}\n'
                                    .format(ref_coord, alt, mate+1))
    dedup_log and dedup_log.write('{:5d} {:1s}:  '.format(ref_coord, alt))
    for mate in (0, 1):
      if dedup_log:
        error_type = these_error_types[mate]
        if error_type is None:
          dedup_log.write('             ')
        else:
          dedup_log.write('  {:2d} errors  '.format(len(error_type)))
      new_errors_lists[mate].append(these_error_types[mate])
    dedup_log and dedup_log.write('\n')
  return new_errors_lists


def log_nonref_errors(nonref_errors, dedup_log):
  dedup_log.write('NO REF COORD:\n')
  for mate in (0, 1):
    for error_type in nonref_errors[mate]:
      error = error_type[0]
      dedup_log.write('{:5d} {:1s}:  '.format(error['coord'], error['alt']))
      if mate == 1:
        dedup_log.write('             ')
      dedup_log.write('  {:2d} errors\n'.format(len(error_type)))


def get_edges(pair):
  edges = [{}, {}]
  for mate in (0, 1):
    edges[mate]['start'] = pair[mate].get_position()
    edges[mate]['end'] = pair[mate].get_end_position(one_based=False)
  overlap = {}
  overlap['start'] = max(edges[0]['start'], edges[1]['start'])
  overlap['end'] = min(edges[0]['end'], edges[1]['end'])
  for mate in (0, 1):
    for edge in ('start', 'end'):
      if edges[mate]['start'] < overlap[edge] < edges[mate]['end']:
        edges[mate]['overlap_'+edge] = overlap[edge]
      else:
        edges[mate]['overlap_'+edge] = None
    if edges[mate]['overlap_start'] is not None and edges[mate]['overlap_end'] is None:
      edges[mate]['overlap_end'] = overlap['end']
    elif edges[mate]['overlap_start'] is None and edges[mate]['overlap_end'] is not None:
      edges[mate]['overlap_start'] = overlap['start']
  return edges


def get_overlap_len(edges):
  """Get the total number of reference base pairs that the two alignments span, as well as the
  number of base pairs in the overlap between the two alignments."""
  overlap_start = edges[0]['overlap_start'] or edges[1]['overlap_start']
  overlap_end = edges[0]['overlap_end'] or edges[1]['overlap_end']
  if overlap_start is None or overlap_end is None:
    overlap_len = 0
  else:
    overlap_len = overlap_end - overlap_start + 1
  non_overlap_lens = [None, None]
  for mate in (0, 1):
    non_overlap_lens[mate] = edges[mate]['end'] - edges[mate]['start'] + 1 - overlap_len
  return overlap_len, non_overlap_lens


def get_revcomp(seq_or_alt):
  if hasattr(seq_or_alt, 'translate'):
    return seq_or_alt.translate(REVCOMP_TABLE)[::-1]
  else:
    # If it's not a string, it might be an integer (deletion length). Return this unchanged.
    return seq_or_alt


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
  try:
    sys.exit(main(sys.argv))
  except IOError as ioe:
    if ioe.errno != errno.EPIPE:
      raise
