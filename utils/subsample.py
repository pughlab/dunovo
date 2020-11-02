#!/usr/bin/env python3
import argparse
import logging
import random
import sys

USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  parser.add_argument('infile', metavar='read-families.tsv', nargs='?',
    help='The input reads, sorted into families.')
  parser.add_argument('-f', '--fraction', type=float, default=0.1,
    help='Fraction of families to output. Default: %(default)s')
  parser.add_argument('-s', '--seed', type=int, default=1,
    help='Random number generator seed. Default: %(default)s')
  options.add_argument('-h', '--help', action='help',
    help='Print this argument help text and exit.')
  logs = parser.add_argument_group('Logging')
  logs.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = logs.add_mutually_exclusive_group()
  volume.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  random.seed(args.seed)

  if args.infile:
    infile = open(args.infile)
  else:
    infile = sys.stdin

  family = []
  last_barcode = None
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    if not fields:
      continue
    barcode = fields[0]
    if barcode != last_barcode:
      if random.random() <= args.fraction:
        sys.stdout.write(''.join(family))
      family = []
    family.append(line)
    last_barcode = barcode
  if random.random() <= args.fraction:
    sys.stdout.write(''.join(family))

  if infile is not sys.stdin:
    infile.close()


def fail(message):
  logging.critical('Error: '+str(message))
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception(message)


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
