#!/usr/bin/env python

# Author: John Hawkins (jsh)

import logging
import optparse
import os
import subprocess
import sys

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')


def main():
  logging.info('Parsing command line.')
  usage = '%prog [options] input_file [...]'
  parser = optparse.OptionParser(usage=usage)
  parser.add_option(
      '--skip_trna',
      action='store_true', dest='skip_trna',
      default=False,
      help='If set, don\'t bother filtering trna.')
  parser.add_option(
      '--bowtie_path', type='string',
      action='store', dest='bowtie_path',
      default='bowtie',
      help='Path to bowtie binary.')
  parser.add_option(
      '--bowtie_parallelism', type='int',
      action='store', dest='bowtie_parallelism',
      default=7,
      help='Tell bowtie to run with this many parallel processes.')
  parser.add_option(
      '--bowtie_error_tolerance', type='int',
      action='store', dest='bowtie_error_tolerance',
      default=3,
      help='Allow at most this many errors in the first 28 nucleotides during bowtie alignments.')
  parser.add_option(
      '--bowtie_phred_quality', type='int',
      action='store', dest='bowtie_phred_quality',
      default=120,
      help='The sum of the Phred quality values at all mismatched positions (not just in the seed) may not exceed this number. Default = 120')
  parser.add_option(
      '--bowtie_max_matches', type='int',
      action='store', dest='bowtie_max_matches',
      default=5,
      help='If there are more matches than this give up and write to .max.')
  parser.add_option(
      '--trna_index', type='string',
      action='store', dest='trna_index',
      default='models/sc_tRNA_snoRNA',
      help='Path to bowtie index of tRNA to be removed.')
  parser.add_option(
      '--rrna_index', type='string',
      action='store', dest='rrna_index',
      default='models/rrna',
      help='Path to bowtie index of rRNA to be removed.')
  (options, args) = parser.parse_args()

  logging.info('Processing input.')
  for input_file in args:
    logging.info('Processing file: {0}'.format(input_file))
    if options.skip_trna:
      trna_free_file = input_file
    else:
      trna_free_file = filter_trna(input_file,
                                   options.trna_index,
                                   options.bowtie_path,
                                   options.bowtie_parallelism,
                                   options.bowtie_error_tolerance,
                                   options.bowtie_phred_quality,
                                   options.bowtie_max_matches)
    rrna_free_file = filter_rrna(trna_free_file,
                                 options.rrna_index,
                                 options.bowtie_path,
                                 options.bowtie_parallelism,
                                 options.bowtie_error_tolerance,
                                 options.bowtie_phred_quality,
                                 options.bowtie_max_matches)
    logging.info('Cleaned sequence file: {0}'.format(rrna_free_file))


def filter_trna(input_file,
                trna_index,
                bowtie_path,
                bowtie_parallelism,
                bowtie_error_tolerance,
                bowtie_phred_quality,
                bowtie_max_matches):
  """Remove tRNA from sequence file.
  Args:
    input_file: Path to fastq file to be processed.
    trna_index: The bowtie index of the RNA to strip out.
    bowtie_path: Path to the bowtie binary.
    bowtie_parallelism: Number of processes bowtie should use.
    bowtie_error_tolerance: Number of alignment errors tolerated per sequence.
    bowtie_max_matches: Maximum allowable matches before fail-fast.
  Returns:
    output_file: Name of the file to which tRNA-free sequences are written.
  """
  logging.info('Removing tRNA from: {0}'.format(input_file))
  matched_trna_file = open(os.path.splitext(input_file)[0] + '.trna', 'w')
  bowtie_error_log = open(
      os.path.splitext(input_file)[0] + '.bowtie_err_trna', 'w')
  output_file = os.path.splitext(input_file)[0] + '.trna_free'
  max_excluded_file = os.path.splitext(input_file)[0] + '.trna_toomany'
  command = [bowtie_path]
  command.extend(['-a'])
  command.extend(['-p', str(bowtie_parallelism)])
  command.extend(['-n', str(bowtie_error_tolerance)])
  command.extend(['-m', str(bowtie_max_matches)])
  command.extend(['-e', str(bowtie_phred_quality)])
  command.extend(['--nofw'])
  command.extend(['--un', output_file])
  command.extend(['--max', max_excluded_file])
  # These two are actual arguments, and need to be last.
  command.extend([trna_index])
  command.extend([input_file])
  logging.info(' '.join(command))
  subprocess.check_call(command,
                        stdout=matched_trna_file,
                        stderr=bowtie_error_log)
  return output_file


def filter_rrna(input_file,
                rrna_index,
                bowtie_path,
                bowtie_parallelism,
                bowtie_error_tolerance,
                bowtie_phred_quality,
                bowtie_max_matches):
  """Remove rRNA from sequence file.
  Args:
    input_file: Path to fastq file to be processed.
    rrna_index: The bowtie index of the RNA to strip out.
    bowtie_path: Path to the bowtie binary.
    bowtie_parallelism: Number of processes bowtie should use.
    bowtie_error_tolerance: Number of alignment errors tolerated per sequence.
    bowtie_max_matches: Maximum allowable matches before fail-fast.
  Returns:
    output_file: Name of the file to which rRNA-free sequences are written.
  """
  logging.info('Removing rRNA from: {0}'.format(input_file))
  bowtie_error_log = open(
      os.path.splitext(input_file)[0] + '.bowtie_err_rrna', 'w')
  matched_rrna_file = open(os.path.splitext(input_file)[0] + '.rrna', 'w')
  output_file = os.path.splitext(input_file)[0] + '.trna_rrna_free'
  max_excluded_file = os.path.splitext(input_file)[0] + '.rrna_toomany'
  command = [bowtie_path]
  command.extend(['-a'])
  command.extend(['-p', str(bowtie_parallelism)])
  command.extend(['-n', str(bowtie_error_tolerance)])
  command.extend(['-m', str(bowtie_max_matches)])
  command.extend(['-e', str(bowtie_phred_quality)])
  command.extend(['--nofw'])
  command.extend(['--un', output_file])
  command.extend(['--max', max_excluded_file])
  # These two are actual arguments, and need to be last.
  command.extend([rrna_index])
  command.extend([input_file])
  logging.info(' '.join(command))
  subprocess.check_call(command,
                        stdout=matched_rrna_file,
                        stderr=bowtie_error_log)
  return output_file


##############################################
if __name__ == "__main__":
    sys.exit(main())
