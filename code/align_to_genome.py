#!/usr/bin/env python

# Author: John Hawkins (jsh)

import logging
import optparse
import os
import subprocess
import shutil
import sys

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')


def main():
  logging.info('Parsing command line.')
  usage = '%prog [options] input_file [...] output_file_base'
  parser = optparse.OptionParser(usage=usage)
  parser.add_option(
      '--genome', type='string',
      action='store', dest='genome',
      default='models/sc_sgd_gff_20091011',
      help='Path to bowtie index of genome against which to align. Default is \'models/sc_sgd_gff_20091011\'')
  parser.add_option(
      '--treat_as_mrna',
      action='store_true', dest='treat_as_mrna',
      default=False,
      help='If set, use tophat to align around introns. Warning, untested code, so it is best to go with the default: False.')
  parser.add_option(
      '--bowtie_path', type='string',
      action='store', dest='bowtie_path',
      default='bowtie',
      help='Path to bowtie binary.Default is \'bowtie\' and probably will not work.')
  parser.add_option(
      '--fna_genome', type='string',
      action='store', dest='fna_genome',
      default='models/sc_sgd_gff_20091011.fna',
      help='Path to reference genome (not index) for samtools fillmd.')
  parser.add_option(
      '--tophat_path', type='string',
      action='store', dest='tophat_path',
      default='tophat',
      help='Path to tophat binary.')
  parser.add_option(
      '--tophat_min_anchor_len', type='int',
      action='store', dest='tophat_min_anchor_len',
      default=12,
      help='Tophat should require anchor coverage at least this wide.')
  parser.add_option(
      '--tophat_max_intron_len', type='int',
      action='store', dest='tophat_max_intron_len',
      default=5000,
      help='Ignore intron candidates larger than this. Default is for yeast.')
  parser.add_option(
      '--bowtie_parallelism', type='int',
      action='store', dest='bowtie_parallelism',
      default=7,
      help='Tell bowtie to run with this many parallel processes. Default is 7.')
  parser.add_option(
      '--bowtie_error_tolerance', type='int',
      action='store', dest='bowtie_error_tolerance',
      default=3,
      help='Allow at most this many errors in bowtie alignments.')
  parser.add_option(
      '--bowtie_phred_quality', type='int',
      action='store', dest='bowtie_phred_quality',
      default=120,
      help='The sum of the Phred quality values at all mismatched positions (not just in the seed) may not exceed this number. Default = 120')
  parser.add_option(
      '--bowtie_max_matches', type='int',
      action='store', dest='bowtie_max_matches',
      default=1,
      help='If there are more matches than this give up and write to .max.')
  (options, args) = parser.parse_args()

  logging.info('Processing input.')
  input_files = args[:-1]
  output_file_base = args[-1]
  if options.treat_as_mrna:
    align_to_genome = align_with_tophat
    binary_path = options.tophat_path
  else:
    align_to_genome = align_with_bowtie
    binary_path = options.bowtie_path
  final_alignment = align_to_genome(input_files,
                                    options.genome,
                                    binary_path,
                                    options.bowtie_parallelism,
                                    options.bowtie_error_tolerance,
                                    options.bowtie_phred_quality,
                                    options.bowtie_max_matches,
                                    options.tophat_min_anchor_len,
                                    options.tophat_max_intron_len,
                                    options.fna_genome,
                                    output_file_base)
  logging.info('Final alignment file: {0}'.format(final_alignment))


def align_with_tophat(input_files,
                      genome,
                      binary_path,
                      bowtie_parallelism,
                      bowtie_error_tolerance,
                      bowtie_phred_quality,
                      bowtie_max_matches,
                      tophat_min_anchor_len,
                      tophat_max_intron_len,
                      fna_genome,
                      output_file_base):
  """Align input sequences in input_file to provided genome using tophat.
  Args:
    input_files: List of fastq files to be processed.
    genome: The bowtie index of the genome to align against.
    binary_path: Path to the tophat binary.
    bowtie_parallelism: Number of processes bowtie should use.
    bowtie_error_tolerance: Number of alignment errors tolerated per sequence.
    bowtie_max_matches: Maximum allowable matches before fail-fast.
    tophat_min_anchor_len: Minimum anchor length for tophat junctions.
    fna_genome: Reference genome for computing MD tags.
    output_file_base: Base name for output files and tophat dump directory.
  Returns:
    output_file: Name of the file to which final alignments are written.
  """
  logging.info('Aligning with tophat.')
  output_file = output_file_base + '.alignment.sam'
  tophat_stdout_log = open(output_file_base + '.tophat_stdout_genome', 'w')
  tophat_error_log = open(output_file_base + '.tophat_stderr_genome', 'w')
  tophat_dir = output_file_base + '.tophat'
  # First run Tophat.
  command = [binary_path]
  command.extend(['--max-multihits', str(bowtie_max_matches)])
  command.extend(['--min-anchor-length', str(tophat_min_anchor_len)])
  command.extend(['--max-intron-length', str(tophat_max_intron_len)])
  command.extend(['--num-threads', str(bowtie_parallelism)])
  command.extend(['--output-dir', tophat_dir])
  command.extend(['--segment-mismatches', str(bowtie_error_tolerance)])
  if not os.path.exists(tophat_dir):
    os.mkdir(tophat_dir)
  # These are positional arguments, and should be last.
  command.append(genome)
  command.append(','.join(input_files))
  logging.info(' '.join(command))
  subprocess.check_call(command,
                        stdout=tophat_stdout_log,
                        stderr=tophat_error_log)
  bam_file = os.path.join(tophat_dir, 'accepted_hits.bam')
  if not os.path.exists(bam_file):
    logging.fatal('No such file: {0}'.format(bam_file))
  if not os.path.exists(fna_genome):
    logging.fatal('No such file: {0}'.format(fna_genome))
  # Run output through samtools fillmd.
  sam_file = open(os.path.splitext(bam_file)[0] + '.sam', 'w')
  command = ['samtools', 'fillmd']
  command.append(bam_file)
  command.append(fna_genome)
  logging.info(' '.join(command) + ' > ' + sam_file.name)
  subprocess.check_call(command,
                        stdout=sam_file)
  sam_file.close()
  # Sort the SAM output back to top level workspace.
  command = ['sort']
  command.extend(['-k', '3d,3'])
  command.extend(['-k', '2g,2'])
  command.extend(['-k', '4g,4'])
  command.append(sam_file.name)
  with open(output_file, 'w') as output:
    logging.info(' '.join(command) + ' > ' + output_file)
    subprocess.check_call(command,
                          stdout=output)
  return output_file


def align_with_bowtie(input_files,
                      genome,
                      binary_path,
                      bowtie_parallelism,
                      bowtie_error_tolerance,
                      bowtie_phred_quality,
                      bowtie_max_matches,
                      tophat_min_anchor_len,
                      tophat_max_intron_len,
                      fna_genome,
                      output_file_base):
  """Align input sequences in input_file to provided genome with just bowtie.
  Args:
    input_files: List of fastq files to be processed.
    genome: The bowtie index of the genome to align against.
    binary_path: Path to the bowtie binary.
    bowtie_parallelism: Number of processes bowtie should use.
    bowtie_error_tolerance: Number of alignment errors tolerated per sequence.
    bowtie_max_matches: Maximum allowable matches before fail-fast.
    tophat_min_anchor_len: Ignored.
    tophat_max_intron_len: Ignored.
    fna_genome: Ignored.
    output_file_base: Base name for output files.
  Returns:
    output_file: Name of the file to which final alignments are written.
  """
  logging.info('Aligning with bowtie.')
  output_file = output_file_base + '.alignment.sam'
  unsorted_output = output_file + '.unsorted'
  bowtie_stdout_log = open(output_file_base + '.bowtie_stdout_genome', 'w')
  bowtie_stderr_log = open(output_file_base + '.bowtie_stderr_genome', 'w')
  unmatchable_seqs = output_file_base + '.unmatchable'
  max_excluded_file = output_file_base + '.toomany'
  command = [binary_path]
  command.extend(['-a'])
  command.extend(['-p', str(bowtie_parallelism)])
  command.extend(['-n', str(bowtie_error_tolerance)])
  command.extend(['-e', str(bowtie_phred_quality)])
  command.extend(['-m', str(bowtie_max_matches)])
  command.extend(['--sam'])
  command.extend(['--sam-nohead'])
  command.extend(['--un', unmatchable_seqs])
  command.extend(['--max', max_excluded_file])
  # These three are actual arguments, and need to be last.
  command.extend([genome])
  command.extend([','.join(input_files)])
  command.extend([unsorted_output])
  logging.info(' '.join(command))
  subprocess.check_call(command,
                        stdout=bowtie_stdout_log,
                        stderr=bowtie_stderr_log)
  command = ['sort']
  command.extend(['-k', '3d,3'])
  command.extend(['-k', '2g,2'])
  command.extend(['-k', '4g,4'])
  command.append(unsorted_output)
  with open(output_file, 'w') as output:
    logging.info(' '.join(command) + ' > ' + output_file)
    subprocess.check_call(command,
                          stdout=output)
  return output_file


##############################################
if __name__ == "__main__":
    sys.exit(main())
