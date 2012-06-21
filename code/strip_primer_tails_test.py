#!/usr/bin/env python

# Author: John Hawkins (jsh)

import strip_primer_tails as strip

import unittest

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def same_sequence(a, b):
  qa = a.letter_annotations['phred_quality']
  qb = b.letter_annotations['phred_quality']
  return str(a) == str(b) and qa == qb

class TestCleaningAndStripping(unittest.TestCase):
  def setUp(self):
    self.inequality_error = (
        '\n{0:fastq-illumina} and \n{1:fastq-illumina} were not equal.')
    self.min_primer_match = 10
    self.min_sequence_len = 18
    self.primer = 'TCGTATGCCGTCTTCTGCTTG'
    self.primer_seq = SeqRecord(self.primer)
    primer_annotations = [10] * len(self.primer_seq)
    self.primer_seq.letter_annotations['phred_quality'] = primer_annotations
    self.baseline = SeqRecord('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT')
    allswell = [33] * 41
    self.baseline.letter_annotations['phred_quality'] = allswell

  def testTrivialCleaning(self):
    cleaned = strip.clean_for_illumina_flag(self.baseline)
    self.assertTrue(same_sequence(cleaned, self.baseline),
        self.inequality_error.format(cleaned, self.baseline))

  def testHappyCleaning(self):
    bad_tail = [2] * 5
    self.baseline.letter_annotations['phred_quality'][-5:] = bad_tail
    cleaned = strip.clean_for_illumina_flag(self.baseline)
    self.assertTrue(same_sequence(cleaned, self.baseline[:-5]),
        self.inequality_error.format(cleaned, self.baseline))

  def testRFindIfNot(self):
    tagger = lambda x: x == 2
    things = [1,1,2,1,2,2]
    self.assertEqual(3, strip.rfind_if_not(tagger, things))
    things = [1,1,1]
    self.assertEqual(2, strip.rfind_if_not(tagger, things))

  def testInternalBadReadSkippedByCleaning(self):
    bad_tail = [2] * 5
    self.baseline.letter_annotations['phred_quality'][-5:] = bad_tail
    self.baseline.letter_annotations['phred_quality'][5:10] = bad_tail
    self.assertEqual(36, len(strip.clean_for_illumina_flag(self.baseline)))

  def testTrivialTrimming(self):
    trimmed = strip.trim_primer(self.primer, self.baseline,
                                self.min_primer_match)
    self.assertTrue(same_sequence(trimmed, self.baseline),
        self.inequality_error.format(trimmed, self.baseline))

  def testSimpleTrimming(self):
    self.baseline = self.baseline[:-len(self.primer)] + self.primer_seq
    trimmed = strip.trim_primer(self.primer, self.baseline,
                                self.min_primer_match)
    self.assertTrue(same_sequence(trimmed, self.baseline[:-len(self.primer)]),
        self.inequality_error.format(trimmed, self.baseline))

  def testRNASeqTrimming(self):
    read = 'CATTCTGTGGAACGGTCCGGTTGGCGCTGTAGGCACCATCAATTCGTATG'
    self.primer_seq = SeqRecord(read)
    self.primer_seq.letter_annotations['phred_quality'] = [10] * len(read)
    primer = 'CTGTAGGCACCATCAAT'
    goal = 'CATTCTGTGGAACGGTCCGGTTGGCG'
    goal_seq = SeqRecord(goal)
    goal_seq.letter_annotations['phred_quality'] = [10] * len(goal)
    trimmed = strip.trim_primer(primer, self.primer_seq,
                                self.min_primer_match)
    self.assertTrue(same_sequence(trimmed, goal_seq),
        self.inequality_error.format(trimmed, self.baseline))

  def testIgnoredShortPrimerTrimming(self):
    self.baseline = self.baseline[:-5] + self.primer_seq[:5]
    trimmed = strip.trim_primer(self.primer, self.baseline,
                                self.min_primer_match)
    self.assertTrue(same_sequence(trimmed, self.baseline),
        self.inequality_error.format(trimmed, self.baseline))

  def testTrivialProcessSequences(self):
    processed = strip.processed_sequences(self.primer, [self.baseline],
                                          self.min_sequence_len,
                                          self.min_primer_match).next()
    self.assertTrue(same_sequence(processed, self.baseline),
        self.inequality_error.format(processed, self.baseline))

  def testTrimAndClean(self):
    s = self.baseline[:28] + self.primer_seq[:12]
    s.letter_annotations['phred_quality'][-20:] = [2] * 20
    processed = strip.processed_sequences(self.primer, [s],
                                          self.min_sequence_len,
                                          self.min_primer_match).next()
    self.assertTrue(same_sequence(processed, s[:20]),
        self.inequality_error.format(processed, s[:20]))

  def testTooSmallAfterTrimming(self):
    l = len(self.primer_seq) + 15
    s = self.baseline[:15] + self.primer_seq + self.baseline[l:]
    s.letter_annotations['phred_quality'][-5:] = [2] * 5
    processed = strip.processed_sequences(self.primer, [s],
                                          self.min_sequence_len,
                                          self.min_primer_match)
    self.assertEqual([], list(processed))

  def testTooSmallAfterCleaning(self):
    bad_tail = [2] * (len(self.baseline) - 16)
    self.baseline.letter_annotations['phred_quality'][16:] = bad_tail
    processed = strip.processed_sequences(self.primer, [self.baseline],
                                          self.min_sequence_len,
                                          self.min_primer_match)
    self.assertEqual([], list(processed))


if __name__ == '__main__':
  unittest.main()
