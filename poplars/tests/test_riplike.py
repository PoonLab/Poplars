import unittest
import os
from poplars.riplike import *

TEST_HIV_GENOME = os.path.join(os.path.dirname(__file__), '../ref_genomes/K03455.fasta')
TEST_REFERENCE = os.path.join(os.path.dirname(__file__), '../ref_genomes/HIV1_Mgroup.fasta')


class testRiplike(unittest.TestCase):

    def setUp(self):
        with open(TEST_HIV_GENOME) as handle, open(TEST_REFERENCE) as ref_handle:
            self.hiv_genome = handle.readlines()[1]
            self.reference = convert_fasta(ref_handle)

    def testHamming(self):
        test = [['query', 'ACGT'], ['CON_OF_CONS', 'ACGT'], ['A', 'ACGG']]
        result = hamming(test)
        expected = {'A': [0, 0, 0, 1]}
        self.assertEqual(expected, result)

    def testUpdateAlignment(self):
        aln = update_alignment(self.hiv_genome, self.reference)
        aln = dict(aln)

        # alignment should have original number plus one
        result = len(aln)
        expected = 13 + 1
        self.assertEqual(expected, result)

        # alignment should contain entry for query sequence
        result = 'query' in aln
        self.assertTrue(result)

        # aligned query should have same number of nucleotides
        result = [aln['query'].count(nt) for nt in 'ACGT']
        expected = [self.hiv_genome.upper().count(nt) for nt in 'ACGT']
        self.assertEqual(expected, result)

        # all sequences should have the same length
        result = len(set([len(s) for h, s in aln.items()]))
        expected = 1
        self.assertEqual(expected, result)

    def testEncode(self):
        s = 'ATGCGC--  T'
        expected = [0B0001, 0B0010, 0B0100, 0B0011, 0B0100, 0B0011, 0B1111, 0B1111, 0B0000, 0B0000, 0B0010]
        result = encode(s)
        self.assertEqual(expected, result)

    def testRiplike(self):
        result = riplike(self.hiv_genome, self.reference)









