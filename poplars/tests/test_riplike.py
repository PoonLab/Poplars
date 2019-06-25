import unittest
import os
from poplars.riplike import *


TEST_HIV_GENOME = os.path.join(os.path.dirname(__file__), '../ref_genomes/K03455.fasta')

class testRiplike(unittest.TestCase):

    def setUp(self):
        with open(TEST_HIV_GENOME) as handle:
            self.hiv_genome = handle.readlines()[1]
        
    def testPdistSimple(self):
        result = pdistance('ACGT', 'ACGC')
        expected = (1, 4)
        self.assertEqual(expected, result)
        
        result = pdistance('ACGT', 'TGCA')
        expected = (4, 4)
        self.assertEqual(expected, result)
    
    def testPdistGapped(self):
        result = pdistance('ACGT', '---T')
        expected = (0, 1)
        self.assertEqual(expected, result)
        
        result = pdistance('ACGT', 'G---')
        expected = (1, 1)
        self.assertEqual(expected, result)
    
    def testBootstrap(self):
        pass  #FIXME: not sure how to test this function just yet
        
    def testUpdateAlignment(self):
        aln = update_alignment(self.hiv_genome)
        aln = dict(aln)
        
        # alignment should have original number plus one
        result = len(aln)
        expected = 11
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
        
    def testRiplike(self):
        result = riplike(self.hiv_genome)
        








