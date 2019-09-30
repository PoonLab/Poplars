import unittest
from io import StringIO
from poplars.common import *


class TestConvertFasta(unittest.TestCase):
    def testSimpleConversion(self):
        # create file stream as test fixture
        handle = StringIO(">a\nACGT\n>b\nGCTA\n")
        result = convert_fasta(handle)
        
        expected = [['a', 'ACGT'], ['b', 'GCTA']]
        self.assertEqual(expected, result)

    def testFailedConversion(self):
        handle = StringIO("a\nACGT\n>b\nGCTA\n")
        
        with self.assertRaises(NameError):
            convert_fasta(handle)
            
    def testMultilineConversion(self):
        handle = StringIO(">a\nACGT\nGGGG\nATGATCGTAA\n>b\nGCTA\n\nAAAAAAAA\n")
        result = convert_fasta(handle)
        
        expected = [
            ['a', 'ACGTGGGGATGATCGTAA'],
            ['b', 'GCTAAAAAAAAA']
        ]
        self.assertEqual(expected, result)


class TestTransposeFasta(unittest.TestCase):

    def testWrongType(self):
        expected = None

        fasta = ('a', {'ATGC'})
        result = transpose_fasta(fasta)
        self.assertEqual(expected, result)

        fasta = ['>header', ['ATGC']]
        result = transpose_fasta(fasta)
        self.assertEqual(expected, result)

        fasta = [['>header', 'ATGC', 'TACG']]
        result = transpose_fasta(fasta)
        self.assertEqual(expected, result)

    def testSimpleTranspose(self):
        fasta = [['a', 'ACGT'], ['b', 'GCTA']]
        expected = [['A', 'G'], ['C', 'C'], ['G', 'T'], ['T', 'A']]
        result = transpose_fasta(fasta)
        self.assertEqual(expected, result)


class TestPluralityConsensus(unittest.TestCase):

    def testConsensusA(self):
        column = ['A', 'A', 'G', 'A', 'T']
        expected = 'A'
        result = plurality_consensus(column)
        self.assertEqual(expected, result)

    def testTiedNucleotides(self):
        column = ['A', 'T', 'G', 'C', 'A', 'T', 'G', 'C']
        expected = '-'
        result = plurality_consensus(column)
        self.assertEqual(expected, result)

        column = ['A', 'A', 'A', 'G', 'G', 'G']
        expected = 'R'
        result = plurality_consensus(column)
        self.assertEqual(expected, result)

        column = ['T', 'C', 'T', '-', '-']
        expected = 'T'
        result = plurality_consensus(column)
        self.assertEqual(expected, result)

        column = ['G', '-', '-', 'C']
        expected = 'S'
        result = plurality_consensus(column)
        self.assertEqual(expected, result)

        column = ['-', '-' '-']
        expected = '-'
        result = plurality_consensus(column)
        self.assertEqual(expected, result)


class TestConsensus(unittest.TestCase):

    def testShortSeq(self):
        fasta = [['a', 'ACGTGGGGATGATCGTAA'],
                 ['b', 'GCTAAAAAAAAATGCGCA'],
                 ['c', 'TGCATGCGAAGCTA--TC']]
        expected = 'DCBADGVGAAGATVSKHA'
        result = consensus(fasta)
        self.assertEqual(expected, result)

        fasta = [['>Seq1', 'TATTACAGGGCTATTATTAACAAGAGATGGTGGTAATAACAATGGGACCGATGGTAATAACAATGGGACCGAAATCTTCAGA'],
                 ['>Seq2', 'TATTACAGGGCTATTATTAACAAGAGATGGTGGTAATAACAATGGGACCGATGGTAATAACAATGGGACCGAAATCTTCAGA'],
                 ['>Seq5', 'TATTACAGGGCTATTATTAACAAGAGATGGTGGTAATAACAATAGGACCGATGGTAATAACAATAGGACCGAAACCTTCAGA'],
                 ['>Seq7', 'TATTACAAGGCTATTATTAACAAGAGATGGTAATAATAACAATAAGACCGATAGTAATAACAATAGGACCGAAATCTTCAGA'],
                 ['>Seq9', 'TATTACAGGGCTATTATTAACAAAAAATAGTAGTAATAACAATGAAACCAATAGTAATAACAATAGAACCGAAACCTTCAAA']]
        expected_consens = 'TATTACAGGGCTATTATTAACAAGAGATGGTGGTAATAACAATGGGACCGATGGTAATAACAATAGGACCGAAATCTTCAGA'
        result = consensus(fasta)
        self.assertEqual(expected_consens, result)


if __name__ == '__main__':
    unittest.main()
