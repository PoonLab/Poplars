import unittest
from io import StringIO
from poplars.common import convert_fasta


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
