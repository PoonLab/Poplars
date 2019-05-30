import unittest
from poplars.sequence_locator import valid_sequence


class TestValidSequence(unittest.TestCase):

    def testNucEmptyString(self):
        expected = False
        result = valid_sequence("nucl", "")
        self.assertEqual(expected, result)

    def testAAEmptyString(self):
        expected = False
        result = valid_sequence("prot", "")
        self.assertEqual(expected, result)

    def testNotDNA(self):
        expected = False
        result = valid_sequence("nucl", "atgcatgcatgcnm")
        self.assertEqual(expected, result)

    def testSimpleUseDNA(self):
        expected = True
        result = valid_sequence("nucl", "atgcgcgatgcagcac")

    def testSimpleUseDNA1(self):
        expected = True
        result = valid_sequence("nucl", "ATGCGATCGATCGAGCTAGC")
        self.assertEqual(expected, result)

    def testNotAA(self):
        expected = False
        result = valid_sequence("prot", "ATGCATGCATGCNMZ")
        self.assertEqual(expected, result)

    def testSimpleUseAA(self):
        expected = True
        result = valid_sequence("prot", "MWAGLSALALPWY")
        self.assertEqual(expected, result)

    def testSimpleUseAA1(self):
        expected = True
        result = valid_sequence("prot", "mwaglsalwgggp")
        self.assertEqual(expected, result)
