import unittest
from poplars.sequence_locator import valid_sequence
from poplars.sequence_locator import get_ref_seq
from poplars.sequence_locator import valid_inputs


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

    def testDNALowerCase(self):
        expected = True
        result = valid_sequence("nucl", "atgcgcgatgcagcac")
        self.assertEqual(expected, result)

    def testDNAUpperCase(self):
        expected = True
        result = valid_sequence("nucl", "ATGCGATCGATCGAGCTAGC")
        self.assertEqual(expected, result)

    def testNotAA(self):
        expected = False
        result = valid_sequence("prot", "ATGCATGCATGCNMZ")
        self.assertEqual(expected, result)

    def testAALowerCase(self):
        expected = True
        result = valid_sequence("prot", "MWAGLSALALPWY")
        self.assertEqual(expected, result)

    def testAAUpperCase(self):
        expected = True
        result = valid_sequence("prot", "mwaglsalwgggp")
        self.assertEqual(expected, result)


class TestGetReferenceSequence(unittest.TestCase):

    def testDefaultHIVGenome(self):
        expected = "tggaagggctaattcactcccaacgaagacaagatatccttgatctgtgg" \
                   "atctaccacacacaaggctacttccctgattagcagaactacacaccagg" \
                   "gccagggatcagatatccactgacctttggatggtgctacaagctagtac" \
                   "cagttgagccagagaagttagaagaagccaacaaaggagagaacaccagc" \
                   "ttgttacaccctgtgagcctgcatggaatggatgacccggagagagaagt" \
                   "gttagagtggaggtttgacagccgcctagcatttcatcacatggcccgag" \
                   "agctgcatccggagtacttcaagaactgctgacatcgagcttgctacaag"
        result = get_ref_seq("hiv")
        self.assertEqual(expected, result)

    def testDefaultSIVGenome(self):
        expected = "gcatgcacattttaaaggcttttgctaaatatagccaaaagtccttctac" \
                   "aaattttctaagagttctgattcaaagcagtaacaggccttgtctcatca" \
                   "tgaactttggcatttcatctacagctaagtttatatcataaatagttctt" \
                   "tacaggcagcaccaacttatacccttatagcatactttactgtgtgaaaa" \
                   "ttgcatctttcattaagcttactgtaaatttactggctgtcttccttgca" \
                   "ggtttctggaagggatttattacagtgcaagaagacatagaatcttagac" \
                   "atatacttagaaaaggaagaaggcatcataccagattggcaggattacac" \
                   "ctcaggaccaggaattagatacccaaagacatttggctggctatggaaat"
        result = get_ref_seq("siv")
        self.assertEqual(expected, result)

    def testInputGenome(self):
        expected = "gcatgcacattttaaaggcttttgctaaatatagccaaaagtccttctac" \
                   "aaattttctaagagttctgattcaaagcagtaacaggccttgtctcatca" \
                   "tgaactttggcatttcatctacagctaagtttatatcataaatagttctt" \
                   "tacaggcagcaccaacttatacccttatagcatactttactgtgtgaaaa" \
                   "ttgcatctttcattaagcttactgtaaatttactggctgtcttccttgca" \
                   "ggtttctggaagggatttattacagtgcaagaagacatagaatcttagac" \
                   "atatacttagaaaaggaagaaggcatcataccagattggcaggattacac" \
                   "ctcaggaccaggaattagatacccaaagacatttggctggctatggaaat"
        result = get_ref_seq("siv", "/home/kwade4/PycharmProjects/Poplars/poplars/ref_genomes/siv-test-genome.fasta")
        self.assertEqual(expected, result)


class TestAlign(unittest.TestCase):

    def testSimpleUse(self):
        pass


class TestValidCoordinates(unittest.TestCase):

    def testInvalidStartHIV(self):
        expected = False
        result = valid_inputs("hiv", -0, 100, "Complete")
        self.assertEqual(expected, result)

    def testInvalidEndSIV(self):
        expected = False
        result = valid_inputs("siv", 0, -90, "Gag")
        self.assertEqual(expected, result)

    def testEndBeforeStartHIV(self):
        expected = False
        result = valid_inputs("hiv", 10, 1, "Pol")
        self.assertEqual(expected, result)

    def testValidHIV(self):
        expected = True
        result = valid_inputs("hiv", 1, 100, "Env")
        self.assertEqual(expected, result)

    def testValidSIV(self):
        expected = True
        result = valid_inputs("siv", 1, 10278, "Nef")
        self.assertEqual(expected, result)

    def testSameStartEndSIV(self):
        expected = False
        result = valid_inputs("siv", 100, 100, "gp160")
        self.assertEqual(expected, result)

    def testSmallRangeHIV(self):
        expected = True
        result = valid_inputs("hiv", 989, 990, "RT")
        self.assertEqual(expected, result)

    def testOutOfBoundsHIV(self):
        expected = False
        result = valid_inputs("hiv", 1, 10278, "p6")
        self.assertEqual(expected, result)

    def testOutofBoundsSIV(self):
        expected = False
        result = valid_inputs("siv", 1, 10279, "gp120")
        self.assertEqual(expected, result)

    def testInvalidRegionHIV(self):
        expected = False
        result = valid_inputs("hiv", 1, 100, "Vpx")
        self.assertEqual(expected, result)

    def testInvalidRegionSIV(self):
        expected = False
        result = valid_inputs("siv", 1, 10278, "Vpu")
        self.assertEqual(expected, result)
