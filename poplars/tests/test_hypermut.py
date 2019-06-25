import unittest
import os
import re
from poplars.common import convert_fasta
from poplars.hypermut import make_results
from poplars.hypermut import hypermut
from poplars.hypermut import rate_ratio
from poplars.hypermut import make_data_file

TEST_HIV = os.path.join(os.path.dirname(__file__), 'fixtures/hiv-test-genome.fasta')
TEST_DATA = os.path.join(os.path.dirname(__file__), 'fixtures/aligned_HIV1_Mgroup.fasta')


class TestMakeResults(unittest.TestCase):
    def setUp(self):
        self.hiv_data = open(TEST_HIV)
        self.data = open(TEST_DATA)

    def testDSMotifs(self):
        mut = re.compile('[AGCT](?=[AG][AGT])')  # Matches potential mutation sites (GRD)
        s = convert_fasta(self.hiv_data)[0][1]
        result = [match.start() for match in mut.finditer(s)]
        expected = [0, 1, 2, 3, 4, 5, 9, 10, 20, 23, 24, 25, 26, 29, 30, 31, 32, 34, 40, 41, 45, 47, 48, 49, 62,
                    63, 64, 76, 77, 80, 83, 84, 85, 96, 97, 98, 102, 103, 104, 105, 106, 109, 110, 111, 113, 120,
                    127, 128, 129, 131, 132, 139, 140, 144, 145, 150, 151, 154, 155, 159, 160, 161, 162, 163, 164,
                    165, 168, 169, 170, 171, 172, 173, 174, 178, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190,
                    196, 201, 211, 213, 214, 221, 223, 224, 225, 226, 228, 229, 230, 232, 237, 238, 239, 240, 241,
                    242, 243, 244, 245, 246, 247, 249, 252, 253, 254, 255, 257, 258, 259, 260, 261, 265, 268, 276,
                    279, 284, 289, 291, 296, 297, 298, 299, 305, 309, 310, 311, 312, 319, 320, 321, 322, 329, 332,
                    335, 336, 346, 347]
        self.assertEqual(expected, result)

    def testCtrlMotifs(self):
        ctrl = re.compile('[AGCT](?=[CT].|[AG]C)')
        s = convert_fasta(self.hiv_data)[0][1]
        result = [match.start() for match in ctrl.finditer(s)]
        expected = [6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 27, 28, 33, 35, 36, 37, 38, 39, 42, 43,
                    44, 46, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 65, 66, 67, 68, 69, 70, 71, 72, 73,
                    74, 75, 78, 79, 81, 82, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 99, 100, 101, 107, 108, 112,
                    114, 115, 116, 117, 118, 119, 121, 122, 123, 124, 125, 126, 130, 133, 134, 135, 136, 137, 138,
                    141, 142, 143, 146, 147, 148, 149, 152, 153, 156, 157, 158, 166, 167, 175, 176, 177, 179, 180,
                    191, 192, 193, 194, 195, 197, 198, 199, 200, 202, 203, 204, 205, 206, 207, 208, 209, 210, 212,
                    215, 216, 217, 218, 219, 220, 222, 227, 231, 233, 234, 235, 236, 248, 250, 251, 256, 262, 263,
                    264, 266, 267, 269, 270, 271, 272, 273, 274, 275, 277, 278, 280, 281, 282, 283, 285, 286, 287,
                    288, 290, 292, 293, 294, 295, 300, 301, 302, 303, 304, 306, 307, 308, 313, 314, 315, 316, 317,
                    318, 323, 324, 325, 326, 327, 328, 330, 331, 333, 334, 337, 338, 339, 340, 341, 342, 343, 344, 345]
        self.assertEqual(expected, result)

    def test_alignmentA1(self):
        fasta = convert_fasta(self.data)
        refseq = fasta[0][1]
        gees = [i for i, nt in enumerate(refseq) if nt.upper() == 'G']
        seq = fasta[1]
        output = make_results(seq, gees)

        ex_num_muts = 47
        self.assertEqual(ex_num_muts, output.num_muts)

        ex_pot_muts = 1068
        self.assertEqual(ex_pot_muts, output.pot_muts)

        ex_ctrl_muts = 53
        self.assertEqual(ex_ctrl_muts, output.ctrl_muts)

        ex_pot_ctrl_muts = 986
        self.assertEqual(ex_pot_ctrl_muts, output.potential_ctrls)

        ex_rate_ratio = 0.82
        self.assertEqual(ex_rate_ratio, round(output.rate_ratio, 2))

        ex_p_value = 0.870295
        self.assertEqual(ex_p_value, output.p_value)

        ex_ctable = [[47, 1068], [53, 986]]
        self.assertEqual(ex_ctable, output.ctable)

        # ex_mut_sites = None
        # ex_ctrl_sites = None

    def tearDown(self):
        self.hiv_data.close()
        self.data.close()


class TestHypermut(unittest.TestCase):

    def setUp(self):
        self.hiv_data = open(TEST_HIV)

    def testGees(self):
        fasta = convert_fasta(self.hiv_data)
        refseq = fasta[0][1]        # Reference sequence is the first entry
        result = [i for i, nt in enumerate(refseq) if nt.upper() == 'G']
        expected = [1, 2, 5, 6, 7, 24, 27, 32, 41, 46, 48, 49, 65, 66, 77, 82, 85, 98, 99, 100, 104, 105, 106, 111,
                    121, 128, 129, 132, 133, 135, 142, 146, 152, 155, 157, 161, 163, 166, 170, 173, 176, 185, 186,
                    188, 190, 198, 202, 212, 214, 216, 220, 224, 225, 229, 230, 233, 238, 239, 241, 243, 245, 248,
                    250, 254, 256, 258, 259, 261, 262, 266, 270, 273, 278, 292, 293, 297, 299, 301, 304, 310, 311,
                    313, 322, 327, 330, 336, 338, 342, 349]
        self.assertEqual(expected, result)

    def tearDown(self):
        self.hiv_data.close()


class TestRateRatio(unittest.TestCase):
    def testRateRatio(self):
        ctable = [[4, 17], [7, 33]]
        expected = 1.1
        output = rate_ratio(ctable)
        result = round(output, 1)
        self.assertEqual(expected, result)

    def test1UndefinedRateRatio(self):
        ctable = [[4, 17], [0, 33]]
        expected = "undef"
        result = rate_ratio(ctable)
        self.assertEqual(expected, result)

    def test2UndefinedRateRatio(self):
        ctable = [[4, 17], [7, 0]]
        expected = "undef"
        result = rate_ratio(ctable)
        self.assertEqual(expected, result)


class testMakeDataFile(unittest.TestCase):
    def test(self):
        pass