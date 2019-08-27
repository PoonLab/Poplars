import unittest
import os
from poplars.common import convert_fasta
from poplars.hypermut import *


class HIVTestCase(unittest.TestCase):

    def setUp(self):

        genome_path = os.path.join(os.path.dirname(__file__), 'fixtures/hiv-test-genome.fasta')
        infile_path = os.path.join(os.path.dirname(__file__), 'fixtures/aligned_HIV1_Mgroup.fasta')
        test_data = os.path.join(os.path.dirname(__file__), 'fixtures/hypermut-test.fasta')

        with open(genome_path) as handle, open(infile_path) as in_handle, open(test_data) as test_handle:
            self.hiv_genome = handle.read()
            self.hiv_data = in_handle.read()
            self.simple_data = convert_fasta(test_handle.read().split())

        self.mut = re.compile('[AGCT](?=[AG][AGT])')  # Matches potential mutation sites (RD)
        self.ctrl = re.compile('[AGCT](?=[CT].|[AG]C)')  # Matches potential control sites (YN or RC)


class Hypermut(HIVTestCase):

    maxDiff = None

    def testSimpleData(self):

        # Test potential mutation sites
        test_matches = ['GAA', 'GAG', 'GAT', 'GGA', 'GGG', 'GGT', 'CGAA']
        for s in test_matches:
            self.assertTrue(self.mut.match(s))
        test_mismatches = ['CTC', 'GTC', 'GGCGC']
        for s in test_mismatches:
            self.assertFalse(self.mut.match(s))

        # Test potential control sites
        test_matches = ['TAC', 'CCC', 'CAC', 'TTGC']
        for s in test_matches:
            self.assertTrue(self.ctrl.match(s))
        test_mismatches = ['TAT', 'CGGC', 'CAA']
        for s in test_mismatches:
            self.assertFalse(self.ctrl.match(s))

        # Test gees
        refseq = self.simple_data[0][1]  # Reference sequence is the first entry
        res_gees = [i for i, nt in enumerate(refseq) if nt.upper() == 'G']
        exp_gees = [5, 8, 15, 16, 19, 24, 27, 30, 33, 35, 36, 39, 46, 51, 62, 69, 84, 89, 92, 96, 102, 111,
                    115, 121, 139, 145, 156, 157, 162, 163, 164, 166, 168, 183, 184, 186, 195, 196, 198, 205,
                    210, 217, 226, 230, 234, 238, 239, 246, 257, 261, 274, 289, 297, 307, 315, 316, 318, 319,
                    320, 321, 327, 333, 338, 343, 352, 354, 355, 357, 358, 359, 360, 373, 389, 397, 403, 404,
                    408, 415, 421, 422, 426, 432, 435, 436, 437, 450, 453, 454, 472, 475, 497, 499, 500, 503,
                    504, 507, 510, 511, 516, 521, 525, 538, 540, 541, 550, 553, 570, 571, 572, 586, 588, 591,
                    592, 594, 595, 606, 607, 608, 612, 615, 616, 627, 628, 629, 633, 643]
        self.assertEqual(exp_gees, res_gees)

        query_seq = self.simple_data[1:]

        outputs = []
        for seq in query_seq:
            outputs.append(make_results(seq, res_gees))

        for i in range(len(outputs)):
            ex_num_muts = [0, 4, 26, 48]
            self.assertEqual(ex_num_muts[i], outputs[i].num_muts)     # Test match sites

            ex_pot_muts = [71, 69, 71, 71]
            self.assertEqual(ex_pot_muts[i], outputs[i].pot_muts)     # Test potential mutation sites

            ex_ctrl_muts = [0, 1, 1, 9]
            self.assertEqual(ex_ctrl_muts[i], outputs[i].ctrl_muts)    # Test control mutations

            ex_pot_ctrl_muts = [54, 52, 54, 54]
            self.assertEqual(ex_pot_ctrl_muts[i], outputs[i].potential_ctrls)  # Test potential controls

            ex_rate_ratio = ['undef', 3.01, 19.77, 4.06]
            if type(ex_rate_ratio) == str:
                self.assertEqual(ex_rate_ratio[i], outputs[i].rate_ratio)   # Test rate ratio
            else:
                self.assertAlmostEqual(ex_rate_ratio[i], outputs[i].rate_ratio, places=2)

            ex_p_value = [1, 0.282669, 5.35061e-07, 8.26961e-09]
            self.assertAlmostEqual(ex_p_value[i], outputs[i].p_value, places=2)    # Test p-value

            ex_ctable = [[[0, 71], [0, 54]],   # Seq2

                         [[4, 65], [1, 51]],   # Seq5

                         [[26, 45], [1, 53]],  # Seq7

                         [[48, 23], [9, 45]]]  # Seq14
            self.assertEqual(ex_ctable[i], outputs[i].ctable)   # Test contingency table

            # Test mutation sites
            ex_mut_sites = [{28: 0, 31: 0, 34: 0, 36: 0, 47: 0, 52: 0, 63: 0, 93: 0, 97: 0, 112: 0, 140: 0, 157: 0,
                            163: 0, 164: 0, 165: 0, 167: 0, 184: 0, 185: 0, 187: 0, 196: 0, 197: 0, 199: 0, 231: 0,
                            235: 0, 239: 0, 240: 0, 258: 0, 290: 0, 308: 0, 316: 0, 317: 0, 319: 0, 320: 0, 321: 0,
                            328: 0, 355: 0, 356: 0, 358: 0, 359: 0, 360: 0, 361: 0, 404: 0, 405: 0, 422: 0, 423: 0,
                            433: 0, 436: 0, 437: 0, 451: 0, 454: 0, 455: 0, 476: 0, 504: 0, 505: 0, 511: 0, 512: 0,
                            539: 0, 541: 0, 551: 0, 571: 0, 587: 0, 589: 0, 592: 0, 595: 0, 607: 0, 608: 0, 613: 0,
                            616: 0, 628: 0, 629: 0, 634: 0},

                            {28: 0, 31: 0, 34: 0, 36: 0, 47: 0, 52: 0, 63: 0, 93: 0, 97: 0, 112: 0, 140: 0, 157: 0,
                             163: 0, 164: 0, 165: 0, 167: 0, 184: 0, 185: 0, 187: 0, 196: 0, 197: 0, 199: 1, 231: 0,
                             235: 0, 239: 0, 240: 0, 258: 0, 290: 1, 316: 0, 317: 0, 319: 0, 320: 0, 321: 0, 328: 0,
                             355: 0, 356: 0, 358: 0, 359: 0, 360: 0, 361: 0, 422: 0, 423: 0, 427: 0, 433: 0, 436: 0,
                             437: 0, 451: 0, 454: 0, 455: 0, 476: 0, 504: 0, 505: 0, 511: 0, 512: 0, 539: 0, 541: 0,
                             551: 0, 571: 0, 587: 0, 589: 0, 592: 0, 595: 0, 607: 1, 608: 0, 613: 0, 616: 0, 628: 1,
                             629: 0, 634: 0},

                            {28: 0, 31: 0, 34: 0, 36: 1, 47: 0, 52: 1, 63: 0, 93: 1, 97: 0, 112: 1, 140: 0, 157: 1,
                             163: 0, 164: 1, 165: 0, 167: 0, 184: 0, 185: 0, 187: 0, 196: 1, 197: 0, 199: 0, 231: 0,
                             235: 0, 239: 1, 240: 1, 258: 0, 290: 0, 308: 0, 316: 0, 317: 0, 319: 1, 320: 1, 321: 0,
                             328: 0, 355: 1, 356: 0, 358: 0, 359: 1, 360: 1, 361: 0, 404: 1, 405: 1, 422: 1, 423: 0,
                             433: 0, 436: 1, 437: 0, 451: 0, 454: 0, 455: 0, 476: 0, 504: 1, 505: 0, 511: 1, 512: 0,
                             539: 0, 541: 0, 551: 0, 571: 1, 587: 0, 589: 0, 592: 0, 595: 1, 607: 1, 608: 1, 613: 0,
                             616: 1, 628: 1, 629: 0, 634: 0},

                            {28: 0, 31: 0, 34: 1, 36: 0, 47: 1, 52: 1, 63: 1, 93: 1, 97: 1, 112: 1, 140: 1, 157: 0,
                             163: 1, 164: 0, 165: 0, 167: 1, 184: 0, 185: 1, 187: 1, 196: 0, 197: 0, 199: 1, 231: 1,
                             235: 1, 239: 1, 240: 1, 258: 0, 290: 0, 308: 1, 316: 0, 317: 1, 319: 1, 320: 0, 321: 1,
                             328: 0, 355: 1, 356: 1, 358: 0, 359: 0, 360: 0, 361: 0, 404: 1, 405: 1, 422: 1, 423: 1,
                             433: 1, 436: 1, 437: 1, 451: 1, 454: 1, 455: 1, 476: 1, 504: 1, 505: 1, 511: 0, 512: 1,
                             539: 1, 541: 1, 551: 1, 571: 0, 587: 1, 589: 1, 592: 1, 595: 1, 607: 0, 608: 1, 613: 1,
                             616: 1, 628: 1, 629: 0, 634: 0}]

            self.assertEqual(ex_mut_sites[i], outputs[i].mut_sites)

            # Test control sites
            ex_ctrl_sites = [{6: 0, 9: 0, 16: 0, 17: 0, 20: 0, 25: 0, 37: 0, 40: 0, 70: 0, 85: 0, 90: 0, 103: 0,
                             116: 0, 122: 0, 146: 0, 158: 0, 169: 0, 206: 0, 211: 0, 218: 0, 227: 0, 247: 0, 262: 0,
                             275: 0, 298: 0, 322: 0, 334: 0, 339: 0, 344: 0, 353: 0, 374: 0, 390: 0, 398: 0, 409: 0,
                             416: 0, 427: 0, 438: 0, 473: 0, 498: 0, 500: 0, 501: 0, 508: 0, 517: 0, 522: 0, 526: 0,
                             542: 0, 554: 0, 572: 0, 573: 0, 593: 0, 596: 0, 609: 0, 617: 0, 630: 0},

                             {6: 0, 9: 0, 16: 0, 17: 0, 20: 0, 25: 0, 37: 0, 40: 0, 70: 0, 85: 0, 90: 0, 103: 0,
                              116: 0, 122: 0, 146: 0, 158: 0, 169: 0, 206: 0, 211: 0, 218: 0, 227: 0, 247: 0, 262: 0,
                              275: 0, 298: 0, 308: 1, 322: 0, 334: 0, 339: 0, 344: 0, 353: 0, 374: 0, 390: 0, 416: 0,
                              438: 0, 473: 0, 498: 0, 500: 0, 501: 0, 508: 0, 517: 0, 522: 0, 526: 0, 542: 0, 554: 0,
                              572: 0, 573: 0, 593: 0, 596: 0, 609: 0, 617: 0, 630: 0},

                             {6: 0, 9: 0, 16: 0, 17: 0, 20: 0, 25: 0, 37: 0, 40: 0, 70: 0, 85: 0, 90: 0, 103: 0,
                              116: 0, 122: 0, 146: 0, 158: 0, 169: 0, 206: 0, 211: 0, 218: 0, 227: 0, 247: 0, 262: 0,
                              275: 0, 298: 0, 322: 0, 334: 0, 339: 0, 344: 0, 353: 0, 374: 0, 390: 0, 398: 0, 409: 0,
                              416: 0, 427: 0, 438: 0, 473: 0, 498: 0, 500: 0, 501: 0, 508: 0, 517: 0, 522: 0, 526: 0,
                              542: 0, 554: 0, 572: 0, 573: 0, 593: 0, 596: 1, 609: 0, 617: 0, 630: 0},

                             {6: 0, 9: 0, 16: 0, 17: 0, 20: 0, 25: 0, 37: 0, 40: 0, 70: 0, 85: 0, 90: 0, 103: 0,
                              116: 0, 122: 1, 146: 0, 158: 1, 169: 0, 206: 1, 211: 0, 218: 0, 227: 0, 247: 0, 262: 0,
                              275: 1, 298: 0, 322: 1, 334: 0, 339: 0, 344: 0, 353: 0, 374: 0, 390: 0, 398: 0, 409: 0,
                              416: 0, 427: 0, 438: 0, 473: 0, 498: 0, 500: 0, 501: 0, 508: 0, 517: 0, 522: 0, 526: 1,
                              542: 1, 554: 0, 572: 0, 573: 0, 593: 0, 596: 0, 609: 1, 617: 0, 630: 1}]
            self.assertEqual(ex_ctrl_sites[i], outputs[i].ctrl_sites)


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


if __name__ == '__main__':
    unittest.main()
