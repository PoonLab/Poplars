import unittest
from poplars.sequence_locator import *

HIV_NSEQ_PATH = os.path.abspath('fixtures/hiv-test-genome.fasta')
HIV_PSEQ_PATH = os.path.abspath('fixtures/hiv-test-proteins.fasta')
HIV_NCOORDS_PATH = os.path.abspath('fixtures/hiv_test_nt_coords.csv')
HIV_PCOORDS_PATH = os.path.abspath('fixtures/hiv_test_prot_coords.csv')

SIV_NSEQ_PATH = os.path.abspath('fixtures/siv-test-genome.fasta')
SIV_PSEQ_PATH = os.path.abspath('fixtures/siv-test-proteins.fasta')
SIV_NCOORDS_PATH = os.path.abspath('fixtures/siv_test_nt_coords.csv')
SIV_PCOORDS_PATH = os.path.abspath('fixtures/siv_test_prot_coords.csv')


class HIV(unittest.TestCase):

    def setUp(self):
        self.hiv_configs = handle_args('hiv', 'nucl')

        self.ref_nt_seq, self.ref_aa_seq = self.hiv_configs[0][0][1], self.hiv_configs[1]
        nt_coords, aa_coords = self.hiv_configs[2], self.hiv_configs[3]
        self.ref_seq = self.hiv_configs[4]

        self.nt_coords_handle = open(nt_coords)
        self.aa_coords_handle = open(aa_coords)

        self.genome = Genome(self.nt_coords_handle, self.ref_nt_seq,
                             self.aa_coords_handle, self.ref_aa_seq, self.ref_seq)

    def tearDown(self):
        self.nt_coords_handle.close()
        self.aa_coords_handle.close()


####################################
# Test cases for Region
####################################
class TestGenome(HIV):

    def testMakeRefRegions(self):
        expected = {'Gag':
                    [[790, 2292], 'ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAA'
                                  'TATAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAGAAGGC'
                                  'TGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACC'
                                  'CTCTATTGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCAAAACAAAAGTAAG'
                                  'AAAAAAGCACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAAAATTACCCTATAGTGCAGAACATCCAGGGG'
                                  'CAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTG'
                                  'ATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACACCATGCTAAACACAGTGGGGGGACATCAAGCA'
                                  'GCCATGCAAATGTTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGAGTGCATCCAGTGCATGCAGGGCCTATTGCA'
                                  'CCAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACAAAT'
                                  'AATCCACCTATCCCAGTAGGAGAAATTTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACC'
                                  'AGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGGTTCTATAAAACTCTAAGAGCCGAGCAA'
                                  'GCTTCACAGGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCA'
                                  'TTGGGACCAGCGGCTACACTAGAAGAAATGATGACAGCATGTCAGGGAGTAGGAGGACCCGGCCATAAGGCAAGAGTTTTGGCT'
                                  'GAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATGATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGT'
                                  'TTCAATTGTGGCAAAGAAGGGCACACAGCCAGAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGA'
                                  'CACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTT'
                                  'CTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAG'
                                  'GAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA',
                     [1, 500], 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYC'
                               'VHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSA'
                               'LSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEI'
                               'YKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMT'
                               'ACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLG'
                               'KIWPSYKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ'],
                    'Matrix(p17/p15)':
                    [[790, 1185], 'ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAA'
                                  'TATAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAGAAGGC'
                                  'TGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACC'
                                  'CTCTATTGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCAAAACAAAAGTAAG'
                                  'AAAAAAGCACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAAAATTAC',
                     [1, 132], 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYC'
                               'VHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY'],
                    'Capsid(p24/p27)':
                    [[1186, 1878], 'CCTATAGTGCAGAACATCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGA'
                                   'AGAGAAGGCTTTCAGCCCAGAAGTGATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACACCATGC'
                                   'TAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATGTTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGAGTG'
                                   'CATCCAGTGCATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCT'
                                   'TCAGGAACAAATAGGATGGATGACAAATAATCCACCTATCCCAGTAGGAGAAATTTATAAAAGATGGATAATCCTGGGATTAA'
                                   'ATAAAATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGAC'
                                   'CGGTTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAGGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGC'
                                   'GAACCCAGATTGTAAGACTATTTTAAAAGCATTGGGACCAGCGGCTACACTAGAAGAAATGATGACAGCATGTCAGGGAGTAG'
                                   'GAGGACCCGGCCATAAGGCAAGAGTTTTG',
                     [133, 363], 'PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHP'
                                 'VHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYK'
                                 'TLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVL'],
                    'p2':
                    [[1879, 1920], 'GCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATG',
                     [364, 377], 'AEAMSQVTNSATIM'],
                    'Nucleocapsid(p7/p8)':
                    [[1921, 2085], 'ATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACACAGCCAGAAATTG'
                                   'CAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAGGCTAAT',
                     [378, 432], 'MQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQAN'],
                    'p1':
                    [[2086, 2133], 'TTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTT',
                     [433, 448], 'FLGKIWPSYKGRPGNF'],
                    'p6':
                    [[2134, 2292], 'CTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCA'
                                   'GGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA',
                     [449, 500], 'LQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ']}
        result = self.genome.make_ref_regions(self.nt_coords_handle, self.aa_coords_handle, self.ref_nt_seq)

        for i, reg in enumerate(result):
            self.assertEqual(list(expected.keys())[i], reg.region_name)
            self.assertEqual(expected[reg.region_name][0], reg.ncoords)
            self.assertEqual(expected[reg.region_name][1], reg.nt_seq)
            self.assertEqual(expected[reg.region_name][2], reg.pcoords)
            self.assertEqual(expected[reg.region_name][3], reg.aa_seq)


####################################
# Test cases for helper functions
####################################
class TestValidSequence(unittest.TestCase):

    def testNucEmptyString(self):
        expected = False
        result = valid_sequence("nucl", [[">", ""]])
        self.assertEqual(expected, result)

    def testAAEmptyString(self):
        expected = False
        result = valid_sequence("prot", [[">", ""]])
        self.assertEqual(expected, result)

    def testNotDNA(self):
        expected = False
        result = valid_sequence("nucl", [[">header", "atgcatgc"], [">header2", "atgcatgcatgcnm"]])
        self.assertEqual(expected, result)

    def testDNAUpperCase(self):
        expected = True
        result = valid_sequence("nucl", [[">header", "ATGCGATCGATCGAGCTAGC"]])
        self.assertEqual(expected, result)

    def testNotAA(self):
        expected = False
        result = valid_sequence("prot", [[">h1", "ATGCATGCATGCNMZ"], [">h2", "GGGALLDMMP"]])
        self.assertEqual(expected, result)

    def testAAUpperCase(self):
        expected = True
        result = valid_sequence("prot", [[">header", "MWAGLSALALPWY"]])
        self.assertEqual(expected, result)

    def testGappySequence(self):
        expected = True
        result = valid_sequence("nucl", [["header", "ATGC-NATGCATGC"]])
        self.assertEqual(expected, result)


class TestValidInputs(unittest.TestCase):

    def testInvalidStartHIV(self):
        expected = False
        result = valid_inputs("hiv", -1, 100, "Complete")
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
        expected = True
        result = valid_inputs("hiv", 1, 10278, "p6")
        self.assertEqual(expected, result)

    def testOutofBoundsSIV(self):
        expected = True
        result = valid_inputs("siv", 1, 10536, "gp120")
        self.assertEqual(expected, result)

    def testInvalidRegionHIV(self):
        expected = False
        result = valid_inputs("hiv", 1, 100, "Vpx")
        self.assertEqual(expected, result)

    def testInvalidStringHIV(self):
        expected = False
        result = valid_inputs("hiv", 78, "kjl", "Complete")
        self.assertEqual(expected, result)

    def testValidStringSIV(self):
        expected = True
        result = valid_inputs("siv", 890, 6000, "end")
        self.assertEqual(expected, result)


class TestGetQuery(unittest.TestCase):

    def testNucleotideQuery(self):
        expected = [["query", "ATGCGCG"]]
        result = get_query("nucl", "ATGCGCG")
        self.assertEqual(expected, result)

    def testProteinQuery(self):
        expected = [["query1", "MPPLMMADLADLGG"]]
        query = ">query1\nMPPLMMADLADLGG"
        result = get_query("prot", query)
        self.assertEqual(expected, result)

    def testLongNucleotideSequence(self):
        expected = [["query", "ATGCGCGAATTAGCGA"]]
        query = "atgcgcg\naattagcga"
        result = get_query("nucl", query)
        self.assertEqual(expected, result)

    def testInvalidNucleotideQuery(self):
        query = ">query\natgcgcg&\n"
        with self.assertRaises(SystemExit) as e:
            get_query("nucl", query)
        self.assertEqual(e.exception.code, 0)

    def testInvalidProteinQuery(self):
        query = ">query\nMPPLMMAD>LADLGG"
        with self.assertRaises(SystemExit) as e:
            get_query("prot", query)
        self.assertEqual(e.exception.code, 0)

    def testRevCompProt(self):
        handle = ">seq2\nMPPLMMADLADLGG\n"
        expected = [["seq2", "MPPLMMADLADLGG"]]
        result = get_query('prot', handle)
        self.assertEqual(expected, result)


class TestGetReferenceSequence(unittest.TestCase):

    def testHIVRefNuclSeq(self):
        expected = [['K03455|HIVHXB2CG',
                     'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACC'
                     'AGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAAC'
                     'ACCAGCTTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACA'
                     'TGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG']]
        res = get_ref_seq(HIV_NSEQ_PATH, 'nucl')
        # Check the first 350 nucleotides in the default reference sequence
        result = [[res[0][0], res[0][1][:350]]]
        self.assertEqual(expected, result)

    def testSIVRefProtSeq(self):
        expected = [['Rev(exon1)|SIVMM239', 'MSNHEREEELRKRLRLIHLLHQT']]
        res = get_ref_seq(SIV_PSEQ_PATH, 'prot')
        result = [[res[0][0], res[0][1]]]
        self.assertEqual(expected, result)


class TestReverseComp(unittest.TestCase):
    def testSimpleUse(self):
        expected = "TCGCTAATTCGCGCATN*"
        result = reverse_comp("*NATGCGCGAATTAGCGA")
        self.assertEqual(expected, result)


class TestHandleNuclArgs(unittest.TestCase):

    def testHIVDefaultArgs(self):
        """
        Tests the scenario when the user selects HIV and nucleotide alignment
        """
        result = handle_args('hiv', 'nucl')

        expected_ref_nt_seq = [
            ['K03455|HIVHXB2CG',
             'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAG'
             'GGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGTTACACC'
             'CTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGT'
             'ACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG']]

        # Check the first 350 nucleotides in the default reference sequence
        seq = result[0][0][1][:350]
        header = result[0][0][0]
        result_ref_nt_seq = [[header, seq]]
        self.assertEqual(expected_ref_nt_seq, result_ref_nt_seq)

        expected_ref_aa_seq = [
            ['Gag|HIVHXB2CG',
             'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIE'
             'EEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAA'
             'EWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKN'
             'WMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGH']]

        # Check the first 400 amino acids in the reference protein sequence
        aa_seq = result[1][0][1][:400]
        header = result[1][0][0]
        result_ref_aa_seq = [[header, aa_seq]]
        self.assertEqual(expected_ref_aa_seq, result_ref_aa_seq)

    def testHIVArgs(self):
        """
        Tests the scenario when the user selects HIV, specifies the reference protein and nucleotide sequences
        """
        result = handle_args('hiv', 'nucl')

        expected_ref_nt = [
            ['K03455|HIVHXB2CG',
             'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAG'
             'GGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGTTACACC'
             'CTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGT'
             'ACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG']]

        # Check the first 350 nucleotides
        seq = result[0][0][1][:350]
        header = result[0][0][0]
        result_ref_nt = [[header, seq]]
        self.assertEqual(expected_ref_nt, result_ref_nt)

        expected_ref_prot = [
            ['Gag|HIVHXB2CG',             'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELR'
                                          'SLYNTVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTL'
                                          'NAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQM'
                                          'REPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTL'
                                          'RAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQR'
                                          'GNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPT'
                                          'APPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ'],
            ['Matrix(p17/p15)|HIVHXB2CG', 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELR'
                                          'SLYNTVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY'],
            ['Capsid(p24/p27)|HIVHXB2CG', 'PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEE'
                                          'AAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILD'
                                          'IRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKA'
                                          'RVL']]
        result_ref_prot = result[1][:3]
        self.assertEqual(expected_ref_prot, result_ref_prot)

    def testSIVDefaultArgs(self):
        """
        Tests the scenario when the user selects SIV and protein alignment
        """
        result = handle_args('siv', 'prot')

        expected_ref_nt_seq = [
            ['M33262|SIVMM239', 'GCATGCACATTTTAAAGGCTTTTGCTAAATATAGCCAAAAGTCCTTCTACAAATTTTCTAAGAGTTCTGATTCAAAGCAGTAACAG'
                                'GCCTTGTCTCATCATGAACTTTGGCATTTCATCTACAGCTAAGTTTATATCATAAATAGTTCTTTACAGGCAGCACCAACTTATAC'
                                'CCTTATAGCATACTTTACTGTGTGAAAATTGCATCTTTCATTAAGCTTACTGTAAATTTACTGGCTGTCTTCCTTGCAGGTTTCTG'
                                'GAAGGGATTTATTACAGTGCAAGAAGACATAGAATCTTAGACATATACTTAGAAAAGGAAGAAGGCATCATACCAGATTGGCAGGA'
                                'TTACACCTCAGGACCAGGAATTAGATACCCAAAGACATTTGGCTGGCTATGGAAAT']]

        # Check the first 400 nucleotides in the default reference sequence
        seq = result[0][0][1][:400]
        header = result[0][0][0]
        result_ref_nt_seq = [[header, seq]]
        self.assertEqual(expected_ref_nt_seq, result_ref_nt_seq)

        expected_ref_aa_seq = [
            ['Capsid(p24/p27)|SIVMM239',     'PVQQIGGNYVHLPLSPRTLNAWVKLIEEKKFGAEVVPGFQALSEGCTPYDINQMLNCVGDHQAAMQIIRDIIN'
                                             'EEAADWDLQHPQPAPQQGQLREPSGSDIAGTTSSVDEQIQWMYRQQNPIPVGNIYRRWIQLGLQKCVRMYNPT'
                                             'NILDVKQGPKEPFQSYVDRFYKSLRAEQTDAAVKNWMTQTLLIQNANPDCKLVLKGLGVNPTLEEMLTACQGV'
                                             'GGPGQKARLM'],
            ['p2|SIVMM239',                  'AEALKEALAPVPIPFAA'],
            ['Nucleocapsid(p7/p8)|SIVMM239', 'AQQRGPRKPIKCWNCGKEGHSARQCRAPRRQGCWKCGKMDHVMAKCPDRQAG']]
        result_ref_aa_seq = result[1][2:5]
        self.assertEqual(expected_ref_aa_seq, result_ref_aa_seq)

    def testSIVProtArgs(self):
        """
        Tests the scenario when the user specifies only the protein sequences
        """
        result = handle_args('siv', 'prot')

        expected_ref_nt = [
            ['M33262|SIVMM239', 'GCATGCACATTTTAAAGGCTTTTGCTAAATATAGCCAAAAGTCCTTCTACAAATTTTCTAAGAGTTCTGATTCAAAGCAGTAACAG'
                                'GCCTTGTCTCATCATGAACTTTGGCATTTCATCTACAGCTAAGTTTATATCATAAATAGTTCTTTACAGGCAGCACCAACTTATAC'
                                'CCTTATAGCATACTTTACTGTGTGAAAATTGCATCTTTCATTAAGCTTACTGTAAATTTACTGGCTGTCTTCCTTGCAGGTTTCTG'
                                'GAAGGGATTTATTACAGTGCAAGAAGACATAGAATCTTAGACATATACTTAGAAAAGGAAGAAGGCATCATACCAGATTGGCAGGA'
                                'TTACAC']]
        seq = result[0][0][1][:350]
        header = result[0][0][0]
        result_ref_nt = [[header, seq]]
        self.assertEqual(expected_ref_nt, result_ref_nt)

        expected_ref_prot = [['Rev(exon1)|SIVMM239', 'MSNHEREEELRKRLRLIHLLHQT'],
                             ['Rev(exon2)|SIVMM239', 'PYPTGPGTANQRRQRKRRWRRRWQQLLALADRIYSFPDPPTDTPLDLAIQQLQNLAIESIPDPP'
                                                     'TNTPEALCDPTEDSRSPQD']]
        result_ref_prot = result[1]
        self.assertEqual(expected_ref_prot, result_ref_prot)


if __name__ == '__main__':
    unittest.main()
