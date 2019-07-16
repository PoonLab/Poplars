import unittest
import os
from io import StringIO
from poplars.sequence_locator import *

NUCL_QUERY = os.path.join(os.path.dirname(__file__), 'fixtures/nucl-query.fasta')
PROT_QUERY = os.path.join(os.path.dirname(__file__), 'fixtures/protein-query.fasta')

DEFAULT_HIV_GENOME = os.path.join(os.path.dirname(__file__), '../ref_genomes/K03455.fasta')
DEFAULT_HIV_PROTS = os.path.join(os.path.dirname(__file__), '../ref_genomes/K03455-protein.fasta')
TEST_HIV_GENOME = os.path.join(os.path.dirname(__file__), 'fixtures/hiv-test-genome.fasta')
TEST_HIV_PROTS = os.path.join(os.path.dirname(__file__), 'fixtures/hiv-test-proteins.fasta')

DEFAULT_SIV_GENOME = os.path.join(os.path.dirname(__file__), '../ref_genomes/M33262.fasta')
DEFAULT_SIV_PROTS = os.path.join(os.path.dirname(__file__), '../ref_genomes/M33262-protein.fasta')
TEST_SIV_GENOME = os.path.join(os.path.dirname(__file__), 'fixtures/siv-test-genome.fasta')
TEST_SIV_PROTS = os.path.join(os.path.dirname(__file__), 'fixtures/siv-test-proteins.fasta')

DEFAULT_HIV_NT_COORDS = os.path.abspath('../ref_genomes/K03455_genome_coordinates.csv')
DEFAULT_SIV_NT_COORDS = os.path.abspath('../ref_genomes/M33262_genome_coordinates.csv')
TEST_HIV_COORDS = os.path.join(os.path.dirname(__file__), 'fixtures/hiv_test_coords.csv')
TEST_SIV_COORDS = os.path.join(os.path.dirname(__file__), 'fixtures/siv_test_coords.csv')


class TestReadCoordinates(unittest.TestCase):

    def setUp(self):
        self.test_siv_coords = open(TEST_SIV_COORDS)

    def testSIVInputCoords(self):

        region_names = ['Rev(with intron)', 'Rev(exon1)', 'Rev(exon2)', 'Env', 'V1', 'V2', 'V3', 'V4', 'V5',
                        'RRE', 'gp120', 'gp41', 'Nef', '3\'LTR', '3\'LTR-R']

        region_coordinates = [[6784, 9315], [6784, 6853], [9062, 9315], [6860, 9499], [7196, 7360], [7364, 7492],
                              [7791, 7892], [8063, 8153], [8273, 8290], [8380, 8735], [6860, 8434], [8435, 9499],
                              [9333, 10124], [9719, 10535], [10235, 10411]]

        result = read_coordinates('siv', self.test_siv_coords)
        for i in range(len(result)):
            self.assertEqual(region_names[i], result[i].region_name)
            self.assertEqual(region_coordinates[i], result[i].nt_coords)

    def tearDown(self):
        self.test_siv_coords.close()


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

    def testDNALowerCase(self):
        expected = True
        result = valid_sequence("nucl", [[">h", "atgcgcgatgcagcac"]])
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

    def testAALowerCase(self):
        expected = True
        result = valid_sequence("prot", [[">header", "mwaglsalwgggp"]])
        self.assertEqual(expected, result)

    def testGappySequence(self):
        expected = True
        result = valid_sequence("nucl", [["header", "atgc-natgcgacgac"]])
        self.assertEqual(expected, result)


class TestValidInputs(unittest.TestCase):

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

    def setUp(self):
        self.hiv_genome_file = open(TEST_HIV_GENOME)
        self.siv_genome_file = open(TEST_SIV_GENOME)

    def testNucleotideQuery(self):
        expected = [["query", "ATGCGCG"]]
        handle = StringIO(">query\natgcgcg\n")
        result = get_query("nucl", handle)
        self.assertEqual(expected, result)

    def testProteinQuery(self):
        expected = [["query", "MPPLMMADLADLGG"]]
        handle = StringIO(">query\nMPPLMMADLADLGG\n")
        result = get_query("prot", handle)
        self.assertEqual(expected, result)

    def testLongNucleotideSequence(self):
        expected = [["query", "ATGCGCGAATTAGCGA"]]
        handle = StringIO(">query\natgcgcg\naattagcga\n")
        result = get_query("nucl", handle)
        self.assertEqual(expected, result)

    def testDefaultHIVGenome(self):
        expected = [["K03455|HIVHXB2CG", "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGG" 
                                         "ATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG" 
                                         "GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTAC" 
                                         "CAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC" 
                                         "TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGT"
                                         "GTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG"
                                         "AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG"]]
        result = get_query("nucl", self.hiv_genome_file)
        self.assertEqual(expected, result)

    def testInvalidNucleotideQuery(self):
        handle = StringIO(">query\natgcgcg&\n")
        with self.assertRaises(SystemExit) as e:
            get_query("nucl", handle)
        self.assertEqual(e.exception.code, 0)

    def testInvalidProteinQuery(self):
        handle = StringIO(">query\nMPPLMMAD>LADLGG\n")
        with self.assertRaises(SystemExit) as e:
            get_query("prot", handle)
        self.assertEqual(e.exception.code, 0)

    def testPlainText(self):
        handle = StringIO("atgatcg\n")
        expected = [["Sequence1", "ATGATCG"]]
        result = get_query("nucl", handle)
        self.assertEqual(expected, result)

    def testMultipleQueries(self):
        handle = StringIO("atgct--agc\natgca---ga\n")
        expected = [["Sequence1", "ATGCT--AGC"], ["Sequence2", "ATGCA---GA"]]
        result = get_query("nucl", handle)
        self.assertEqual(expected, result)

    def testMultipleFasta(self):
        handle = StringIO(">q1\natgct--agc\n>q2\natgca---ga\n")
        expected = [["q1", "ATGCT--AGC"], ["q2", "ATGCA---GA"]]
        result = get_query("nucl", handle)
        self.assertEqual(expected, result)

    def tearDown(self):
        self.hiv_genome_file.close()
        self.siv_genome_file.close()


class TestReverseComp(unittest.TestCase):

    def testSimpleUse(self):
        expected = "TCGCTAATTCGCGCATN*"
        result = reverse_comp("*NATGCGCGAATTAGCGA")
        self.assertEqual(expected, result)


class TestGetReferenceNucl(unittest.TestCase):

    def testDefaultHIVGenome(self):
        expected = [["K03455|HIVHXB2CG", "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGG"
                                         "ATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG"
                                         "GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTAC" 
                                         "CAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC" 
                                         "TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGT" 
                                         "GTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG" 
                                         "AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG"]]
        res = get_ref_nt_seq(DEFAULT_HIV_GENOME)
        # Check the first 350 nucleotides in the default reference sequence
        result = [[res[0][0], res[0][1][:350]]]
        self.assertEqual(expected, result)

    def testDefaultSIVGenome(self):
        expected = [["M33262|SIVMM239", "GCATGCACATTTTAAAGGCTTTTGCTAAATATAGCCAAAAGTCCTTCTAC"
                                        "AAATTTTCTAAGAGTTCTGATTCAAAGCAGTAACAGGCCTTGTCTCATCA"
                                        "TGAACTTTGGCATTTCATCTACAGCTAAGTTTATATCATAAATAGTTCTT"
                                        "TACAGGCAGCACCAACTTATACCCTTATAGCATACTTTACTGTGTGAAAA"
                                        "TTGCATCTTTCATTAAGCTTACTGTAAATTTACTGGCTGTCTTCCTTGCA"
                                        "GGTTTCTGGAAGGGATTTATTACAGTGCAAGAAGACATAGAATCTTAGAC"
                                        "ATATACTTAGAAAAGGAAGAAGGCATCATACCAGATTGGCAGGATTACAC"
                                        "CTCAGGACCAGGAATTAGATACCCAAAGACATTTGGCTGGCTATGGAAAT"]]
        res = get_ref_nt_seq(DEFAULT_SIV_GENOME)
        # Check the first 400 nucleotides in the default reference sequence
        result = [[res[0][0], res[0][1][:400]]]
        self.assertEqual(expected, result)


class TestGetReferenceProt(unittest.TestCase):

    def testDefaultSIVProteins(self):
        expected = [["Capsid|SIVMM239", "PVQQIGGNYVHLPLSPRTLNAWVKLIEEKKFGAEVVPGFQALSEGCTPYD"
                                        "INQMLNCVGDHQAAMQIIRDIINEEAADWDLQHPQPAPQQGQLREPSGSD"
                                        "IAGTTSSVDEQIQWMYRQQNPIPVGNIYRRWIQLGLQKCVRMYNPTNILD"
                                        "VKQGPKEPFQSYVDRFYKSLRAEQTDAAVKNWMTQTLLIQNANPDCKLVL"
                                        "KGLGVNPTLEEMLTACQGVGGPGQKARLM"],
                    ["p2|SIVMM239", "AEALKEALAPVPIPFAA"],
                    ["Nucleocapsid|SIVMM239", "AQQRGPRKPIKCWNCGKEGHSARQCRAPRRQGCWKCGKMDHVMAKCPDRQAG"]]
        result = get_ref_aa_seq(DEFAULT_SIV_PROTS)[2:5]
        self.assertEqual(expected, result)

    def testInputHIVProteins(self):
        expected = [["Gag|HIVHXB2CG", "MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPG"
                                      "LLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTK"
                                      "EALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAI"
                                      "SPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQA"
                                      "AMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQE"
                                      "QIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFR"
                                      "DYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATL"
                                      "EEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKC"
                                      "FNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSY"
                                      "KGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLR"
                                      "SLFGNDPSSQ"],
                    ["Matrix|HIVHXB2CG", "MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGL"
                                         "LETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEA"
                                         "LDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY"],
                    ["Capsid|HIVHXB2CG", "PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQ"
                                         "DLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPR"
                                         "GSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSI"
                                         "LDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKT"
                                         "ILKALGPAATLEEMMTACQGVGGPGHKARVL"]]
        result = get_ref_aa_seq(TEST_HIV_PROTS)
        self.assertEqual(expected, result)

    def testInputSIVProteins(self):
        expected = [["RNase|SIVMM239", "YTDGSCNKQSKEGKAGYITDRGKDKVKVLEQTTNQQAELEAFLMALTDSG"
                                       "PKANIIVDSQYVMGIITGCPTESESRLVNQIIEEMIKKSEIYVAWVPAHK"
                                       "GIGGNQEIDHLVSQGIRQVL"],
                    ["Integrase|SIVMM239", "FLEKIEPAQEEHDKYHSNVKELVFKFGLPRIVARQIVDTCDKCHQKGEAI"
                                           "HGQANSDLGTWQMDCTHLEGKIIIVAVHVASGFIEAEVIPQETGRQTALF"
                                           "LLKLAGRWPITHLHTDNGANFASQEVKMVAWWAGIEHTFGVPYNPQSQGV"
                                           "VEAMNHHLKNQIDRIREQANSVETIVLMAVHCMNFKRRGGIGDMTPAERL"
                                           "INMITTEQEIQFQQSKNSKFKNFRVYYREGRDQLWKGPGELLWKGEGAVI"
                                           "LKVGTDIKVVPRRKAKIIKDYGGGKEVDSSSHMEDTGEAREVA"],
                    ["Vif|SIVMM239", "MEEEKRWIAVPTWRIPERLERWHSLIKYLKYKTKDLQKVCYVPHFKVGWA"
                                     "WWTCSRVIFPLQEGSHLEVQGYWHLTPEKGWLSTYAVRITWYSKNFWTDV"
                                     "TPNYADILLHSTYFPCFTAGEVRRAIRGEQLLSCCRFPRAHKYQVPSLQY"
                                     "LALKVVSDVRSQGENPTWKQWRRDNRRGLRMAKQNSRGDKQRGGKPPTKG"
                                     "ANFPGLAKVLGILA"],
                    ["Vpx|SIVMM239", "MSDPRERIPPGNSGEETIGEAFEWLNRTVEEINREAVNHLPRELIFQVWQ"
                                     "RSWEYWHDEQGMSPSYVKYRYLCLIQKALFMHCKKGCRCLGEGHGAGGWR"
                                     "PGPPPPPPPGLA"],
                    ["Vpr|SIVMM239", "MEERPPENEGPQREPWDEWVVEVLEELKEEALKHFDPRLLTALGNHIYNR"
                                     "HGDTLEGAGELIRILQRALFMHFRGGCIHSRIGQPGGGNPLSAIPPSRSML"]]
        result = get_ref_aa_seq(TEST_SIV_PROTS)
        self.assertEqual(expected, result)


class TestFindGenomicRegions(unittest.TestCase):

    def setUp(self):
        self.hiv_genome_file = open(TEST_HIV_GENOME)
        self.siv_genome_file = open(TEST_SIV_GENOME)
        self.hiv_prot_file = open(TEST_HIV_PROTS)
        self.siv_prot_file = open(TEST_SIV_PROTS)

    def testSimpleUse(self):
        query = get_query('nucl', self.hiv_genome_file)
        reference_sequence = get_ref_nt_seq(DEFAULT_HIV_GENOME)
        sequence_alignment = sequence_align(query, reference_sequence)
        coordinates = get_region_coordinates(sequence_alignment[-1])

        result = find_genomic_regions('hiv', reference_sequence, coordinates)
        expected = ['5\'LTR', '5\'LTR-U3', 'Nef', '3\'LTR']
        self.assertEqual(expected, result)

    def testSIVAlignment(self):
        query = get_query('nucl', self.siv_genome_file)
        reference_sequence = get_ref_seq('hiv', 'nucl')
        sequence_alignment = sequence_align(query, reference_sequence)
        coordinates = get_region_coordinates(sequence_alignment[-1])
        result = find_genomic_regions('hiv', reference_sequence, coordinates)
        print(result)
        # self.assertEqual(expected, result)

    def tearDown(self):
        self.hiv_genome_file.close()
        self.siv_genome_file.close()
        self.hiv_prot_file.close()
        self.siv_prot_file.close()


class TestGetRegionCoordinates(unittest.TestCase):

    def setUp(self):
        self.hiv_genome_file = open(TEST_HIV_GENOME)
        self.siv_genome_file = open(TEST_SIV_GENOME)
        self.hiv_prot_file = open(TEST_HIV_PROTS)
        self.siv_prot_file = open(TEST_SIV_PROTS)

    def testHIVWithHIVGenome(self):
        query = get_query('nucl', self.hiv_genome_file)
        reference_sequence = get_ref_nt_seq('hiv', 'nucl')     # Reference genome is K03455.fasta
        sequence_alignment = sequence_align(query, reference_sequence)
        expected = [[0, 349]]
        result = get_region_coordinates(sequence_alignment[-1])
        self.assertEqual(expected, result)

    def testHIVWithSIVGenome(self):
        query = get_query('nucl', self.hiv_genome_file)
        reference_sequence = get_ref_nt_seq('siv', 'nucl')
        sequence_alignment = sequence_align(query, reference_sequence)
        expected = [[256, 402],   [414, 429],   [595, 611],   [2453, 2466], [2855, 2869], [3390, 3406],
                    [3560, 3573], [3856, 3867], [4501, 4518], [4556, 4566], [4750, 4758], [5912, 5927],
                    [6005, 6014], [6634, 6641], [8649, 8667], [8740, 8748], [10132, 10146]]
        result = get_region_coordinates(sequence_alignment[-1])
        self.assertEqual(expected, result)

    def tearDown(self):
        self.hiv_genome_file.close()
        self.siv_genome_file.close()
        self.hiv_prot_file.close()
        self.siv_prot_file.close()


class TestRetrieve(unittest.TestCase):

    def setUp(self):
        self.hiv_genome_file = open(TEST_HIV_GENOME)
        self.siv_genome_file = open(TEST_SIV_GENOME)
        self.hiv_prot_file = open(TEST_HIV_PROTS)
        self.siv_prot_file = open(TEST_SIV_PROTS)

    def testDefaultInput(self):
        reference_sequence = get_ref_seq('hiv', 'nucl')
        expected_seq = 'TTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTT'
        result_seq = retrieve('hiv', reference_sequence, 'p1')
        self.assertEqual(expected_seq, result_seq)

    def testSIVInput(self):
        reference_sequence = get_ref_seq('siv', 'nucl')
        expected_seq = 'GTATTCAAATTTGGATTACCCAGAATAGTGGCCAGACAGATAGTAGACACCTGTGATAAATGTCATCAGAAAGGAGAGG' \
                       'CTATACATGGGCAGGCAAATTCAGATCTAGGGACTTGGCAAAT'
        result_seq = retrieve('siv', reference_sequence, 'Integrase', 67, 188)
        self.assertEqual(expected_seq, result_seq)

    def testBeyondEnd(self):
        reference_sequence = get_ref_seq('hiv', 'nucl')
        expected_seq = 'ATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCACACAATGAATGGACACTAGAGCTTTTAGAGGAGCTTAA' \
                       'GAATGAAGCTGTTAGACATTTTCCTAGGATTTGGCTCCATGGCTTAGGGCAACATATCTATGAAACTTATGGGGATACTT' \
                       'GGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTGTTTATCCATTTTCAGAATTGGGTGTCGACATAGCAG' \
                       'AATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCCTAG'
        result_seq = retrieve('hiv', reference_sequence, 'Vpr', 1, 9000)
        self.assertEqual(expected_seq, result_seq)

    def tearDown(self):
        self.hiv_genome_file.close()
        self.siv_genome_file.close()
        self.hiv_prot_file.close()
        self.siv_prot_file.close()


class TestHandleArgs(unittest.TestCase):

    maxDiff = None

    def setUp(self):
        self.nucl_query = open(NUCL_QUERY)
        self.prot_query = open(PROT_QUERY)

        self.hiv_default_genome = open(DEFAULT_HIV_GENOME)
        self.hiv_default_prot = open(DEFAULT_HIV_PROTS)
        self.hiv_test_genome = open(TEST_HIV_GENOME)
        self.hiv_test_prot = open(TEST_HIV_PROTS)

        self.siv_default_genome = open(DEFAULT_SIV_GENOME)
        self.siv_default_prot = open(DEFAULT_SIV_PROTS)
        self.siv_test_genome = open(TEST_SIV_GENOME)
        self.siv_test_prot = open(TEST_SIV_PROTS)

    def testDefaultHIVNucl(self):
        """
        Tests the scenario when the user selects HIV and nucleotide alignment
        """

        result = handle_args('hiv', 'nucl', self.nucl_query, revcomp='n', ref_nt=None, ref_aa=None)

        expected_query = [['query', 'GACTCGAAAGCGAAAGTTCCAGAGAAGTTCTCTCGACGCAGGACTCGGCTTGCTGAGGTG'
                                    'CACACAGCAAGAGGCGAGAGCGGCGACTGGTGAGTACGCCAAATTTTGACTAGCGGAGGC'
                                    'TAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGTGGGGGAAAATTAGATGCATG'
                                    'GGAAAAAATTCGGTTACGGCCAGGGGGAAAGAAAAAATATAGAATGAAACATTTAGTATG'
                                    'GGCAAGCAGAGAGTTAGAAAGATTCGCACTTAACCC']]
        result_query = result[0]
        self.assertEqual(expected_query, result_query)

        expected_ref_nt_seq = [["K03455|HIVHXB2CG", "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGG"
                                                    "ATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG"
                                                    "GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTAC"
                                                    "CAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC"
                                                    "TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGT"
                                                    "GTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG"
                                                    "AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG"]]

        # Check the first 350 nucleotides in the default reference sequence
        seq = result[1][0][1][:350]
        header = result[1][0][0]
        result_ref_nt_seq = [[header, seq]]
        self.assertEqual(expected_ref_nt_seq, result_ref_nt_seq)

        expected_ref_aa_seq = [
            ['Gag|HIVHXB2CG',    'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRS'
                                 'LYNTVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNA'
                                 'WVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREP'
                                 'RGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQ'
                                 'ASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRN'
                                 'QRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPTAPPEES'
                                 'FRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ'],
            ['Matrix|HIVHXB2CG', 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRS'
                                 'LYNTVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY'],
            ['Capsid|HIVHXB2CG', 'PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEA'
                                 'AEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDI'
                                 'RQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVL'],
            ['p2|HIVHXB2CG',     'AEAMSQVTNSATIM']]

        # Check the first 4 proteins in the default reference protein sequence
        result_ref_aa_seq = result[2][:4]
        self.assertEqual(expected_ref_aa_seq, result_ref_aa_seq)

    def testDefaultSIVProt(self):
        """
        Tests the scenario when the user selects SIV and protein alignment
        """

        result = handle_args('siv', 'prot', self.nucl_query, revcomp='n', ref_nt=None, ref_aa=None)

        expected_query = [['query', 'GACTCGAAAGCGAAAGTTCCAGAGAAGTTCTCTCGACGCAGGACTCGGCTTGCTGAGGTG'
                                    'CACACAGCAAGAGGCGAGAGCGGCGACTGGTGAGTACGCCAAATTTTGACTAGCGGAGGC'
                                    'TAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGTGGGGGAAAATTAGATGCATG'
                                    'GGAAAAAATTCGGTTACGGCCAGGGGGAAAGAAAAAATATAGAATGAAACATTTAGTATG'
                                    'GGCAAGCAGAGAGTTAGAAAGATTCGCACTTAACCC']]
        result_query = result[0]
        self.assertEqual(expected_query, result_query)

        expected_ref_nt_seq = [["M33262|SIVMM239", "GCATGCACATTTTAAAGGCTTTTGCTAAATATAGCCAAAAGTCCTTCTAC"
                                                   "AAATTTTCTAAGAGTTCTGATTCAAAGCAGTAACAGGCCTTGTCTCATCA"
                                                   "TGAACTTTGGCATTTCATCTACAGCTAAGTTTATATCATAAATAGTTCTT"
                                                   "TACAGGCAGCACCAACTTATACCCTTATAGCATACTTTACTGTGTGAAAA"
                                                   "TTGCATCTTTCATTAAGCTTACTGTAAATTTACTGGCTGTCTTCCTTGCA"
                                                   "GGTTTCTGGAAGGGATTTATTACAGTGCAAGAAGACATAGAATCTTAGAC"
                                                   "ATATACTTAGAAAAGGAAGAAGGCATCATACCAGATTGGCAGGATTACAC"
                                                   "CTCAGGACCAGGAATTAGATACCCAAAGACATTTGGCTGGCTATGGAAAT"]]

        # Check the first 400 nucleotides in the default reference sequence
        seq = result[1][0][1][:400]
        header = result[1][0][0]
        result_ref_nt_seq = [[header, seq]]
        self.assertEqual(expected_ref_nt_seq, result_ref_nt_seq)

        expected_ref_aa_seq = [["Capsid|SIVMM239",      "PVQQIGGNYVHLPLSPRTLNAWVKLIEEKKFGAEVVPGFQALSEGCTPYD"
                                                        "INQMLNCVGDHQAAMQIIRDIINEEAADWDLQHPQPAPQQGQLREPSGSD"
                                                        "IAGTTSSVDEQIQWMYRQQNPIPVGNIYRRWIQLGLQKCVRMYNPTNILD"
                                                        "VKQGPKEPFQSYVDRFYKSLRAEQTDAAVKNWMTQTLLIQNANPDCKLVL"
                                                        "KGLGVNPTLEEMLTACQGVGGPGQKARLM"],
                              ["p2|SIVMM239",           "AEALKEALAPVPIPFAA"],
                              ["Nucleocapsid|SIVMM239", "AQQRGPRKPIKCWNCGKEGHSARQCRAPRRQGCWKCGKMDHVMAKCPDRQAG"]]
        result_ref_aa_seq = result[2][2:5]
        self.assertEqual(expected_ref_aa_seq, result_ref_aa_seq)

    def testDefaultRevCompHIVProt(self):
        """
        Tests the scenario when the user selects HIV and nucleotide alignment with the reverse complement of the query
        """

        result = handle_args('hiv', 'prot', self.nucl_query, revcomp='y', ref_nt=None, ref_aa=None)

        expected_query = [['query', 'GACTCGAAAGCGAAAGTTCCAGAGAAGTTCTCTCGACGCAGGACTCGGCTTGCTGAGGTG'
                                    'CACACAGCAAGAGGCGAGAGCGGCGACTGGTGAGTACGCCAAATTTTGACTAGCGGAGGC'
                                    'TAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGTGGGGGAAAATTAGATGCATG'
                                    'GGAAAAAATTCGGTTACGGCCAGGGGGAAAGAAAAAATATAGAATGAAACATTTAGTATG'
                                    'GGCAAGCAGAGAGTTAGAAAGATTCGCACTTAACCC']]
        result_query = result[0]
        self.assertEqual(expected_query, result_query)

        expected_ref_nt_seq = [["K03455|HIVHXB2CG", "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGG"
                                                    "ATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG"
                                                    "GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTAC"
                                                    "CAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC"
                                                    "TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGT"
                                                    "GTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG"
                                                    "AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG"]]

        # Check the first 350 nucleotides in the default reference sequence
        seq = result[1][0][1][:350]
        header = result[1][0][0]
        result_ref_nt_seq = [[header, seq]]
        self.assertEqual(expected_ref_nt_seq, result_ref_nt_seq)

        expected_ref_aa_seq = [
            ['Gag|HIVHXB2CG', 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRS'
                              'LYNTVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNA'
                              'WVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREP'
                              'RGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQ'
                              'ASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRN'
                              'QRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPTAPPEES'
                              'FRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ'],
            ['Matrix|HIVHXB2CG', 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRS'
                                 'LYNTVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY'],
            ['Capsid|HIVHXB2CG', 'PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEA'
                                 'AEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDI'
                                 'RQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVL'],
            ['p2|HIVHXB2CG', 'AEAMSQVTNSATIM']]

        # Check the first 4 proteins in the default reference protein sequence
        result_ref_aa_seq = result[2][:4]
        self.assertEqual(expected_ref_aa_seq, result_ref_aa_seq)

    def testDefaultRevCompSIVNucl(self):
        """
        Tests the scenario when the user selects SIV and nucleotide alignment with the reverse complement of the query
        """

        result = handle_args('siv', 'nucl', self.nucl_query, revcomp='y', ref_nt=None, ref_aa=None)

        expected_query = [['query', 'GGGTTAAGTGCGAATCTTTCTAACTCTCTGCTTGCCCATACTAAATGTTTCATTCTATATTTTTTC'
                                    'TTTCCCCCTGGCCGTAACCGAATTTTTTCCCATGCATCTAATTTTCCCCCACTTAATACTGACGCT'
                                    'CTCGCACCCATCTCTCTCCTTCTAGCCTCCGCTAGTCAAAATTTGGCGTACTCACCAGTCGCCGCT'
                                    'CTCGCCTCTTGCTGTGTGCACCTCAGCAAGCCGAGTCCTGCGTCGAGAGAACTTCTCTGGAACTTT'
                                    'CGCTTTCGAGTC']]
        result_query = result[0]
        self.assertEqual(expected_query, result_query)

        expected_ref_nt_seq = [["M33262|SIVMM239", "GCATGCACATTTTAAAGGCTTTTGCTAAATATAGCCAAAAGTCCTTCTAC"
                                                   "AAATTTTCTAAGAGTTCTGATTCAAAGCAGTAACAGGCCTTGTCTCATCA"
                                                   "TGAACTTTGGCATTTCATCTACAGCTAAGTTTATATCATAAATAGTTCTT"
                                                   "TACAGGCAGCACCAACTTATACCCTTATAGCATACTTTACTGTGTGAAAA"
                                                   "TTGCATCTTTCATTAAGCTTACTGTAAATTTACTGGCTGTCTTCCTTGCA"
                                                   "GGTTTCTGGAAGGGATTTATTACAGTGCAAGAAGACATAGAATCTTAGAC"
                                                   "ATATACTTAGAAAAGGAAGAAGGCATCATACCAGATTGGCAGGATTACAC"
                                                   "CTCAGGACCAGGAATTAGATACCCAAAGACATTTGGCTGGCTATGGAAAT"]]

        # Check the first 400 nucleotides in the default reference sequence
        seq = result[1][0][1][:400]
        header = result[1][0][0]
        result_ref_nt_seq = [[header, seq]]
        self.assertEqual(expected_ref_nt_seq, result_ref_nt_seq)

        expected_ref_aa_seq = [["Capsid|SIVMM239", "PVQQIGGNYVHLPLSPRTLNAWVKLIEEKKFGAEVVPGFQALSEGCTPYD"
                                                   "INQMLNCVGDHQAAMQIIRDIINEEAADWDLQHPQPAPQQGQLREPSGSD"
                                                   "IAGTTSSVDEQIQWMYRQQNPIPVGNIYRRWIQLGLQKCVRMYNPTNILD"
                                                   "VKQGPKEPFQSYVDRFYKSLRAEQTDAAVKNWMTQTLLIQNANPDCKLVL"
                                                   "KGLGVNPTLEEMLTACQGVGGPGQKARLM"],
                               ["p2|SIVMM239", "AEALKEALAPVPIPFAA"],
                               ["Nucleocapsid|SIVMM239", "AQQRGPRKPIKCWNCGKEGHSARQCRAPRRQGCWKCGKMDHVMAKCPDRQAG"]]
        result_ref_aa_seq = result[2][2:5]
        self.assertEqual(expected_ref_aa_seq, result_ref_aa_seq)

    def testInputNuclProtHIV(self):
        """
        Tests the scenario when the user selects HIV, specifies the reference sequences,
        and selects nucleotide alignment
        """
        result = handle_args('hiv', 'nucl', self.nucl_query, revcomp='n',
                             ref_nt=TEST_HIV_GENOME, ref_aa=TEST_HIV_PROTS)

        expected_query = [['query', 'GACTCGAAAGCGAAAGTTCCAGAGAAGTTCTCTCGACGCAGGACTCGGCTTGCTGAGGTG'
                                    'CACACAGCAAGAGGCGAGAGCGGCGACTGGTGAGTACGCCAAATTTTGACTAGCGGAGGC'
                                    'TAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGTGGGGGAAAATTAGATGCATG'
                                    'GGAAAAAATTCGGTTACGGCCAGGGGGAAAGAAAAAATATAGAATGAAACATTTAGTATG'
                                    'GGCAAGCAGAGAGTTAGAAAGATTCGCACTTAACCC']]
        result_query = result[0]
        self.assertEqual(expected_query, result_query)

        expected_ref_nt = [["K03455|HIVHXB2CG", "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGG"
                                                "ATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG"
                                                "GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTAC"
                                                "CAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC"
                                                "TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGT"
                                                "GTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG"
                                                "AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG"]]
        result_ref_nt = result[1]
        self.assertEqual(expected_ref_nt, result_ref_nt)

        expected_ref_prot = [["Gag|HIVHXB2CG",    "MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPG"
                                                  "LLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTK"
                                                  "EALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAI"
                                                  "SPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQA"
                                                  "AMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQE"
                                                  "QIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFR"
                                                  "DYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATL"
                                                  "EEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKC"
                                                  "FNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSY"
                                                  "KGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLR"
                                                  "SLFGNDPSSQ"],
                            ["Matrix|HIVHXB2CG",  "MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGL"
                                                  "LETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEA"
                                                  "LDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY"],
                            ["Capsid|HIVHXB2CG",  "PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQ"
                                                  "DLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPR"
                                                  "GSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSI"
                                                  "LDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKT"
                                                  "ILKALGPAATLEEMMTACQGVGGPGHKARVL"]]
        result_ref_prot = result[2]
        self.assertEqual(expected_ref_prot, result_ref_prot)

    def testInputAllSIV(self):
        """
        Tests the scenario when the user selects SIV, specifies the reference sequences, and selects protein alignment
        """
        result = handle_args('siv', 'prot', self.prot_query, revcomp='n',
                             ref_nt=TEST_SIV_GENOME, ref_aa=TEST_SIV_PROTS)

        expected_query = [['query', 'MAYTQTATTSALLDTVRGNNSLVNDLAKRRLYDTAVEEFNARDRRPKVNFSKVISEEQTL'
                                    'IATRAYPEFQITFYNTQNAVHSLAGGLRSLELEYLMMQIPYGSLTYDIGGNFASHLFKGR'
                                    'AYVHCCMPNLDVRDIMRHEGQKDSIELYLSRLERGGKTVPNFQKEAFDRYAEIPEDAVCH'
                                    'NTFQTMRHQPMQQSGRVYAIALHSIYDIPADEFGAALLRKNVHTCYAAFHFSENLLLEDS'
                                    'YVNLDEINACFSRDGDKLTFSFASESTLNYCHSYSNILKYVCKTYFPASNREVYMKEFLV']]
        result_query = result[0]
        self.assertEqual(expected_query, result_query)

        expected_ref_nt = [["M33262|SIVMM239", "GCATGCACATTTTAAAGGCTTTTGCTAAATATAGCCAAAAGTCCTTCTAC"
                                               "AAATTTTCTAAGAGTTCTGATTCAAAGCAGTAACAGGCCTTGTCTCATCA"
                                               "TGAACTTTGGCATTTCATCTACAGCTAAGTTTATATCATAAATAGTTCTT"
                                               "TACAGGCAGCACCAACTTATACCCTTATAGCATACTTTACTGTGTGAAAA"
                                               "TTGCATCTTTCATTAAGCTTACTGTAAATTTACTGGCTGTCTTCCTTGCA"
                                               "GGTTTCTGGAAGGGATTTATTACAGTGCAAGAAGACATAGAATCTTAGAC"
                                               "ATATACTTAGAAAAGGAAGAAGGCATCATACCAGATTGGCAGGATTACAC"
                                               "CTCAGGACCAGGAATTAGATACCCAAAGACATTTGGCTGGCTATGGAAAT"]]
        result_ref_nt = result[1]
        self.assertEqual(expected_ref_nt, result_ref_nt)

        expected_ref_prot = [["RNase|SIVMM239",     "YTDGSCNKQSKEGKAGYITDRGKDKVKVLEQTTNQQAELEAFLMALTDSG"
                                                    "PKANIIVDSQYVMGIITGCPTESESRLVNQIIEEMIKKSEIYVAWVPAHK"
                                                    "GIGGNQEIDHLVSQGIRQVL"],
                            ["Integrase|SIVMM239",  "FLEKIEPAQEEHDKYHSNVKELVFKFGLPRIVARQIVDTCDKCHQKGEAI"
                                                    "HGQANSDLGTWQMDCTHLEGKIIIVAVHVASGFIEAEVIPQETGRQTALF"
                                                    "LLKLAGRWPITHLHTDNGANFASQEVKMVAWWAGIEHTFGVPYNPQSQGV"
                                                    "VEAMNHHLKNQIDRIREQANSVETIVLMAVHCMNFKRRGGIGDMTPAERL"
                                                    "INMITTEQEIQFQQSKNSKFKNFRVYYREGRDQLWKGPGELLWKGEGAVI"
                                                    "LKVGTDIKVVPRRKAKIIKDYGGGKEVDSSSHMEDTGEAREVA"],
                            ["Vif|SIVMM239",        "MEEEKRWIAVPTWRIPERLERWHSLIKYLKYKTKDLQKVCYVPHFKVGWA"
                                                    "WWTCSRVIFPLQEGSHLEVQGYWHLTPEKGWLSTYAVRITWYSKNFWTDV"
                                                    "TPNYADILLHSTYFPCFTAGEVRRAIRGEQLLSCCRFPRAHKYQVPSLQY"
                                                    "LALKVVSDVRSQGENPTWKQWRRDNRRGLRMAKQNSRGDKQRGGKPPTKG"
                                                    "ANFPGLAKVLGILA"],
                            ["Vpx|SIVMM239",        "MSDPRERIPPGNSGEETIGEAFEWLNRTVEEINREAVNHLPRELIFQVWQ"
                                                    "RSWEYWHDEQGMSPSYVKYRYLCLIQKALFMHCKKGCRCLGEGHGAGGWR"
                                                    "PGPPPPPPPGLA"],
                            ["Vpr|SIVMM239",        "MEERPPENEGPQREPWDEWVVEVLEELKEEALKHFDPRLLTALGNHIYNR"
                                                    "HGDTLEGAGELIRILQRALFMHFRGGCIHSRIGQPGGGNPLSAIPPSRSML"]]
        result_ref_prot = result[2]
        self.assertEqual(expected_ref_prot, result_ref_prot)

    def tearDown(self):
        self.nucl_query.close()
        self.prot_query.close()

        self.hiv_default_genome.close()
        self.hiv_default_prot.close()
        self.hiv_test_genome.close()
        self.hiv_test_prot.close()

        self.siv_default_genome.close()
        self.siv_default_prot.close()
        self.siv_test_genome.close()
        self.siv_test_prot.close()

