# TODO: remove redundant setUp and tearDown functions

import unittest
import os
from poplars.sequence_locator import valid_sequence
from poplars.sequence_locator import valid_inputs
from poplars.sequence_locator import get_query
from poplars.sequence_locator import get_ref_seq
from poplars.sequence_locator import sequence_align
from poplars.sequence_locator import make_aa_dict
from poplars.sequence_locator import find_genomic_regions
from poplars.sequence_locator import find_aa_regions
from poplars.sequence_locator import get_region_coordinates
from poplars.sequence_locator import get_matches
from poplars.sequence_locator import retrieve

TEST_HIV_GENOME = os.path.join(os.path.dirname(__file__), 'fixtures/hiv-test-genome.fasta')
TEST_SIV_GENOME = os.path.join(os.path.dirname(__file__), 'fixtures/siv-test-genome.fasta')
TEST_HIV_PROTS = os.path.join(os.path.dirname(__file__), 'fixtures/hiv-test-proteins.fasta')
TEST_SIV_PROTS = os.path.join(os.path.dirname(__file__), 'fixtures/siv-test-proteins.fasta')


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
        expected = "ATGCGCG"
        result = get_query("nucl", [">query", "atgcgcg"])
        self.assertEqual(expected, result)

    def testProteinQuery(self):
        expected = "MPPLMMADLADLGG"
        result = get_query("prot", [">query", "MPPLMMADLADLGG"])
        self.assertEqual(expected, result)

    def testLongNucleotideSequence(self):
        expected = "ATGCGCGAATTAGCGA"
        result = get_query("nucl", [">query", "atgcgcg", "aattagcga"])
        self.assertEqual(expected, result)

    def testDefaultHIVGenome(self):
        expected = ("TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGG" 
                    "ATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG" 
                    "GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTAC" 
                    "CAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC" 
                    "TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGT"
                    "GTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG"
                    "AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG")
        result = get_query("nucl", self.hiv_genome_file)
        self.assertEqual(expected, result)

    def testInvalidNucleotideQuery(self):
        with self.assertRaises(SystemExit) as e:
            get_query("nucl", [">query", "atgcgcg*"])
        self.assertEqual(e.exception.code, 0)

    def testInvalidProteinQuery(self):
        with self.assertRaises(SystemExit) as e:
            get_query("prot", [">query", "MPPLMMAD>LADLGG"])
        self.assertEqual(e.exception.code, 0)

    def tearDown(self):
        self.hiv_genome_file.close()
        self.siv_genome_file.close()


class TestGetReferenceSequence(unittest.TestCase):
    """
    Adapted from https://stackoverflow.com/questions/8047736/how-to-load-data-from-a-file-for-a-unit-test-in-python
    """

    def setUp(self):
        self.hiv_genome_file = open(TEST_HIV_GENOME)
        self.siv_genome_file = open(TEST_SIV_GENOME)
        self.hiv_prot_file = open(TEST_HIV_PROTS)
        self.siv_prot_file = open(TEST_SIV_PROTS)

    def testDefaultHIVGenome(self):
        expected = [["K03455|HIVHXB2CG", "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGG"
                                         "ATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG"
                                         "GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTAC" 
                                         "CAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC" 
                                         "TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGT" 
                                         "GTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG" 
                                         "AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG"]]
        res = get_ref_seq("hiv", "nucl")
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
        res = get_ref_seq("siv", "nucl")
        # Check the first 400 nucleotides in the default reference sequence
        result = [[res[0][0], res[0][1][:400]]]
        self.assertEqual(expected, result)

    def testInputHIVGenome(self):
        expected = [["K03455|HIVHXB2CG", "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGG"
                                         "ATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG"
                                         "GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTAC"
                                         "CAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC"
                                         "TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGT"
                                         "GTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG"
                                         "AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG"]]
        result = get_ref_seq("hiv", "nucl", self.hiv_genome_file)
        self.assertEqual(expected, result)

    def testInputSIVGenome(self):
        expected = [["M33262|SIVMM239", "GCATGCACATTTTAAAGGCTTTTGCTAAATATAGCCAAAAGTCCTTCTAC"
                                        "AAATTTTCTAAGAGTTCTGATTCAAAGCAGTAACAGGCCTTGTCTCATCA"
                                        "TGAACTTTGGCATTTCATCTACAGCTAAGTTTATATCATAAATAGTTCTT"
                                        "TACAGGCAGCACCAACTTATACCCTTATAGCATACTTTACTGTGTGAAAA"
                                        "TTGCATCTTTCATTAAGCTTACTGTAAATTTACTGGCTGTCTTCCTTGCA"
                                        "GGTTTCTGGAAGGGATTTATTACAGTGCAAGAAGACATAGAATCTTAGAC"
                                        "ATATACTTAGAAAAGGAAGAAGGCATCATACCAGATTGGCAGGATTACAC"
                                        "CTCAGGACCAGGAATTAGATACCCAAAGACATTTGGCTGGCTATGGAAAT"]]
        result = get_ref_seq("siv", "nucl", self.siv_genome_file)
        self.assertEqual(expected, result)

    def testDefaultHIVProteins(self):
        expected = [["Gag|HIVHXB2CG",    "MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGL"
                                         "LETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEA"
                                         "LDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPR"
                                         "TLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQM"
                                         "LKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWM"
                                         "TNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRF"
                                         "YKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTAC"
                                         "QGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGH"
                                         "TARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSYKGRPGNFLQ"
                                         "SRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ"],
                    ["Matrix|HIVHXB2CG", "MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGL"
                                         "LETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEA"
                                         "LDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY"]]
        result = get_ref_seq("hiv", "prot")[:2]
        self.assertEqual(expected, result)

    def testDefaultSIVProteins(self):
        expected = [["Capsid|SIVMM239",       "PVQQIGGNYVHLPLSPRTLNAWVKLIEEKKFGAEVVPGFQALSEGCTPYD"
                                              "INQMLNCVGDHQAAMQIIRDIINEEAADWDLQHPQPAPQQGQLREPSGSD"
                                              "IAGTTSSVDEQIQWMYRQQNPIPVGNIYRRWIQLGLQKCVRMYNPTNILD"
                                              "VKQGPKEPFQSYVDRFYKSLRAEQTDAAVKNWMTQTLLIQNANPDCKLVL"
                                              "KGLGVNPTLEEMLTACQGVGGPGQKARLM"],
                    ["p2|SIVMM239",           "AEALKEALAPVPIPFAA"],
                    ["Nucleocapsid|SIVMM239", "AQQRGPRKPIKCWNCGKEGHSARQCRAPRRQGCWKCGKMDHVMAKCPDRQAG"]]
        result = get_ref_seq("siv", "prot")[2:5]
        self.assertEqual(expected, result)

    def testInputHIVProteins(self):
        expected = [["Gag|HIVHXB2CG",    "MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPG"
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
        result = get_ref_seq("hiv", "prot", self.hiv_prot_file)
        self.assertEqual(expected, result)

    def testInputSIVProteins(self):
        expected = [["RNase|SIVMM239",     "YTDGSCNKQSKEGKAGYITDRGKDKVKVLEQTTNQQAELEAFLMALTDSG"
                                           "PKANIIVDSQYVMGIITGCPTESESRLVNQIIEEMIKKSEIYVAWVPAHK"
                                           "GIGGNQEIDHLVSQGIRQVL"],
                    ["Integrase|SIVMM239", "FLEKIEPAQEEHDKYHSNVKELVFKFGLPRIVARQIVDTCDKCHQKGEAI"
                                           "HGQANSDLGTWQMDCTHLEGKIIIVAVHVASGFIEAEVIPQETGRQTALF"
                                           "LLKLAGRWPITHLHTDNGANFASQEVKMVAWWAGIEHTFGVPYNPQSQGV"
                                           "VEAMNHHLKNQIDRIREQANSVETIVLMAVHCMNFKRRGGIGDMTPAERL"
                                           "INMITTEQEIQFQQSKNSKFKNFRVYYREGRDQLWKGPGELLWKGEGAVI"
                                           "LKVGTDIKVVPRRKAKIIKDYGGGKEVDSSSHMEDTGEAREVA"],
                    ["Vif|SIVMM239",       "MEEEKRWIAVPTWRIPERLERWHSLIKYLKYKTKDLQKVCYVPHFKVGWA"
                                           "WWTCSRVIFPLQEGSHLEVQGYWHLTPEKGWLSTYAVRITWYSKNFWTDV"
                                           "TPNYADILLHSTYFPCFTAGEVRRAIRGEQLLSCCRFPRAHKYQVPSLQY"
                                           "LALKVVSDVRSQGENPTWKQWRRDNRRGLRMAKQNSRGDKQRGGKPPTKG"
                                           "ANFPGLAKVLGILA"],
                    ["Vpx|SIVMM239",       "MSDPRERIPPGNSGEETIGEAFEWLNRTVEEINREAVNHLPRELIFQVWQ"
                                           "RSWEYWHDEQGMSPSYVKYRYLCLIQKALFMHCKKGCRCLGEGHGAGGWR"
                                           "PGPPPPPPPGLA"],
                    ["Vpr|SIVMM239",       "MEERPPENEGPQREPWDEWVVEVLEELKEEALKHFDPRLLTALGNHIYNR"
                                           "HGDTLEGAGELIRILQRALFMHFRGGCIHSRIGQPGGGNPLSAIPPSRSML"]]
        result = get_ref_seq("siv", "prot", self.siv_prot_file)
        self.assertEqual(expected, result)

    def tearDown(self):
        self.hiv_genome_file.close()
        self.siv_genome_file.close()
        self.hiv_prot_file.close()
        self.siv_prot_file.close()


class TestMakeAADictionary(unittest.TestCase):

    def setUp(self):
        self.hiv_prot_file = open(TEST_HIV_PROTS)
        self.siv_prot_file = open(TEST_SIV_PROTS)

    def testDefaultSequence(self):
        expected = {"Matrix|SIVMM239":  "MGVRNSVLSGKKADELEKIRLRPNGKKKYMLKHVVWAANELDRFGLAESL"
                                        "LENKEGCQKILSVLAPLVPTGSENLKSLYNTVCVIWCIHAEEKVKHTEEA"
                                        "KQIVQRHLVVETGTTETMPKTSRPTAPSSGRGGNY"}
        ref_aa_seq = get_ref_seq('siv', 'prot')
        result = make_aa_dict(ref_aa_seq[1:2])      # second sequence
        self.assertEqual(expected, result)

    def testInputSequence(self):
        expected = {"Gag|HIVHXB2CG":    "MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGL"
                                        "LETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEA"
                                        "LDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPR"
                                        "TLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQM"
                                        "LKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWM"
                                        "TNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRF"
                                        "YKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTAC"
                                        "QGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGH"
                                        "TARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSYKGRPGNFLQ"
                                        "SRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ"}
        ref_aa_seq = get_ref_seq('hiv', 'prot', self.hiv_prot_file)
        result = make_aa_dict(ref_aa_seq[:1])       # first sequence
        self.assertEqual(expected, result)

    def tearDown(self):
        self.hiv_prot_file.close()
        self.siv_prot_file.close()


class TestFindGenomicRegions(unittest.TestCase):

    def setUp(self):
        self.hiv_genome_file = open(TEST_HIV_GENOME)
        self.siv_genome_file = open(TEST_SIV_GENOME)
        self.hiv_prot_file = open(TEST_HIV_PROTS)
        self.siv_prot_file = open(TEST_SIV_PROTS)

    def testSimpleUse(self):
        query = get_query('nucl', self.hiv_genome_file)
        reference_sequence = get_ref_seq('hiv', 'nucl')
        sequence_alignment = sequence_align(query, reference_sequence)
        coordinates = get_region_coordinates(sequence_alignment[-1])
        result = find_genomic_regions('hiv', reference_sequence, coordinates)
        expected = [('5\'LTR',    "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGG"
                                  "ATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG"
                                  "GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTAC"
                                  "CAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC"
                                  "TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGT"
                                  "GTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG"
                                  "AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG"),

                    ('5\'LTR-U3', "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGG"
                                  "ATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG"
                                  "GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTAC"
                                  "CAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC"
                                  "TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGT"
                                  "GTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG"
                                  "AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG")]
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
        reference_sequence = get_ref_seq('hiv', 'nucl')     # Reference genome is K03455.fasta
        sequence_alignment = sequence_align(query, reference_sequence)
        expected = [[0, 349]]
        result = get_region_coordinates(sequence_alignment[-1])
        self.assertEqual(expected, result)

    def testHIVWithSIVGenome(self):
        query = get_query('nucl', self.hiv_genome_file)
        reference_sequence = get_ref_seq('siv', 'nucl')
        sequence_alignment = sequence_align(query, reference_sequence)
        expected = [[256, 402], [414, 429], [595, 611], [2453, 2466], [2855, 2869], [3390, 3406],
                    [3560, 3573], [3856, 3867], [4501, 4518], [4556, 4566], [4750, 4758], [5912, 5927],
                    [6005, 6014], [6634, 6641], [8649, 8667], [8740, 8748], [10132, 10146]]
        result = get_region_coordinates(sequence_alignment[-1])
        self.assertEqual(expected, result)

    def tearDown(self):
        self.hiv_genome_file.close()
        self.siv_genome_file.close()
        self.hiv_prot_file.close()
        self.siv_prot_file.close()
