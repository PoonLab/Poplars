import unittest
from io import StringIO
from poplars.sequence_locator import *

DEFAULT_HIV_GENOME = os.path.join(os.path.dirname(__file__), '../ref_genomes/K03455.fasta')
TEST_HIV_GENOME = os.path.join(os.path.dirname(__file__), 'fixtures/hiv-test-genome.fasta')
DEFAULT_SIV_GENOME = os.path.join(os.path.dirname(__file__), '../ref_genomes/M33262.fasta')
TEST_SIV_GENOME = os.path.join(os.path.dirname(__file__), 'fixtures/siv-test-genome.fasta')

DEFAULT_HIV_NT_COORDS = os.path.abspath('../ref_genomes/K03455_genome_coordinates.csv')
TEST_HIV_NT_COORDS = os.path.join(os.path.dirname(__file__), 'fixtures/hiv_test_nt_coords.csv')
DEFAULT_SIV_NT_COORDS = os.path.abspath('../ref_genomes/M33262_genome_coordinates.csv')
TEST_SIV_NT_COORDS = os.path.join(os.path.dirname(__file__), 'fixtures/siv_test_nt_coords.csv')

DEFAULT_HIV_PROTS = os.path.join(os.path.dirname(__file__), '../ref_genomes/K03455-protein.fasta')
TEST_HIV_PROTS = os.path.join(os.path.dirname(__file__), 'fixtures/hiv-test-proteins.fasta')
DEFAULT_SIV_PROTS = os.path.join(os.path.dirname(__file__), '../ref_genomes/M33262-protein.fasta')
TEST_SIV_PROTS = os.path.join(os.path.dirname(__file__), 'fixtures/siv-test-proteins.fasta')

DEFAULT_HIV_AA_COORDS = os.path.abspath('../ref_genomes/K03455_protein_coordinates.csv')
TEST_HIV_AA_COORDS = os.path.join(os.path.dirname(__file__), 'fixtures/hiv_test_prot_coords.csv')
DEFAULT_SIV_AA_COORDS = os.path.abspath('../ref_genomes/M33262_protein_coordinates.csv')
TEST_SIV_AA_COORDS = os.path.join(os.path.dirname(__file__), 'fixtures/siv_test_prot_coords.csv')


class TestGetCoords(unittest.TestCase):

    def setUp(self):
        self.test_nt_coords_file = open(TEST_HIV_NT_COORDS)
        self.nt_coords = self.test_nt_coords_file.read()

        self.test_aa_coords_file = open(TEST_SIV_AA_COORDS)
        self.aa_coords = self.test_aa_coords_file.read()

    def testNuclCoords(self):
        region = GenomeRegion('test', self.nt_coords)
        result = region.get_coords('nucl')
        expected = '5\'LTR,1,634\nGag,790,2292\nMatrix,790,1185\nCapsid,1186,1878\np2,1879,1920\n' \
                   'Nucleocapsid,1921,2085\np1,2086,2133\np6,2134,2292'
        self.assertEqual(expected, result)

    def testProtCoords(self):
        region = GenomeRegion('test', None, None, self.aa_coords)
        result = region.get_coords('prot')
        expected = 'Rev,2087,2193\nRev(exon1),2087,2110\nRev(exon2),2111,2193\ngp160,2194,3072\n' \
                   'gp120,2194,2718\ngp41,2719,3072\nNef,3073,3335'
        self.assertEqual(expected, result)

    def tearDown(self):
        self.test_nt_coords_file.close()
        self.test_aa_coords_file.close()


class TestGetSequence(unittest.TestCase):

    def setUp(self):
        self.siv_genome_file = open(TEST_SIV_GENOME)
        self.nt_seq = self.siv_genome_file.read()

        self.hiv_protein_file = open(TEST_HIV_PROTS)
        self.aa_seq = self.hiv_protein_file.read()

    def testGetNuclSeq(self):
        region = GenomeRegion('test', None, self.nt_seq)
        result = region.get_sequence('nucl')[0:1218]
        expected = '>M33262|SIVMM239\n' \
                   'AATGGTGATTATTCAGAAGTGGCCCTTAATGTTACAGAAAGCTTTGATGCCTGGAATAATACAGTCACAGAACAGGCAAT\n' \
                   'AGAGGATGTATGGCAACTCTTTGAGACCTCAATAAAGCCTTGTGTAAAATTATCCCCATTATGCATTACTATGAGATGCA\n' \
                   'ATAAAAGTGAGACAGATAGATGGGGATTGACAAAATCAATAACAACAACAGCATCAACAACATCAACGACAGCATCAGCA\n' \
                   'AAAGTAGACATGGTCAATGAGACTAGTTCTTGTATAGCCCAGGATAATTGCACAGGCTTGGAACAAGAGCAAATGATAAG\n' \
                   'CTGTAAATTCAACATGACAGGGTTAAAAAGAGACAAGAAAAAAGAGTACAATGAAACTTGGTACTCTGCAGATTTGGTAT\n' \
                   'GTGAACAAGGGAATAACACTGGTAATGAAAGTAGATGTTACATGAACCACTGTAACACTTCTGTTATCCAAGAGTCTTGT\n' \
                   'GACAAACATTATTGGGATGCTATTAGATTTAGGTATTGTGCACCTCCAGGTTATGCTTTGCTTAGATGTAATGACACAAA\n' \
                   'TTATTCAGGCTTTATGCCTAAATGTTCTAAGGTGGTGGTCTCTTCATGCACAAGGATGATGGAGACACAGACTTCTACTT\n' \
                   'GGTTTGGCTTTAATGGAACTAGAGCAGAAAATAGAACTTATATTTACTGGCATGGTAGGGATAATAGGACTATAATTAGT\n' \
                   'TTAAATAAGTATTATAATCTAACAATGAAATGTAGAAGACCAGGAAATAAGACAGTTTTACCAGTCACCATTATGTCTGG\n' \
                   'ATTGGTTTTCCACTCACAACCAATCAATGATAGGCCAAAGCAGGCATGGTGTTGGTTTGGAGGAAAATGGAAGGATGCAA\n' \
                   'TAAAAGAGGTGAAGCAGACCATTGTCAAACATCCCAGGTATACTGGAACTAACAATACTGATAAAATCAATTTGACGGCT\n' \
                   'CCTGGAGGAGGAGATCCGGAAGTTACCTTCATGTGGACAAATTGCAGAGGAGAGTTCCTCTACTGTAAAATGAATTGGTT\n' \
                   'TCTAAATTGGGTAGAAGATAGGAATACAGCTAACCAGAAGCCAAAGGAACAGCATAAAAGGAATTACGTGCCATGTCATA\n' \
                   'TTAGACAAATAATCAACACTTGGCATAAAGTAGGCAAAAATGTTTATTTGCCTCCAAGAGAGGGAGA'
        self.assertEqual(expected, result)

    def testGetProtSeq(self):
        region = GenomeRegion('test', None, None, None, self.aa_seq)
        result = region.get_sequence('prot')
        expected = '>test|HIVHXB2CG\n' \
                   'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYN\n' \
                   'TVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVE\n' \
                   'EKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTT\n' \
                   'STLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTET\n' \
                   'LLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGH\n' \
                   'TARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPID\n' \
                   'KELYPLTSLRSLFGNDPSSQ'
        self.assertEqual(expected, result)

    def tearDown(self):
        self.siv_genome_file.close()
        self.hiv_protein_file.close()


class TestSetCoords(unittest.TestCase):

    def testSetNuclCoords(self):
        region = GenomeRegion('test', 'ATGCGACGACGACGACGTTAGCAGTCGATCATCATGCTGATC', None, None, None)
        region.set_coords([1, 10], 'nucl')
        expected = [1, 10]
        result = region.nt_coords
        self.assertEqual(expected, result)

    def testSetProtCoords(self):
        region = GenomeRegion('test', None, None, 'ATVLHGLKLKAMAIMTIM')
        region.set_coords([1, 900], 'prot')
        expected = [1, 900]
        result = region.aa_coords
        self.assertEqual(expected, result)


class TestSetSeqFromRef(unittest.TestCase):

    def testSetNuclSeq(self):
        region = GenomeRegion('test', [1, 8])
        region.set_seq_from_ref('ATGCGACGACGACGACGTTAGCAGTCGATCATCATGCTGATC', 'nucl')
        expected = 'ATGCGACG'
        result = region.nt_seq
        self.assertEqual(expected, result)

    def testSetProtSeq(self):
        region = GenomeRegion('test', None, None, [1, 9])
        region.set_seq_from_ref('ATVLHGLKLKAMAIMTIM', 'prot')
        expected = 'ATVLHGLKL'
        result = region.aa_seq
        self.assertEqual(expected, result)


class TestSetSequence(unittest.TestCase):

    def testNuclSequence(self):
        region = GenomeRegion('test')
        region.set_sequence('ATGCGACGACGACGACGTTAGCAGTCGATCATCATGCTGATC', 'nucl')
        expected = 'ATGCGACGACGACGACGTTAGCAGTCGATCATCATGCTGATC'
        result = region.nt_seq
        self.assertEqual(expected, result)

    def testProtSequence(self):
        region = GenomeRegion('test')
        region.set_sequence('ATVLHGLKLKAMAIMTIM', 'prot')
        expected = 'ATVLHGLKLKAMAIMTIM'
        result = region.aa_seq
        self.assertEqual(expected, result)


class TestSetPosFromCDS(unittest.TestCase):

    def testHIVLTR5Start(self):
        region = GenomeRegion('5\'LTR')
        region.set_pos_from_cds('hiv')
        expected = ['N/A']
        result = region.pos_from_cds
        self.assertEqual(expected, result)

    def testFromSIVLTR5Start(self):
        region = GenomeRegion('5\'LTR')
        region.set_pos_from_cds('siv')
        expected = ['N/A']
        result = region.pos_from_cds
        self.assertEqual(expected, result)

    def testFromHIVStart(self):
        region = GenomeRegion('Gag', [790, 2292])
        region.set_pos_from_cds('hiv')
        expected = [1, 1503]
        result = region.pos_from_cds
        self.assertEqual(expected, result)

    def testFromSIVStart(self):
        region = GenomeRegion('Integrase', [4785, 5666])
        region.set_pos_from_cds('siv')
        expected = [3477, 4358]
        result = region.pos_from_cds
        self.assertEqual(expected, result)


class TestSetPosFromAAStart(unittest.TestCase):

    def testCapsidFromAAStart(self):
        region = GenomeRegion('Capsid', [1186, 1878], None, None,
                              'PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEW'
                              'DRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEP'
                              'FRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVL')
        region.set_pos_from_aa_start('hiv')
        expected = [1, 231]
        result = region.pos_from_aa_start
        self.assertEqual(expected, result)

    def testLTR5FromAAStart(self):
        region = GenomeRegion('5\'LTR')
        region.set_pos_from_aa_start('hiv')
        expected = None
        result = region.pos_from_aa_start
        self.assertEqual(expected, result)

    def testGagFromAAStart(self):
        region = GenomeRegion('Gag', [1186, 1878], None, None,
                              'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYN'
                              'TVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVE'
                              'EKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTT'
                              'STLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTET'
                              'LLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGH'
                              'TARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPID'
                              'KELYPLTSLRSLFGNDPSSQ')
        region.set_pos_from_aa_start('hiv')
        expected = [133, 363]
        result = region.pos_from_aa_start
        self.assertEqual(expected, result)


class TestMakeCodonAln(unittest.TestCase):

    def testSimpleUse(self):
        region = GenomeRegion('p6', [2134, 2298], 'CTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAAC'
                                                  'TCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTG'
                                                  'GCAACGACCCCTCGTCACAATAA',
                                    [449, 500],   'LQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ')
        region.make_codon_aln()
        expected = '-L--Q--S--R--P--E--P--T--A--P--P--E--E--S--F--R--S--G--V--E--T--T--T--P--P--Q-' \
                   '-K--Q--E--P--I--D--K--E--L--Y--P--L--T--S--L--R--S--L--F--G--N--D--P--S--S--Q--*-'
        result = region.codon_aln
        self.assertEqual(expected, result)

        # Check that codon_aln aligns with the nucleotide sequence
        len_nt_seq = len(region.nt_seq)
        len_codon_aln = len(region.codon_aln)
        self.assertEqual(len_nt_seq, len_codon_aln)


class TestGlobalToLocalIndex(unittest.TestCase):

    def testNuclCoords(self):
        region = GenomeRegion('p1', [2086, 2133], 'TTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTT',
                                    [433, 448],   'FLGKIWPSYKGRPGNF')
        expected = [1, 48]
        result = region.global_to_local_index([2086, 2133], 'nucl')
        self.assertEqual(expected, result)

    def testProtCoords(self):
        region = GenomeRegion('p2', [1879, 1920], 'GCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATG',
                                    [364, 377],   'AEAMSQVTNSATIM')
        expected = [1, 14]
        result = region.global_to_local_index([364, 377], 'prot')
        self.assertEqual(expected, result)


class TestLocalToGlobalIndex(unittest.TestCase):

    def testNuclCoords(self):
        region = GenomeRegion('Nucleocapsid', [1921, 2085], 'ATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGTTTCAATTGTGGCA'
                                                            'AAGAAGGGCACACAGCCAGAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATG'
                                                            'TGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAGGCTAAT',
                                              [1, 55],      'MQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQAN')
        expected = [1921, 2085]
        result = region.local_to_global_index([1, 165], 'nucl')
        self.assertEqual(expected, result)

    def testProtCoords(self):
        region = GenomeRegion('RNase', [3870, 4229], 'TATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTAATAGAGG'
                                                     'AAGACAAAAAGTTGTCACCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATC'
                                                     'TAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATC'
                                                     'ATTCAAGCACAACCAGATCAAAGTGAATCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAAA'
                                                     'AAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAG'
                                                     'ATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTA',
                                       [1096, 1215], 'YVDGAANRETKLGKAGYVTNRGRQKVVTLTDTTNQKTELQAIYLALQDSGLEVNIVTDSQYALGI'
                                                     'IQAQPDQSESELVNQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKVL')
        expected = [1096, 1215]
        result = region.local_to_global_index([1, 120], 'prot')
        self.assertEqual(expected, result)


class TestGetOverlap(unittest.TestCase):

    def testNuclOverlap(self):
        region = GenomeRegion('gp41', [7758, 8795], 'GCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCC'
                                                    'TCAATGACGCTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAGCAGCAGAACAATTTG'
                                                    'CTGAGGGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATCAAGCAGCTCCAG'
                                                    'GCAAGAATCCTGGCTGTGGAAAGATACCTAAAGGATCAACAGCTCCTGGGGATTTGGGGTTGCTCT'
                                                    'GGAAAACTCATTTGCACCACTGCTGTGCCTTGGAATGCTAGTTGGAGTAATAAATCTCTGGAACAG'
                                                    'ATTTGGAATCACACGACCTGGATGGAGTGGGACAGAGAAATTAACAATTACACAAGCTTAATACAC'
                                                    'TCCTTAATTGAAGAATCGCAAAACCAGCAAGAAAAGAATGAACAAGAATTATTGGAATTAGATAAA'
                                                    'TGGGCAAGTTTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAAATTATTCATAATG'
                                                    'ATAGTAGGAGGCTTGGTAGGTTTAAGAATAGTTTTTGCTGTACTTTCTATAGTGAATAGAGTTAGG'
                                                    'CAGGGATATTCACCATTATCGTTTCAGACCCACCTCCCAACCCCGAGGGGACCCGACAGGCCCGAA'
                                                    'GGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGATTAGTGAACGGATCCTTG'
                                                    'GCACTTATCTGGGACGATCTGCGGAGCCTGTGCCTCTTCAGCTACCACCGCTTGAGAGACTTACTC'
                                                    'TTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAGCCCTCAAATATTGGTGG'
                                                    'AATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAGCTTGCTCAATGCCACAGCC'
                                                    'ATAGCAGTAGCTGAGGGGACAGATAGGGTTATAGAAGTAGTACAAGGAGCTTGTAGAGCTATTCGC'
                                                    'CACATACCTAGAAGAATAAGACAGGGCTTGGAAAGGATTTTGCTATAA',
                                      [2602, 2946], 'AVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQ'
                                                    'ARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSNKSLEQIWNHTTWMEWDREINNYTSLIH'
                                                    'SLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKLFIMIVGGLVGLRIVFAVLSIVNRVR'
                                                    'QGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLL'
                                                    'LIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQGACRAIR'
                                                    'HIPRRIRQGLERILL')
        expected = ('CAGATAGGGTTATAGAAGTAGTACAAGGAGCTTGTAGAGCTATTCGCCACATACCTAGAAGAATAAGACAGGGCTTGGAAAGGATTTTGCTATAA',
                    [8700, 8795])
        result = region.get_overlap([8700, 8795], 'nucl')
        self.assertEqual(expected, result)

    def testProtOverlap(self):
        region = GenomeRegion('Matrix', [790, 1185], 'ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAG'
                                                     'GCCAGGGGGAAAGAAAAAATATAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGAT'
                                                     'TCGCAGTTAATCCTGGCCTGTTAGAAACATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAA'
                                                     'CCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACCCTCTATTG'
                                                     'TGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCAAA'
                                                     'ACAAAAGTAAGAAAAAAGCACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAA'
                                                     'AATTAC',
                                        [1, 132], 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSL'
                                                  'QTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY')
        expected = ('MGARASVLSGGELDR', [1, 15])
        result = region.get_overlap([1, 15], 'prot')
        self.assertEqual(expected, result)


class TestSetRegions(unittest.TestCase):

    def testSIVInputCoords(self):
        region_names = ['Rev(with intron)', 'Rev(exon1)', 'Rev(exon2)', 'gp160', 'gp120', 'gp41', 'Nef', '3\'LTR']
        region_coordinates = [[6784, 9315], [6784, 6853], [9062, 9315], [6860, 9499],
                              [6860, 8434], [8435, 9499], [9333, 10124], [9719, 10535]]
        result = set_regions('siv', TEST_SIV_GENOME, TEST_SIV_NT_COORDS, TEST_SIV_PROTS, TEST_SIV_AA_COORDS)
        for i in range(len(result)):
            self.assertEqual(region_names[i], result[i].region_name)
            self.assertEqual(region_coordinates[i], result[i].nt_coords)


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
        handle = StringIO(">query\natgcgcg\n")
        result = get_query("nucl", handle, 'n')
        self.assertEqual(expected, result)

    def testProteinQuery(self):
        expected = [["query", "MPPLMMADLADLGG"]]
        handle = StringIO(">query\nMPPLMMADLADLGG\n")
        result = get_query("prot", handle, 'n')
        self.assertEqual(expected, result)

    def testLongNucleotideSequence(self):
        expected = [["query", "ATGCGCGAATTAGCGA"]]
        handle = StringIO(">query\natgcgcg\naattagcga\n")
        result = get_query("nucl", handle, 'n')
        self.assertEqual(expected, result)

    def testRevCompNucl(self):
        handle = StringIO(">query\nTCGCTAATTCGCGCATN*")
        expected = [["query", "*NATGCGCGAATTAGCGA"]]
        result = get_query("nucl", handle, 'y')
        self.assertEqual(expected, result)

    def testInvalidNucleotideQuery(self):
        handle = StringIO(">query\natgcgcg&\n")
        with self.assertRaises(SystemExit) as e:
            get_query("nucl", handle, 'n')
        self.assertEqual(e.exception.code, 0)

    def testInvalidProteinQuery(self):
        handle = StringIO(">query\nMPPLMMAD>LADLGG\n")
        with self.assertRaises(SystemExit) as e:
            get_query("prot", handle, 'n')
        self.assertEqual(e.exception.code, 0)

    def testRevCompProt(self):
        handle = StringIO(">query\nMPPLMMADLADLGG\n")
        expected = [["query", "MPPLMMADLADLGG"]]
        result = get_query('prot', handle, 'y')
        self.assertEqual(expected, result)

    def testPlainText(self):
        handle = StringIO("atgatcg\n")
        expected = [["Sequence1", "ATGATCG"]]
        result = get_query("nucl", handle, 'n')
        self.assertEqual(expected, result)

    def testMultipleQueries(self):
        handle = StringIO("atgct--agc\natgca---ga\n")
        expected = [["Sequence1", "ATGCT--AGC"], ["Sequence2", "ATGCA---GA"]]
        result = get_query("nucl", handle, 'n')
        self.assertEqual(expected, result)

    def testMultipleFasta(self):
        handle = StringIO(">q1\natgct--agc\n>q2\natgca---ga\n")
        expected = [["q1", "ATGCT--AGC"], ["q2", "ATGCA---GA"]]
        result = get_query("nucl", handle, 'n')
        self.assertEqual(expected, result)


class TestReverseComp(unittest.TestCase):

    def testSimpleUse(self):
        expected = "TCGCTAATTCGCGCATN*"
        result = reverse_comp("*NATGCGCGAATTAGCGA")
        self.assertEqual(expected, result)


class TestGetReferenceNucl(unittest.TestCase):

    def testDefaultHIVGenome(self):
        expected = [['K03455|HIVHXB2CG',
                     'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACC'
                     'AGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAAC'
                     'ACCAGCTTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACA'
                     'TGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG']]
        res = get_ref_nt_seq(DEFAULT_HIV_GENOME)
        # Check the first 350 nucleotides in the default reference sequence
        result = [[res[0][0], res[0][1][:350]]]
        self.assertEqual(expected, result)

    def testDefaultSIVGenome(self):
        expected = [['M33262|SIVMM239',
                     'GCATGCACATTTTAAAGGCTTTTGCTAAATATAGCCAAAAGTCCTTCTACAAATTTTCTAAGAGTTCTGATTCAAAGCAGTAACAGGCCTTGTCTCA'
                     'TCATGAACTTTGGCATTTCATCTACAGCTAAGTTTATATCATAAATAGTTCTTTACAGGCAGCACCAACTTATACCCTTATAGCATACTTTACTGTG'
                     'TGAAAATTGCATCTTTCATTAAGCTTACTGTAAATTTACTGGCTGTCTTCCTTGCAGGTTTCTGGAAGGGATTTATTACAGTGCAAGAAGACATAGA'
                     'ATCTTAGACATATACTTAGAAAAGGAAGAAGGCATCATACCAGATTGGCAGGATTACACCTCAGGACCAGGAATTAGATACCCAAAGACATTTGGCT'
                     'GGCTATGGAAAT']]
        res = get_ref_nt_seq(DEFAULT_SIV_GENOME)
        # Check the first 400 nucleotides in the default reference sequence
        result = [[res[0][0], res[0][1][:400]]]
        self.assertEqual(expected, result)


class TestGetReferenceProt(unittest.TestCase):

    def testInputHIVProteins(self):
        expected = [['test|HIVHXB2CG',
                     'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDT'
                     'KEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGH'
                     'QAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKE'
                     'PFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRK'
                     'IVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYP'
                     'LTSLRSLFGNDPSSQ']]
        result = get_ref_aa_seq(TEST_HIV_PROTS)
        self.assertEqual(expected, result)

    def testInputSIVProtein(self):
        expected = [['test|SIVMM239',
                     'MSNHEREEELRKRLRLIHLLHQTNPYPTGPGTANQRRQRKRRWRRRWQQLLALADRIYSFPDPPTDTPLDLAIQQLQNLAIESIPDPPTNTPEALCD'
                     'PTEDSRSPQDMGCLGNQLLIAILLLSVYGIYCTLYVTVFYGVPAWRNATIPLFCATKNRDTWGTTQCLPDNGDYSEVALNVTESFDAWNNTVTEQAI'
                     'EDVWQLFETSIKPCVKLSPLCITMRCNKSETDRWGLTKSITTTASTTSTTASAKVDMVNETSSCIAQDNCTGLEQEQMISCKFNMTGLKRDKKKEYN'
                     'ETWYSADLVCEQGNNTGNESRCYMNHCNTSVIQESCDKHYWDAIRFRYCAPPGYALLRCNDTNYSGFMPKCSKVVVSSCTRMMETQTSTWFGFNGTR'
                     'AENRTYIYWHGRDNRTIISLNKYYNLTMKCRRPGNKTVLPVTIMSGLVFHSQPINDRPKQAWCWFGGKWKDAIKEVKQTIVKHPRYTGTNNTDKINL'
                     'TAPGGGDPEVTFMWTNCRGEFLYCKMNWFLNWVEDRNTANQKPKEQHKRNYVPCHIRQIINTWHKVGKNVYLPPREGDLTCNSTVTSLIANIDWIDG'
                     'NQTNITMSAEVAELYRLELGDYKLVEITPIGLAPTDVKRYTTGGTSRNKRGVFVLGFLGFLATAGSAMGAASLTLTAQSRTLLAGIVQQQQQLLDVV'
                     'KRQQELLRLTVWGTKNLQTRVTAIEKYLKDQAQLNAWGCAFRQVCHTTVPWPNASLTPKWNNETWQEWERKVDFLEENITALLEEAQIQQEKNMYEL'
                     'QKLNSWDVFGNWFDLASWIKYIQYGVYIVVGVILLRIVIYIVQMLAKLRQGYRPVFSSPPSYFQQTHIQQDPALPTREGKERDGGEGGGNSSWPWQI'
                     'EYIHFLIRQLIRLLTWLFSNCRTLLSRVYQILQPILQRLSATLQRIREVLRTELTYLQYGWSYFHEAVQAVWRSATETLAGAWGDLWETLRRGGRWI'
                     'LAIPRRIRQGLELTLLMGGAISMRRSRPSGDLRQRLLRARGETYGRLLGEVEDGYSQSPGGLDKGLSSLSCEGQKYNQGQYMNTPWRNPAEEREKLA'
                     'YRKQNMDDIDEEDDDLVGVSVRPKVPLRTMSYKLAIDMSHFIKEKGGLEGIYYSARRHRILDIYLEKEEGIIPDWQDYTSGPGIRYPKTFGWLWKLV'
                     'PVNVSDEAQEDEEHYLMHPAQTSQWDDPWGEVLAWKFDPTLAYTYEAYVRYPEEFGSKSGLSEEEVRRRLTARGLLNMADKKETR']]
        result = get_ref_aa_seq(TEST_SIV_PROTS)
        self.assertEqual(expected, result)


class TestSequenceAlign(unittest.TestCase):

    def testNuclAlign(self):
        expected = \
            ['query',
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------------------------------------------------------------------------------'
             '---------------------------------CTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCC'
             'CCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA']

        query = [['query', 'CTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGA'
                           'AGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA']]

        reference_sequence = \
            [['>K03455|HIVHXB2CG',
              'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCC'
              'AGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGTTA'
              'CACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATC'
              'CGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGC'
              'CCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACT'
              'GCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTG'
              'TGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGCACGGC'
              'AAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGA'
              'ATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTT'
              'AATCCTGGCCTGTTAGAAACATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATA'
              'CAGTAGCAACCCTCTATTGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAAGC'
              'ACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAAAATTACCCTATAGTGCAGAACATCCAGGGGCAAATGGTACATCAGGCCATATCACCT'
              'AGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTGATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATT'
              'TAAACACCATGCTAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATGTTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGAGTGCATCCAGT'
              'GCATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACAAAT'
              'AATCCACCTATCCCAGTAGGAGAAATTTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGAC'
              'AAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGGTTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAGGAGGTAAAAAATTGGATGACAGAAACCTT'
              'GTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCATTGGGACCAGCGGCTACACTAGAAGAAATGATGACAGCATGTCAGGGAGTAGGAGGA'
              'CCCGGCCATAAGGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATGATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTG'
              'TTAAGTGTTTCAATTGTGGCAAAGAAGGGCACACAGCCAGAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAA'
              'AGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCA'
              'CCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCT'
              'TTGGCAACGACCCCTCGTCACAATAA']]
        # Aligned query sequence
        result = sequence_align(query, reference_sequence, None)[1]
        self.assertEqual(expected, result)

    def testProtAlign(self):

        query = [['query', 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVA'
                           'TLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY']]

        ref_seq = [['reference',
                    'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTK'
                    'EALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQA'
                    'AMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFR'
                    'DYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKC'
                    'FNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLR'
                    'SLFGNDPSSQ']]

        expected = ['query',
                    'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTK'
                    'EALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY----------------------------------------------------------------'
                    '--------------------------------------------------------------------------------------------------'
                    '--------------------------------------------------------------------------------------------------'
                    '--------------------------------------------------------------------------------------------------'
                    '----------']

        result = sequence_align(query, ref_seq, None)[1]
        self.assertEqual(expected, result)


class TestGetRegionCoordinates(unittest.TestCase):
    """
    Region coordinates use inclusive 0-based indexing
    """
    def testNuclAlignment(self):
        aln = ['query',
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------CTTCAGAGCAGACCAGAGCCAACAGCCCCA'
               'CCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCT'
               'TTGGCAACGACCCCTCGTCACAATAA']
        expected = [(2133, 2292)]
        result = get_region_coordinates(aln)
        self.assertEqual(expected, result)

    def testProtAlignment(self):
        aln = ['query',
               'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDK'
               'IEEEQNKSKKKAQQAAADTGHSNQVSQNY--------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '----------------------------------------------------------------------------------------']
        expected = [(0, 132)]
        result = get_region_coordinates(aln)
        self.assertEqual(expected, result)


class TestFindMatches(unittest.TestCase):

    def testNuclMatches(self):
        ref_regions = [
            GenomeRegion(
                '5\'LTR',
                [1, 634],
                'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGC'
                'CAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGT'
                'TACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGC'
                'ATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGC'
                'GAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAAC'
                'CCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAG'
                'TCAGTGTGGAAAATCTCTAGCA',
                None, None),
            GenomeRegion(
                'Gag',
                [790, 2292],
                'ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATATA'
                'GTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTT'
                'CAGACAGGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACCCTCTATTGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGAC'
                'AAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAAGCACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAAAATTACCCTATAGTGCAG'
                'AACATCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTGATACCC'
                'ATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACACCATGCTAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATGTTAAAAGAGACC'
                'ATCAATGAGGAAGCTGCAGAATGGGATAGAGTGCATCCAGTGCATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGA'
                'ACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACAAATAATCCACCTATCCCAGTAGGAGAAATTTATAAAAGATGGATAATCCTGGGATTAAATAAA'
                'ATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGGTTCTATAAAACTCTAAGAGCC'
                'GAGCAAGCTTCACAGGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCATTGGGACCAGCG'
                'GCTACACTAGAAGAAATGATGACAGCATGTCAGGGAGTAGGAGGACCCGGCCATAAGGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACAAATTCAGCT'
                'ACCATAATGATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACACAGCCAGAAATTGCAGGGCCCCT'
                'AGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCTACAAG'
                'GGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAG'
                'GAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA',
                [1, 500],
                'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALD'
                'KIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKET'
                'INEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRA'
                'EQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGHTARNCRAP'
                'RKKGCWKCGKEGHQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ'),
            GenomeRegion(
                'Matrix',
                [790, 1185],
                'ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATATA'
                'GTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTT'
                'CAGACAGGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACCCTCTATTGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGAC'
                'AAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAAGCACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAAAATTAC',
                [1, 132],
                'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALD'
                'KIEEEQNKSKKKAQQAAADTGHSNQVSQNY'),
            GenomeRegion(
                'Capsid',
                [1186, 1878],
                'CCTATAGTGCAGAACATCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCA'
                'GAAGTGATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACACCATGCTAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATG'
                'TTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGAGTGCATCCAGTGCATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACCAAGGGGAAGT'
                'GACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACAAATAATCCACCTATCCCAGTAGGAGAAATTTATAAAAGATGGATAATCCTG'
                'GGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGGTTCTATAAA'
                'ACTCTAAGAGCCGAGCAAGCTTCACAGGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCA'
                'TTGGGACCAGCGGCTACACTAGAAGAAATGATGACAGCATGTCAGGGAGTAGGAGGACCCGGCCATAAGGCAAGAGTTTTG',
                [133, 363],
                'PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGS'
                'DIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKA'
                'LGPAATLEEMMTACQGVGGPGHKARVL'),
            GenomeRegion(
                'p2', [1879, 1920], 'GCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATG', [364, 377], 'AEAMSQVTNSATIM'),
            GenomeRegion(
                'Nucleocapsid',
                [1921, 2085], 'ATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACACAGCCAGAAATTGCAGGG'
                              'CCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAGGCTAAT',
                [378, 342], 'MQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQAN'),
            GenomeRegion(
                'p1', [2086, 2133], 'TTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTT', [433, 448], 'FLGKIWPSYKGRPGNF'),
            GenomeRegion(
                'p6', [2134, 2292], 'CTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGC'
                                    'AGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA',
                [449, 500],   'LQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ')]

        coordinates = [(2133, 2292)]
        result = find_matches('hiv', 'nucl', ref_regions, coordinates)
        region_names = ['Gag', 'p6']
        for i in range(len(region_names)):
            self.assertEqual(region_names[i], result[i].region_name)


class TestRetrieve(unittest.TestCase):

    def setUp(self):
        self.hiv_nt_seq_file = open(TEST_HIV_GENOME)
        self.hiv_aa_seq_file = open(TEST_HIV_PROTS)
        self.hiv_nt_seq = convert_fasta(self.hiv_nt_seq_file)[0][1]
        self.hiv_aa_seq = convert_fasta(self.hiv_aa_seq_file)[0][1]

        self.hiv_nt_coords = open(TEST_HIV_NT_COORDS)
        self.hiv_aa_coords = open(TEST_HIV_AA_COORDS)

        self.siv_genome_file = open(TEST_SIV_GENOME)
        self.siv_prot_file = open(TEST_SIV_PROTS)

    def testDefaultInput(self):
        ref_regions = set_regions('hiv', self.hiv_nt_seq, TEST_HIV_NT_COORDS, self.hiv_aa_seq, TEST_HIV_AA_COORDS)
        for ref in ref_regions:
            print(ref.region_name)
            print(ref.nt_seq)
            print(ref.nt_coords)
            print(ref.aa_seq)
            print(ref.aa_coords)
            print()
        expected_seq = 'TTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTT'
        result = retrieve('hiv', 'nucl', ref_regions, 'p1')
        result_seq = result.nt_seq
        self.assertEqual(expected_seq, result_seq)

    def testSIVInput(self):
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
        self.hiv_nt_seq_file.close()
        self.hiv_aa_seq_file.close()
        self.hiv_nt_coords.close()
        self.hiv_aa_coords.close()

        self.siv_genome_file.close()
        self.siv_prot_file.close()


class TestHandleArgs(unittest.TestCase):

    maxDiff = None

    def setUp(self):

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

        expected_query = [
            ['query',
             'GACTCGAAAGCGAAAGTTCCAGAGAAGTTCTCTCGACGCAGGACTCGGCTTGCTGAGGTGCACACAGCAAGAGGCGAGAGCGGCGACTGGTGAGTACGCCAAATT'
             'TTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGTGGGGGAAAATTAGATGCATGGGAAAAAATTCGGTTACGGCCAGGGGGAAA'
             'GAAAAAATATAGAATGAAACATTTAGTATGGGCAAGCAGAGAGTTAGAAAGATTCGCACTTAACCC']]
        result_query = result[0]
        self.assertEqual(expected_query, result_query)

        expected_ref_nt_seq = [
            ['K03455|HIVHXB2CG',
             'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAG'
             'GGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGTTACACC'
             'CTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGT'
             'ACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG']]

        # Check the first 350 nucleotides in the default reference sequence
        seq = result[1][0][1][:350]
        header = result[1][0][0]
        result_ref_nt_seq = [[header, seq]]
        self.assertEqual(expected_ref_nt_seq, result_ref_nt_seq)

        expected_ref_aa_seq = [
            ['Gag|HIVHXB2CG',
             'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIE'
             'EEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAA'
             'EWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKN'
             'WMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEG'
             'HQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ'],
            ['Matrix|HIVHXB2CG',
             'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIE'
             'EEQNKSKKKAQQAAADTGHSNQVSQNY'],
            ['Capsid|HIVHXB2CG',
             'PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIA'
             'GTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAAT'
             'LEEMMTACQGVGGPGHKARVL'],
            ['p2|HIVHXB2CG', 'AEAMSQVTNSATIM']]

        # Check the first 4 proteins in the default reference protein sequence
        result_ref_aa_seq = result[2][:4]
        self.assertEqual(expected_ref_aa_seq, result_ref_aa_seq)

    def testDefaultSIVProt(self):
        """
        Tests the scenario when the user selects SIV and protein alignment
        """

        result = handle_args('siv', 'prot', self.nucl_query, revcomp='n', ref_nt=None, ref_aa=None)

        expected_query = [
            ['query', 'GACTCGAAAGCGAAAGTTCCAGAGAAGTTCTCTCGACGCAGGACTCGGCTTGCTGAGGTGCACACAGCAAGAGGCGAGAGCGGCGACTGGTGAGTA'
                      'CGCCAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGTGGGGGAAAATTAGATGCATGGGAAAAAATTCG'
                      'GTTACGGCCAGGGGGAAAGAAAAAATATAGAATGAAACATTTAGTATGGGCAAGCAGAGAGTTAGAAAGATTCGCACTTAACCC']]
        result_query = result[0]
        self.assertEqual(expected_query, result_query)

        expected_ref_nt_seq = [
            ['M33262|SIVMM239', 'GCATGCACATTTTAAAGGCTTTTGCTAAATATAGCCAAAAGTCCTTCTACAAATTTTCTAAGAGTTCTGATTCAAAGCAGTAACAG'
                                'GCCTTGTCTCATCATGAACTTTGGCATTTCATCTACAGCTAAGTTTATATCATAAATAGTTCTTTACAGGCAGCACCAACTTATAC'
                                'CCTTATAGCATACTTTACTGTGTGAAAATTGCATCTTTCATTAAGCTTACTGTAAATTTACTGGCTGTCTTCCTTGCAGGTTTCTG'
                                'GAAGGGATTTATTACAGTGCAAGAAGACATAGAATCTTAGACATATACTTAGAAAAGGAAGAAGGCATCATACCAGATTGGCAGGA'
                                'TTACACCTCAGGACCAGGAATTAGATACCCAAAGACATTTGGCTGGCTATGGAAAT']]

        # Check the first 400 nucleotides in the default reference sequence
        seq = result[1][0][1][:400]
        header = result[1][0][0]
        result_ref_nt_seq = [[header, seq]]
        self.assertEqual(expected_ref_nt_seq, result_ref_nt_seq)

        expected_ref_aa_seq = [
            ['Capsid|SIVMM239', 'PVQQIGGNYVHLPLSPRTLNAWVKLIEEKKFGAEVVPGFQALSEGCTPYDINQMLNCVGDHQAAMQIIRDIINEEAADWDLQHPQP'
                                'APQQGQLREPSGSDIAGTTSSVDEQIQWMYRQQNPIPVGNIYRRWIQLGLQKCVRMYNPTNILDVKQGPKEPFQSYVDRFYKSLRA'
                                'EQTDAAVKNWMTQTLLIQNANPDCKLVLKGLGVNPTLEEMLTACQGVGGPGQKARLM'],
            ['p2|SIVMM239', 'AEALKEALAPVPIPFAA'],
            ['Nucleocapsid|SIVMM239', 'AQQRGPRKPIKCWNCGKEGHSARQCRAPRRQGCWKCGKMDHVMAKCPDRQAG']]
        result_ref_aa_seq = result[2][2:5]
        self.assertEqual(expected_ref_aa_seq, result_ref_aa_seq)

    def testDefaultRevCompHIVProt(self):
        """
        Tests the scenario when the user selects HIV and nucleotide alignment with the reverse complement of the query
        """

        result = handle_args('hiv', 'prot', self.nucl_query, revcomp='y', ref_nt=None, ref_aa=None)

        expected_query = [
            ['query', 'GACTCGAAAGCGAAAGTTCCAGAGAAGTTCTCTCGACGCAGGACTCGGCTTGCTGAGGTGCACACAGCAAGAGGCGAGAGCGGCGACTGGTGAGTA'
                      'CGCCAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGTGGGGGAAAATTAGATGCATGGGAAAAAATTCG'
                      'GTTACGGCCAGGGGGAAAGAAAAAATATAGAATGAAACATTTAGTATGGGCAAGCAGAGAGTTAGAAAGATTCGCACTTAACCC']]
        result_query = result[0]
        self.assertEqual(expected_query, result_query)

        expected_ref_nt_seq = [
            ['K03455|HIVHXB2CG',
             'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAG'
             'GGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGTTACACC'
             'CTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGT'
             'ACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG']]

        # Check the first 350 nucleotides in the default reference sequence
        seq = result[1][0][1][:350]
        header = result[1][0][0]
        result_ref_nt_seq = [[header, seq]]
        self.assertEqual(expected_ref_nt_seq, result_ref_nt_seq)

        expected_ref_aa_seq = [
            ['Gag|HIVHXB2CG',
             'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIE'
             'EEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAA'
             'EWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKN'
             'WMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEG'
             'HQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ'],
            ['Matrix|HIVHXB2CG', 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRS'
                                 'LYNTVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY'],
            ['Capsid|HIVHXB2CG', 'PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHP'
                                 'VHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYK'
                                 'TLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVL'],
            ['p2|HIVHXB2CG', 'AEAMSQVTNSATIM']]

        # Check the first 4 proteins in the default reference protein sequence
        result_ref_aa_seq = result[2][:4]
        self.assertEqual(expected_ref_aa_seq, result_ref_aa_seq)

    def testDefaultRevCompSIVNucl(self):
        """
        Tests the scenario when the user selects SIV and nucleotide alignment with the reverse complement of the query
        """
        result = handle_args('siv', 'nucl', self.nucl_query, revcomp='y', ref_nt=None, ref_aa=None)

        expected_query = [
            ['query', 'GGGTTAAGTGCGAATCTTTCTAACTCTCTGCTTGCCCATACTAAATGTTTCATTCTATATTTTTTCTTTCCCCCTGGCCGTAACCGAATTTTTTCC'
                      'CATGCATCTAATTTTCCCCCACTTAATACTGACGCTCTCGCACCCATCTCTCTCCTTCTAGCCTCCGCTAGTCAAAATTTGGCGTACTCACCAGTC'
                      'GCCGCTCTCGCCTCTTGCTGTGTGCACCTCAGCAAGCCGAGTCCTGCGTCGAGAGAACTTCTCTGGAACTTTCGCTTTCGAGTC']]
        result_query = result[0]
        self.assertEqual(expected_query, result_query)

        expected_ref_nt_seq = [
            ['M33262|SIVMM239', 'GCATGCACATTTTAAAGGCTTTTGCTAAATATAGCCAAAAGTCCTTCTACAAATTTTCTAAGAGTTCTGATTCAAAGCAGTAACAG'
                                'GCCTTGTCTCATCATGAACTTTGGCATTTCATCTACAGCTAAGTTTATATCATAAATAGTTCTTTACAGGCAGCACCAACTTATAC'
                                'CCTTATAGCATACTTTACTGTGTGAAAATTGCATCTTTCATTAAGCTTACTGTAAATTTACTGGCTGTCTTCCTTGCAGGTTTCTG'
                                'GAAGGGATTTATTACAGTGCAAGAAGACATAGAATCTTAGACATATACTTAGAAAAGGAAGAAGGCATCATACCAGATTGGCAGGA'
                                'TTACACCTCAGGACCAGGAATTAGATACCCAAAGACATTTGGCTGGCTATGGAAAT']]

        # Check the first 400 nucleotides in the default reference sequence
        seq = result[1][0][1][:400]
        header = result[1][0][0]
        result_ref_nt_seq = [[header, seq]]
        self.assertEqual(expected_ref_nt_seq, result_ref_nt_seq)

        expected_ref_aa_seq = [
            ['Capsid|SIVMM239', 'PVQQIGGNYVHLPLSPRTLNAWVKLIEEKKFGAEVVPGFQALSEGCTPYDINQMLNCVGDHQAAMQIIRDIINEEAADWDLQHPQP'
                                'APQQGQLREPSGSDIAGTTSSVDEQIQWMYRQQNPIPVGNIYRRWIQLGLQKCVRMYNPTNILDVKQGPKEPFQSYVDRFYKSLRA'
                                'EQTDAAVKNWMTQTLLIQNANPDCKLVLKGLGVNPTLEEMLTACQGVGGPGQKARLM'],
            ['p2|SIVMM239', 'AEALKEALAPVPIPFAA'],
            ['Nucleocapsid|SIVMM239', 'AQQRGPRKPIKCWNCGKEGHSARQCRAPRRQGCWKCGKMDHVMAKCPDRQAG']]
        result_ref_aa_seq = result[2][2:5]
        self.assertEqual(expected_ref_aa_seq, result_ref_aa_seq)

    def testInputNuclProtHIV(self):
        """
        Tests the scenario when the user selects HIV, specifies the reference sequences,
        and selects nucleotide alignment
        """
        result = handle_args('hiv', 'nucl', self.nucl_query, revcomp='n',
                             ref_nt=TEST_HIV_GENOME, ref_aa=TEST_HIV_PROTS)

        expected_query = [
            ['query', 'GACTCGAAAGCGAAAGTTCCAGAGAAGTTCTCTCGACGCAGGACTCGGCTTGCTGAGGTGCACACAGCAAGAGGCGAGAGCGGCGACTGGTGAGTA'
                      'CGCCAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGTGGGGGAAAATTAGATGCATGGGAAAAAATTCG'
                      'GTTACGGCCAGGGGGAAAGAAAAAATATAGAATGAAACATTTAGTATGGGCAAGCAGAGAGTTAGAAAGATTCGCACTTAACCC']]
        result_query = result[0]
        self.assertEqual(expected_query, result_query)

        expected_ref_nt = [
            ['K03455|HIVHXB2CG',
             'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAG'
             'GGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGTTACACC'
             'CTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGT'
             'ACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG']]
        result_ref_nt = result[1]
        self.assertEqual(expected_ref_nt, result_ref_nt)

        expected_ref_prot = [
            ['Gag|HIVHXB2CG', 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCV'
                              'HQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALS'
                              'EGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKR'
                              'WIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQG'
                              'VGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPS'
                              'YKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ'],
            ['Matrix|HIVHXB2CG', 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATL'
                                 'YCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY'],
            ['Capsid|HIVHXB2CG', 'PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHP'
                                 'VHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYK'
                                 'TLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVL']]
        result_ref_prot = result[2]
        self.assertEqual(expected_ref_prot, result_ref_prot)

    def testInputAllSIV(self):
        """
        Tests the scenario when the user selects SIV, specifies the reference sequences, and selects protein alignment
        """
        result = handle_args('siv', 'prot', self.prot_query, revcomp='n',
                             ref_nt=TEST_SIV_GENOME, ref_aa=TEST_SIV_PROTS)

        expected_query = [
            ['query', 'MAYTQTATTSALLDTVRGNNSLVNDLAKRRLYDTAVEEFNARDRRPKVNFSKVISEEQTLIATRAYPEFQITFYNTQNAVHSLAGGLRSLELEYLM'
                      'MQIPYGSLTYDIGGNFASHLFKGRAYVHCCMPNLDVRDIMRHEGQKDSIELYLSRLERGGKTVPNFQKEAFDRYAEIPEDAVCHNTFQTMRHQPMQ'
                      'QSGRVYAIALHSIYDIPADEFGAALLRKNVHTCYAAFHFSENLLLEDSYVNLDEINACFSRDGDKLTFSFASESTLNYCHSYSNILKYVCKTYFPA'
                      'SNREVYMKEFLV']]
        result_query = result[0]
        self.assertEqual(expected_query, result_query)

        expected_ref_nt = [
            ['M33262|SIVMM239', 'GCATGCACATTTTAAAGGCTTTTGCTAAATATAGCCAAAAGTCCTTCTACAAATTTTCTAAGAGTTCTGATTCAAAGCAGTAACAG'
                                'GCCTTGTCTCATCATGAACTTTGGCATTTCATCTACAGCTAAGTTTATATCATAAATAGTTCTTTACAGGCAGCACCAACTTATAC'
                                'CCTTATAGCATACTTTACTGTGTGAAAATTGCATCTTTCATTAAGCTTACTGTAAATTTACTGGCTGTCTTCCTTGCAGGTTTCTG'
                                'GAAGGGATTTATTACAGTGCAAGAAGACATAGAATCTTAGACATATACTTAGAAAAGGAAGAAGGCATCATACCAGATTGGCAGGA'
                                'TTACACCTCAGGACCAGGAATTAGATACCCAAAGACATTTGGCTGGCTATGGAAAT']]
        result_ref_nt = result[1]
        self.assertEqual(expected_ref_nt, result_ref_nt)

        expected_ref_prot = [
            ['RNase|SIVMM239', 'YTDGSCNKQSKEGKAGYITDRGKDKVKVLEQTTNQQAELEAFLMALTDSGPKANIIVDSQYVMGIITGCPTESESRLVNQIIEEMIK'
                               'KSEIYVAWVPAHKGIGGNQEIDHLVSQGIRQVL'],
            ['Integrase|SIVMM239', 'FLEKIEPAQEEHDKYHSNVKELVFKFGLPRIVARQIVDTCDKCHQKGEAIHGQANSDLGTWQMDCTHLEGKIIIVAVHVASGF'
                                   'IEAEVIPQETGRQTALFLLKLAGRWPITHLHTDNGANFASQEVKMVAWWAGIEHTFGVPYNPQSQGVVEAMNHHLKNQIDRIR'
                                   'EQANSVETIVLMAVHCMNFKRRGGIGDMTPAERLINMITTEQEIQFQQSKNSKFKNFRVYYREGRDQLWKGPGELLWKGEGAV'
                                   'ILKVGTDIKVVPRRKAKIIKDYGGGKEVDSSSHMEDTGEAREVA'],
            ['Vif|SIVMM239', 'MEEEKRWIAVPTWRIPERLERWHSLIKYLKYKTKDLQKVCYVPHFKVGWAWWTCSRVIFPLQEGSHLEVQGYWHLTPEKGWLSTYAVRI'
                             'TWYSKNFWTDVTPNYADILLHSTYFPCFTAGEVRRAIRGEQLLSCCRFPRAHKYQVPSLQYLALKVVSDVRSQGENPTWKQWRRDNRRG'
                             'LRMAKQNSRGDKQRGGKPPTKGANFPGLAKVLGILA'],
            ['Vpx|SIVMM239', 'MSDPRERIPPGNSGEETIGEAFEWLNRTVEEINREAVNHLPRELIFQVWQRSWEYWHDEQGMSPSYVKYRYLCLIQKALFMHCKKGCRC'
                             'LGEGHGAGGWRPGPPPPPPPGLA'],
            ['Vpr|SIVMM239', 'MEERPPENEGPQREPWDEWVVEVLEELKEEALKHFDPRLLTALGNHIYNRHGDTLEGAGELIRILQRALFMHFRGGCIHSRIGQPGGGN'
                             'PLSAIPPSRSML']]
        result_ref_prot = result[2]
        self.assertEqual(expected_ref_prot, result_ref_prot)

    def tearDown(self):

        self.hiv_default_genome.close()
        self.hiv_default_prot.close()
        self.hiv_test_genome.close()
        self.hiv_test_prot.close()

        self.siv_default_genome.close()
        self.siv_default_prot.close()
        self.siv_test_genome.close()
        self.siv_test_prot.close()
