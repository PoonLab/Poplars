import unittest

from sequence_locator import *

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

        self.hiv_genome = Genome(self.nt_coords_handle, self.ref_nt_seq,
                             self.aa_coords_handle, self.ref_aa_seq, self.ref_seq)

        self.p2_ref = RefRegion('p2', self.hiv_genome, [1879, 1920], [364, 377])
        self.p6_ref = RefRegion('p6', self.hiv_genome, [2134, 2292], [449, 500])

    def tearDown(self):
        self.nt_coords_handle.close()
        self.aa_coords_handle.close()


class SIV(unittest.TestCase):

    def setUp(self):
        self.siv_configs = handle_args('siv', 'nucl')
        self.ref_nt_seq, self.ref_aa_seq = self.siv_configs[0][0][1], self.siv_configs[1]
        nt_coords, aa_coords = self.siv_configs[2], self.siv_configs[3]
        self.ref_seq = self.siv_configs[4]

        self.nt_coords_handle = open(nt_coords)
        self.aa_coords_handle = open(aa_coords)
        self.siv_genome = Genome(self.nt_coords, self.ref_nt_seq,
                             self.aa_coords, self.ref_aa_seq, self.ref_seq)

    def tearDown(self):
        self.nt_coords.close()
        self.aa_coords.close()


####################################
# Test cases for Region
####################################
class TestRegion(HIV):

    def testCoords(self):
        region = Region('Gag', self.hiv_genome)

        region.set_coords([790, 2292], 'nucl')
        expected_ncoords = [790, 2292]
        result_ncoords = region.get_coords('nucl')
        self.assertEqual(expected_ncoords, result_ncoords)

        region.set_coords([1, 500], 'prot')
        expected_pcoords = [1, 500]
        result_pccords = region.get_coords('prot')
        self.assertEqual(expected_pcoords, result_pccords)

    def testSequenceFromGenome(self):
        region = Region('Vpu', self.hiv_genome, [6062, 6310], [2009, 2090])

        region.set_nt_seq_from_genome()
        expected_nt = 'ACGCAACCTATACCAATAGTAGCAATAGTAGCATTAGTAGTAGCAATAATAATAGCAATAGTTGTGTGGTCCATAGTAAT' \
                      'CATAGAATATAGGAAAATATTAAGACAAAGAAAAATAGACAGGTTAATTGATAGACTAATAGAAAGAGCAGAAGACAGTG' \
                      'GCAATGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGCTCCTTGGGATGTTGATGATCTGTAG'
        result_nt = region.get_sequence_from_genome('nucl')
        self.assertEqual(expected_nt, result_nt)

        region.set_aa_seq_from_genome()
        expected_aa = 'TQPIPIVAIVALVVAIIIAIVVWSIVIIEYRKILRQRKIDRLIDRLIERAEDSGNESEGEISALVEMGVEMGHHAPWDVDDL'
        result_aa = region.get_sequence_from_genome('prot')
        self.assertEqual(expected_aa, result_aa)


####################################
# Test cases for RefRegion
####################################
class TestRefRegion(HIV):

    def testMakeCodonAln(self):
        self.p2_ref.set_nt_seq_from_genome()
        self.p2_ref.set_aa_seq_from_genome()
        expected = '-A--E--A--M--S--Q--V--T--N--S--A--T--I--M-'
        self.p2_ref.make_codon_aln()
        self.assertEqual(expected, self.p2_ref.codon_aln)

        self.p6_ref.set_nt_seq_from_genome()
        self.p6_ref.set_aa_seq_from_genome()
        expected = '-L--Q--S--R--P--E--P--T--A--P--P--E--E--S--F--R--S--G--V--E--T--T--T--P--P--Q--K-' \
                   '-Q--E--P--I--D--K--E--L--Y--P--L--T--S--L--R--S--L--F--G--N--D--P--S--S--Q--*-'
        self.p6_ref.make_codon_aln()
        self.assertEqual(expected, self.p6_ref.codon_aln)

    def testSetProteinEquivalent(self):
        self.p2_ref.set_nt_seq_from_genome()
        self.p2_ref.set_aa_seq_from_genome()
        self.p2_ref.make_codon_aln()
        expected_prot_seq = 'AEAMSQVTNSATIM'
        expected_prot_coords = [1, 14]

        # Local coordinates: 0, 42 (relative to the start of the p2 coding region)
        # Global coordinates: 1879, 1920 (relative to the start of the genome)
        overlap = QueryRegion('p2', self.p2_ref, 'nucl', self.hiv_genome, [0, 42], [1879, 1920])
        overlap.set_nt_seq_from_genome()
        overlap.set_aa_seq_from_genome()
        result = self.p2_ref.set_protein_equivalents(overlap)
        self.assertEqual(expected_prot_seq, result[0])
        self.assertEqual(expected_prot_coords, result[1])

        self.p6_ref.set_nt_seq_from_genome()
        self.p6_ref.set_aa_seq_from_genome()
        self.p6_ref.make_codon_aln()
        expected_prot_seq = 'QSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ'    # excludes the first and last amino acids
        expected_prot_coords = [2, 52]

        # Local coordinates: 1, 52 (relative to the start of the p6 coding region)
        # Global coordinates: 2131, 2289 (relative to the start of the genome)
        overlap = QueryRegion('p6', self.p6_ref, 'nucl', self.hiv_genome, [2, 155], [2131, 2289])
        result = self.p6_ref.set_protein_equivalents(overlap)
        self.assertEqual(expected_prot_seq, result[0])
        self.assertEqual(expected_prot_coords, result[1])


####################################
# Test cases for QueryRegion
####################################
class TestQueryRegion(HIV):

    def test_set_relative_positions(self):

        # Test p2
        self.p2_ref.set_nt_seq_from_genome()
        self.p2_ref.set_aa_seq_from_genome()
        self.p2_ref.codon_aln = self.p2_ref.make_codon_aln()

        p2_q = QueryRegion('p2', self.p2_ref, 'nucl', self.hiv_genome, [0, 42], [1879, 1920])
        p2_q.set_nt_seq_from_genome()
        p2_q.set_aa_seq_from_genome()
        self.p2_ref.set_protein_equivalents(p2_q)

        self.assertEqual([1, 42], p2_q.set_pos_from_cds())
        self.assertEqual([1879, 1920], p2_q.set_pos_from_gstart())
        self.assertEqual([1, 14], p2_q.set_pos_from_pstart())
        self.assertEqual([1, 42], p2_q.set_pos_from_qstart('nucl'))
        self.assertEqual([1, 42], p2_q.set_pos_from_rstart('nucl'))

        # Test p6
        self.p6_ref.set_nt_seq_from_genome()
        self.p6_ref.set_aa_seq_from_genome()
        self.p6_ref.codon_aln = self.p6_ref.make_codon_aln()

        p6_q = QueryRegion('p6', self.p6_ref, 'nucl', self.hiv_genome, [2, 155],  [2131, 2289])
        p6_q.set_nt_seq_from_genome()
        p6_q.set_aa_seq_from_genome()
        self.p6_ref.set_protein_equivalents(p6_q)

        self.assertEqual([2, 158], p6_q.set_pos_from_cds())
        self.assertEqual([2135, 2291], p6_q.set_pos_from_gstart())
        self.assertEqual([2, 52], p6_q.set_pos_from_pstart())
        self.assertEqual([2, 158], p6_q.set_pos_from_qstart('nucl'))
        self.assertEqual([2, 158], p6_q.set_pos_from_rstart('nucl'))


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

        expected_ref_aa_seq = 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYC'\
                              'VHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSA'\
                              'LSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEI'\
                              'YKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMT'\
                              'ACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGH'

        aa_seq = result[1]['Gag'][:400]
        self.assertEqual(expected_ref_aa_seq, aa_seq)

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

        expected_ref_prot = 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCV' \
                            'HQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY'
        result_ref_prot = result[1]['Matrix(p17/p15)']
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

        expected_ref_aa_seq = 'PVQQIGGNYVHLPLSPRTLNAWVKLIEEKKFGAEVVPGFQALSEGCTPYDINQMLNCVGDHQAAMQIIRDIINEEAADWDLQHPQPA'\
                              'PQQGQLREPSGSDIAGTTSSVDEQIQWMYRQQNPIPVGNIYRRWIQLGLQKCVRMYNPTNILDVKQGPKEPFQSYVDRFYKSLRAEQ'\
                              'TDAAVKNWMTQTLLIQNANPDCKLVLKGLGVNPTLEEMLTACQGVGGPGQKARLM'
        result_ref_aa_seq = result[1]['Capsid(p24/p27)']
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

        expected_ref_prot = 'MSDPRERIPPGNSGEETIGEAFEWLNRTVEEINREAVNHLPRELIFQVWQRSWEYWHDEQGMSPSYVKYRYLCLI' \
                            'QKALFMHCKKGCRCLGEGHGAGGWRPGPPPPPPPGLA'
        result_ref_prot = result[1]['Vpx']
        self.assertEqual(expected_ref_prot, result_ref_prot)


if __name__ == '__main__':
    unittest.main()
