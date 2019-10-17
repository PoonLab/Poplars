import unittest
from unittest.mock import patch
from poplars.sequence_locator import *


class HIV(unittest.TestCase):

    def setUp(self):
        self.hiv_nt_seq_path = os.path.join(os.path.dirname(__file__), 'fixtures/hiv-test-genome.fasta')
        self.hiv_aa_seq_path = os.path.join(os.path.dirname(__file__), 'fixtures/hiv-test-proteins.fasta')
        self.hiv_ncoords_path = os.path.join(os.path.dirname(__file__), 'fixtures/hiv_test_nt_coords.csv')
        self.hiv_pcoords_path = os.path.join(os.path.dirname(__file__), 'fixtures/hiv_test_prot_coords.csv')

        with open(self.hiv_nt_seq_path) as hiv_nt, open(self.hiv_aa_seq_path) as hiv_aa, \
                open(self.hiv_ncoords_path) as hiv_ncoords, open(self.hiv_pcoords_path) as hiv_pcoords:
            self.hiv_nt_seq = convert_fasta(hiv_nt.read().split())[0][1]
            self.hiv_aa_seq = convert_fasta(hiv_aa.read().split())[0][1]
            self.hiv_ncoords = hiv_ncoords
            self.hiv_pcoords = hiv_pcoords


class TestHIV(HIV):

    maxDiff = None

    def setUp(self):
        self.hiv_genome_regions = make_regions(self.hiv_ncoords, self.hiv_nt_seq, self.hiv_pcoords, self.hiv_aa_seq)

    def testGetCoords(self):
        region = GenomeRegion('test', self.hiv_ncoords)
        result = region.get_coords('nucl')
        expected = '5\'LTR,1,634\nGag,790,2292\nMatrix(p17/p15),790,1185\nCapsid(p24/p27),1186,1878\np2,1879,1920\n' \
                   'Nucleocapsid(p7/p8),1921,2085\np1,2086,2133\np6,2134,2292'
        self.assertEqual(expected, result)

    def testGetSequence(self):
        region = GenomeRegion('test', None, None, None, self.hiv_aa_seq)
        result = region.get_sequence('prot')
        expected = 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTK'\
                   'EALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQA'\
                   'AMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFR'\
                   'DYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKC'\
                   'FNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLR'\
                   'SLFGNDPSSQ'
        self.assertEqual(expected, result)

    def testSetRegions(self):

        region_names = {
            'Gag':
                [[790, 2292], 'ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATA'
                              'TAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAGAAGGCTGTA'
                              'GACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACCCTCTAT'
                              'TGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAAGC'
                              'ACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAAAATTACCCTATAGTGCAGAACATCCAGGGGCAAATGGTAC'
                              'ATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTGATACCCATGTTT'
                              'TCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACACCATGCTAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATGTT'
                              'AAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGAGTGCATCCAGTGCATGCAGGGCCTATTGCACCAGGCCAGATGAGAG'
                              'AACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACAAATAATCCACCTATCCCAGTA'
                              'GGAGAAATTTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGACA'
                              'AGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGGTTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAGGAGGTAAAAAATT'
                              'GGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCATTGGGACCAGCGGCTACACTAGAA'
                              'GAAATGATGACAGCATGTCAGGGAGTAGGAGGACCCGGCCATAAGGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACAAATTC'
                              'AGCTACCATAATGATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACACAG'
                              'CCAGAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAG'
                              'GCTAATTTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACC'
                              'AGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTT'
                              'CCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA',
                 [1, 500], 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLY'
                           'CVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMF'
                           'SALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPV'
                           'GEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLE'
                           'EMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQ'
                           'ANFLGKIWPSYKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ'],
            'Matrix(p17/p15)':
                [[790, 1185], 'ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATA'
                              'TAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAGAAGGCTGTA'
                              'GACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACCCTCTAT'
                              'TGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAAGC'
                              'ACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAAAATTAC',
                 [1, 132], 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLY'
                           'CVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY'],
            'Capsid(p24/p27)':
                [[1186, 1878], 'CCTATAGTGCAGAACATCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGA'
                               'GAAGGCTTTCAGCCCAGAAGTGATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACACCATGCTAAACA'
                               'CAGTGGGGGGACATCAAGCAGCCATGCAAATGTTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGAGTGCATCCAGTG'
                               'CATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAAT'
                               'AGGATGGATGACAAATAATCCACCTATCCCAGTAGGAGAAATTTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAA'
                               'TGTATAGCCCTACCAGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGGTTCTATAAAACTCTA'
                               'AGAGCCGAGCAAGCTTCACAGGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTAT'
                               'TTTAAAAGCATTGGGACCAGCGGCTACACTAGAAGAAATGATGACAGCATGTCAGGGAGTAGGAGGACCCGGCCATAAGGCAAGAG'
                               'TTTTG',
                 [133, 363], 'PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPV'
                             'HAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTL'
                             'RAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVL'],
            'p2':
                [[1879, 1920], 'GCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATG',
                 [364, 377], 'AEAMSQVTNSATIM'],
            'Nucleocapsid(p7/p8)':
                [[1921, 2085], 'ATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACACAGCCAGAAATTGCAGG'
                               'GCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAGGCTAAT',
                 [378, 432], 'MQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQAN'],
            'p1':
                [[2086, 2133], 'TTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTT',
                 [433, 448], 'FLGKIWPSYKGRPGNF'],
            'p6':
                [[2134, 2292], 'CTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAG'
                               'CCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA',
                 [449, 500], 'LQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ']}

        configs = handle_args('hiv', 'prot', self.hiv_nt_seq_path, self.hiv_ncoords_path,
                              self.hiv_aa_seq_path, self.hiv_pcoords_path)

        ref_nt_seq = configs[0][0][1]
        ref_aa_seq = configs[1]
        nt_coords = configs[2]
        aa_coords = configs[3]

        with open(nt_coords) as ncoords, open(aa_coords) as pcoords:
            result = make_regions(ncoords, ref_nt_seq, pcoords, ref_aa_seq)

            for i, reg in enumerate(result):
                self.assertEqual(list(region_names.keys())[i], reg.region_name)
                self.assertEqual(region_names[reg.region_name][0], reg.ncoords)
                self.assertEqual(region_names[reg.region_name][1], reg.nt_seq)
                self.assertEqual(region_names[reg.region_name][2], reg.pcoords)
                self.assertEqual(region_names[reg.region_name][3], reg.aa_seq)

    def testGetRefSeq(self):
        expected = [['K03455|HIVHXB2CG',
                     'TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACC'
                     'AGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAAC'
                     'ACCAGCTTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACA'
                     'TGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAG']]
        res = get_ref_seq(self.hiv_nt_seq_path, 'nucl')
        # Check the first 350 nucleotides in the default reference sequence
        result = [[res[0][0], res[0][1][:350]]]
        self.assertEqual(expected, result)

    def testSequenceAlign(self):
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

        # Aligned query sequence
        result = sequence_align(query, self.hiv_nt_seq, None)[1]
        self.assertEqual(expected, result)

    def testRetrieve(self):

        configs = handle_args('hiv', 'nucl', self.hiv_nt_seq_path, self.hiv_ncoords_path,
                              self.hiv_aa_seq_path, self.hiv_pcoords_path)
        ref_nt_seq = configs[0][0][1]
        ref_aa_seq = configs[1]
        nt_coords = configs[2]
        aa_coords = configs[3]

        with open(nt_coords) as ncoords, open(aa_coords) as pcoords:
            ref_regions = make_regions(ncoords, ref_nt_seq, pcoords, ref_aa_seq)

        result = retrieve('hiv', 'nucl', ref_regions, 'p2')
        query_region = result[0]
        expected_region = 'p2'
        expected_seq = 'GCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATG'
        self.assertEqual(expected_region, query_region.region_name)
        self.assertEqual(expected_seq, query_region.get_sequence('nucl'))

        overlap_regions = result[1]
        exp_region_names = ['Gag', 'p1']
        exp_pos_from_cds = [[1090, 1131], [1, 42]]
        exp_pos_from_qstart = [[1, 42], [1, 42]]
        exp_pos_from_gstart = [[1879, 1920], [1879, 1920]]
        exp_pos_from_pstart = [[364, 377], [1, 14]]
        expected_proteins = ['AEAMSQVTNSATIM', 'AEAMSQVTNSATIM']

        self.assertListEqual(list(overlap_regions.keys()), exp_region_names)
        for reg in result:
            if reg == 'Gag':
                self.assertEqual(exp_pos_from_cds[0], result[reg].rel_pos['CDS'], )
                self.assertEqual(exp_pos_from_qstart[0], result[reg].rel_pos['qstart'])
                self.assertEqual(exp_pos_from_gstart[0], result[reg].rel_pos['gstart'])
                self.assertEqual(exp_pos_from_pstart[0], result[reg].rel_pos['pstart'])
                self.assertEqual(expected_proteins[0], result[reg].get_sequence('prot'))

            elif reg == 'p1':
                self.assertEqual(exp_pos_from_cds[1], result[reg].rel_pos['CDS'], )
                self.assertEqual(exp_pos_from_qstart[1], result[reg].rel_pos['qstart'])
                self.assertEqual(exp_pos_from_gstart[1], result[reg].rel_pos['gstart'])
                self.assertEqual(exp_pos_from_pstart[1], result[reg].rel_pos['pstart'])
                self.assertEqual(expected_proteins[1], result[reg].get_sequence('prot'))

    def testHandleNuclArgs(self):
        """
        Tests the scenario when the user selects HIV and nucleotide alignment
        """
        result = handle_args('hiv', 'nucl', None, None, None, None)

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

    def testHandleNuclProtArgs(self):
        """
        Tests the scenario when the user selects HIV, specifies the reference protein and nucleotide sequences
        """

        result = handle_args('hiv', 'nucl', self.hiv_nt_seq_path, self.hiv_ncoords_path,
                             self.hiv_aa_seq_path, self.hiv_pcoords_path)

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

    def testNuclEquiv(self):
        configs = handle_args('hiv', 'prot', self.hiv_nt_seq_path, self.hiv_ncoords_path,
                              self.hiv_aa_seq_path, self.hiv_pcoords_path)

        ref_nt_seq = configs[0][0][1]
        ref_aa_seq = configs[1]
        nt_coords = configs[2]
        aa_coords = configs[3]

        with open(nt_coords) as ncoords, open(aa_coords) as pcoords:
            ref_regions = make_regions(ncoords, ref_nt_seq, pcoords, ref_aa_seq)
        for region in ref_regions:
            region.make_codon_aln()

        query_region = GenomeRegion('Gag', [2133, 2292], 'GCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATG',
                                    [364, 377], 'AEAMSQVTNSATIM')

        expected = 'GCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATG'
        result = set_nucleotide_equivalents(query_region)
        self.assertEqual(expected, result)

    def testNuclEquivWithGaps(self):
        ref_regions = [GenomeRegion('test1', [1, 31], 'GGGGGCCCGGGTTTAAACCCGGGTTTAAATTTC', [1, 11], 'GGPGLNPGLNF')]
        for region in ref_regions:
            region.make_codon_aln()
        q_region = GenomeRegion('test1', [12, 27],  'TTAAACCCGGGTTTA', [1, 5], 'LNPGL')
        expected = 'TTAAACCCGGGTTTA'
        result = set_nucleotide_equivalents(q_region)
        self.assertEqual(expected, result)

    def testProtEquiv(self):
        configs = handle_args('hiv', 'nucl', self.hiv_nt_seq_path, self.hiv_ncoords_path,
                              self.hiv_aa_seq_path, self.hiv_pcoords_path)

        ref_nt_seq = configs[0][0][1]
        ref_aa_seq = configs[1]
        nt_coords = configs[2]
        aa_coords = configs[3]

        with open(nt_coords) as ncoords, open(aa_coords) as pcoords:
            ref_regions = make_regions(ncoords, ref_nt_seq, pcoords, ref_aa_seq)

        for reg in ref_regions:
            reg.make_codon_aln()

        query_region = GenomeRegion('Gag', [2133, 2292], 'GCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATG',
                                    [364, 377], 'AEAMSQVTNSATIM')
        expected = 'AEAMSQVTNSATIM'
        result = set_protein_equivalents(query_region)
        self.assertEqual(expected, result)


class SIV(unittest.TestCase):
    def setUp(self):
        self.siv_nt_seq_path = os.path.join(
            os.path.join(os.path.dirname(__file__), os.pardir), 'ref_genomes/M33262.fasta')
        self.siv_aa_seq_path = os.path.join(os.path.dirname(__file__), 'fixtures/siv-test-proteins.fasta')
        self.siv_ncoords_path = os.path.join(os.path.dirname(__file__), 'fixtures/siv_test_nt_coords.csv')
        self.siv_pcoords_path = os.path.join(os.path.dirname(__file__), 'fixtures/siv_test_prot_coords.csv')

        with open(self.siv_nt_seq_path) as siv_nt, open(self.siv_aa_seq_path) as siv_aa, \
                open(self.siv_ncoords_path) as siv_ncoords, open(self.siv_pcoords_path) as siv_pcoords:
            self.siv_nt_seq = convert_fasta(siv_nt.read().split())[0][1]
            self.siv_aa_seq = convert_fasta(siv_aa.read().split())[0][1]
            self.siv_ncoords = siv_ncoords
            self.siv_pcoords = siv_pcoords


class TestSIV(SIV):

    maxDiff = None

    def setUp(self):
        self.siv_genome_regions = make_regions(self.siv_ncoords, self.siv_nt_seq, self.siv_pcoords, self.siv_aa_seq)

    def testGetCoords(self):
        region = GenomeRegion('test', None, None, self.siv_pcoords)
        result = region.get_coords('prot')
        expected = 'Rev(exon1),2087,2110\nRev(exon2),2111,2193'
        self.assertEqual(expected, result)

    def testGetSequence(self):
        region = GenomeRegion('test', None, self.siv_nt_seq)
        result = region.get_sequence('nucl')[0:1218]
        expected = 'GCATGCACATTTTAAAGGCTTTTGCTAAATATAGCCAAAAGTCCTTCTACAAATTTTCTAAGAGTTCTGATTCAAAGCAGTAACAGGCCTTGTCTCAT'\
                   'CATGAACTTTGGCATTTCATCTACAGCTAAGTTTATATCATAAATAGTTCTTTACAGGCAGCACCAACTTATACCCTTATAGCATACTTTACTGTGTG'\
                   'AAAATTGCATCTTTCATTAAGCTTACTGTAAATTTACTGGCTGTCTTCCTTGCAGGTTTCTGGAAGGGATTTATTACAGTGCAAGAAGACATAGAATC'\
                   'TTAGACATATACTTAGAAAAGGAAGAAGGCATCATACCAGATTGGCAGGATTACACCTCAGGACCAGGAATTAGATACCCAAAGACATTTGGCTGGCT'\
                   'ATGGAAATTAGTCCCTGTAAATGTATCAGATGAGGCACAGGAGGATGAGGAGCATTATTTAATGCATCCAGCTCAAACTTCCCAGTGGGATGACCCTT'\
                   'GGGGAGAGGTTCTAGCATGGAAGTTTGATCCAACTCTGGCCTACACTTATGAGGCATATGTTAGATACCCAGAAGAGTTTGGAAGCAAGTCAGGCCTG'\
                   'TCAGAGGAAGAGGTTAGAAGAAGGCTAACCGCAAGAGGCCTTCTTAACATGGCTGACAAGAAGGAAACTCGCTGAAACAGCAGGGACTTTCCACAAGG'\
                   'GGATGTTACGGGGAGGTACTGGGGAGGAGCCGGTCGGGAACGCCCACTTTCTTGATGTATAAATATCACTGCATTTCGCTCTGTATTCAGTCGCTCTG'\
                   'CGGAGAGGCTGGCAGATTGAGCCCTGGGAGGTTCTCTCCAGCACTAGCAGGTAGAGCCTGGGTGTTCCCTGCTAGACTCTCACCAGCACTTGGCCGGT'\
                   'GCTGGGCAGAGTGACTCCACGCTTGCTTGCTTAAAGCCCTCTTCAATAAAGCTGCCATTTTAGAAGTAAGCTAGTGTGTGTTCCCATCTCTCCTAGCC'\
                   'GCCGCCTGGTCAACTCGGTACTCAATAATAAGAAGACCCTGGTCTGTTAGGACCCTTTCTGCTTTGGGAAACCGAAGCAGGAAAATCCCTAGCAGATT'\
                   'GGCGCCTGAACAGGGACTTGAAGGAGAGTGAGAGACTCCTGAGTACGGCTGAGTGAAGGCAGTAAGGGCGGCAGGAACCAACCACGACGGAGTGCTCC'\
                   'TATAAAGGCGCGGGTCGGTACCAGACGGCGTGAGGAGCGGGA'
        self.assertEqual(expected, result)

    def testMakeRegions(self):
        region_names = {
            'Rev(exon1)':
                [[6784, 6853], 'ATGAGCAATCACGAAAGAGAAGAAGAACTCCGAAAAAGGCTAAGGCTAATACATCTTCTGCATCAAACAA',
                 [2087, 2110], 'MSNHEREEELRKRLRLIHLLHQT'],
            'Rev(exon2)':
                [[9062, 9315], 'ACCCATATCCAACAGGACCCGGCACTGCCAACCAGAGAAGGCAAAGAAAGAGACGGTGGAGAAGGCGGTGGCAACAGCTCCTGGCCT'
                               'TGGCAGATAGAATATATTCATTTCCTGATCCGCCAACTGATACGCCTCTTGACTTGGCTATTCAGCAACTGCAGAACCTTGCTATCG'
                               'AGAGTATACCAGATCCTCCAACCAATACTCCAGAGGCTCTCTGCGACCCTACAGAGGATTCGAGAAGTCCTCAGGACTGA',
                 [2111, 2193], 'PYPTGPGTANQRRQRKRRWRRRWQQLLALADRIYSFPDPPTDTPLDLAIQQLQNLAIESIPDPPTNTPEALCDPTEDSRSPQD']}

        configs = handle_args('siv', 'nucl', self.siv_nt_seq_path, self.siv_ncoords_path,
                              self.siv_aa_seq_path, self.siv_pcoords_path)

        ref_nt_seq = configs[0][0][1]
        ref_aa_seq = configs[1]
        nt_coords = configs[2]
        aa_coords = configs[3]

        with open(nt_coords) as ncoords, open(aa_coords) as pcoords:
            result = make_regions(ncoords, ref_nt_seq, pcoords, ref_aa_seq)

            for i, reg in enumerate(result):
                self.assertEqual(list(region_names.keys())[i], reg.region_name)
                self.assertEqual(region_names[reg.region_name][0], reg.ncoords)
                self.assertEqual(region_names[reg.region_name][1], reg.nt_seq)
                self.assertEqual(region_names[reg.region_name][2], reg.pcoords)
                self.assertEqual(region_names[reg.region_name][3], reg.aa_seq)

    def testGetRefSeq(self):
        expected = [['Rev(exon1)|SIVMM239', 'MSNHEREEELRKRLRLIHLLHQT']]
        res = get_ref_seq(self.siv_aa_seq_path, 'prot')
        result = [[res[0][0], res[0][1]]]
        self.assertEqual(expected, result)

    def testSequenceAlign(self):
        query = [['query', 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVA'
                           'TLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY']]

        expected = ['query',
                    'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTK'
                    'EALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY----------------------------------------------------------------'
                    '--------------------------------------------------------------------------------------------------'
                    '--------------------------------------------------------------------------------------------------'
                    '--------------------------------------------------------------------------------------------------'
                    '----------']

        result = sequence_align(query, self.siv_aa_seq, None)[1]
        self.assertEqual(expected, result)

    def testRetrieve(self):

        configs = handle_args('siv', 'nucl', self.siv_nt_seq_path, self.siv_ncoords_path,
                              self.siv_aa_seq_path, self.siv_pcoords_path)

        ref_nt_seq = configs[0][0][1]
        ref_aa_seq = configs[1]
        nt_coords = configs[2]
        aa_coords = configs[3]

        with open(nt_coords) as ncoords, open(aa_coords) as pcoords:
            ref_regions = make_regions(ncoords, ref_nt_seq, pcoords, ref_aa_seq)

        result = retrieve('nucl', 'Nef', 20, 80)
        result_region = result[0]
        expected_region = 'Nef'
        expected_seq = 'TGAGGCGGTCCAGGCCGTCTGGAGATCTGCGACAGAGACTCTTGCGGGCGCGTGGGGAGAC'
        self.assertEqual(expected_region, result_region.region_name)
        self.assertEqual(expected_seq, result_region.get_sequence('nucl'))

        overlap_regions = result[1]
        exp_region_names = ['Env(gp160)', 'gp41', 'Nef']
        exp_pos_from_cds = [[2493, 2553], [918, 978], [20, 80]]
        exp_pos_from_qstart = [[1, 61], [1, 61], [1, 61]]
        exp_pos_from_gstart = [[9096, 9156], [9096, 9156], [9096, 9156]]
        exp_pos_from_pstart = [[832, 851], [307, 326], [7, 27]]
        expected_proteins = ['VWRSATETLAGAWGD', 'VWRSATETLAGAWGD', 'SGDLRQRLLRARGE']

        self.assertListEqual(list(overlap_regions.keys()), exp_region_names)
        for reg in result:
            if reg == 'Env(gp160)':
                self.assertEqual(exp_pos_from_cds[0], result[reg].rel_pos['CDS'], )
                self.assertEqual(exp_pos_from_qstart[0], result[reg].rel_pos['qstart'])
                self.assertEqual(exp_pos_from_gstart[0], result[reg].rel_pos['gstart'])
                self.assertEqual(exp_pos_from_pstart[0], result[reg].rel_pos['pstart'])
                self.assertEqual(expected_proteins[0], result[reg].get_sequence('prot'))

            elif reg == 'gp41':
                self.assertEqual(exp_pos_from_cds[1], result[reg].rel_pos['CDS'], )
                self.assertEqual(exp_pos_from_qstart[1], result[reg].rel_pos['qstart'])
                self.assertEqual(exp_pos_from_gstart[1], result[reg].rel_pos['gstart'])
                self.assertEqual(exp_pos_from_pstart[1], result[reg].rel_pos['pstart'])
                self.assertEqual(expected_proteins[1], result[reg].get_sequence('prot'))

            elif reg == 'Nef':
                self.assertEqual(exp_pos_from_cds[2], result[reg].rel_pos['CDS'], )
                self.assertEqual(exp_pos_from_qstart[2], result[reg].rel_pos['qstart'])
                self.assertEqual(exp_pos_from_gstart[2], result[reg].rel_pos['gstart'])
                self.assertEqual(exp_pos_from_pstart[2], result[reg].rel_pos['pstart'])
                self.assertEqual(expected_proteins[2], result[reg].get_sequence('prot'))

    def testHandleProtArgs(self):
        """
        Tests the scenario when the user selects SIV and protein alignment
        """

        result = handle_args('siv', 'prot', None, None, None, None)

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

    def testHandleArgs(self):
        """
        Tests the scenario when the user specifies only the protein sequences
        """
        result = handle_args('siv', 'prot', None, None, self.siv_aa_seq_path, self.siv_pcoords_path)

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


class TestSetCoords(TestHIV):
    def testSetNuclCoords(self):
        region = GenomeRegion('test', None, None, 'ATGCGACGACGACGACGTTAGCAGTCGATCATCATGCTGATC')
        region.set_coords([1, 10], 'nucl')
        expected = [1, 10]
        result = region.ncoords
        self.assertEqual(expected, result)

    def testSetProtCoords(self):
        region = GenomeRegion('test', None, None, None, 'ATVLHGLKLKAMAIMTIM')
        region.set_coords([1, 900], 'prot')
        expected = [1, 900]
        result = region.pcoords
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
        region.set_coords([1, 634], 'nucl')
        region.set_pos_from_cds([1, 634])
        expected = []
        result = region.cds_offset
        self.assertEqual(expected, result)

    def testFromSIVLTR5Start(self):
        region = GenomeRegion('5\'LTR')
        region.set_pos_from_cds([1, 1110])
        expected = []
        result = region.cds_offset
        self.assertEqual(expected, result)

    def testFromHIVStart(self):
        ref_region = GenomeRegion('Gag', [790, 2292])
        query_region = GenomeRegion('Gag', [790, 2292])
        query_region.set_pos_from_cds(ref_region)
        expected = [1, 1503]
        result = query_region.cds_offset
        self.assertEqual(expected, result)

    def testFromSIVStart(self):
        ref_region = GenomeRegion('Integrase', [4785, 5666])
        query_region = GenomeRegion('Integrase', [4785, 5666])
        query_region.set_pos_from_cds(ref_region)
        expected = [1, 882]
        result = query_region.cds_offset
        self.assertEqual(expected, result)

    def testGagfromP1Start(self):
        region = GenomeRegion('Gag', [2086, 2133])
        region.set_pos_from_cds([2086, 2133])
        expected = [1297, 1344]
        result = region.cds_offset
        self.assertEqual(expected, result)

    def testPolFromP1Start(self):
        region = GenomeRegion('Pol', [2086, 2133])
        region.set_pos_from_cds([2086, 2133])
        expected = [2, 49]
        result = region.cds_offset
        self.assertEqual(expected, result)


class TestSetPosFromAAStart(unittest.TestCase):

    def testCapsidFromAAStart(self):
        region = GenomeRegion('Capsid',
                              [1186, 1878], 'CCTATAGTGCAGAACATCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAA'
                                            'AGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTGATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCAC'
                                            'AAGATTTAAACACCATGCTAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATGTTAAAAGAGACCATCAAT'
                                            'GAGGAAGCTGCAGAATGGGATAGAGTGCATCCAGTGCATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACC'
                                            'AAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACAAATAATCCACCTA'
                                            'TCCCAGTAGGAGAAATTTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACC'
                                            'AGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGGTTCTATAAAACTCTAAG'
                                            'AGCCGAGCAAGCTTCACAGGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATT'
                                            'GTAAGACTATTTTAAAAGCATTGGGACCAGCGGCTACACTAGAAGAAATGATGACAGCATGTCAGGGAGTAGGA'
                                            'GGACCCGGCCATAAGGCAAGAGTTTTG',
                              [1, 231], 'PIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAA'
                                        'EWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQG'
                                        'PKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVL')
        region.set_pos_from_pstart()
        expected = [1, 231]
        result = region.pstart
        self.assertEqual(expected, result)

    def testLTR5FromAAStart(self):
        region = GenomeRegion('5\'LTR')
        region.set_pos_from_pstart()
        expected = []
        result = region.pstart
        self.assertEqual(expected, result)

    def testGagFromAAStart(self):
        region = GenomeRegion('Gag',
                              [2137, 2289], 'CAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTC'
                                            'AGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTC'
                                            'GTCACAA',
                              [450, 500],   'QSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ')
        region.set_pos_from_pstart()
        expected = [450, 500]
        result = region.pstart
        self.assertEqual(expected, result)


class TestMakeCodonAln(unittest.TestCase):

    def testSimpleUse(self):
        region = GenomeRegion('p6', [2134, 2292], 'CTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAAC'
                                                  'TCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTG'
                                                  'GCAACGACCCCTCGTCACAATAA',
                                    [449, 500],   'LQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ*')
        region.make_codon_aln()
        expected = '-L--Q--S--R--P--E--P--T--A--P--P--E--E--S--F--R--S--G--V--E--T--T--T--P--P--Q-' \
                   '-K--Q--E--P--I--D--K--E--L--Y--P--L--T--S--L--R--S--L--F--G--N--D--P--S--S--Q--*-'
        result = region.codon_aln
        self.assertEqual(expected, result)


class TestGlobalToLocalIndex(unittest.TestCase):

    def testNuclCoords(self):
        region = GenomeRegion('p1', [2086, 2133], 'TTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTT',
                              [433, 448], 'FLGKIWPSYKGRPGNF')
        expected = [1, 48]
        result = region.global_to_local_index([2086, 2133], 'nucl')
        self.assertEqual(expected, result)

    def testProtCoords(self):
        region = GenomeRegion('p2', [1879, 1920], 'GCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATG',
                              [364, 377], 'AEAMSQVTNSATIM')
        expected = [1, 14]
        result = region.global_to_local_index([364, 377], 'prot')
        self.assertEqual(expected, result)


class TestLocalToGlobalIndex(unittest.TestCase):

    def testNuclCoords(self):
        region = GenomeRegion('Nucleocapsid',
                              [1921, 2085], 'ATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACACAGC'
                                            'CAGAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTA'
                                            'CTGAGAGACAGGCTAAT',
                              [1, 55], 'MQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQAN')
        expected = [1921, 2085]
        result = region.local_to_global_index(region, [1, 165], 'nucl')
        self.assertEqual(expected, result)

    def testProtCoords(self):
        region = GenomeRegion('RNase',
                              [3870, 4229], 'TATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTAATAGAGGAAGACAAAA'
                                            'AGTTGTCACCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATCTAGCTTTGCAGGATTCGG'
                                            'GATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATCATTCAAGCACAACCAGATCAAAGTGAA'
                                            'TCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAAAAAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACA'
                                            'CAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTA',
                              [1096, 1215], 'YVDGAANRETKLGKAGYVTNRGRQKVVTLTDTTNQKTELQAIYLALQDSGLEVNIVTDSQYALGIIQAQPDQSE'
                                            'SELVNQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKVL')
        expected = [1096, 1215]
        result = region.local_to_global_index(region, [1, 120], 'prot')
        self.assertEqual(expected, result)


class TestProtCoordsFromNuclCoords(unittest.TestCase):

    def testSimpleUse(self):
        region = GenomeRegion('Gag', [2086, 2133], 'TTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTT', None, None)
        region.set_pos_from_cds([2086, 2133])
        region.set_pcoords_from_ncoords()
        expected = [433, 448]
        self.assertEqual(expected, region.pcoords)


class TestGetOverlap(unittest.TestCase):

    def testNuclOverlap(self):
        region = GenomeRegion(
            'gp41',
            [7758, 8795], 'GCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCCTCAATGACGCTGACGGTACAGGCCAG'
                          'ACAATTATTGTCTGGTATAGTGCAGCAGCAGAACAATTTGCTGAGGGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCA'
                          'TCAAGCAGCTCCAGGCAAGAATCCTGGCTGTGGAAAGATACCTAAAGGATCAACAGCTCCTGGGGATTTGGGGTTGCTCTGGAAAACTCATT'
                          'TGCACCACTGCTGTGCCTTGGAATGCTAGTTGGAGTAATAAATCTCTGGAACAGATTTGGAATCACACGACCTGGATGGAGTGGGACAGAGA'
                          'AATTAACAATTACACAAGCTTAATACACTCCTTAATTGAAGAATCGCAAAACCAGCAAGAAAAGAATGAACAAGAATTATTGGAATTAGATA'
                          'AATGGGCAAGTTTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAATTATTCATAATGATAGTAGGAGGCTTGGTAGGTTTAA'
                          'GAATAGTTTTTGCTGTACTTTCTATAGTGAATAGAGTTAGGCAGGGATATTCACCATTATCGTTTCAGACCCACCTCCCAACCCCGAGGGGA'
                          'CCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGATTAGTGAACGGATCCTTGGCACTTATCTG'
                          'GGACGATCTGCGGAGCCTGTGCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCA'
                          'GGGGGTGGGAAGCCCTCAAATATTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAGCTTGCTCAATGCCACA'
                          'GCCATAGCAGTAGCTGAGGGGACAGATAGGGTTATAGAAGTAGTACAAGGAGCTTGTAGAGCTATTCGCCACATACCTAGAAGAATAAGACA'
                          'GGGCTTGGAAAGGATTTTGCTATAA',
            [2602, 2946], 'AVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLI'
                          'CTTAVPWNASWSNKSLEQIWNHTTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKLFIMIVGGLVGL'
                          'RIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGR'
                          'RGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL'
        )
        expected = ['CAGATAGGGTTATAGAAGTAGTACAAGGAGCTTGTAGAGCTATTCGCCACATACCTAGAAGAATAAGACAGGGCTTGGAAAGGATTTTGCTATAA',
                    [8700, 8795]]
        result = region.get_overlap(region, [8700, 8795], 'nucl')
        self.assertEqual(expected, result)

    def testProtOverlap(self):
        region = GenomeRegion(
            'Matrix',
            [790, 1185], 'ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTA'
                         'AAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAGAAGGCTGTAGACAAATACTGGGA'
                         'CAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACCCTCTATTGTGTGCATCAAAGGATAGAG'
                         'ATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAAGCACAGCAAGCAGCAGCTGACACAGGACAC'
                         'AGCAATCAGGTCAGCCAAAATTAC',
            [1, 132], 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIE'
                      'IKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNY'
        )
        expected = ['MGARASVLSGGELDR', [1, 15]]
        result = region.get_overlap(region, [1, 15], 'prot')
        self.assertEqual(expected, result)


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
        result = get_query("nucl", "ATGCGCG", False)
        self.assertEqual(expected, result)

    def testProteinQuery(self):
        expected = [["query1", "MPPLMMADLADLGG"]]
        query = ">query1\nMPPLMMADLADLGG"
        result = get_query("prot", query, False)
        self.assertEqual(expected, result)

    def testLongNucleotideSequence(self):
        expected = [["query", "ATGCGCGAATTAGCGA"]]
        query = "atgcgcg\naattagcga"
        result = get_query("nucl", query, False)
        self.assertEqual(expected, result)

    def testRevCompNucl(self):
        query = ">seq1\nTCGCTAATTCGCGCATN*"
        expected = [["seq1", "*NATGCGCGAATTAGCGA"]]
        result = get_query("nucl", query, True)
        self.assertEqual(expected, result)

    def testInvalidNucleotideQuery(self):
        query = ">query\natgcgcg&\n"
        with self.assertRaises(SystemExit) as e:
            get_query("nucl", query, 'n')
        self.assertEqual(e.exception.code, 0)

    def testInvalidProteinQuery(self):
        query = ">query\nMPPLMMAD>LADLGG"
        with self.assertRaises(SystemExit) as e:
            get_query("prot", query, 'n')
        self.assertEqual(e.exception.code, 0)

    def testRevCompProt(self):
        handle = ">seq2\nMPPLMMADLADLGG\n"
        expected = [["seq2", "MPPLMMADLADLGG"]]
        result = get_query('prot', handle, True)
        self.assertEqual(expected, result)


class TestReverseComp(unittest.TestCase):

    def testSimpleUse(self):
        expected = "TCGCTAATTCGCGCATN*"
        result = reverse_comp("*NATGCGCGAATTAGCGA")
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
        result = query_region_coordinates(aln)
        self.assertEqual(expected, result)

    def testProtAlignment(self):
        aln = ['query',
               'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDK'
               'IEEEQNKSKKKAQQAAADTGHSNQVSQNY--------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '-------------------------------------------------------------------------------------------------------'
               '----------------------------------------------------------------------------------------']
        expected = [(0, 132)]
        result = query_region_coordinates(aln)
        self.assertEqual(expected, result)


class TestFindMatches(unittest.TestCase):

    def testNuclMatches(self):

        coordinates = [(2133, 2292)]
        result = find_matches('nucl', coordinates)

        exp_region_names = ['Gag',           'p6']
        exp_pos_from_cds = [[1344, 1503],    [1, 159]]
        exp_pos_from_qstart = [[1, 160],     [2, 160]]
        exp_pos_from_gstart = [[2133, 2292], [2134, 2292]]
        exp_pos_from_pstart = [[449, 501],   [1, 53]]
        exp_aa_seq = ['XLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ',
                      'LQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ']
        exp_nt_seq = ['TCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGA'
                      'CAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA',
                      'TCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGA'
                      'CAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA']

        self.assertListEqual(exp_region_names, list(result.keys()))
        for reg in result:
            if reg == 'Gag':
                self.assertEqual(exp_pos_from_cds[0], result[reg].rel_pos['CDS'], )
                self.assertEqual(exp_pos_from_qstart[0], result[reg].rel_pos['qstart'])
                self.assertEqual(exp_pos_from_gstart[0], result[reg].rel_pos['gstart'])
                self.assertEqual(exp_pos_from_pstart[0], result[reg].rel_pos['pstart'])
                self.assertEqual(exp_aa_seq[0], result[reg].get_sequence('prot'))
                self.assertEqual(exp_nt_seq[0], result[reg].get_sequence('nucl'))

            elif reg == 'p6':
                self.assertEqual(exp_pos_from_cds[1], result[reg].rel_pos['CDS'], )
                self.assertEqual(exp_pos_from_qstart[1], result[reg].rel_pos['qstart'])
                self.assertEqual(exp_pos_from_gstart[1], result[reg].rel_pos['gstart'])
                self.assertEqual(exp_pos_from_pstart[1], result[reg].rel_pos['pstart'])
                self.assertEqual(exp_aa_seq[1], result[reg].get_sequence('prot'))
                self.assertEqual(exp_nt_seq[1], result[reg].get_sequence('nucl'))


def create_suite():
    test_suite = unittest.TestSuite()
    test_suite.addTest(unittest.TestCase())
    test_suite.addTest(TestHIV())
    test_suite.addTest(TestSIV())
    return test_suite


if __name__ == '__main__':
    suite = create_suite()

    runner = unittest.TextTestRunner()
    runner.run(suite)
