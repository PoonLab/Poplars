"""
Python implementation of HIV Sequence Locator from https://www.hiv.lanl.gov
HIV and SIV genomic region coordinates are based on the HXB2 and Mac239
annotation resources from https://www.hiv.lanl.gov/content/sequence/HIV/MAP/annotation.html

Note: The first 256 nucleotides of SIVMM239 correspond to the flanking sequence,
and are included in the complete SIV genome (https://www.ncbi.nlm.nih.gov/nucleotide/M33262)
"""

import re
from poplars.mafft import *
from poplars.common import convert_fasta
import textwrap


class GenomeRegion:
    """
    Represents information about a genomic region
    """

    def __init__(self, region_name, ncoords=None, nt_seq=None, pcoords=None, aa_seq=None):
        """
        Stores information about each genomic region
        :param region_name: The name of the genomic region
        :param ncoords: A list containing the global start and end indices of the nucleotide region
        :param nt_seq:  The nucleotide sequence
        :param pcoords: A list containing the global start and end indices of the protein region
        :param aa_seq: The amino acid sequence
        """
        self.region_name = region_name
        self.ncoords = ncoords
        self.nt_seq = nt_seq
        self.pcoords = pcoords
        self.aa_seq = aa_seq
        self.rel_pos = {'CDS': [], 'gstart': [], 'qstart': [], 'pstart': []}
        self.codon_aln = ''

    def get_coords(self, base):
        if base == 'nucl':
            coords = self.ncoords
        else:
            coords = self.pcoords
        return coords

    def get_sequence(self, base):
        if base == 'nucl':
            sequence = self.nt_seq
        else:
            sequence = self.aa_seq
        return sequence

    def set_coords(self, coord_pair, base):
        if base == 'nucl':
            self.ncoords = coord_pair
        else:
            self.pcoords = coord_pair

    def set_seq_from_ref(self, sequence, base):
        if base == 'nucl':
            self.nt_seq = sequence[self.ncoords[0] - 1: self.ncoords[1]]
        else:
            self.aa_seq = sequence[self.pcoords[0] - 1: self.pcoords[1]]

    def set_sequence(self, seq, base):
        if base == 'nucl':
            self.nt_seq = seq
        else:
            self.aa_seq = seq

    def set_pos_from_cds(self, region_coords):
        """
        Gives the position of a sequence relative to the start of the coding sequence
        """
        if self.ncoords is not None:
            local_ncoords = self.global_to_local_index(self.get_coords('nucl'), 'nucl')
            print('local ncoords {}'.format(local_ncoords))

            if self.region_name != '5\'LTR':
                len_region = region_coords[1] - region_coords[0]
                cds_start = region_coords[0] - local_ncoords[0] + 1
                cds_end = cds_start + len_region
                self.rel_pos['CDS'].append(cds_start)
                self.rel_pos['CDS'].append(cds_end)

    def set_pos_from_gstart(self):
        self.rel_pos['gstart'] = self.ncoords

    def set_pos_from_pstart(self):
        """
        Gives the position of the sequence relative to the start of the protein sequence
        """
        if self.pcoords is not None:
            if 'LTR' in self.region_name:
                self.rel_pos['pstart'] = None
            else:
                self.rel_pos['pstart'] = self.pcoords

    def set_pos_from_qstart(self, q_coords, base):
        """
        Gives the position of the sequence relative to the start of the region of interest
        :param q_coords: The coodinates of the query region
        :param base: The base of the sequence (nucleotide or protein)
        :return: The position relative to the start of the region of interest
        """
        r_coords = self.get_coords(base)
        r_seq = self.get_sequence(base)
        if r_coords is not None and r_seq is not None:
            start_offset = r_coords[0] - q_coords[0] + 1
            end_offset = start_offset + len(r_seq) - 1
            self.rel_pos['qstart'] = [start_offset, end_offset]

    def make_codon_aln(self):
        """
        Aligns the protein sequence with its associated nucleotide sequence
        """
        if self.nt_seq is not None and self.aa_seq is not None:
            codon_aln = []
            codons = [''.join(t) for t in zip(*[iter(self.nt_seq)] * 3)]
            for aa, codon in zip(self.aa_seq, codons):
                # Check if stop codon
                if codon == 'TAA' or codon == 'TGA' or codon == 'TAG':
                    codon_aln.append('-*-')
                else:
                    codon_aln.append('-{}-'.format(aa))
            self.codon_aln = ''.join(codon_aln)
        return self.codon_aln

    def global_to_local_index(self, coord_pair, base):
        """
        Converts a pair of global indices to local indices, relative to the sequence region of interest
        """
        start = self.get_coords(base)[0] - 1  # 0-based inclusive indexing
        start_offset = coord_pair[0] - start
        end_offset = coord_pair[1] - start
        local_pair = [start_offset, end_offset]
        return local_pair

    @staticmethod
    def local_to_global_index(region, local_pair, base):
        """
        Converts a pair of local indices (relative to the region of interest) to global indices
        """
        r_coords = region.get_coords(base)
        if r_coords is not None:
            global_start = r_coords[0] + local_pair[0] - 1
            global_end = r_coords[0] + local_pair[1] - 1
            global_pair = [global_start, global_end]
            return global_pair

    def set_pcoords_from_ncoords(self):
        """
        Sets protein coordinates relative to the protein start, given the nucleotide coordinates
        """
        if self.rel_pos['CDS'] is not None:
            prot_start = self.rel_pos['CDS'][0] // 3 + 1
            prot_end = self.rel_pos['CDS'][1] // 3
            self.pcoords = [prot_start, prot_end]

    def set_ncoords_from_pcoords(self):
        """
        Sets nucleotide coordinates given the protein coordinates relative to the protein start
        :return ncoords: the nucleotide coordinates
        """
        if self.pcoords is not None:
            nucl_start = self.pcoords[0] - 1 * 3
            nucl_end = self.pcoords[1] * 3
            return [nucl_start, nucl_end]

    def get_overlap(self, region, coord_pair, base):
        """
        Gets the sequence regions that overlap with the region of interest
        :param region: the region of interest
        :param coord_pair: the coordinates
        :param base: The base of the sequence
        :return: the sequence of the overlap, and the indices of the overlap
        """
        overlap = []
        local_pair = self.global_to_local_index(coord_pair, base)
        seq = self.get_sequence(base)
        coords = region.get_coords(base)

        if coords[1] > coord_pair[0]:
            # Check if local coordinates are in the range of the sequence
            start = max(local_pair[0], 1)
            end = min(local_pair[1] + 1, len(seq))
            if base == 'nucl':
                overlap.append(seq[start - 1: end + 1])
                overlap.append(self.local_to_global_index(region, [start, end], base))
            else:
                overlap.append(seq[start - 1: end - 1])
                overlap.append(self.local_to_global_index(region, [start, end - 1], base))

        return overlap


def set_regions(base, nt_coords, nt_seq, aa_coords, aa_seq):
    """
    Reads in the start and end coordinates of the genomic regions and associates the region with its coordinates.
    If no coordinate files are specified, set default nucleotide and protein coordinate files.
    :param base: The base of the sequence
    :param nt_coords: Path to the csv file containing the global coordinates of the nucleotide region.
            The file stream has one genomic entry per line and has the following format: region_name,start,end
    :param nt_seq: The nucleotide sequence
    :param aa_coords: Path to the csv file containing the global coordinates of the protein region.
            The file stream has one genomic entry per line and has the following format: region_name,start,end
    :param aa_seq: A list of lists containing the amino acid sequences
    :return: A list of GenomeRegions
    """

    genome_regions = []

    if base == 'nucl':

        # Parse nucleotide region coordinates file
        for nt_line in nt_coords:
            nt_line = nt_line.strip()
            nt_line = nt_line.split(',')
            nucl_coords = [int(nt_line[1]), int(nt_line[2])]

            seq_region = GenomeRegion(nt_line[0])

            # Set global and local nucleotide coordinates
            seq_region.set_coords(nucl_coords, 'nucl')
            seq_region.set_seq_from_ref(nt_seq, 'nucl')

            # Set relative positions
            seq_region.set_pos_from_cds(nucl_coords)
            seq_region.set_pos_from_gstart()
            seq_region.set_pos_from_pstart()

            genome_regions.append(seq_region)

        # Parse protein coordinates file
        prot_names = []
        prot_coords = []
        for aa_line in aa_coords:
            aa_line = aa_line.strip()
            aa_line = aa_line.split(',')
            prot_names.append(aa_line[0])
            prot_coords.append([int(aa_line[1]), int(aa_line[2])])

        for i, coords in enumerate(prot_coords):
            for seq_region in genome_regions:
                if prot_names[i].startswith(seq_region.region_name):
                    # Set global and local protein coordinates
                    seq_region.set_coords(coords, 'prot')
                    seq_region.set_sequence(aa_seq[i][1], 'prot')
                    seq_region.set_pos_from_pstart()
                    seq_region.make_codon_aln()

    else:
        # Parse protein region coordinates file
        for i, aa_line in enumerate(aa_coords):
            aa_line = aa_line.strip()
            aa_line = aa_line.split(',')
            prot_coords = [int(aa_line[1]), int(aa_line[2])]

            seq_region = GenomeRegion(aa_line[0])

            # Set global and local nucleotide coordinates
            seq_region.set_coords(prot_coords, 'prot')
            seq_region.set_sequence(aa_seq[i][1], 'prot')

            # Set relative positions
            seq_region.set_pos_from_cds(prot_coords)
            seq_region.set_pos_from_gstart()
            seq_region.set_pos_from_pstart()

            genome_regions.append(seq_region)

        # Parse nucleotide coordinates file
        nucl_names = []
        nucl_coords = []
        for nt_line in nt_coords:
            nt_line = nt_line.strip()
            nt_line = nt_line.split(',')
            nucl_names.append(nt_line[0])
            nucl_coords.append([int(nt_line[1]), int(nt_line[2])])

        for i, name in enumerate(nucl_names):
            for seq_region in genome_regions:
                if name.startswith(seq_region.region_name):
                    # Set global and local protein coordinates
                    seq_region.set_coords(nucl_coords[i], 'nucl')
                    seq_region.set_seq_from_ref(nt_seq, 'nucl')
                    seq_region.set_pos_from_pstart()
                    seq_region.make_codon_aln()

    return genome_regions


def valid_sequence(base, sequence):
    """
    Verifies that input sequence is valid
    :param base: The base of the sequence (nucl or prot)
    :param sequence: A list of lists containing header and sequence pairs
    :raises ValueError: If the sequence is empty or if it contains invalid characters
    :return: <True> If the input sequence uses the correct alphabet, <False> otherwise
    """
    dna_alphabet = 'ATGC-*XN'
    aa_alphabet = 'ARDNCEQGHILKMFPSTWYV-*X'

    if not sequence:
        print("Invalid sequence: sequence length is 0\n")
        return False

    for h, s in sequence:
        if not s:
            print("Invalid sequence: sequence length is 0\n")
            return False
        elif base == 'nucl':
            if not all(pos in dna_alphabet for pos in s):
                print("Invalid nucleotide sequence:\n{}\n{}\n".format(h, s))
                return False
        else:
            if not all(pos in aa_alphabet for pos in s):
                print("Invalid amino acid sequence:\n{}\n{}\n".format(h, s))
                return False

    return True


def valid_inputs(virus, start_coord, end_coord, region):
    """
    Checks that the start and end coordinates of the genomic region are valid
    :param virus: The reference virus
    :param start_coord: The starting coordinate of the genomic region
    :param end_coord: The ending coordinate of the genomic region
    :param region: The genomic region
    :return: <True> If the coordinates are valid, <False> otherwise
    """

    if type(start_coord) == str:
        print("Invalid start coordinate type: {}".format(type(start_coord)))

    if start_coord < 0:
        print("Invalid start coordinate: {}".format(start_coord))
        return False

    if type(end_coord) != str and type(end_coord) != int:
        print("Invalid end coordinate type: {}".format(type(end_coord)))
        return False

    if type(end_coord) == str:
        if end_coord != "end":
            print("Invalid end coordinate: {}".format(end_coord))
            return False
    else:
        if end_coord <= 0 or start_coord >= end_coord:
            print("Invalid range: {} to {}".format(start_coord, end_coord))
            return False

    if virus == 'hiv' and region == 'Vpx' or virus == 'siv' and region == 'Vpu':
        print("Invalid region: {} in {}".format(region, virus))
        return False

    return True


def get_query(base, query, revcomp):
    """
    Gets the query sequence and checks that it is valid
    :param base: The base (nucleotide or protein)
    :param query: The query sequence as a string or the file path to the query sequence
    :param revcomp: Option to align the reverse complement of the nucleotide sequence
    :return: A list of lists containing the sequence identifiers and the query sequences
    """

    # If the query is a string
    if not os.path.exists(query):
        query_lines = query.split('\n')
        [line.strip('\n') for line in query_lines]

        # Convert string into a list of lists
        if query_lines[0].startswith('>'):
            header = query_lines[0].strip('>')
            seq = ''.join(query_lines[1:])
        else:
            header = "query"
            seq = ''.join(query_lines)
        query_seq = [[header, seq.upper()]]

    # If the query is a file path
    else:
        with open(query, 'r') as query_handle:
            line = query_handle.readline()
            query_handle.seek(0)  # Reset pointer to beginning

            # Parse query sequence
            if line.startswith('>'):  # Fasta file
                query_seq = convert_fasta(query_handle)

            else:
                query_seq = []
                count = 1
                for line in query_handle:
                    line = line.strip('\n')
                    if len(line) > 0:
                        query_seq.extend([["Sequence{}".format(count), "{}".format(line.upper())]])
                        count += 1

    if not valid_sequence(base, query_seq):
        sys.exit(0)

    # At this point, the sequence is valid
    if revcomp:
        if base == 'prot':
            print("Invalid option: reverse complement is not available for proteins.")
        else:
            rc_query = reverse_comp(query_seq[0][1])
            header = query_seq[0][0]
            query_seq = [[header, rc_query]]

    return query_seq


def reverse_comp(query_sequence):
    """
    Reverses and complements the query sequence
    :param query_sequence: The query sequence
    :return: The reverse complement of the query sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '*': '*', 'N': 'N', '-': '-'}
    rev_comp = "".join(complement.get(nt, nt) for nt in reversed(query_sequence))
    return rev_comp


def get_ref_seq(ref_seq, base):
    """
    Converts the reference sequence to a string and checks if the sequence is valid
    :param ref_seq: The path to the reference sequence
    :param base: The base of the sequence
    :return reference_sequence: The reference sequence as a list where the first element is the header
                    and the second element is the sequence as a string
    """

    with open(ref_seq, 'r') as ref_handle:
        reference_sequence = convert_fasta(ref_handle)

    if valid_sequence(base, reference_sequence):
        return reference_sequence
    else:
        sys.exit(0)


def sequence_align(query_sequence, reference_sequence, outfile=None):
    """
    Aligns the query sequence to the reference genome
    :param query_sequence: The query sequence
    :param reference_sequence: The reference sequence
    :param outfile: <option> The file stream of the output file in write mode
    """

    result = align(query_sequence[0][1], reference_sequence)

    if outfile is not None:
        outfile.write("Alignment:\n")
    else:
        print("Alignment:")

    for h, s in result:
        if outfile is not None:
            outfile.write(">{}\n".format(h))
            outfile.write(textwrap.fill(s, 60))
        else:
            print(">{}".format(h))
            print(textwrap.fill(s, 60))

    return result


def get_region_coordinates(alignment):
    """
    Gets the indices of regions where the query aligns with the reference sequence
    :param alignment: A list containing the header of the query sequence, and the aligned query sequence
    :return: The positions where the query sequence aligns with the reference sequences(s) with no gaps
    """
    pat = re.compile('[A-Z]{2,}')  # Match the aligned region (minimum alignment length is 2)
    coordinates = [(match.start(), match.end()) for match in pat.finditer(alignment[1])]
    return coordinates


def find_matches(base, ref_regions, match_coordinates):
    """
    Finds the genomic regions where the query sequence aligns with the reference sequence
    :param base: The base of the query sequence
    :param ref_regions: A list of GenomeRegion objects
    :param match_coordinates: A list of indices where the query sequence aligns with the reference sequence
    :return query_regions: A dictionary of matches, with keys as the region name and values as the GenomeRegion object
    """
    query_regions = {}

    for coord in match_coordinates:
        start_aln = coord[0]
        end_aln = coord[1]

        for ref_region in ref_regions:
            if ref_region.region_name != "Complete":
                overlap = ref_region.get_overlap(ref_region, [start_aln, end_aln], base)
                if len(overlap) == 2:
                    ov_seq = overlap[0]
                    ov_coord = overlap[1]

                    if ov_seq:
                        query_region = GenomeRegion(ref_region.region_name)
                        query_region.set_sequence(ov_seq, base)

                        # Set protein and nucleotide coordinates
                        query_region.set_coords(ov_coord, base)

                        if base == 'nucl':
                            prot_start = ov_coord[0] // 3
                            prot_end = ov_coord[1] // 3
                            query_region.set_coords([prot_start, prot_end], 'prot')
                            set_protein_equivalents(query_region, ref_regions)
                        else:
                            nucl_start = ov_coord[0] * 3
                            nucl_end = ov_coord[1] * 3
                            query_region.set_coords([nucl_start, nucl_end], 'nucl')
                            set_nucleotide_equivalents(query_region, ref_regions)

                        # Set relative positions
                        query_region.set_pos_from_cds(ref_region.get_coords(base))
                        query_region.set_pos_from_gstart()
                        query_region.set_pos_from_qstart(coord, base)
                        query_region.set_pos_from_pstart()

                        query_regions[query_region.region_name] = query_region

    return query_regions


def set_protein_equivalents(query_reg, ref_regions):
    """
    Finds the protein equivalent of the nucleotide sequence
    :param query_reg: the region that aligns with the query sequence
    :param ref_regions: a list of GenomeRegion objects
    :return prot_equiv: the protein equivalent of the nucleotide sequence
    """
    prot_seq = ''
    prot_coords = []
    non_coding = ["5'LTR", "TAR", "3'LTR"]
    for ref_reg in ref_regions:
        if ref_reg.region_name == query_reg.region_name and ref_reg.region_name not in non_coding:
            if ref_reg.codon_aln is not None and query_reg.pcoords is not None:
                prot_coords = query_reg.get_overlap(query_reg.region, query_reg.get_coords('prot'), 'prot')
                prot_seq = ref_reg.codon_aln[query_reg.pcoords[0]: query_reg.pcoords[1]]
                prot_seq = re.sub('[-]', '', prot_seq)
                query_reg.set_sequence('prot', prot_seq)
                query_reg.set_coords(prot_coords, 'prot')

    return prot_seq, prot_coords


def set_nucleotide_equivalents(query_reg, ref_regions):
    """
    Finds the nucleotide equivalent of the protein sequence
    :param query_reg: the region that aligns with the query sequence
    :param ref_regions: a list of GenomeRegion objects
    :return nt_equiv: the nucleotide equivalent of the protein sequence
    """
    nt_seq = ''
    nt_coords = []
    for ref_reg in ref_regions:
        if ref_reg.region_name == query_reg.region_name and ref_reg.codon_aln is not None:
            if query_reg.ncoords is not None:
                query_reg.make_codon_aln()
                regex = re.compile(query_reg.codon_aln)
                coords = regex.search(ref_reg.codon_aln).span()
                nt_coords = query_reg.get_overlap(query_reg.region_name, query_reg.get_coords('nucl'), 'nucl')
                nt_seq = ref_reg.nt_seq[coords[0]: coords[1]]
                query_reg.set_sequence('nucl', nt_seq)
                query_reg.set_coords(nt_coords, 'nucl')

    return nt_seq, nt_coords


def output_retrieved_region(region, outfile=None):
    """
    Outputs the retrieved region
    :param region: A list of GenomeRegions where the query sequence aligns with the reference sequence
    :param outfile: The file stream for the output file in write mode
    """

    if outfile is None:
        print("\nRetrieved Region:\t{}".format(region.region_name))
        print("\tNucleotide Sequence:")
        seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
        for line in seq_lines:
            print('\t\t{}'.format(line))

        if region.aa_seq:
            print("\tProtein Sequence:")
            seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
            for line in seq_lines:
                print('\t\t{}\n'.format(line))

        print("\n\tRelative Positions:")

        if region.rel_pos['CDS']:
            print("\t\tNucleotide position relative to CDS start:\t{} --> {}"
                  .format(region.rel_pos['CDS'][0], region.rel_pos['CDS'][1]))
        else:
            print("\t\tNucleotide position relative to CDS start:\tN/A")

        if region.rel_pos['gstart']:
            print("\t\tNucleotide position relative to genome start:\t{} --> {}"
                  .format(region.rel_pos['gstart'][0], region.rel_pos['gstart'][1]))

        if region.rel_pos['pstart']:
            print("\t\tAmino acid position relative to protein start:\t{} --> {}"
                  .format(region.rel_pos['pstart'][0], region.rel_pos['pstart'][1]))

        if region.rel_pos['qstart']:
            print("\t\tPosition relative to query start:\t{} --> {}"
                  .format(region.rel_pos['qstart'][0], region.rel_pos['qstart'][1]))

    else:
        outfile.write("\nRetrieved Region:\t{}".format(region.region_name))
        outfile.write("\tNucleotide Sequence:\n")
        seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
        for line in seq_lines:
            outfile.write('\t\t{}\n'.format(line))

        if region.aa_seq:
            outfile.write("\tProtein Sequence:\n")
            seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
            for line in seq_lines:
                outfile.write('\t\t{}\n'.format(line))

        outfile.write("\n\tRelative Positions:")

        if region.rel_pos['CDS']:
            outfile.write("\t\tNucleotide position relative to CDS start:\t{} --> {}"
                          .format(region.rel_pos['CDS'][0], region.rel_pos['CDS'][1]))
        else:
            outfile.write("\t\tNucleotide position relative to CDS start:\tN/A")

        if region.rel_pos['gstart']:
            outfile.write("\t\tNucleotide position relative to genome start:\t{} --> {}"
                          .format(region.rel_pos['gstart'][0] + 1, region.rel_pos['gstart'][1] + 1))

        if region.rel_pos['pstart']:
            outfile.write("\t\tAmino acid position relative to protein start:\t{} --> {}"
                          .format(region.rel_pos['pstart'][0], region.rel_pos['pstart'][1]))

        if region.rel_pos['qstart']:
            outfile.write("\t\tPosition relative to query start:\t{} --> {}"
                          .format(region.rel_pos['qstart'][0], region.rel_pos['qstart'][1]))


def output_overlap(overlap_regions, outfile=None):
    """
    Outputs the regions that overlap with the query
    :param overlap_regions: A dictionary of GenomeRegions where the keys are the region name and
                                the values are the query sequence aligns with the reference sequence
    :param outfile: The file stream for the output file in write mode
    """

    if outfile is None:

        for key in overlap_regions:
            region = overlap_regions[key]
            if region.region_name.startswith('5\'LTR'):
                print("\t3'LTR\n")

            print("\nRegion:\t{}".format(region.region_name))
            print("\n\tNucleotide Sequence:")
            seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
            for line in seq_lines:
                print('\t\t{}'.format(line))

            if region.aa_seq is not None:
                print("\n\tProtein Sequence:")
                seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
                for line in seq_lines:
                    print('\t\t{}\n'.format(line))

            print("\n\tRelative Positions:")

            if region.rel_pos['CDS']:
                print("\t\tNucleotide position relative to CDS start:\t{} --> {}"
                      .format(region.rel_pos['CDS'][0], region.rel_pos['CDS'][1]))
            else:
                print("\t\tNucleotide position relative to CDS start: N/A")

            if region.rel_pos['gstart']:
                print("\t\tNucleotide position relative to genome start:\t{} --> {}"
                      .format(region.rel_pos['gstart'][0], region.rel_pos['gstart'][1]))

            if region.rel_pos['pstart']:
                print("\t\tAmino acid position relative to protein start:\t{} --> {}"
                      .format(region.rel_pos['pstart'][0], region.rel_pos['pstart'][1]))

            if region.rel_pos['qstart']:
                print("\t\tPosition relative to query start:\t{} --> {}"
                      .format(region.rel_pos['qstart'][0], region.rel_pos['qstart'][1]))

    else:
        for key in overlap_regions:
            region = overlap_regions[key]
            if region.region_name.startswith('5\'LTR'):
                outfile.write("\t3'LTR\n")

            outfile.write("\nRegion:\t{}".format(region.region_name))
            outfile.write("\n\tNucleotide Sequence:\n")
            seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
            for line in seq_lines:
                outfile.write('\t\t{}'.format(line))

            if region.aa_seq is not None:
                outfile.write("\n\tProtein Sequence:\n")
                seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
                for line in seq_lines:
                    print('\t\t{}'.format(line))

            outfile.write("\n\tRelative Positions:")

            if region.rel_pos['CDS']:
                outfile.write("\t\tNucleotide position relative to CDS start:\t{} --> {}"
                              .format(region.rel_pos['CDS'][0], region.rel_pos['CDS'][1]))
            else:
                outfile.write("\t\tNucleotide position relative to CDS start: N/A")

            if region.rel_pos['gstart']:
                print("\t\tNucleotide position relative to genome start:\t{} --> {}"
                      .format(region.rel_pos['gstart'][0] + 1, region.rel_pos['gstart'][1] + 1))

            if region.rel_pos['pstart']:
                outfile.write("\t\tAmino acid position relative to protein start:\t{} --> {}"
                              .format(region.rel_pos['pstart'][0], region.rel_pos['pstart'][1]))

            if region.rel_pos['qstart']:
                outfile.write("\t\tPosition relative to query start:\t{} --> {}"
                              .format(region.rel_pos['qstart'][0], region.rel_pos['qstart'][1]))


def retrieve(base, ref_regions, region, qstart=1, qend='end'):
    """
    Retrieves a sequence given its coordinates
    :param base: The base of the sequence (nucleotide or protein)
    :param ref_regions: A list of GenomeRegion objects
    :param region: The genomic region
    :param qstart: <option> The start coordinate of the query region (given as local coordinate)
    :param qend: <option> The end coordinate of the query region (given as local coordinate)
    :return: The genomic region defined by the starting and ending coordinates
    """

    query_region = None
    overlap_regions = {}
    for ref_region in ref_regions:

        if ref_region.region_name == region:
            global_range = ref_region.get_coords(base)

            # Convert global region coordinates to local coordinates
            local_range = ref_region.global_to_local_index(global_range, base)

            region_start, region_end = local_range[0], local_range[1]

            if qend == 'end':
                qend = region_end - region_start

            if qstart < region_start or qend > region_end:
                print("Invalid {} coordinates: {} to {}.\nValid range: 1 to {}"
                      .format(region, qstart, qend, (region_end - region_start)))
                sys.exit(0)

            query_region = GenomeRegion(region)

            # Set local and global coordinates
            global_coords = GenomeRegion.local_to_global_index(ref_region, [qstart, qend], base)
            query_region.set_coords(global_coords, base)

            # Set sequences protein and nucleotide sequences
            seq = ref_region.get_sequence(base)[qstart - 1: qend + 1]
            query_region.set_sequence(seq, base)

            # Set equivalent sequence
            if base == 'nucl':
                equiv_seq = set_nucleotide_equivalents(query_region, ref_regions)
                query_region.set_sequence(equiv_seq, 'prot')
            else:
                equiv_seq = set_protein_equivalents(query_region, ref_regions)
                query_region.set_sequence(equiv_seq, 'nucl')

        if query_region is not None:
            retrieved_regions = find_matches(base, ref_regions, [query_region.get_coords(base)])

            # Remove duplicated retrieved region
            for key in retrieved_regions:
                if retrieved_regions[key].region_name != region:
                    overlap_regions[key] = retrieved_regions[key]

    return query_region, overlap_regions


def handle_args(virus, base, ref_nt, nt_coords, ref_aa, aa_coords):
    """
    Handles the possible execution paths for the program
    :param virus: The reference virus
    :param base: The base of the reference sequence
    :param ref_nt: Path to the file containing the reference nucleotide sequence
    :param nt_coords: Path to the csv file containing the coordinates of the genomic regions
    :param ref_aa: Path to the file stream containing the reference amino acid sequence
    :param aa_coords: Path to the csv file containing the coordinates of the protein regions
    :return configs: a list containing the query sequence and paths to the reference nucleotide and protein sequences
    """
    if virus == 'hiv':
        if ref_nt is None:
            ref_nt = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455.fasta")
        if nt_coords is None:
            nt_coords = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455_genome_coordinates.csv")
        if ref_aa is None:
            ref_aa = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455-protein.fasta")
        if aa_coords is None:
            aa_coords = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455_protein_coordinates.csv")

    else:
        if ref_nt is None:
            ref_nt = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262.fasta")
        if nt_coords is None:
            nt_coords = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262_genome_coordinates.csv")
        if ref_aa is None:
            ref_aa = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262-protein.fasta")
        if aa_coords is None:
            aa_coords = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262_protein_coordinates.csv")

    ref_nt_seq = get_ref_seq(ref_nt, 'nucl')
    ref_aa_seq = get_ref_seq(ref_aa, 'prot')

    # Set the reference sequence to be used in the sequence alignment
    if base == 'nucl':
        reference_sequence = ref_nt_seq
    else:
        reference_sequence = ref_aa_seq

    configs = [ref_nt_seq, ref_aa_seq, nt_coords, aa_coords, reference_sequence]

    return configs


def parse_args():
    """
    Parses command line arguments
    """
    regions = ['5\'LTR', '5\'LTR-R', '5\'LTR-U3', '5\'LTR-U5', 'TAR', 'Gag-Pol', 'Gag', 'Matrix',
               'Capsid', 'p2', 'Nucleocapsid', 'p1', 'p6', 'Pol', 'GagPolTF', 'Protease', 'RT', 'RNase',
               'Integrase', 'Vif', 'Vpr', 'Tat(with intron)', 'Tat(exon1)', 'Tat(exon2)', 'Rev(with intron)',
               'Rev(exon1)', 'Rev(exon2)', 'Vpu', 'Vpx', 'Env', 'gp160', 'V1', 'V2', 'V3', 'V4', 'V5',
               'RRE', 'gp120', 'gp41', 'Nef', '3\'LTR', '3\'LTR-R', '3\'LTR-U3', '3\'LTR-U5', 'Complete']

    parser = argparse.ArgumentParser(
        description='An implementation of the HIV Sequence Locator tool by the Los Alamos National Laboratory.'
                    'This tool aligns a nucleotide or protein sequence relative to HIV or SIV reference genomes; '
                    'or retrieves a sequence in the HXB2 or SIVmm239 reference genome from its coordinates.',
    )
    parser.add_argument('virus', metavar='virus', choices=['hiv', 'siv'],
                        help='The reference virus')
    parser.add_argument('base', metavar='base', choices=['nucl', 'prot'],
                        help='Sequence base type. Allowed bases are \'nucl\' and \'prot\'')
    parser.add_argument('-nt_coords', default=None,
                        help='Path to the csv file containing the coordinates of the nucleotide region.'
                             'The file must be contain the region name, start coordinate, and end coordinate.'
                             'coordinates separated by a comma (,). Ex: region_name    start,end')
    parser.add_argument('-aa_coords', default=None,
                        help='Path to the csv file containing the coordinates of the amino acid region.'
                             'The file must contain the region name, start coordinate, and end coordinate')
    parser.add_argument('-ref_nt', metavar='',
                        help='Path to the file containing the reference nucleotide sequence')
    parser.add_argument('-ref_aa', metavar='',
                        help='Path to the file containing the reference amino acid sequence')
    subparsers = parser.add_subparsers(dest='subcommand')

    # Create subparser for 'align' mode
    parser_align = subparsers.add_parser('align', formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                         description='Align a nucleotide or protein sequence '
                                                     'relative to the HIV or SIV reference genome')
    parser_align.add_argument('query', help='The query sequence. This can be the a string or a file path.')
    parser_align.add_argument('--revcomp', action="store_true", default=False,
                              help='Align the reverse complement of the query sequence with the reference sequence')
    parser_align.add_argument('-outfile', type=argparse.FileType('w'),
                              help='Path to the file where results will be written. '
                                   'If no file is specified, the results will be printed to the console')

    # Create subparser for 'retrieve' mode
    parser_retrieve = subparsers.add_parser('retrieve', formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                            description='Retrieve a sequence in HXB2 or SIVmm239 from its coordinates')
    parser_retrieve.add_argument('-region', metavar='', default="Complete", choices=regions,
                                 help='List of case sensitive genomic regions. '
                                      'Allowed regions are: ' + ', '.join(regions))
    parser_retrieve.add_argument('-start', default=1, type=int,
                                 help='Coordinate relative to the start of the region.')
    parser_retrieve.add_argument('-end', default='end',
                                 help='Coordinate relative to the end of the region. Enter an integer or \'end\'.')
    parser_retrieve.add_argument('-outfile', metavar='', type=argparse.FileType('w'),
                                 help='File where results will be written. If no file is specified'
                                      'the results will be printed to the console')
    return parser.parse_args()


def main():
    args = parse_args()

    # Ensure proper configuration files are set
    configs = handle_args(args.virus, args.base, args.ref_nt, args.nt_coords, args.ref_aa, args.aa_coords)
    ref_nt_seq = configs[0][0][1]
    ref_aa_seq = configs[1]
    nt_coords = configs[2]
    aa_coords = configs[3]
    reference_sequence = configs[4]

    with open(nt_coords) as nt_coords_handle, open(aa_coords) as aa_coords_handle:

        # Create genomic region objects based on configuration files
        ref_regions = set_regions(args.base, nt_coords_handle, ref_nt_seq, aa_coords_handle, ref_aa_seq)

        if args.subcommand == "align":
            query = get_query(args.base, args.query, args.revcomp)
            alignment = sequence_align(query, reference_sequence, args.outfile)

            # Find indices where the query sequence aligns with the reference sequence
            match_coords = get_region_coordinates(alignment[-1])  # Query sequence will be the last item in the list
            query_regions = find_matches(args.base, ref_regions, match_coords)

            output_overlap(query_regions, args.outfile)

        else:
            valid_in = valid_inputs(args.virus, args.start, args.end, args.region)

            if not valid_in:
                sys.exit(0)
            else:
                result = retrieve(args.base, ref_regions, args.region, args.start, args.end)
                query_region = result[0]
                overlap_regions = result[1]

                # Print retrieved region first
                output_retrieved_region(query_region, args.outfile)

                if args.outfile is None:
                    print("\n\nRegions touched by the query sequence:")
                else:
                    args.outfile.write("\n\nRegions touched by the query sequence:\n\n")

                output_overlap(overlap_regions, args.outfile)


if __name__ == '__main__':
    main()
