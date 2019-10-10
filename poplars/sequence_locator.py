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
        self.cds_offset, self.gstart, self.qstart, self.pstart = [], [], [], []
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

    def set_sequence(self, seq, base):
        if base == 'nucl':
            self.nt_seq = seq
        else:
            self.aa_seq = seq

    def set_seq_from_ref(self, sequence, base):
        if base == 'nucl':
            self.nt_seq = sequence[self.ncoords[0] - 1: self.ncoords[1]]
        else:
            self.aa_seq = sequence[self.pcoords[0] - 1: self.pcoords[1]]

    def set_pos_from_cds(self, ref_reg_name):
        """
        Gives the position of a the query sequence relative to the start of the coding sequence
        """

        ref_region = GENOME_REGIONS[ref_reg_name]
        if self.ncoords and ref_region.ncoords:

            query_start = self.ncoords[0]
            query_end = self.ncoords[1]
            ref_start = ref_region.ncoords[0]
            ref_end = ref_region.ncoords[1]

            if query_start == ref_start:
                start = 1
            else:
                start = query_start - ref_start

            if query_end == ref_end:
                end = query_end - query_start
            else:
                end = query_end - start

            self.cds_offset = [start, end]

    def set_pos_from_gstart(self):
        self.gstart = self.ncoords

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
            self.qstart = [start_offset, end_offset]

    def set_pos_from_pstart(self):
        """
        Gives the position of the sequence relative to the start of the protein sequence
        """
        if self.pcoords is not None:
            if 'LTR' in self.region_name:
                self.pstart = None
            else:
                self.pstart = self.pcoords

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

    @staticmethod
    def global_to_local_index(coord_pair):
        """
        Converts a pair of global indices to local indices, relative to the sequence region of interest
        """
        start = coord_pair[0]  # 0-based inclusive indexing
        start_offset = coord_pair[0] - start
        end_offset = coord_pair[1] - start
        local_pair = [start_offset, end_offset]
        return local_pair

    @staticmethod
    def local_to_global_index(r_coords, local_pair):
        """
        Converts a pair of local indices (relative to the region of interest) to global indices
        """
        if r_coords is not None:
            global_start = r_coords[0] + local_pair[0] - 1
            global_end = r_coords[0] + local_pair[1] - 1
            global_pair = [global_start, global_end]
            return global_pair

    def set_pcoords_from_ncoords(self):
        """
        Sets protein coordinates relative to the protein start, given the nucleotide coordinates
        """
        if self.cds_offset is not None:
            prot_start = self.cds_offset[0] // 3 + 1
            prot_end = self.cds_offset[1] // 3
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


def make_regions(nt_coords, nt_seq, aa_coords, aa_seq):
    """
    Reads in the start and end coordinates of the genomic regions and associates the region with its coordinates.
    If no coordinate files are specified, set default nucleotide and protein coordinate files.
    :param nt_coords: Path to the csv file containing the global coordinates of the nucleotide region.
            The file stream has one genomic entry per line and has the following format: region_name,start,end
    :param nt_seq: The nucleotide sequence
    :param aa_coords: Path to the csv file containing the global coordinates of the protein region.
            The file stream has one genomic entry per line and has the following format: region_name,start,end
    :param aa_seq: A list of lists containing the protein sequences
    :return genome_regions: A dictionary with keys as the region names and the values are GenomeRegion objects
    """

    genome_regions = {}
    # Parse nucleotide region coordinates file
    for nt_line in nt_coords:
        nt_line = nt_line.strip()
        nt_line = nt_line.split(',')
        nucl_coords = [int(nt_line[1]), int(nt_line[2])]
        seq_region = GenomeRegion(nt_line[0])

        # Set nucleotide coordinates and sequences
        seq_region.set_coords(nucl_coords, 'nucl')
        seq_region.set_seq_from_ref(nt_seq, 'nucl')
        genome_regions[nt_line[0]] = seq_region

    # Parse protein coordinates file
    prot_names, prot_coords = [], []
    for aa_line in aa_coords:
        aa_line = aa_line.strip()
        aa_line = aa_line.split(',')
        prot_names.append(aa_line[0])
        prot_coords.append([int(aa_line[1]), int(aa_line[2])])

    # Match protein regions to nucleotide regions
    for i, coords in enumerate(prot_coords):
        for region_name in genome_regions:
            if prot_names[i].startswith(genome_regions[region_name].region_name):
                # Set protein coordinates and sequences
                genome_regions[region_name].set_coords(coords, 'prot')
                genome_regions[region_name].set_sequence(aa_seq[i][1], 'prot')

    # Align nucleotide sequence with protein sequence
    for name in genome_regions:
        genome_regions[name].make_codon_aln()

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


def get_ref_seq(ref_seq, base):
    """
    Converts the reference sequence to a string and checks if the sequence is valid
    :param ref_seq: The path to the reference sequence
    :param base: The base of the sequence
    :return reference_sequence: A list of the form [header, sequence]
    """
    with open(ref_seq, 'r') as ref_handle:
        reference_sequence = convert_fasta(ref_handle)

    if valid_sequence(base, reference_sequence):
        return reference_sequence
    else:
        sys.exit(0)


def reverse_comp(query_sequence):
    """
    Reverses and complements the query sequence
    :param query_sequence: The query sequence
    :return: The reverse complement of the query sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '*': '*', 'N': 'N', '-': '-'}
    rev_comp = "".join(complement.get(nt, nt) for nt in reversed(query_sequence))
    return rev_comp


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


def query_region_coordinates(alignment):
    """
    Gets the indices of regions where the query aligns with the reference sequence
    :param alignment: A list containing the header of the query sequence, and the aligned query sequence
    :return query_matches: The positions where the query sequence aligns with the reference sequences(s) with no gaps
    """
    pat = re.compile('[A-Z]{2,}')  # Match the aligned region (minimum alignment length is 2)
    query_matches = [(match.start(), match.end()) for match in pat.finditer(alignment[1])]
    return query_matches


def find_matches(base, query_matches):
    """
    Finds the genomic regions where the query sequence aligns with the reference sequence
    :param base: The base of the query sequence
    :param query_matches: A list of indices where the query sequence aligns with the reference sequence
    :return query_regions: A dictionary of matches, with keys as the region name and values as the GenomeRegion object
    """
    query_regions = {}

    for q_coords in query_matches:
        for name in GENOME_REGIONS:
            if GENOME_REGIONS[name].region_name != 'Complete':

                # Use coordinates to find which regions overlap with the query region
                overlap = find_overlap(base, GENOME_REGIONS[name].region_name, q_coords)
                if overlap['region_name'] is not None:
                    query_region = GenomeRegion(overlap['region_name'])
                    query_region.set_sequence(overlap['nt_seq'], 'nucl')
                    query_region.set_coords(overlap['nt_coords'], 'nucl')
                    query_region.set_sequence(overlap['aa_seq'], 'prot')
                    query_region.set_coords(overlap['aa_coords'], 'prot')

    return query_regions


def find_overlap(base, reg, q_coords):
    """
    Gets the sequence regions that overlap with the region of interest
    :param base: the base of the sequence (nucleotide or protein)
    :param reg: the genome region
    :param q_coords: the coordinates of the query region (global coordinates)
    :return: the sequence of the overlap, and the indices of the overlap
    """

    overlap = {'region_name': None, 'nt_seq': None, 'nt_coords': None, 'aa_seq': None, 'aa_coords': None}
    start, end = 0, 0
    ref_coords = GENOME_REGIONS[reg].get_coords(base)       # Coordinates are global coordinates

    # If the ref_coords are in the range of the query region and the q_coords are in the range of the ref region
    if not ((ref_coords[0] < ref_coords[1] < q_coords[0]) or (q_coords[0] < q_coords[1] < ref_coords[0])):

        # If the end of the query region exceeds the end of the reference region
        if q_coords[1] > ref_coords[1]:
            end = ref_coords[1]
        else:
            end = q_coords[1]

        # If the query region starts before the reference region starts
        if q_coords[0] < ref_coords[0]:
            start = ref_coords[0]
        else:
            start = q_coords[0]

    overlap['region_name'] = GENOME_REGIONS[reg].region_name

    if base == 'nucl':
        q_nt_seq = GENOME_REGIONS[reg].get_sequence('nucl')[start: end]
        overlap['nt_seq'] = q_nt_seq
        overlap['nt_coords'] = q_coords

        prot_overlap = set_protein_equivalents(overlap)
        overlap['aa_seq'] = prot_overlap[0]
        overlap['aa_coords'] = prot_overlap[1]

    else:
        q_aa_seq = reg.get_sequence('prot')[start: end]
        overlap['aa_seq'] = q_aa_seq
        overlap['aa_coords'] = q_coords

    return overlap


def set_protein_equivalents(overlap):
    """
    Finds the protein equivalent of the nucleotide sequence
    :param overlap: the region that aligns with the query sequence
    :return aa_seq, aa_coords: the sequence and coordinates of the overlapping protein sequence
    """
    non_coding = ["5'LTR", "TAR", "3'LTR"]
    aa_seq, aa_coords = '', []

    for ref_reg in GENOME_REGIONS:
        if GENOME_REGIONS[ref_reg].region_name == overlap['region_name'] \
                and GENOME_REGIONS[ref_reg].region_name not in non_coding:
            if GENOME_REGIONS[ref_reg].codon_aln is not None:

                # Convert query nt_coordinates to local coordinates
                local_nt_coords = GenomeRegion.global_to_local_index(overlap['nt_coords'])

                # Slice aligned nucleotide sequence at the nucleotide position
                aligned_prot = GENOME_REGIONS[ref_reg].codon_aln[local_nt_coords[0]: local_nt_coords[1]]

                aligned_prot_seq = re.sub('[-]', '', aligned_prot)        # Get corresponding protein sequence
                aa_seq = aligned_prot_seq
                aa_coords = [1, len(aligned_prot_seq)]

    return aa_seq, aa_coords


def set_nucleotide_equivalents(overlap):
    """
    Finds the nucleotide equivalent of the protein sequence
    :param overlap: the region that aligns with the query sequence
    :return nt_seq, nt_coords: the sequence and coordinates of the overlapping nucleotide sequence
    """
    nt_seq, nt_coords = '', []

    for ref_reg in GENOME_REGIONS:
        if ref_reg.region_name == overlap['region_name'] and ref_reg.codon_aln is not None:
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


def retrieve(base, region, qstart=1, qend='end'):
    """
    Retrieves a sequence given its coordinates
    :param base: The base of the sequence (nucleotide or protein)
    :param region: The genomic region
    :param qstart: <option> The start coordinate of the query region (given as local coordinate)
    :param qend: <option> The end coordinate of the query region (given as local coordinate)
    :return: The genomic region defined by the starting and ending coordinates
    """
    query_region = None
    overlap_regions = {}
    for ref_region in GENOME_REGIONS:

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
            global_coords = GenomeRegion.local_to_global_index(ref_region, [qstart, qend])
            query_region.set_coords(global_coords, base)

            # Set sequences protein and nucleotide sequences
            seq = ref_region.get_sequence(base)[qstart - 1: qend + 1]
            query_region.set_sequence(seq, base)

            # Set equivalent sequence
            if base == 'nucl':
                equiv_seq = set_nucleotide_equivalents(query_region)
                query_region.set_sequence(equiv_seq, 'prot')
            else:
                equiv_seq = set_protein_equivalents(query_region)
                query_region.set_sequence(equiv_seq, 'nucl')

        if query_region is not None:
            retrieved_regions = find_matches(base, [query_region.get_coords(base)])

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
                             'Ex: region_name,start,end')
    parser.add_argument('-aa_coords', default=None,
                        help='Path to the csv file containing the coordinates of the amino acid region.'
                             'The file must contain the region name, start coordinate, and end coordinate.'
                             'Ex: region_name,start,end')
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
    ref_nt_seq, ref_aa_seq = configs[0][0][1], configs[1]
    nt_coords, aa_coords = configs[2], configs[3]
    reference_sequence = configs[4]

    with open(nt_coords) as nt_coords_handle, open(aa_coords) as aa_coords_handle:

        # Create genomic region objects based on configuration files
        global GENOME_REGIONS
        GENOME_REGIONS = make_regions(nt_coords_handle, ref_nt_seq, aa_coords_handle, ref_aa_seq)

        if args.subcommand == "align":
            query = get_query(args.base, args.query, args.revcomp)
            alignment = sequence_align(query, reference_sequence, args.outfile)

            # Find indices where the query sequence aligns with the reference sequence
            match_coords = query_region_coordinates(alignment[-1])  # Query sequence will be the last item in the list
            query_regions = find_matches(args.base, match_coords)
            output_overlap(query_regions, args.outfile)

        else:
            valid_in = valid_inputs(args.virus, args.start, args.end, args.region)

            if not valid_in:
                sys.exit(0)
            else:
                result = retrieve(args.base, args.region, args.start, args.end)
                query_region = result[0]
                overlap_regions = result[1]

                output_retrieved_region(query_region, args.outfile)        # Ouptput retrieved region first
                if args.outfile is None:
                    print("\n\nRegions touched by the query sequence:")
                else:
                    args.outfile.write("\n\nRegions touched by the query sequence:\n\n")

                output_overlap(overlap_regions, args.outfile)


if __name__ == '__main__':
    main()
