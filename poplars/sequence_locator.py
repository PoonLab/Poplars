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

    def __init__(self, region_name, nt_coords=None, nt_seq=None, aa_coords=None, aa_seq=None):
        """
        Stores information about each genomic region
        :param region_name: The name of the genomic region
        :param nt_coords: A list containing the start and end coordinates of the nucleotide region. Ex: [[1, 890]]
        :param nt_seq: <option> The nucleotide sequence of the genomic region
        :param aa_coords: <option> A list containing the start and end coordinates of the protein region. Ex: [[8,78]]
        :param aa_seq: <option> The amino acid sequence of the genomic region
        """
        self.region_name = region_name
        self.nt_coords = nt_coords
        self.nt_seq = nt_seq
        self.aa_coords = aa_coords
        self.aa_seq = aa_seq

        self.pos_from_cds = []         # Nucleotide position relative to CDS start
        self.pos_from_gstart = None    # Nucleotide position relative to the reference genome start
        self.pos_from_qstart = None    # Nucleotide position relative to the start of the query sequence
        self.pos_from_rstart = None
        self.pos_from_aa_start = None  # Amino acid position relative to protein start
        self.codon_aln = ''            # Amino acid sequence aligned with the nucleotide sequence

    def get_coords(self, base):
        if base == 'nucl':
            coords = self.nt_coords
        else:
            coords = self.aa_coords
        return coords

    def get_sequence(self, base):
        if base == 'nucl':
            sequence = self.nt_seq
        else:
            sequence = self.aa_seq
        return sequence

    def set_coords(self, coord_pair, base):
        if base == 'nucl':
            self.nt_coords = coord_pair
        else:
            self.aa_coords = coord_pair

    def set_seq_from_ref(self, sequence, base):
        if base == 'nucl':
            self.nt_seq = sequence[self.nt_coords[0] - 1: self.nt_coords[1]]
        else:
            self.aa_seq = sequence[self.aa_coords[0] - 1: self.aa_coords[1]]

    def set_sequence(self, seq, base):
        if base == 'nucl':
            self.nt_seq = seq
        else:
            self.aa_seq = seq

    def set_pos_from_cds(self, virus):
        """
        Gives the position of a sequence relative to the start of the coding sequence
        """
        if '5\'LTR' in self.region_name:
            self.pos_from_cds = ['N/A']
        else:
            if virus == 'hiv':
                cds_start = 790
            else:
                cds_start = 1309
            self.pos_from_cds.append((self.nt_coords[0] + 1 - cds_start))
            self.pos_from_cds.append((self.nt_coords[1] + 1 - cds_start))

    def set_pos_from_aa_start(self, virus):
        """
        Gives the position of the sequence relative to the start of the protein sequence
        """
        if self.aa_seq is not None:
            if self.pos_from_cds is not None:
                self.set_pos_from_cds(virus)

            # If the region is within the CDS
            if len(self.pos_from_cds) == 2:
                # If the whole protein sequence is encompassed in the CDS range
                if len(self.aa_seq) == (((self.pos_from_cds[1] - self.pos_from_cds[0]) // 3) + 1):
                    self.pos_from_aa_start = [1, (((self.pos_from_cds[1] - self.pos_from_cds[0]) // 3) + 1)]
                else:
                    self.pos_from_aa_start = [(self.pos_from_cds[0] // 3) + 1, self.pos_from_cds[1] // 3]

            # If the region is outside the CDS
            else:
                self.pos_from_aa_start = None

    def make_codon_aln(self):
        """
        Aligns the protein sequence with its associated nucleotide sequence
        """
        if self.nt_seq is not None and self.aa_seq is not None:
            codon_aln = []
            codons = [''.join(t) for t in zip(*[iter(self.nt_seq)] * 3)]
            for i in range(len(codons)):
                # Check if stop codon
                if codons[i] == 'TAA' or codons[i] == 'TGA' or codons[i] == 'TAG':
                    codon_aln.append('-*-')
                else:
                    codon_aln.append('-{}-'.format(self.aa_seq[i]))
            self.codon_aln = ''.join(codon_aln)
        return self.codon_aln

    def global_to_local_index(self, coord_pair, base):
        """
        Converts a pair of global indices to local indices, relative to the sequence region of interest
        """
        start = self.get_coords(base)[0] - 1         # 0-based inclusive indexing
        start_offset = coord_pair[0] - start
        end_offset = coord_pair[1] - start
        local_pair = [start_offset, end_offset]
        return local_pair

    def local_to_global_index(self, local_pair, base):
        """
        Converts a pair of local indices (relative to the region of interest) to global indices
        """
        start = local_pair[0] + self.get_coords(base)[0] - 1
        end = self.get_coords(base)[0] + local_pair[1] - 1
        global_pair = [start, end]
        return global_pair

    def get_overlap(self, coord_pair, base):
        """
        Gets the sequence regions that overlap with the region of interest
        :param coord_pair: the coordinates
        :param base: The base of the sequence
        :return: the sequence of the overlap, and the indices of the overlap
        """
        local_pair = self.global_to_local_index(coord_pair, base)
        seq = self.get_sequence(base)

        # Check if local coordinates are in the range of the sequence
        start = max(local_pair[0], 0)
        end = min(local_pair[1] + 1, len(seq))
        if base == 'nucl':
            overlap = (seq[start: end], self.local_to_global_index([start, end], base))
        else:
            overlap = (seq[start - 1: end - 1], self.local_to_global_index([start, end - 1], base))
        return overlap


def set_regions(virus, nt_reference, nt_coords, aa_reference, aa_coords):
    """
    Reads in the start and end coordinates of the genomic regions and associates the region with its coordinates.
    If no coordinate files are specified, set default nucleotide and protein coordinate files.
    :param virus: The organism (HIV or SIV)
    :param nt_reference: The file stream containing the reference nucleotide sequence in read mode
    :param nt_coords: Path to the csv file containing the coordinates of the nucleotide region.
            The file stream has one genomic entry per line and has the following format: region_name,start,end
    :param aa_reference: The file stream containing the reference nucleotide sequence in read mode
    :param aa_coords: Path to the csv file containing the coordinates of the protein region.
            The file stream has one genomic entry per line and has the following format: region_name,start,end
    :return: A list of GenomeRegions
    """

    genome_regions = []

    # Parse nucleotide region coordinates file
    with open(nt_coords, 'r') as nt_handle:
        for nt_line in nt_handle:
            nt_line = nt_line.strip()
            nt_line = nt_line.split(',')
            nucl_coords = [int(nt_line[1]), int(nt_line[2])]
            seq_region = GenomeRegion(nt_line[0])
            seq_region.set_coords(nucl_coords, 'nucl')
            seq_region.set_seq_from_ref(nt_reference, 'nucl')
            seq_region.set_pos_from_cds(virus)
            seq_region.pos_from_gstart = nucl_coords
            seq_region.set_pos_from_aa_start(virus)
            genome_regions.append(seq_region)

    # Parse protein coordinates file
    with open(aa_coords, 'r') as aa_handle:
        prot_names = []
        prot_coords_list = []
        for aa_line in aa_handle:
            aa_line = aa_line.strip()
            aa_line = aa_line.split(',')
            prot_coords_list.append([int(aa_line[1]), int(aa_line[2])])
            prot_names.append(aa_line[0])

        for i in range(len(prot_names)):
            for genome_region in genome_regions:
                if prot_names[i] in genome_region.region_name:
                    genome_region.set_coords(prot_coords_list[i], 'prot')
                    genome_region.set_seq_from_ref(aa_reference, 'prot')

    return genome_regions


def valid_sequence(base, sequence):
    """
    Verifies that input sequence is valid
    :param base: The base of the sequence (nucl or prot)
    :param sequence: A list of lists containing header and sequence pairs
    :raises ValueError: If the sequence is empty or if it contains invalid characters
    :return: <True> If the input sequence uses the correct alphabet, <False> otherwise
    """
    dna_alphabet = 'atgc-*xn'
    aa_alphabet = 'ARDNCEQGHILKMFPSTWYV-*X'

    if not sequence:
        print("Invalid sequence: sequence length is 0\n")
        return False

    for h, s in sequence:
        if not s:
            print("Invalid sequence: sequence length is 0\n")
            return False
        elif base == 'nucl':
            # Nucleotide sequences are converted to lowercase
            s = s.lower()
            if not all(pos in dna_alphabet for pos in s):
                print("Invalid nucleotide sequence:\n{}\n{}\n".format(h, s))
                return False
        else:
            # Amino acid sequences are converted to uppercase
            s = s.upper()
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

    if start_coord <= 0:
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


def get_query(base, query_file, revcomp):
    """
    Gets the query sequence and checks that it is valid
    :param base: The base (nucleotide or protein)
    :param query_file: The file stream containing the query sequence in read mode
    :param revcomp: Option to align the reverse complement of the nucleotide sequence ('y' or 'n')
    :return: A list of lists containing the sequence identifiers and the query sequences
    """
    line = query_file.readline()
    query_file.seek(0)              # Reset pointer to beginning

    # Parse query sequence
    if line.startswith('>'):        # Fasta file
        query = convert_fasta(query_file)

    else:
        query = []
        count = 1
        for line in query_file:
            line = line.strip('\n')
            if len(line) > 0:
                query.extend([["Sequence{}".format(count), "{}".format(line.upper())]])
                count += 1

    if not valid_sequence(base, query):
        sys.exit(0)

    else:
        if revcomp == 'y':
            if base == 'prot':
                print("Invalid option: reverse complement is not available for proteins.")
            else:
                rc_query = reverse_comp(query[0][1])
                header = query[0][0]
                query = [[header, rc_query]]

        return query


def reverse_comp(query_sequence):
    """
    Reverses and complements the query sequence
    :param query_sequence: The query sequence
    :return: The reverse complement of the query sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '*': '*', 'N': 'N', '-': '-'}
    rev_comp = "".join(complement.get(nt, nt) for nt in reversed(query_sequence))
    return rev_comp


def get_ref_nt_seq(ref_seq):
    """
    Converts the reference sequence to a string and checks if the sequence is valid
    :param ref_seq: The path to the reference sequence
    :return reference_sequence: The reference sequence as a list where the first element is the header
                    and the second element is the sequence as a string
    """
    with open(ref_seq, 'r') as ref_handle:
        reference_sequence = convert_fasta(ref_handle)

    if valid_sequence('nucl', reference_sequence):
        return reference_sequence
    else:
        sys.exit(0)


def get_ref_aa_seq(ref_seq):
    """
    Converts the reference sequence to a string and checks if the sequence is valid
    :param ref_seq: The path to the reference sequence
    :return reference_sequence: The reference sequence as a list where the first element is the header
                    and the second element is the sequence as a string
    """
    with open(ref_seq, 'r') as ref_handle:
        reference_sequence = convert_fasta(ref_handle)

    if valid_sequence('prot', reference_sequence):
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
        print("\033[1mAlignment:\033[0m")

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
    pat = re.compile('[A-Z]{2,}')     # Match the aligned region (minimum alignment length is 2)
    coordinates = [(match.start(), match.end()) for match in pat.finditer(alignment[1])]
    return coordinates


def find_matches(virus, base, ref_regions, match_coordinates):
    """
    Finds the genomic regions where the query sequence aligns with the reference sequence
    :param virus: The organism (HIV or SIV)
    :param base: The base of the query sequence
    :param ref_regions: A list of GenomeRegion objects
    :param match_coordinates: A list of indices where the query sequence aligns with the reference sequence
    :return matches: A dictionary of matches, where the key is the region name and the value is the coordinates
    """
    query_regions = []

    for coord in match_coordinates:
        start_aln = coord[0]
        end_aln = coord[1] - 1  # -1 to account for match.end() in regex

        for ref_region in ref_regions:
            if ref_region.region_name != "Complete":
                ov_seq, ov_coord = ref_region.get_overlap([start_aln, end_aln], base)

                if ov_seq:
                    query_region = GenomeRegion(ref_region.region_name)
                    query_region.set_coords(ov_coord, base)
                    query_region.set_sequence(ov_seq, base)
                    query_region.set_pos_from_cds(virus)
                    query_region.pos_from_gstart = [start_aln, end_aln]
                    query_region.set_pos_from_aa_start(virus)

                    if base == 'nucl':
                        set_protein_equivalents(query_region, ref_regions)
                    else:
                        set_nucleotide_equivalents(query_region, ref_regions)
                    query_regions.append(query_region)

    return query_regions


def set_protein_equivalents(query_reg, ref_regions):
    """
    Finds the protein equivalent of the nucleotide sequence
    :param query_reg: the region that aligns with the query sequence
    :param ref_regions: a list of GenomeRegion objects
    :return prot_equiv: the protein equivalent of the nucleotide sequence
    """
    prot_equiv = None
    non_coding = ["5'LTR", "TAR", "3'LTR"]
    for ref_reg in ref_regions:
        if ref_reg.region_name == query_reg.region_name and ref_reg.region_name not in non_coding:
            if ref_reg.codon_aln is not None and query_reg.aa_coords is not None:
                prot_equiv = ref_reg.codon_aln[query_reg.nt_coords[0]: query_reg.nt_coords[1]]
                prot_equiv = re.sub('[-]', '', prot_equiv)
                query_reg.set_aa_seq(prot_equiv)

    return prot_equiv


def set_nucleotide_equivalents(query_reg, ref_regions):
    """
    Finds the nucleotide equivalent of the protein sequence
    :param query_reg: the region that aligns with the query sequence
    :param ref_regions: a list of GenomeRegion objects
    :return nt_equiv: the nucleotide equivalent of the protein sequence
    """
    nt_equiv = None
    for ref_reg in ref_regions:
        if ref_reg.region_name == query_reg.region_name and ref_reg.codon_aln is not None:
            if query_reg.nt_coords is not None:
                query_reg.make_codon_aln()
                regex = re.compile(query_reg.codon_aln)
                coords = regex.search(ref_reg.codon_aln).span()
                nt_equiv = ref_reg.nt_seq[coords[0]: coords[1]]
                query_reg.set_nt_seq(nt_equiv)

    return nt_equiv


def output(query_regions, outfile=None):
    """
    Prints the output
    :param query_regions: A list of GenomeRegions where the query sequence aligns with the reference sequence
    :param outfile: The file stream for the output file in write mode
    """
    if outfile is None:
        print("\n\n\033[1mRegions touched by the query sequence:\033[0m")

        for reg in query_regions:
            if reg.region_name.startswith('5\'LTR'):
                print("\t3'LTR\n")

            print("\nRegion:\t{}".format(reg.region_name))
            print(textwrap.fill(reg.nt_seq, 60))

            if reg.aa_seq is not None:
                print(textwrap.fill(reg.aa_seq, 60))

            print("\n\033[1mRelative Positions:\033[0m\n")

            if len(reg.pos_from_cds) == 2:
                print("\tNucleotide position relative to CDS start: {} --> {}"
                      .format(reg.pos_from_cds[0], reg.pos_from_cds[1]))
            else:
                print("\tNucleotide position relative to CDS start: N/A")

            if reg.pos_from_gstart is not None:
                print("\tNucleotide position relative to genome start: {} --> {}"
                      .format(reg.pos_from_gstart[0] + 1, reg.pos_from_gstart[1] + 1))

            if reg.pos_from_aa_start is not None:
                print("\tAmino acid position relative to protein start: {} --> {}"
                      .format(reg.pos_from_aa_start[0], reg.pos_from_aa_start[1]))

    else:
        outfile.write("\n\nRegions touched by the query sequence:")
        for reg in query_regions:
            if reg.region_name.startswith('5\'LTR'):
                outfile.write("\t3'LTR\n")

            outfile.write("\nRegion:\t{}".format(reg.region_name))
            outfile.write(textwrap.fill(reg.nt_seq, 60))
            outfile.write('\n')

            if reg.aa_seq is not None:
                outfile.write(textwrap.fill(reg.aa_seq, 60))
                outfile.write('\n')

            outfile.write("\n\nRelative Positions:\n")
            if len(reg.pos_from_cds) == 2:
                outfile.write("\tNucleotide position relative to CDS: {} to {}\n"
                              .format(reg.pos_from_cds[0], reg.pos_from_cds[1]))
            else:
                outfile.write("\tNucleotide position relative to CDS start: N/A\n")

            if reg.pos_from_gstart is not None:
                outfile.write("\tNucleotide position relative to start of genome: {} --> {}\n"
                              .format(reg.pos_from_gstart[0] + 1, reg.pos_from_gstart[1]))

            if reg.pos_from_aa_start is not None:
                outfile.write("\tAmino acid position relative to protein start: {} --> {}\n"
                              .format(reg.pos_from_aa_start[0], reg.pos_from_aa_start[1]))


def retrieve(virus, base, ref_regions, region, outfile=None, start_offset=1, end_offset='end'):
    """
    Retrieves a sequence given its coordinates
    :param virus: The organism (HIV or SIV)
    :param base: The base of the sequence (nucleotide or protein)
    :param ref_regions: A list of GenomeRegion objects
    :param region: The genomic region
    :param outfile: The file stream of the output file
    :param start_offset: <option> The start coordinate
    :param end_offset: <option> The end coordinate
    :return: The genomic region defined by the starting and ending coordinates
    """
    for ref_region in ref_regions:
        sequence_range = ref_region.get_coords(base)
        region_start = sequence_range[0]
        region_end = sequence_range[1]

        if start_offset <= region_start:
            start = region_start
        else:
            start = region_start + (start_offset - region_start)

        # If end_coord is greater than the region's end coordinate, set end_coord to region's end coordinate
        if end_offset == 'end' or end_offset > region_end:
            end = region_end
        else:
            end = region_end + (region_end - end_offset)

        retrieved_region = None
        if ref_region.region_name == region:
            s = ref_region.get_sequence(base)
            region_to_retrieve = s[start - 1: end]

            if base == 'nucl':
                retrieved_region = GenomeRegion(region, [start, end], region_to_retrieve,
                                                ref_region.aa_coords, ref_region.aa_seq)
                retrieved_region.set_sequence(region_to_retrieve, 'nucl')
            else:
                retrieved_region = GenomeRegion(region, ref_region.nt_coords, ref_region.nt_seq,
                                                [start, end], region_to_retrieve)
                retrieved_region.set_sequence(region_to_retrieve, 'prot')


            retrieved_region.set_pos_from_cds(virus)
            retrieved_region.pos_from_gstart = retrieved_region.local_to_global_index([start, end], base)
            retrieved_region.set_pos_from_aa_start(virus)

        if retrieved_region:
            if outfile is None:
                print("\033[1mRetrieved sequence: \033[0m\n")
                print("Region:\t{}".format(retrieved_region.region_name))
                print(textwrap.fill(retrieved_region.get_sequence(base)))

                print("\n\033[1mRelative Positions: \033[0m")
                if len(retrieved_region.pos_from_cds) == 2:
                    print("\tNucleotide position relative to CDS start: {} --> {}"
                          .format(retrieved_region.pos_from_cds[0], retrieved_region.pos_from_cds[1]))
                else:
                    print("\tNucleotide position relative to CDS start: N/A")

                if retrieved_region.pos_from_gstart is not None:
                    print("\tNucleotide position relative to genome start: {} --> {}"
                          .format(retrieved_region.pos_from_gstart[0] + 1, retrieved_region.pos_from_gstart[1]))

                if retrieved_region.pos_from_aa_start is not None:
                    print("\tAmino acid position relative to protein start: {} --> {}"
                          .format(retrieved_region.pos_from_aa_start[0], retrieved_region.pos_from_aa_start[1]))

            else:
                outfile.write("\n\nRetrieved sequence:\n")
                outfile.write("\nRegion:\t{}".format(retrieved_region.region_name))
                outfile.write("\n" + textwrap.fill(retrieved_region.get_sequence(base)))

                outfile.write("\n\nRelative Positions: \n")
                if len(retrieved_region.pos_from_cds) == 2:
                    outfile.write("\tNucleotide position relative to CDS: {} to {}\n"
                                  .format(retrieved_region.pos_from_cds[0], retrieved_region.pos_from_cds[1]))
                else:
                    outfile.write("\tNucleotide position relative to CDS start: N/A\n")

                if retrieved_region.pos_from_gstart is not None:
                    outfile.write("\tNucleotide position relative to start of genome: {} --> {}\n"
                                  .format(retrieved_region.pos_from_gstart[0] + 1, retrieved_region.pos_from_gstart[1]))

                if retrieved_region.pos_from_aa_start is not None:
                    outfile.write("\tAmino acid position relative to protein start: {} --> {}\n"
                                  .format(retrieved_region.pos_from_aa_start[0], retrieved_region.pos_from_aa_start[1]))

    return retrieved_region


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

    # Get the reference sequences
    ref_nt_seq = get_ref_nt_seq(ref_nt)
    ref_aa_seq = get_ref_aa_seq(ref_aa)

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
    parser_align.add_argument('query', type=argparse.FileType('r'),
                              help='Path to the file containing the query sequence.')
    parser_align.add_argument('-revcomp', default='n', choices=['y', 'n'],
                              help='Align the reverse complement of the query sequence with the reference sequence')
    parser_align.add_argument('-outfile', type=argparse.FileType('w'),
                              help='Path to the file where results will be written. '
                                   'If no file is specified, the results will be printed to the console')

    # Create subparser for 'retrieve' mode
    parser_retrieve = subparsers.add_parser('retrieve', formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                            description='Retrieve a sequence in HXB2 or SIVmm239 from its coordinates')
    parser_retrieve.add_argument('-start', default=1, type=int,
                                 help='Starting coordinate of the genomic region')
    parser_retrieve.add_argument('-end', default="end",
                                 help='Ending coordinate of the genomic region. Enter an integer or \'end\'.')
    parser_retrieve.add_argument('-region', metavar='', default="Complete", choices=regions,
                                 help='List of case sensitive genomic regions. '
                                      'Allowed regions are: ' + ', '.join(regions))
    parser_retrieve.add_argument('-outfile', metavar='', type=argparse.FileType('w'),
                                 help='File where results will be written. If no file is specified'
                                      'the results will be printed to the console')
    return parser.parse_args()


def main():
    args = parse_args()

    # Ensure proper configuration files are set
    configs = handle_args(args.virus, args.base, args.ref_nt, args.nt_coords, args.ref_aa, args.aa_coords)
    ref_nt_seq = configs[0][0][1]
    ref_aa_seq = configs[1][0][1]
    nt_coords = configs[2]
    aa_coords = configs[3]
    reference_sequence = configs[4]

    # Create genomic region objects based on configuration files
    ref_regions = set_regions(args.virus, ref_nt_seq, nt_coords, ref_aa_seq, aa_coords)

    if args.subcommand == "align":
        query = get_query(args.base, args.query, args.revcomp)
        alignment = sequence_align(query, reference_sequence, args.outfile)

        # Find indices where the query sequence aligns with the reference sequence
        match_coords = get_region_coordinates(alignment[-1])  # Query sequence will be the last item in the list
        query_regions = find_matches(args.virus, args.base, ref_regions, match_coords)

        output(query_regions, args.outfile)

    else:
        valid_in = valid_inputs(args.virus, args.start, args.end, args.region)

        if not valid_in:
            sys.exit(0)
        else:
            retrieve(args.virus, args.base, ref_regions, args.region, args.outfile, args.start, args.end)


if __name__ == '__main__':
    main()
