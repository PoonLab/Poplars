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
from math import ceil
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
        :param aa_coords: <option> A list containing the start and end coordinates of the protein region. Ex: [(8,78)]
        :param aa_seq: <option> The amino acid sequence of the genomic region
        """
        self.region_name = region_name
        self.nt_coords = nt_coords
        self.nt_seq = nt_seq
        self.aa_coords = aa_coords
        self.aa_seq = aa_seq

    def set_nt_seq(self, nt_reference):
        """
        Sets the sequence for the nucleotide region
        :param nt_reference: The reference nucleotide sequence
        """
        self.nt_seq = nt_reference[self.nt_coords[0] - 1: self.nt_coords[1]]

    def set_aa_seq(self, aa_reference):
        """
        Sets the sequence for the amino acid region
        :param aa_reference: The reference amino acid sequence
        """
        self.aa_seq = aa_reference[self.aa_coords[0] - 1: self.aa_coords[1]]

    def set_aa_coords(self, aa_coords):
        """
        Sets the coordinates of the protein region
        :param aa_coords:  the coordinates
        """
        self.aa_coords = aa_coords

    def set_nt_coords(self, nt_coords):
        """
        Sets the coordinates of the genome region
        :param nt_coords: the coordinates
        :return:
        """
        self.nt_coords = nt_coords


def set_regions(nt_reference, nt_coords, aa_reference, aa_coords):
    """
    Reads in the start and end coordinates of the genomic regions and associates the region with its coordinates.
    If no coordinate files are specified, set default nucleotide and protein coordinate files.
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
            seq_region.set_nt_coords(nucl_coords)
            seq_region.set_nt_seq(nt_reference)
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
                    genome_region.set_aa_coords(prot_coords_list[i])
                    genome_region.set_aa_seq(aa_reference)

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

    if type(start_coord) != str:
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


# TODO: Check for LTR3 if LTR5 match


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
        else:
            print(">{}".format(h))

        # Ouput 60 characters per line
        seq_lines = [s[i:i + 60] for i in range(0, len(s), 60)]
        for line in seq_lines:
            if outfile is not None:
                outfile.write('{}\n'.format(line))
            else:
                print("{}".format(line))

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


def get_matches(alignment):
    """
    Retrieves sequences where the query sequence matches the reference sequence
    :param alignment: A list containing the aligned regions
    :return coordinates: The sequences where the query sequence aligns with the reference sequences(s) with no gaps
    """
    pat = re.compile('[A-Z]{2,}')    # Match the aligned region (minimum alignment length is 2)
    matches = [match for match in pat.finditer(alignment[1])]
    print("Matches (get_matches): ", matches)
    return matches


def find_nt_regions(genome_regions, match_coordinates):
    """
    Finds the genomic regions where the query sequence aligns with the reference sequence
    :param genome_regions: A list of GenomeRegion objects
    :param match_coordinates: A list of indices where the query sequence aligns with the reference sequence
    :return matches: A dictionary of matches, where the key is the region name and the value is the coordinates
    """
    matches = {}

    for coord in match_coordinates:
        start_aln = coord[0]
        end_aln = coord[1]

        for genome_region in genome_regions:
            if genome_region.nt_coords is not None and genome_region.nt_seq is not None:

                start_nt = genome_region.nt_coords[0] - 1       # Start coordinate of the genomic region
                end_nt = genome_region.nt_coords[1] - 1         # End coordinate of the genomic region
                start_offset = start_aln - start_nt
                end_offset = end_aln - end_nt

                # Check if aligned region is contained in a genomic region
                if start_nt <= start_aln < end_aln <= end_nt:
                    relative_coords = []
                    if genome_region.region_name != "Complete":
                        relative_coords.append(genome_region.nt_seq[start_aln: end_aln])
                        relative_coords.append([start_nt + start_offset, end_nt + end_offset])
                        matches[genome_region.region_name] = relative_coords

    return matches


def find_aa_regions(genome_regions, match_coordinates):

    matches = {}
    for coord in match_coordinates:
        start_aln = coord[0]
        end_aln = coord[1]

        for genome_region in genome_regions:
            if genome_region.aa_coords is not None and genome_region.aa_seq is not None:

                start_aa = genome_region.aa_coords[0] - 1  # Start coordinate of the genomic region
                end_aa = genome_region.aa_coords[1] - 1  # End coordinate of the genomic region
                start_offset = start_aln - start_aa
                end_offset = end_aln - end_aa

                # Check if aligned region is contained in a genomic region
                if start_aa <= start_aln < end_aln <= end_aa:
                    relative_coords = []
                    if genome_region.region_name != "Complete":
                        relative_coords.append(genome_region.aa_seq[start_aln: end_aln])
                        relative_coords.append([start_aa + start_offset, end_aa + end_offset])
                        matches[genome_region.region_name] = relative_coords

    return matches


def find_adjacent_prots(matches, genome_regions):
    """
    Finds amino acid regions where the query aligns with the references genome
    :param matches: A list of match objects
    :param genome_regions: A list of GenomeRegion objects
    :return regions: A dictionary containing the region name and the sequence of the aligned region
    """
    regions = {}
    for match in matches:
        seq = match.group()
        sub_region = []
        for genome_region in genome_regions:
            if genome_region.aa_seq is not None:
                if seq in genome_region.aa_coords:
                    regions[genome_region.region_name] = sub_region.append(
                        [genome_region.aa_coords, genome_region.aa_seq, [match.start(), match.end() - 1]])
    return regions


def find_adjacent_nucls(matches, genome_regions):
    """
    Finds nculeotide regions where the query aligns with the references genome
    :param matches: A list of match objects
    :param genome_regions: A list of GenomeRegion objects
    :return regions: A dictionary containing the region name and the sequence of the aligned region
    """
    regions = {}
    for match in matches:
        seq = match.group()
        sub_region = []
        for genome_region in genome_regions:
            if genome_region.nt_seq is not None:
                if seq in genome_region.nt_coords:
                    regions[genome_region.region_name] = sub_region.append(
                        [genome_region.nt_coords, genome_region.nt_seq, [match.start(), match.end() - 1]])
    return regions


def find_relative_pos(genome_regions, virus, base, ref_seq, result, match_coordinates, regions):
    genomic_region_matches = []
    if virus == 'hiv':
        if base == 'nucl':
            location_in_reference = []
            for reg in regions:
                for seq in reg[regions]:
                    # Find position of match relative to the start of the region
                    regex = re.compile(seq)
                    match_pos = [(match.start(), match.end()) for match in regex.finditer(ref_seq)]
                    for coord in match_coordinates:
                        aln_match = result[coord[0]:coord[1]]
                        regex = re.compile(aln_match)
                        location_in_reference.append([match.start(), match.end()] for match in regex.finditer(ref_seq))

            for genome_region in genome_regions:
                g_region = find_nt_regions(genome_regions, match_coordinates)
                if g_region:
                    genomic_region_matches.append([genome_region.region_name, g_region])

    else:
        if base == 'nucl':
            location_in_reference = []
            for reg in regions:
                for seq in reg[regions]:
                    # Find position of match relative to the start of the region
                    regex = re.compile(seq)
                    match_pos = [(match.start(), match.end()) for match in regex.finditer(ref_seq)]
                    for coord in match_coordinates:
                        aln_match = result[coord[0]:coord[1]]
                        regex = re.compile(aln_match)
                        location_in_reference.append([match.start(), match.end()] for match in regex.finditer(ref_seq))

            for genome_region in genome_regions:
                g_region = find_aa_regions(genome_regions, match_coordinates)
                if g_region:
                    genomic_region_matches.append([genome_region.region_name, g_region])

    return genomic_region_matches


def print_nt_output(outfile, matches, associated_prots=None):
    """
    Prints the output of 'align' mode for a nucleotide query
    :param outfile: The file stream for the output file in write mode
    :param matches: The nucleotide regions where the query aligns with the reference sequence
    :param associated_prots: The protein regions where the query sequence aligned with the reference sequence
    """

    if len(matches) > 0:
        if outfile is not None:
            outfile.write("\n\nNucleotide regions touched by the query sequence:")
            for key in matches:
                outfile.write("\nRegion:\t{}".format(key))
                outfile.write(matches[key][0].textwrap(60))
                outfile.write('\n')

        else:
            print("\n\nNucleotide regions touched by the query sequence:")
            for key in matches:
                print("\nRegion:\t{}".format(key))
                print(textwrap.fill(matches[key][0], 60))

    if associated_prots is not None:
        print_aa_output(outfile, matches)


def print_aa_output(outfile, matches, associated_nucls=None):
    """
    Prints the output of 'align' mode for a protein query
    :param outfile: The file stream for the output file in write mode
    :param matches: The protein regions where the query aligns with the reference sequence
    :param associated_nucls: The nculeotide regions where the query sequence aligned with the reference sequence
    """

    if len(matches) > 0:
        if outfile is not None:
            outfile.write("\n\nProtein regions touched by the query sequence:")
            for key in matches:
                outfile.write("\nRegion:\t{}\n".format(key))
                outfile.write(textwrap.fill(matches[key][0], 60))
                outfile.write('\n')
                outfile.write("\nCoordinates: [{}, {}]\n".format(matches[key][1][0] + 1, matches[key][1][1]))

        else:
            print("\n\nProtein regions touched by the query sequence:")
            for key in matches:
                print("\nRegion:\t{}".format(key))
                print(textwrap.fill(matches[key][0], 60))
                print("Coordinates: [{}, {}]\n".format(matches[key][1][0] + 1, matches[key][1][1]))

    if associated_nucls is not None:
        print(outfile, matches)


def retrieve(base, genome_regions, reference_sequence, region, outfile, start_offset=1, end_offset='end'):
    """
    Retrieves a sequence given its coordinates
    :param base: The base of the sequence (nucleotide or protein)
    :param genome_regions: A list of GenomeRegion objects
    :param reference_sequence: The reference genome sequence
    :param region: The genomic region
    :param outfile: The file stream of the output file
    :param start_offset: <option> The start coordinate
    :param end_offset: <option> The end coordinate
    :return: The genomic region defined by the starting and ending coordinates
    """

    region_to_retrieve = ""

    for genome_region in genome_regions:

        if base == 'nucl':
            sequence_range = genome_region.nt_coords
        else:
            sequence_range = genome_region.aa_coords

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

        region_to_retrieve = reference_sequence[0][1][start-1:end]          # -1 to account for 0-based indexing

        if outfile is None:
            print("\033[1mRetrieved sequence: \033[0m")
            # Print 60 characters per line
            seq_lines = [region_to_retrieve[i:i + 60] for i in range(0, len(region_to_retrieve), 60)]
            for line in seq_lines:
                print("{}".format(line))

        else:
            outfile.write("Retrieved sequence: {}\n".format(region_to_retrieve))

        output_relative_positions(start, end, sequence_range, region, outfile)

    return region_to_retrieve


def output_relative_positions(start, end, sequence_range, region, outfile):

    if outfile is None:
        print("\n\033[1mNucleotide position relative to CDS start: \033[0m{} --> {}\n"
              .format(start - 790, end + 1 - 790))      # +1 to account for 0-based indexing
        print("\033[1mNucleotide position relative to query sequence start: \033[0m{} --> {}\n"
              .format(1, (end-start + 1)))
        print("\033[1mNucleotide position relative to genome start: \033[0m{} --> {}\n"
              .format(start, end + 1))
        print("\033[1mAmino acid position relative to protein start: \033[0m{} --> {}"
              .format(ceil((start + 1 - sequence_range[0])/3), (ceil((end - start)/3))))

    else:
        outfile.write("\nNucleotide position relative to CDS: {} to {}\n"
                      .format(start, end))
        outfile.write("\nNucleotide position relative to query sequence start: {} --> {}\n"
                      .format(1, (start - end + 1)))
        outfile.write("\nNucleotide position relative to start of genome: {} --> {}\n"
                      .format((region[0] + start), (region[1] + end)))
        outfile.write("\nAmino acid position relative to protein start: {} --> {}"
                      .format(ceil(start / 3), ceil(end / 3)))


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
    subparsers = parser.add_subparsers(dest='subcommand')

    # Create subparser for 'align' mode
    parser_align = subparsers.add_parser('align', formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                         description='Align a nucleotide or protein sequence '
                                                     'relative to the HIV or SIV reference genome')

    parser_align.add_argument('query', type=argparse.FileType('r'),
                              help='Path to the file containing the query sequence.')
    parser_align.add_argument('-revcomp', default='n', choices=['y', 'n'],
                              help='Align the reverse complement of the query sequence with the reference sequence')
    parser_align.add_argument('-ref_nt', metavar='',
                              help='Path to the file containing the reference nucleotide sequence')
    parser_align.add_argument('-ref_aa', metavar='',
                              help='Path to the file containing the reference amino acid sequence')
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
    genome_regions = set_regions(ref_nt_seq, nt_coords, ref_aa_seq, aa_coords)

    if args.subcommand == "align":

        query = get_query(args.base, args.query, args.revcomp)
        alignment = sequence_align(query, reference_sequence, args.outfile)

        match_coords = get_region_coordinates(alignment[-1])  # Query sequence will be the last item in the list

        if args.base == 'nucl':
            matches = find_nt_regions(genome_regions, match_coords)
            associated_prots = find_adjacent_prots(matches, genome_regions)     # Find associated protein sequences
            print_nt_output(args.outfile, matches, associated_prots)

        else:
            matches = find_aa_regions(genome_regions, match_coords)
            associated_nucls = find_adjacent_nucls(matches, genome_regions)    # Find associated nucleotide sequences
            print_aa_output(args.output, matches, associated_nucls)

    else:
        valid_in = valid_inputs(args.virus, args.start, args.end, args.region)

        if not valid_in:
            sys.exit(0)
        else:
            retrieve(args.virus, reference_sequence, args.region, args.start, args.end, args.outfile)


if __name__ == '__main__':
    main()
