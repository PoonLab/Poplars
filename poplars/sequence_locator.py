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


class GenomeRegion:
    """
    Represents information about a genomic region
    """

    def __init__(self, region_name, nt_coords, nt_seq=None, aa_coords=None, aa_seq=None):
        """
        Stores information about each genomic region
        :param region_name: The name of the genomic region
        :param nt_coords: A list containing the start and end coordinates of the nucleotide region. Ex: [(1, 890)]
        :param nt_seq: <option> The nucleotide sequence of the genomic region
        :param aa_coords: <option> A list containing the start and end coordinates of the protein region. Ex: [(8,78)]
        :param aa_seq: <option> The amino acid sequence of the genomic region
        """
        self.region_name = region_name
        self.nt_coords = nt_coords
        self.nt_seq = nt_seq
        self.aa_coords = aa_coords
        self.aa_seq = aa_seq

    def set_nt_coords(self, nt_reference):
        """
        Sets the coordinates of the nucleotide regions
        :param nt_reference: The reference nucleotide sequence
        """
        self.nt_seq = nt_reference[self.nt_coords[0] - 1: self.nt_coords[1]]

    def set_aa_coords(self, aa_reference):
        """
        Sets the coordinates of the amino acid regions
        :param aa_reference: The reference amino acid sequence
        """
        self.aa_seq = aa_reference[self.aa_coords[0] - 1: self.aa_coords[1]]

    def set_aa_seq(self, ref_aa_seq, genome_regions):
        """
        Sets the amino acid sequence of the genomic region
        :param ref_aa_seq: The reference amino acid sequence
        :param genome_regions: A list of GenomeRegion
        """
        for genome_reg in genome_regions:
            for h, seq in ref_aa_seq:
                if h == genome_reg.region_name:
                    genome_reg.aa_seq = seq
                else:
                    genome_reg.aa_seq = None


def read_coordinates(virus, coord_infile=None):
    """
    Reads in the start and end coordinates of the genomic regions and associates the region with its coordinates
    :param virus: The virus (hiv or siv)
    :param coord_infile: The file stream containing the coordinates of the genomic regions.
                            The file stream has one genomic entry per line and has the following format:
                            <region_name>   start,end
    :return: A list of GenomeRegions
    """
    if coord_infile is None:
        # If no coordinate sequence is specified, set default
        if virus == 'hiv':
            coord_infile = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455_genome_coordinates.txt")
        if virus == 'siv':
            coord_infile = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262_genome_coordinates.txt")

    genome_regions = []
    for line in coord_infile:
        line = line.split('\t')
        coord = line[1].split(',')
        coords = [int(coord[0]), int(coord[1])]
        seq_region = GenomeRegion(line[0], coords)
        genome_regions.append(seq_region)

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
    if start_coord <= 0:
        return False

    if type(end_coord) != str and type(end_coord) != int:
        return False

    if type(end_coord) == str:
        if end_coord != "end":
            return False
    else:
        if end_coord <= 0 or start_coord >= end_coord:
            return False

    if virus == 'hiv' and region == 'Vpx' or virus == 'siv' and region == 'Vpu':
        print("Invalid region: {} in {}".format(region, virus))
        return False

    return True


def get_query(base, query_file):
    """
    Gets the query sequence and checks that it is valid
    :param base: The base (nucleotide or protein)
    :param query_file: The file stream containing the query sequence in read mode
    :return: A list of lists containing the sequence identifiers and the query sequences
    """
    line = query_file.readline()
    query_file.seek(0)              # reset pointer to beginning

    # Fasta file
    if line.startswith('>'):
        query = convert_fasta(query_file)

    else:
        count = 1
        query = []
        for line in query_file:
            line = line.strip('\n')
            if len(line) > 0:
                query.append(["Sequence{}".format(count), line.upper()])
                count += 1

    if not valid_sequence(base, query):
        sys.exit(0)
    else:
        return query


def reverse_comp(query_sequence):
    """
    Reverses and complements the query sequence
    :param query_sequence: The query sequence
    :return: The reverse complement of the queyr sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '*': '*', 'N': 'N', '-': '-'}
    rev_comp = "".join(complement.get(nt, nt) for nt in reversed(query_sequence))
    return rev_comp


def get_ref_seq(base, ref_seq):
    """
    Converts the reference sequence to a string and checks if the sequence is valid
    :param base: The base (nucleotide or protein)
    :param ref_seq: The path to the reference sequence
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
    :param reference_sequence: The reference sequence.
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
    print("Coordinates (get_region_coordinates): ", coordinates)
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


def find_genomic_regions(virus, reference_nt_sequence, coordinates):
    """
    Finds the genomic regions where the query sequence aligns with the reference sequence
    :param virus: The reference virus
    :param reference_nt_sequence: The nucleotide reference sequence
    :param coordinates: A list of indices where the query sequence aligns with the reference sequence
    :return regions: A list of lists containing the name of the region and positions of the alignments in this region
    """
    regions = {}
    if virus == 'hiv':

        for align_coord in coordinates:
            start_aln = align_coord[0]
            end_aln = align_coord[1]

            for key in HIV_NT_REGIONS:
                # -1 to account for 1-based indexing
                start_region_coord = HIV_NT_REGIONS[key][0] - 1     # Start coordinate of the genomic feature
                end_region_coord = HIV_NT_REGIONS[key][1] - 1          # End coordinate of the genomic feature

                start_offset = start_aln - start_region_coord
                end_offset = end_aln - end_region_coord

                # Check if aligned region is contained in a genomic region
                if start_region_coord <= start_aln < end_aln <= end_region_coord:
                    common_seq_coords = []
                    if key != "Complete":
                        common_seq_coords.append(reference_nt_sequence[0][1][start_aln + 1: end_aln])
                        common_seq_coords.append([start_region_coord + start_offset, end_region_coord + end_offset])
                        regions[key] = common_seq_coords

                # # Check if aligned region overlaps with genomic regions
                # if start_aln < end_region_coord and end_aln > start_region_coord:
                #     common_seq_coords = []
                #     if key != "Complete":
                #         common_seq_coords.append(reference_nt_sequence[0][1][start_aln + 1: end_aln])
                #         common_seq_coords.append([start_region_coord + start_offset, end_region_coord + end_offset])
                #         regions[key] = common_seq_coords

    else:
        common_seq_coords = []
        for align_coord in coordinates:
            start_aln = align_coord[0]
            end_aln = align_coord[1]

            for key in SIV_NT_REGIONS:
                # -1 to account for 1-based indexing
                start_region_coord = SIV_NT_REGIONS[key][0] - 1  # Start coordinate of the genomic feature
                end_region_coord = SIV_NT_REGIONS[key][1] - 1  # End coordinate of the genomic feature

                start_offset = start_aln - start_region_coord
                end_offset = end_aln - end_region_coord

                # Check if aligned region is contained in a genomic region
                if start_region_coord <= start_aln < end_aln <= end_region_coord:
                    if key != "Complete":
                        common_seq_coords.append(reference_nt_sequence[0][1][start_aln + 1: end_aln])
                        common_seq_coords.append([start_region_coord + start_offset, end_region_coord + end_offset])
                        regions[key] = common_seq_coords

                # Check if aligned region overlaps with genomic regions
                if start_aln < end_region_coord and end_aln > start_region_coord:
                    if key != "Complete":
                        common_seq_coords.append(reference_nt_sequence[0][1][start_aln + 1: end_aln])
                        common_seq_coords.append([start_region_coord + start_offset, end_region_coord + end_offset])
                        regions[key] = common_seq_coords

    for key in regions:
        print(key, regions[key])
        print()
    return regions


def find_aa_regions(matches, viral_prots):
    """
    Finds amino acid regions where the query aligns with the references genome
    :param matches: A list of match objects
    :param viral_prots: A dictionary of the viral protein sequence
    :return regions: A list of lists containing the region and the sequence of the aligned region
    """
    regions = {}
    for match in matches:
        seq = match.group()
        sub_region = []
        for key in viral_prots:
            if seq in viral_prots[key]:
                regions[key] = sub_region.append([key, viral_prots[key], [match.start(), match.end()-1]])
    return regions


def find_relative_pos(virus, base, ref_seq, result, coordinates, regions):
    genomic_region_matches = []
    if virus == 'hiv':
        if base == 'nucl':
            location_in_reference = []
            for reg in regions:
                for seq in reg[regions]:
                    # Find position of match relative to the start of the region
                    regex = re.compile(seq)
                    match_pos = [(match.start(), match.end()) for match in regex.finditer(ref_seq)]
                    for coord in coordinates:
                        aln_match = result[coord[0]:coord[1]]
                        regex = re.compile(aln_match)
                        location_in_reference.append([match.start(), match.end()] for match in regex.finditer(ref_seq))

            for region_name in HIV_NT_REGIONS:
                g_region = find_genomic_regions(virus, ref_seq, coordinates)
                if g_region:
                    genomic_region_matches.append([region_name, g_region])

    else:
        if base == 'nucl':
            location_in_reference = []
            for reg in regions:
                for seq in reg[regions]:
                    # Find position of match relative to the start of the region
                    regex = re.compile(seq)
                    match_pos = [(match.start(), match.end()) for match in regex.finditer(ref_seq)]
                    for coord in coordinates:
                        aln_match = result[coord[0]:coord[1]]
                        regex = re.compile(aln_match)
                        location_in_reference.append([match.start(), match.end()] for match in regex.finditer(ref_seq))

            for region_name in SIV_NT_REGIONS:
                g_region = find_genomic_regions(virus, ref_seq, coordinates)
                if g_region:
                    genomic_region_matches.append([region_name, g_region])

    return genomic_region_matches


def pretty_output(outfile, nt_regions=None, aa_regions=None):
    """
    Prints the output of 'align' mode
    :param outfile: <option> The file stream for the output file in write mode
    :param nt_regions: <option> The nucleotide regions where the query aligned with the reference sequence
    :param aa_regions: <option> The protein regions where the query sequence aligned with th reference sequence
    """
    if outfile is not None:
        if nt_regions is not None:
            outfile.write("Nucleotide regions touched by the query sequence:")
            for nt_region in nt_regions:
                outfile.write("\tRegion:\t{}"
                              "\n\tSequence:\t{}\n".format(nt_region, nt_region[1]))

        if aa_regions is not None:
            outfile.write("Protein regions touched by the query sequence:")
            for aa_region in aa_regions:
                outfile.write("\tRegion:\t{}"
                              "\n\tSequence:\t{}"
                              "\n\tCoordinates:\t{}\n".format(aa_region, aa_region[1], aa_region[2]))

    else:
        if nt_regions is not None:
            print("Nucleotide regions touched by the query sequence:")
            for nt_region in nt_regions:
                print("\tRegion:\t{}"
                      "\n\tSequence:\t{}\n".format(nt_region, nt_region[1]))

        if aa_regions is not None:
            print("Protein regions touched by the query sequence:")
            for aa_region in aa_regions:
                print("\tRegion:\t{}"
                      "\n\tSequence:\t{}"
                      "\n\tCoordinates:\t{}\n".format(aa_region[0], aa_region[1], aa_region[2]))


def retrieve(virus, reference_sequence, region, outfile, start_offset=1, end_offset='end'):
    """
    Retrieves a sequence given its coordinates
    :param virus: The reference virus
    :param reference_sequence: The reference genome sequence
    :param region: The genomic region
    :param outfile: The file stream of the output file
    :param start_offset: <option> The start coordinate
    :param end_offset: <option> The end coordinate
    :return: The genomic region defined by the starting and ending coordinates
    """
    if virus == 'hiv':
        sequence_range = HIV_NT_REGIONS[region]
    else:
        sequence_range = SIV_NT_REGIONS[region]

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


def parse_args():
    """
    Parses command line arguments
    """
    # regions = list(HIV_NT_REGIONS)  # The same genomic region options are provided for HIV and SIV

    parser = argparse.ArgumentParser(
        description='An implementation of the HIV Sequence Locator tool by the Los Alamos National Laboratory.'
                    'This tool aligns a nucleotide or protein sequence relative to HIV or SIV reference genomes; '
                    'or retrieves a sequence in the HXB2 or SIVmm239 reference genome from its coordinates.',
    )
    parser.add_argument('virus', metavar='virus', choices=['hiv', 'siv'],
                        help='The reference virus')
    parser.add_argument('base', metavar='base', choices=['nucl', 'prot'],
                        help='Sequence base type. Allowed bases are \'nucl\' and \'prot\'')
    parser.add_argument('-nt_region_coords', type=argparse.FileType('r'),
                        help='Path to the file containing the coordinates of the nucleotide region.'
                             'The file must be tab-delimited and contain the region name and the start and end '
                             'coordinates separated by a comma (,). Ex: region_name    start,end')
    parser.add_argument('-aa_region_coords', type=argparse.FileType('r'),
                        help='Path to the file containing the coordinates of the amino acid region.'
                             'The file must be tab-delimited and contain the region name and the start and end '
                             'coordinates separated by a comma (,). Ex: region_name    start,end')
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


def handle_args(virus, base, query, nt_coords, aa_coords, revcomp, ref_nt, ref_aa):
    """
    Handles the possible execution paths for the program
    :param virus: The reference virus
    :param base: The base of the reference sequence
    :param query: The file stream containing the query sequence
    :param nt_coords: The file stream containing the coordinates of the nucleotide region in read mode
    :param aa_coords: The file stream containing the coordinates of the protein region in read mode
    :param revcomp: Option to align the reverse complement of the nucleotide sequence ('y' or 'n')
    :param ref_nt: Path to the file containing the reference nucleotide sequence
    :param ref_aa: Path to the file stream containing the reference amino acid sequence
    """

    # Get the query sequence
    query = get_query(base, query)

    # Handle reverse complement option
    if revcomp == 'y':
        if base == 'prot':
            print("Invalid option: reverse complement is not available for proteins.")
    else:
        query = reverse_comp(query[0][1])

    # Set the reference nucleotide and/or protein sequence(s).
    # If user specifies both references files or if neither file is specified, the program
    # will give both amino acid and nucleotide regions that are touched by the query sequence.
    if (ref_nt is None) == (ref_aa is None):

        # If no nucleotide reference sequence is specified, set default reference sequence file
        if ref_nt is None:
            if virus == 'hiv':
                ref_nt = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455.fasta")
            else:
                ref_nt = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262.fasta")

            # If no nucleotide sequence coordinates are specified, set default
            if nt_coords is None:
                if virus == 'hiv':
                    nt_coords = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455_genome_coordinates.txt")
                else:
                    nt_coords = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262_genome_coordinates.txt")

        # If no protein reference sequence is specified, set default reference sequence file
        if ref_aa is None:
            if virus == 'hiv':
                ref_aa = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455-protein.fasta")
            else:
                ref_aa = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262-protein.fasta")

            # If no protein sequence coordinates are specified, set default
            if aa_coords is None:
                if virus == 'hiv':
                    aa_coords = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455_protein_coordinates.txt")
                else:
                    aa_coords = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262_protein_coordinates.txt")

        # Get the reference nucleotide and protein sequences
        ref_nt_seq = get_ref_seq(base, ref_nt)
        ref_aa_seq = get_ref_seq(base, ref_aa)

    elif (ref_nt is None) and (ref_aa is not None):
        ref_aa_seq = get_ref_seq(base, ref_aa)
        ref_nt_seq = ref_nt

    else:
        ref_aa_seq = ref_aa
        ref_nt_seq = get_ref_seq(base, ref_nt)

    # Set the reference sequence to be used in the sequence alignment
    if base == 'nucl':
        reference_sequence = ref_nt_seq
    else:
        reference_sequence = ref_aa_seq

    return query, reference_sequence, ref_nt_seq, ref_aa_seq


def main():
    args = parse_args()

    if args.subcommand == "align":

        sequences = handle_args(args.virus, args.base, args.query, args.nt_region_coords, args.aa_region_coords,
                                args.revcomp, args.ref_nt, args.ref_aa)

        query = sequences[0]
        reference_sequence = sequences[1]
        ref_nt_seq = sequences[2]
        ref_aa_seq = sequences[3]

        # Create genomic region objects based on configuration files
        genome_regions = read_coordinates(args.virus, args.region_coords)
        for genome_region in genome_regions:
            genome_region.set_nt_coords(ref_nt_seq)
            genome_region.set_nt_coords(ref_aa_seq)

        alignment = sequence_align(query, reference_sequence, args.outfile)

        nt_regions, aa_regions = None, None

        if args.base == 'nucl':
            coordinates = get_region_coordinates(alignment[-1])  # Query sequence will be the last item in the list
            nt_regions = find_genomic_regions(args.virus, ref_nt_seq, coordinates)

            if (args.ref_nt is None) == (args.ref_aa is None):
                matches = get_matches(alignment[-1])
                aa_regions = find_aa_regions(matches, viral_prots)

        else:
            matches = get_matches(alignment[-1])
            aa_regions = find_aa_regions(matches, viral_prots)

            if (args.ref_nt is None) == (args.ref_aa is None):
                coordinates = get_region_coordinates(alignment[-1])
                nt_regions = find_genomic_regions(args.virus, ref_nt_seq, coordinates)

        pretty_output(args.outfile, nt_regions, aa_regions)

    else:
        # Read reference_sequence from file
        reference_sequence = get_ref_seq(args.virus, args.base)
        valid_in = valid_inputs(args.virus, args.start, args.end, args.region)

        if not valid_in:
            print("Invalid input")
            sys.exit(0)

        if valid_in:
            retrieve(args.virus, reference_sequence, args.region, args.start, args.end, args.outfile)


if __name__ == '__main__':
    main()
