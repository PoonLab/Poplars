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

HIV_NT_REGIONS = {"5'LTR":               (1, 634),      "5'LTR-R":          (456, 551),
                  "5'LTR-U3":            (1, 455),      "5'LTR-U5":         (552, 634),
                  "TAR":                 (453, 513),    "Gag-Pol":          (790, 5096),
                  "Gag":                 (790, 2292),   "Matrix(p17/p15)":  (790, 1185),
                  "Capsid(p24/p27)":     (1186, 1878),  "p2":               (1879, 1920),
                  "Nucleocapsid(p7/p8)": (1921, 2085),  "p1":               (2086, 2133),
                  "p6":                  (2134, 2292),  "Pol":              (2085, 5096),
                  "GagPolTF":            (2085, 2252),  "Protease":         (2253, 2549),
                  "RT":                  (2250, 3869),  "RNase":            (3870, 4229),
                  "Integrase":           (4230, 5096),  "Vif":              (5041, 5619),
                  "Vpr":                 (5559, 5850),  "Tat(with intron)": (5831, 8469),
                  "Tat(exon1)":          (5831, 6045),  "Tat(exon2)":       (8379, 8469),
                  "Rev(with intron)":    (5970, 8653),  "Rev(exon1)":       (5970, 6045),
                  "Rev(exon2)":          (8739, 8653),  "Vpu":              (6062, 6310),
                  "Env(gp160)":          (6225, 8795),  "V1":               (6615, 6692),
                  "V2":                  (6693, 6812),  "V3":               (7110, 7217),
                  "V4":                  (7377, 7478),  "V5":               (7602, 7634),
                  "RRE":                 (7710, 8061),  "gp120":            (6225, 7757),
                  "gp41":                (7758, 8795),  "Nef":              (8797, 9417),
                  "3'LTR":               (9086, 9719),  "3'LTR-R":          (9541, 9636),
                  "3'LTR-U3":            (9086, 9540),  "3'LTR-U5":         (9637, 9719)}

SIV_NT_REGIONS = {"5'LTR":               (257, 1074),   "5'LTR-R":          (777, 950),
                  "5'LTR-U3":            (257, 776),    "5'LTR-U5":         (951, 1074),
                  "TAR":                 (774, 898),    "Gag-Pol":          (1309, 5666),
                  "Gag":                 (1309, 2842),  "Matrix(p17/p15)":  (1309, 1713),
                  "Capsid(p24/p27)":     (1714, 2400),  "p2":               (2401, 2451),
                  "Nucleocapsid(p7/p8)": (2452, 2607),  "p1":               (2608, 2649),
                  "p6":                  (2650, 2842),  "Pol":              (2607, 5666),
                  "GagPolTF":            (2607, 2810),  "Protease":         (2811, 3107),
                  "RT":                  (3108, 4424),  "RNase":            (4425, 4784),
                  "Integrase":           (4785, 5666),  "Vif":              (5596, 6240),
                  "Vpx":                 (6068, 6406),  "Vpr":              (6407, 6712),
                  "Tat(with intron)":    (6558, 9158),  "Tat(exon1)":       (6558, 6853),
                  "Tat(exon2)":          (9062, 9158),  "Rev(with intron)": (6784, 9315),
                  "Rev(exon1)":          (6784, 6853),  "Rev(exon2)":       (9062, 9315),
                  "Env(gp160)":          (6860, 9499),  "V1":               (7196, 7360),
                  "V2":                  (7364, 7492),  "V3":               (7791, 7892),
                  "V4":                  (8063, 8153),  "V5":               (8273, 8290),
                  "RRE":                 (8380, 8735),  "gp120":            (6860, 8434),
                  "gp41":                (8435, 9499),  "Nef":              (9333, 10124),
                  "3'LTR":               (9719, 10535), "3'LTR-R":          (10235, 10411),
                  "3'LTR-U3":            (9719, 10234), "3'LTR-U5":         (10412, 10535)}


def valid_sequence(base, sequence):
    """
    Verifies that input sequence is valid
    :param base: the base of the sequence (nucl or prot)
    :param sequence: a list of lists containing header and sequence pairs
    :raises ValueError: if the sequence is empty or if it contains invalid characters
    :return: <true> if the input sequence uses the correct alphabet
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
    :param virus: the reference virus
    :param start_coord: the starting coordinate of the genomic region
    :param end_coord: the ending coordinate of the genomic region
    :param region: the genomic region
    :return: true if the coordinates are valid, false otherwise
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
    :param base: the base (nucleotide or protein)
    :param query_file: the file stream containing the query sequence in read mode
    :return: a list of lists containing the sequence identifiers and the query sequences
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
    :param query_sequence: the query sequence
    :return: the reverse complement of the queyr sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '*': '*', 'N': 'N', '-': '-'}
    rev_comp = "".join(complement.get(nt, nt) for nt in reversed(query_sequence))
    return rev_comp


def get_ref_seq(virus, base, ref_seq=None):
    """
    Converts the reference sequence to a string and checks if the sequence is valid
    :param virus: The reference virus (HIV or SIV)
    :param base: The base (nucleotide or protein)
    :param ref_seq: <option> The file stream containing the reference sequence in read mode
                    If no file is specified, the reference genomes for HIV and SIV will be used
    :return reference_sequence: The reference sequence as a string
    """

    if ref_seq is None:
        if base == 'nucl':
            # If no reference sequence is specified, set default reference sequence
            if virus == 'hiv':
                ref_seq = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455.fasta")
            if virus == 'siv':
                ref_seq = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262.fasta")
        else:
            if virus == 'hiv':
                ref_seq = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455-protein.fasta")
            if virus == 'siv':
                ref_seq = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262-protein.fasta")

        with open(ref_seq, 'r') as ref_handle:
            reference_sequence = convert_fasta(ref_handle)
    else:
        reference_sequence = convert_fasta(ref_seq)

    if valid_sequence(base, reference_sequence):
        return reference_sequence
    else:
        sys.exit(0)


def make_aa_dict(ref_aa_seq):
    """
    Creates a dictionary where the keys are the name of the protein and the values are the protein sequence
    :param ref_aa_seq: the file stream containing the reference sequence in read mode
    :return viral_prots: a dictionary of the viral proteins
    """
    viral_prots = {}
    for h, seq in ref_aa_seq:
        viral_prots[h] = seq
    return viral_prots


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
    :param alignment: a list containing the header of the query sequence, and the aligned query sequence
    :return: the positions where the query sequence aligns with the reference sequences(s) with no gaps
    """
    pat = re.compile('[A-Z]{2,}')     # Match the aligned region (minimum alignment length is 2)
    coordinates = [(match.start(), match.end()) for match in pat.finditer(alignment[1])]
    print("Coordinates (get_region_coordinates): ", coordinates)
    return coordinates


def get_matches(alignment):
    """
    Retrieves sequences where the query sequence matches the reference sequence
    :param alignment: a list containing the aligned regions
    :return coordinates: the sequences where the query sequence aligns with the reference sequences(s) with no gaps
    """
    pat = re.compile('[A-Z]{2,}')    # Match the aligned region (minimum alignment length is 2)
    matches = [match for match in pat.finditer(alignment[1])]
    print("Matches (get_matches): ", matches)
    return matches


def find_genomic_regions(virus, reference_nt_sequence, coordinates):
    """
    Finds the genomic regions where the query sequence aligns with the reference sequence
    :param virus: the virus (HIV or SIV)
    :param reference_nt_sequence: the nucleotide reference sequence
    :param coordinates: a list of indices where the query sequence aligns with the reference sequence
    :return regions: a list of lists containing the name of the region and positions of the alignments in this region
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
    :param matches: a list of match objects
    :param viral_prots: a dictionary of the viral protein sequence
    :return regions: a list of lists containing the region and the sequence of the aligned region
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
    :param outfile: <option> the file stream for the output file in write mode
    :param nt_regions: <option> the nucleotide regions where the query aligned with the reference sequence
    :param aa_regions: <option> the protein regions where the query sequence aligned with th reference sequence
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
    :param virus: the reference virus
    :param reference_sequence: reference genome sequence
    :param region: the genomic region
    :param outfile: file stream of the output file
    :param start_offset: the starting coordinate
    :param end_offset: the end coordinate
    :return: return the genomic region defined by the starting and ending coordinates
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
    regions = list(HIV_NT_REGIONS)  # The same genomic region options are provided for HIV and SIV

    parser = argparse.ArgumentParser(
        description='An implementation of the HIV Sequence Locator tool by the Los Alamos National Laboratory.'
                    'This tool aligns a nucleotide or protein sequence relative to HIV or SIV reference genomes; '
                    'or retrieves a sequence in the HXB2 or SIVmm239 reference genome from its coordinates.',
    )
    parser.add_argument('virus', metavar='virus', choices=['hiv', 'siv'],
                        help='The reference virus')
    parser.add_argument('base', metavar='base', choices=['nucl', 'prot'],
                        help='Sequence base type. Allowed bases are \'nucl\' and \'prot\'')
    subparsers = parser.add_subparsers(dest='subcommand')

    # Create subparser for 'align' mode
    parser_align = subparsers.add_parser('align', formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                         description='Align a nucleotide or protein sequence '
                                                     'relative to the HIV or SIV reference genome')

    parser_align.add_argument('query', type=argparse.FileType('r'),
                              help='File containing the query sequence.')
    parser_align.add_argument('-revcomp', default='n', choices=['y', 'n'],
                              help='Align the reverse complement of the query sequence with the reference sequence')
    parser_align.add_argument('-ref_nt', metavar='', type=argparse.FileType('r'),
                              help='File containing the reference nucleotide sequence')
    parser_align.add_argument('-ref_aa', metavar='', type=argparse.FileType('r'),
                              help='File containing the reference amino acid sequence')
    parser_align.add_argument('-outfile', type=argparse.FileType('w'),
                              help='File where results will be written. '
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

    if args.subcommand == "align":
        ref_nt_seq = get_ref_seq(args.virus, args.base, args.ref_nt)
        ref_aa_seq = get_ref_seq(args.virus, args.base, args.ref_aa)
        query = get_query(args.base, args.query)

        # Handle reverse complement option
        if args.revcomp == 'y':
            if args.base == 'prot':
                print("Invalid option: reverse complement is not available for proteins.")
            else:
                query = reverse_comp(query[0][1])

        # Set the reference sequence
        if args.base == 'nucl':
            reference_sequence = ref_nt_seq
        else:
            reference_sequence = ref_aa_seq

        viral_prots = make_aa_dict(ref_aa_seq)
        alignment = sequence_align(query, reference_sequence, args.outfile)

        nt_regions, aa_regions = None, None

        if args.base == 'nucl':
            coordinates = get_region_coordinates(alignment[-1])  # Query sequence will be the last item in the list
            nt_regions = find_genomic_regions(args.virus, ref_nt_seq, coordinates)

            # If user specifies both references files or if neither file is specified, the program
            # will give both amino acid and nucleotide regions that are touched by the query sequence.
            # If only one reference file is specified (eg: a reference nucleotide sequence), the program
            # will give only nucleotide regions that are touched by the query sequence.
            if (args.ref_nt is None) == (args.ref_aa is None):
                matches = get_matches(alignment[-1])
                aa_regions = find_aa_regions(matches, viral_prots)

        else:
            matches = get_matches(alignment[-1])
            aa_regions = find_aa_regions(matches, viral_prots)

            # If user specifies both references files or if neither file is specified, the program
            # will give both amino acid and nucleotide regions that are touched by the query sequence.
            # If only one reference file is specified (eg: a reference protein sequence), the program
            # will give only protein regions that are touched by the query sequence.
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
