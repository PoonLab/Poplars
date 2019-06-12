"""
Python implementation of HIV Sequence Locator from https://www.hiv.lanl.gov
HIV and SIV genomic region coordinates are based on the HXB2 and Mac239
annotation resources from https://www.hiv.lanl.gov/content/sequence/HIV/MAP/annotation.html

Note: The first 256 nucleotides of SIVMM239 correspond to the flanking sequence,
and are included in the complete SIV genome (https://www.ncbi.nlm.nih.gov/nucleotide/M33262)
"""

import argparse
import textwrap
import tempfile
from poplars.mafft import *

HIV_REGIONS = {"Complete":              (1, 9719),
               "5'LTR": 	            (1, 634),
               "5'LTR-R": 	            (456, 551),
               "5'LTR-U3": 	            (1, 455),
               "5'LTR-U5": 	            (552, 634),
               "TAR": 	                (453, 513),
               "Gag-Pol": 	            (790, 5096),
               "Gag":                   (790, 2292),
               "Matrix(p17/p15)":       (790, 1185),
               "Capsid(p24/p27)":       (1186, 1878),
               "p2":                    (1879, 1920),
               "Nucleocapsid(p7/p8)":   (1921, 2085),
               "p1":                    (2086, 2133),
               "p6":                    (2134, 2292),
               "Pol":                   (2085, 5096),
               "GagPolTF":              (2085, 2252),
               "Protease":              (2253, 2549),
               "RT":                    (2250, 3869),
               "RNase":                 (3870, 4229),
               "Integrase":             (4230, 5096),
               "Vif":                   (5041, 5619),
               "Vpx":                   (-1, -1),          # Not a region in HIV
               "Vpr":                   (5559, 5850),
               "Tat(with intron)":      (5831, 8469),
               "Tat(exon1)":            (5831, 6045),
               "Tat(exon2)":            (8379, 8469),
               "Rev(with intron)":      (5970, 8653),
               "Rev(exon1)":            (5970, 6045),
               "Rev(exon2)":            (8739, 8653),
               "Vpu":                   (6062, 6310),
               "Env(gp160)":            (6225, 8795),
               "V1": 	                (6615, 6692),
               "V2": 	                (6693, 6812),
               "V3": 	                (7110, 7217),
               "V4": 	                (7377, 7478),
               "V5": 	                (7602, 7634),
               "RRE": 	                (7710, 8061),
               "gp120":                 (6225, 7757),
               "gp41":                  (7758, 8795),
               "Nef":                   (8797, 9417),
               "3'LTR":                 (9086, 9719),
               "3'LTR-R": 	            (9541, 9636),
               "3'LTR-U3": 	            (9086, 9540),
               "3'LTR-U5": 	            (9637, 9719)}


SIV_REGIONS = {"Complete":              (1, 10535),
               "5'LTR":                 (257, 1074),
               "5'LTR-R":               (777, 950),
               "5'LTR-U3":              (257, 776),
               "5'LTR-U5": 	            (951, 1074),
               "TAR": 	                (774, 898),
               "Gag-Pol": 	            (1309, 5666),
               "Gag":                   (1309, 2842),
               "Matrix(p17/p15)":       (1309, 1713),
               "Capsid(p24/p27)":       (1714, 2400),
               "p2":                    (2401, 2451),
               "Nucleocapsid(p7/p8)":   (2452, 2607),
               "p1":                    (2608, 2649),
               "p6":                    (2650, 2842),
               "Pol":                   (2607, 5666),
               "GagPolTF":              (2607, 2810),
               "Protease":              (2811, 3107),
               "RT":                    (3108, 4424),
               "RNase":                 (4425, 4784),
               "Integrase":             (4785, 5666),
               "Vif":                   (5596, 6240),
               "Vpx":                   (6068, 6406),
               "Vpr":                   (6407, 6712),
               "Tat(with intron)":      (6558, 9158),
               "Tat(exon1)":            (6558, 6853),
               "Tat(exon2)":            (9062, 9158),
               "Rev(with intron)":      (6784, 9315),
               "Rev(exon1)":            (6784, 6853),
               "Rev(exon2)":            (9062, 9315),
               "Vpu":                   (-1, -1),           # Not a region in SIV
               "Env(gp160)":            (6860, 9499),
               "V1": 	                (7196, 7360),
               "V2": 	                (7364, 7492),
               "V3": 	                (7791, 7892),
               "V4": 	                (8063, 8153),
               "V5": 	                (8273, 8290),
               "RRE": 	                (8380, 8735),
               "gp120":                 (6860, 8434),
               "gp41":                  (8435, 9499),
               "Nef":                   (9333, 10124),
               "3'LTR":                 (9719, 10535),
               "3'LTR-R": 	            (10235, 10411),
               "3'LTR-U3": 	            (9719, 10234),
               "3'LTR-U5": 	            (10412, 10535)}

# TODO: Make a 'region object to store nucleotide, amino acid, and position information for HIV and SIV.


def parse_args():
    """
    Parses command line arguments
    """

    regions = list(HIV_REGIONS)        # The same genomic region options are provided for HIV and SIV

    parser = argparse.ArgumentParser(
        description='Aligns a nucleotide or protein sequence relative to the HIV or SIV reference genomes; '
                    'or retrieves a sequence in the HXB2 or SIVmm239 reference genome from its coordinates',
    )
    subparsers = parser.add_subparsers(dest='subcommand')

    # Create subparser for 'align' mode
    parser_align = subparsers.add_parser('align',
                                         description='Align a nucleotide or protein sequence '
                                                     'relative to the HIV or SIV reference genome',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_align.add_argument('virus',
                              help='The reference virus',
                              metavar='',
                              choices=['hiv', 'siv']
                              )
    parser_align.add_argument('base',
                              help='Sequence base type. '
                                   'Allowed bases are: ' + ', '.join("\'nucl\' and \'prot\'"),
                              metavar='base',
                              choices=['nucl', 'prot']
                              )
    parser_align.add_argument('query',
                              help='Path to the query file.',
                              )
    parser_align.add_argument('-ref',
                              help='Path to the FASTA file with the reference genome',
                              metavar=''
                              )
    parser_align.add_argument('-outfile',
                              help='File where results will be written. If no file is specified, '
                                   'the results will be written to \'sequence-locator-output.fasta\'',
                              metavar='',
                              default=os.path.abspath("sequence-locator-output.fasta")
                              )

    # Create subparser for 'retrieve' mode
    parser_retrieve = subparsers.add_parser('retrieve',
                                            description='Retrieve a sequence in HXB2 '
                                                        'or SIVmm239 from its coordinates',
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                            )
    parser_retrieve.add_argument('virus',
                                 help='The reference virus',
                                 metavar='',
                                 choices=['hiv', 'siv']
                                 )
    parser_retrieve.add_argument('base',
                                 help='Sequence base type. '
                                      'Allowed bases are: ' + ', '.join("\'nucl\' and \'prot\'"),
                                 metavar='base',
                                 choices=['nucl', 'prot']
                                 )
    parser_retrieve.add_argument('-start',
                                 help='Starting coordinate of the genomic region',
                                 type=int
                                 )
    parser_retrieve.add_argument('-end',
                                 help='Ending coordinate of the genomic region. Enter an integer or \'end\'.',
                                 default="end"
                                 )
    parser_retrieve.add_argument('-region',
                                 help='List of case sensitive genomic regions. '
                                      'Allowed regions are: ' + ', '.join(regions),
                                 metavar='',
                                 choices=regions
                                 )
    parser_retrieve.add_argument('-outfile',
                                 help='File where results will be written. If no file is specified, '
                                      'the results will be written to \'sequence-locator-output.fasta\'',
                                 metavar='',
                                 default=os.path.abspath("sequence-locator-output.fasta")
                                 )
    return parser.parse_args()


def valid_sequence(base, query):
    """
    Verifies that input sequence uses the correct output
    :param base: the base of the sequence (nucl or prot)
    :param query: input sequence
    :raises ValueError: if the sequence is empty or if it contains invalid characters
    :return: <true> if the input sequence uses the correct alphabet
    """
    dna_alphabet = 'atgc'
    aa_alphabet = 'ARDNCEQGHILKMFPSTWYV'

    if not query:
        print("Invalid sequence: sequence length is 0\n")
        return False

    if base == 'nucl':
        # Nucleotide sequences are converted to lowercase
        query = query.lower()
        valid = all(pos in dna_alphabet for pos in query)
        if not valid:
            print("Invalid nucleotide sequence: {}\n".format(query))
            return False

    else:
        # Amino acid sequences are converted to uppercase
        query = query.upper()
        valid = all(pos in aa_alphabet for pos in query)
        if not valid:
            print("Invalid amino acid sequence: {}\n".format(query))
            return False

    return True


def get_ref_seq(virus, base, ref_seq_path=None):
    """
    Converts the reference sequence to a string
    :param virus: the reference virus (HIV or SIV)
    :param ref_seq_path: <option> path to the reference sequence
    :return reference_sequence: the reference sequence as a string
    """

    if base == 'nucl':
        # If no reference sequence is specified, set default reference sequence
        if virus == 'hiv' and ref_seq_path is None:
            seq_path = os.path.dirname(os.path.abspath(__file__))
            ref_seq_path = os.path.join(seq_path, "ref_genomes/K03455.fasta")

        if virus == 'siv' and ref_seq_path is None:
            seq_path = os.path.dirname(os.path.abspath(__file__))
            ref_seq_path = os.path.join(seq_path, "ref_genomes/M33262.fasta")

    # else:
    #     if virus == 'hiv' and ref_seq_path is None:
    #         seq_path = os.path.dirname(os.path.abspath(__file__))
    #         ref_seq_path = os.path.join(seq_path, "ref_genomes/")

    reference_sequence = ""
    with open(ref_seq_path) as ref_handle:
        for line in ref_handle:
            # Skip fist line if the file is a FASTA file
            if not reference_sequence.startswith(">"):
                reference_sequence = ref_handle.read().replace("\n", "")
                reference_sequence.join(line)

    return reference_sequence


def align(query_sequence, reference_sequence, outfile=None):
    """
    Sequence locator for 'align mode'
    :param query_sequence: The query sequence
    :param reference_sequence: The reference sequence.
    :param outfile: <option> The file where the output will be written.
                    If no output file is specified, the output will be printed.
    """

    with tempfile.NamedTemporaryFile('w+', delete=False) as handle:
        handle.write('>reference sequence\n')
        handle.write(reference_sequence)
        handle.write('\n>query sequence\n{}\n'.format(query_sequence))
        handle.seek(0)        # Move position back to allow subprocess to use file

        # Path to the temporary query file for MAFFT
        file_path = os.path.join(tempfile.gettempdir(), handle.name)

        raw_output = run_mafft(file_path)

    output = raw_output.decode('utf-8')

    if outfile is not None:
        with open(outfile, "w+") as out_handle:
            out_handle.write(output)

    # TODO: Print to console if no outfile is specified
    else:
        print(output)

# TODO: write a function that reads output of subprocess and returns the corresponding gene
#  - sequence_retrieve(start_coord=start of longest alignment, end_coord=end of longest alignment))
#  - helper function that does not need ref_seq and needs a new outfile


def retrieve(reference_sequence, virus, base, region=None, start_coord=None, end_coord=None):
    """
    Sequence locator for 'retrieve' mode
    :param reference_sequence: reference genome sequence
    :param virus: the reference virus
    :param base: the base (nucleotide or protein)
    :param region: genomic region
    :param start_coord: starting coordinate
    :param end_coord: ending coordinate
    :return: return the genomic region defined by the starting and ending coordinates
    """

    if region is None:
        region = "Complete"

    if virus == 'hiv':
        sequence_range = HIV_REGIONS[region]
    else:
        sequence_range = SIV_REGIONS[region]

    if start_coord is None:
        start_coord = 1

    # If start_coord is smaller than the region's start coordinate, set start_coord to region's start coordinate
    if start_coord < sequence_range[0]:
        start_coord = sequence_range[0]

    # If end_coord is greater than the region's end coordinate, set end_coord to region's end coordinate
    if end_coord == "end" or end_coord > sequence_range[1]:
        end_coord = sequence_range[1]

    region_to_retrieve = reference_sequence[(start_coord-1):end_coord]      # -1 to account for 0-based indexing

    print("\033[1mRetrieved sequence: \033[0m")
    print(textwrap.fill(region_to_retrieve, 50))


def valid_inputs(virus, start_coord=None, end_coord=None, region=None):
    """
    Checks that the start and end coordinates of the genomic region are valid
    :param virus: the reference virus
    :param start_coord: <option> the starting coordinate of the genomic region
    :param end_coord: <option> the ending coordinate of the genomic region
    :param region: <option> the genomic region
    :return: true if the coordinates are valid, false otherwise
    """

    if type(end_coord) != str and type(end_coord) != int:
        return False

    if type(end_coord) == str:
        if end_coord != "end":
            return False
        if start_coord <= 0:
            return False

    else:
        if start_coord <= 0 or end_coord <= 0:
            return False
        if start_coord >= end_coord:
            return False

    if virus == 'hiv' and region == 'Vpx' or virus == 'siv' and region == 'Vpu':
        print("Invalid region: {} in {}".format(region, virus))
        return False

    return True


def main():
    args = parse_args()

    if args.subcommand == "align":
        virus = args.virus
        base = args.base
        query_file = args.query
        ref_seq_path = args.ref
        outfile = args.outfile

        reference_sequence = get_ref_seq(virus, ref_seq_path)

        # Get query sequence
        with open(query_file) as query_handle:
            if base == 'nucl':
                query_sequence = query_handle.read().lower()
            else:
                query_sequence = query_handle.read().upper()

        if valid_sequence(base, query_sequence) and valid_sequence(base, reference_sequence):
            align(query_sequence, reference_sequence, outfile)

    else:
        virus = args.virus
        base = args.base
        start_coord = args.start
        end_coord = args.end
        region = args.region
        outfile = args.outfile

        # Read reference_sequence from file
        reference_sequence = get_ref_seq(virus, base)

        valid_seq = valid_sequence(base, reference_sequence)
        valid_in = valid_inputs(virus, start_coord, end_coord, region)

        if not valid_seq:
            print("Invalid sequence: {}".format(reference_sequence))

        if not valid_in:
            print("Invalid input")

        if valid_seq and valid_in:
            retrieve(reference_sequence, virus, region, start_coord, end_coord)


if __name__ == '__main__':
    main()
