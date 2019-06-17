"""
Python implementation of HIV Sequence Locator from https://www.hiv.lanl.gov
"""

import argparse
import tempfile
from poplars.mafft import *

HIV_REGIONS = {"Complete":              [(1, 9719)],
               "LTR5":                  [(1, 634)],
               "Gag":                   [(790, 2292)],
               "Matrix(p17/p15)":       [(790, 1185)],
               "Capsid(p24/p27)":       [(1186, 1878)],
               "p2":                    [(1879, 1920)],
               "Nucleocapsid(p7/p8)":   [(1921, 2085)],
               "p1":                    [(2086, 2133)],
               "p6":                    [(2134, 2292)],
               "Pol":                   [(2085, 5096)],
               "GagPolTF":              [(2085, 2252)],
               "Protease":              [(2253, 2549)],
               "RT":                    [(2250, 3869)],
               "RNase":                 [(3870, 4229)],
               "Integrase":             [(4230, 5096)],
               "Vif":                   [(5041, 5619)],
               "Vpx":                   [()],          # Not a region in HIV
               "Vpr":                   [(5559, 5850)],
               "Tat(exons)":            [(5831, 6045), (8379, 8469)],
               "Tat(exon1)":            [(5831, 6045)],
               "Tat(exon2)":            [(8379, 8469)],
               "Rev(exons)":            [(5970, 6045), (8739, 8653)],
               "Rev(exon1)":            [(5970, 6045)],
               "Rev(exon2)":            [(8739, 8653)],
               "Vpu":                   [(6062, 6310)],
               "Env(gp160":             [(6225, 8795)],
               "gp120":                 [(6225, 7757)],
               "gp41":                  [(7758, 8795)],
               "Nef":                   [(8797, 9417)],
               "LTR3":                  [(9086, 9719)]}

SIV_REGIONS = {"Complete":              [(1, 10278)],
               "LTR5":                  [(1, 817)],
               "Gag":                   [(1053, 2585)],
               "Matrix(p17/p15)":       [(1053, 1457)],
               "Capsid(p24/p27)":       [(1458, 2144)],
               "p2":                    [(2145, 2195)],
               "Nucleocapsid(p7/p8)":   [(2196, 2351)],
               "p1":                    [(2352, 2393)],
               "p6":                    [(2394, 2585)],
               "Pol":                   [(2351, 5410)],
               "GagPolTF":              [(2351, 2554)],
               "Protease":              [(2555, 2851)],
               "RT":                    [(2852, 4168)],
               "RNase":                 [(4169, 4528)],
               "Integrase":             [(4529, 5410)],
               "Vif":                   [(5340, 5984)],
               "Vpx":                   [(5812, 6150)],
               "Vpr":                   [(6151, 6456)],
               "Tat(exons)":            [(6302, 6597), (8806, 8902)],
               "Tat(exon1)":            [(6302, 6597)],
               "Tat(exon2)":            [(8806, 8902)],
               "Rev(exons)":            [(6528, 6597), (8806, 9059)],
               "Rev(exon1)":            [(6528, 6597)],
               "Rev(exon2)":            [(8806, 9059)],
               "Vpu":                   [()],           # Not a region in SIV
               "Env(gp160":             [(6604, 9243)],
               "gp120":                 [(6604, 8178)],
               "gp41":                  [(8179, 9243)],
               "Nef":                   [(9077, 9868)],
               "LTR3":                  [(9462, 10278)]}


def parse_args():
    """
    Parses command line arguments
    """

    regions = list(HIV_REGIONS)         # The same genomic region options are provided for HIV and SIV

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
                                 help='Ending coordinate of the genomic region',
                                 type=int
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


def get_ref_seq(virus, ref_seq_path=None):
    """
    Converts the reference sequence to a string
    :param virus: the reference virus (HIV or SIV)
    :param ref_seq_path: <option> path to the reference sequence
    :return reference_sequence: the reference sequence as a string
    """

    # If no reference sequence is specified, set default reference sequence
    if virus == 'hiv' and ref_seq_path is None:
        seq_path = os.path.dirname(os.path.abspath(__file__))
        ref_seq_path = os.path.join(seq_path, "ref_genomes/K03455.fasta")

    if virus == 'siv' and ref_seq_path is None:
        seq_path = os.path.dirname(os.path.abspath(__file__))
        ref_seq_path = os.path.join(seq_path, "ref_genomes/M33262.fasta")

    reference_sequence = ""
    with open(ref_seq_path) as ref_handle:
        for line in ref_handle:
            # Skip fist line if the file is a FASTA file
            if not reference_sequence.startswith(">"):
                reference_sequence = ref_handle.read().replace("\n", "")
                reference_sequence.join(line)

    return reference_sequence




# TODO: write a function that reads output of subprocess and returns the corresponding gene
#  - sequence_retrieve(start_coord=start of longest alignment, end_coord=end of longest alignment))
#  - helper function that does not need ref_seq and needs a new outfile


def retrieve(reference_sequence, virus, region, start_coord, end_coord, outfile=None):
    """
    Sequence locator for 'retrieve' mode
    :param reference_sequence: reference genome sequence
    :param virus: the reference virus
    :param region: genomic region
    :param start_coord: starting coordinate
    :param end_coord: ending coordinate
    :param outfile: <option> The file where the output will be written
                    If no output file is specified, the output will be printed.
    :return: return the genomic region defined by the starting and ending coordinates
    """

    if virus == 'hiv':
        sequence_range = HIV_REGIONS[region]
    else:
        sequence_range = SIV_REGIONS[region]

    # If end_coord is greater than the region's end coordinate, set end_coord to region's end coordinate
    if end_coord > sequence_range[0][1]:
        end_coord = sequence_range[0][1]

    # If start_coord is smaller than the region's start coordinate, set start_coord to region's start coordinate
    if start_coord < sequence_range[0][0]:
        start_coord = sequence_range[0][0]

    region_to_retrieve = reference_sequence[(start_coord-1):(end_coord-1)]      # -1 to account for 0-based indexing

    print(reference_sequence)
    print(region_to_retrieve)


def valid_inputs(virus, start_coord=None, end_coord=None, region=None):
    """
    Checks that the start and end coordinates of the genomic region are valid
    :param virus: the reference virus
    :param start_coord: <option> the starting coordinate of the genomic region
    :param end_coord: <option> the ending coordinate of the genomic region
    :param region: <option> the genomic region
    :return: true if the coordinates are valid, false otherwise
    """

    if start_coord <= 0 or end_coord <= 0:
        print("Invalid coordinate: {} \n"
              "Valid HIV coordinates: 1 to 9719\n"
              "Valid SIV coordinates: 1 to 10278\n".format(start_coord))
        return False

    if start_coord >= end_coord:
        print("Invalid start and end coordinates: ({}, {})\n"
              "The start coordinate must be less than the end coordinate\n".format(start_coord, end_coord))
        return False

    if virus == 'hiv' and end_coord > 9719:
        print("Invalid end coordinate: {}\n"
              "Valid HIV end coordinates: 2 to 9719\n".format(end_coord))
        return False

    if virus == 'siv' and end_coord > 10278:
        print("Invalid end coordinate: {}\n"
              "Valid SIV end coordinates: 2 to 10278\n".format(end_coord))
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
        reference_sequence = get_ref_seq(virus)

        if valid_sequence(base, reference_sequence) and valid_inputs(virus, start_coord, end_coord, region):
            retrieve(reference_sequence, virus, region, start_coord, end_coord, outfile)


if __name__ == '__main__':
    main()
