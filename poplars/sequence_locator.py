"""
Python implementation of HIV Sequence Locator from https://www.hiv.lanl.gov
"""

import argparse
import tempfile
import os
from mafft import *


def parse_args():
    """
    Parses command line arguments
    """
    bases = ["nucl", "prot"]

    regions = ["LTR5", "Gag", "Matrix", "p17/p15", "Capsid",
               "p24/p27", "p2", "Nucleocapsid", "p7/p8", "p1",
               "p6", "GagPolTF", "Protease", "RT", "RNase",
               "Integrase", "Vif", "Vpx", "Vpr", "Tat(exons)",
               "Tat(exon1)", "Tat(exon2)", "Rev(exons)",
               "Rev(exon1)", "Rev(exon2)", "Vpu", "Env",
               "gp160", "gp120", "gp41", "Nef", "LTR3"
               ]

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
    parser_align.add_argument('query',
                              help='Path to the query file.',
                              )
    parser_align.add_argument('base',
                              help='Sequence base type. Allowed bases are: ' + ', '.join(bases),
                              metavar='base',
                              choices=bases
                              )
    parser_align.add_argument('virus',
                              help='The reference virus',
                              metavar='',
                              choices=['hiv', 'siv']
                              )
    parser_align.add_argument('--refseq',
                              help='Path to the FASTA file with the reference genome',
                              metavar=''
                              )
    parser_align.add_argument('--outfile',
                              help='File where results will be written. If no file is specified, '
                                   'the results will be written to \'seqeunce-locator-output.fasta\'',
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
                                 description='Sequence base type. Allowed bases are: ' + ', '.join(bases),
                                 metavar='base',
                                 choices=bases
                                 )
    parser_retrieve.add_argument('--ref',
                                 help='FASTA file with the reference genome',
                                 metavar=''
                                 )
    parser_retrieve.add_argument('--first',
                                 default="start",
                                 help='Starting coordinate of the genomic region'
                                 )
    parser_retrieve.add_argument('--last',
                                 default="end",
                                 help='Ending coordinate of the genomic region'
                                 )
    parser_retrieve.add_argument('--region',
                                 help='List of case sensitive genomic regions. '
                                      'Allowed regions are: ' + ', '.join(regions),
                                 metavar='',
                                 choices=regions
                                 )
    parser_retrieve.add_argument('-outfile',
                                 help='File where results will be written. If no file is specified, '
                                      'the results will be written to \'seqeunce-locator-output.fasta\'',
                                 metavar='',
                                 default=os.path.abspath("sequence-locator-output.fasta")
                                 )
    return parser.parse_args()


def valid_sequence(query, base):
    """
    Verifies that input sequence uses the correct output
    :param query: input sequence
    :param base: the base of the sequence (nucl or prot)
    :return valid: <true> if the input sequence uses the correct alphabet; <false> otherwise
    """
    dna_alphabet = 'atgc'
    aa_alphabet = 'ARDNCEQGHILKMFPSTWYV'

    if not query:
        raise ValueError("Invalid sequence: sequence length is 0")

    if base == "nucl":
        # Nucleotide sequences are converted to lowercase
        query = query.lower()
        valid = all(pos in dna_alphabet for pos in query)
        if not valid:
            raise ValueError("Invalid amino acid sequence: {}".format(query))

    else:
        # Amino acid sequences are converted to uppercase
        query = query.upper()
        valid = all(pos in aa_alphabet for pos in query)
        if not valid:
            raise ValueError("Invalid amino acid sequence: {}".format(query))

    return valid


def get_ref_seq(ref_seq_path, virus):
    """
    Converts the reference sequence to a string
    :param ref_seq_path: path to the reference sequence
    :param virus: the reference virus (HIV or SIV)
    :return reference_sequence: the reference sequence as a string
    """

    # TODO: check that it is not a fasta file

    # If no reference sequence is specified, set default reference sequence
    if virus is 'hiv' and ref_seq_path is None:
        ref_seq = os.path.abspath("ref_genomes/K03455.fasta")
        print("HIV reference sequence: {}".format(ref_seq))

    if virus is 'siv' and ref_seq_path is None:
        ref_seq = os.path.abspath("ref_genomes/M33262.fasta")
        print("SIV reference sequence: {}".format(ref_seq))

    # Get reference sequence
    with open(ref_seq_path) as ref_handle:
        reference_sequence = ref_handle.read()

    return reference_sequence


def align(query, ref_seq, outfile=None):
    """
    Sequence locator for 'align mode'
    :param query: The query sequence
    :param ref_seq: <option> The reference sequence. If no sequence is specified and HIV is the reference virus,
                    the reference sequence will be ref_genomes/K03455.fasta. If no sequence is specified and SIV
                    is the reference virus, the reference sequence will be ref_genomes/M33262.fasta
    :param outfile: <option> The file where the output will be written.
                    If no output file is specified, the output will be printed.
    """

    # Read in reference sequence
    with open(ref_seq, "r") as ref_handle:
        reference_sequence = ref_handle.read()
        print("Reference sequence: {}".format(reference_sequence))

    with tempfile.NamedTemporaryFile('w+', delete=True) as handle:
        handle.write(reference_sequence)
        handle.write('\n>query sequence\n{}\n'.format(query))
        print("Inside tempfile block... \n "
              "Reference sequence: {} \n "
              "Query seqeunce: {} \n".format(reference_sequence, query))

        # Path to the temporary query file
        file_path = os.path.join(tempfile.gettempdir(), handle.name)

        raw_output = run_mafft(file_path)

    output = raw_output.decode('utf-8')
    print("Outfile: {}".format(outfile))

    if outfile is not None:
        with open(outfile, "w+") as out_handle:
            out_handle.write(output)

        # Print to console if no outfile is specified
    else:
        # TODO: print output of alignment to console
        print(output)

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
    :param outfile: the file where the output will be written
    :return: return the genomic region defined by the starting and ending coordinates
    """

    return None


def main():
    args = parse_args()

    if args.subcommand == "align":
        query = args.query
        base = args.base
        virus = args.virus
        ref_seq = args.refseq
        outfile = args.outfile

        # Read reference sequence from file
        reference_sequence = get_ref_seq(ref_seq, virus)

        # Get query sequence
        with open(query) as query_handle:
            if base == "nucl":
                query_seq = query_handle.read().lower()
            else:
                query_seq = query_handle.read().upper()
        print("Query file: {}".format(query))
        print("Query sequence: {}".format(query_seq))

        # Check that query and reference sequences use the correct alphabet
        if valid_sequence(base, query) and valid_sequence(base, reference_sequence):
            print("Valid sequences")
            align(query, reference_sequence, outfile)

    else:
        virus = args.virus
        base = args.virus
        ref_seq = args.refseq
        start_coord = args.first
        end_coord = args.last
        region = args.region
        outfile = args.outfile

        # Read reference_sequence from file
        reference_sequence = get_ref_seq(ref_seq, virus)

        # Check that reference sequence uses the correct alphabet
        if valid_sequence(base, reference_sequence):
            print("Valid sequence")
            retrieve(reference_sequence, virus, region, start_coord, end_coord, outfile)


if __name__ == '__main__':
    main()
