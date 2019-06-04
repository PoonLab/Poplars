"""
Python implementation of HIV Sequence Locator from https://www.hiv.lanl.gov
"""

import argparse
import subprocess
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
    parser_align.add_argument('infile',
                              help='Path to the input file.',
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

    valid = False
    if not query:
        return valid

    if base == "nucl":
        # Nucleotide sequences are converted to lowercase
        query = query.lower()
        valid = all(pos in dna_alphabet for pos in query)
        if not valid:
            print("Invalid nucleotide sequence: {}".format(query))

    else:
        # Amino acid sequences are converted to uppercase
        query = query.upper()
        valid = all(pos in aa_alphabet for pos in query)
        if not valid:
            print("Invalid amino acid sequence: {}".format(query))

    return valid


def align(query, virus, ref_seq=None, outfile=None):
    """
    Pairwise alignment of input sequence relative to HIV reference genome
    :param query: query sequence
    :param virus: the reference virus
    :param ref_seq: reference sequence
    :param outfile: the file where the output will be written
    """

    if virus is 'hiv' and ref_seq is None:
        ref_seq = os.path.abspath("ref_genomes/K03455.fasta")

    if virus is 'siv' and ref_seq is None:
        ref_seq = os.path.abspath("ref_genomes/M33262.fasta")

    # Create a temporary file containing the reference genome and query sequence
    with open(ref_seq, "r") as ref_handle, tempfile.NamedTemporaryFile('w+', delete=True) as handle:
        reference_sequence = ref_handle.read()
        handle.write(reference_sequence)
        handle.write('\n>query sequence\n{}\n'.format(query))

        # Path to the temporary file
        file_path = os.path.join(tempfile.gettempdir(), handle.name)

        raw_output = run_mafft(file_path)

    output = raw_output.decode('utf-8')

    if outfile is not None:
        with open(outfile, "w+") as out_handle:
            out_handle.write(output)
    # Print to console if no outfile is specified
    else:
        print(output)

# test_path = os.path.abspath("siv-test.fasta")
# raw_output = run_mafft(test_path)
# output = raw_output.decode('utf-8')
# with open("sl-siv-test.txt", "w+") as out_handle:
#     out_handle.write(output)


def sequence_align(infile, base, virus, ref_seq=None, outfile=None):
    """
    Sequence locator for 'align mode'
    :param infile: input file containing query sequence
    :param base: base of the sequence (nucleotide or protein)
    :param virus: the reference virus
    :param ref_seq: the reference sequence
    :param outfile: the file where the output will be written
    """

    with open(infile) as handle:
        if base == "nucl":
            query = handle.read().lower()
        else:
            query = handle.read().upper()

        # Check if query uses the correct alphabet
        if valid_sequence(base, query):
            align(query, virus, ref_seq, outfile)
        else:
            raise ValueError("Invalid input")


def sequence_retrieve(ref_seq, virus, region, start_coord, end_coord, outfile=None):
    """
    Sequence locator for 'retrieve' mode
    :param ref_seq: reference genome sequence
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
        base = args.base
        infile = args.infile
        virus = args.virus
        ref_seq = args.refseq
        outfile = args.outfile

        sequence_align(ref_seq, virus, base, infile, outfile)

    # else:
    #     region = args.region
    #     virus = args.virus
    #     ref_seq = args.refseq
    #     start_coord = args.first
    #     end_coord = args.last
    #     outfile = args.outfile

        # sequence_retrieve(ref_seq, virus, region, start_coord, end_coord, outfile)


if __name__ == '__main__':
    main()
