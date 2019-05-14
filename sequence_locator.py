"""
Python implementation of HIV Sequence Locator from https://www.hiv.lanl.gov
"""

import argparse
import sys

from Bio import SeqIO
from Bio import AlignIO

DNA_ALPHABET = "atgc"
AA_ALPHABET = "ardnceqghilkmfpstwyv"
# HIV_GENOME = ""
# SIV_GENOME = ""

class ArgParse(object):
    """
    Parses command line arguments
    """

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='Align a nucleotide or protein sequence relative to the HXB2 or SIVmm239 reference strains; '
                        'or retrieve a sequence in HXB2 or SIVmm239 from its coordinates',
        )

        parser.add_argument('mode',
                            help='Align input with reference sequence; or retrieve a sequence from its coordinates',
                            choices=["align", "retrieve"]
                            )
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.mode):
            print('Invalid mode')
            parser.print_help()
            exit(1)
        getattr(self, args.mode)()

    def align(self):
        """
        Parse arguments related to the 'align' mode
        :return: arguments for 'align' mode
        """

        organisms = ["HIV", "SIV"]

        parser = argparse.ArgumentParser(
            description='Align a nucleotide or protein sequence relative to the HXB2 or SIVmm239 reference strains'
        )
        parser.add_argument('in',
                            help='Path to the input file'
                            )
        parser.add_argument('base',
                            help='Sequence base type',
                            choices=["nucl", "prot"]
                            )
        parser.add_argument('-o, --org',
                            help='Reference organism. Allowed organisms are: ' + ', '.join(organisms),
                            metavar='',
                            choices=organisms
                            )
        args = parser.parse_args(sys.argv[2:])
        return args

    def retrieve(self):
        """
        Parse arguments related to the 'retrieve' mode
        :return: arguments for 'retrieve' mode
        """

        regions = ["LTR5", "Gag", "Matrix", "p17/p15", "Capsid",
                   "p24/p27", "p2", "Nucleocapsid", "p7/p8", "p1",
                   "p6", "GagPolTF", "Protease", "RT", "RNase",
                   "Integrase", "Vif", "Vpx", "Vpr", "Tat(exons)",
                   "Tat(exon1)", "Tat(exon2)", "Rev(exons)",
                   "Rev(exon1)", "Rev(exon2)", "Vpu", "Env",
                   "gp160", "gp120", "gp41", "Nef", "LTR3"
                   ]

        parser = argparse.ArgumentParser(
            description='Retrieve a sequence in HXB2 or SIVmm239 from its coordinates',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        parser.add_argument('strain',
                            help='Reference strain',
                            choices=["HXB2", "SIVmm239"]
                            )
        parser.add_argument('-f, --first',
                            default="start",
                            help='Starting coordinate of the genomic region',
                            metavar=''
                            )
        parser.add_argument('-l, --last',
                            default="end",
                            help='Ending coordinate of the genomic region',
                            metavar=''
                            )
        parser.add_argument('-r, --region',
                            choices=regions,
                            help='List of case sensitive genomic regions.'
                                 'Allowed regions are ' + ', '.join(regions),
                            metavar=''
                            )
        args = parser.parse_args(sys.argv[2:])
        return args


def is_fasta(infile):
    """
    Checks whether an input file is a fasta file
    :param infile: the input file
    :return: true if the file is a fasta file; false otherwise
    """

    fasta = SeqIO.parse(infile, "fasta")
    return any(fasta)


def valid_sequence(base, query):
    """
    Verifies that input sequence uses the correct output
    :param base: the base of the sequence (nucl or prot)
    :param query: input sequence
    :return: true if the input sequence uses the correct alphabet; false otherwise
    """

    if base == "nucl":
        for pos in query:
            if pos not in DNA_ALPHABET:
                print("Invalid nucleotide sequence input")
                return False
        return True
    else:
        for pos in query:
            if pos not in AA_ALPHABET:
                print("Invalid amino acid sequence")
                return False
        return True


def seq_locator():
    if sys.argv[1] == "align":
        cmd_args = vars(ArgParse().align())
        organism = cmd_args.values()[0]
        base = cmd_args.values()[1]
        infile = cmd_args.values()[2]

        with open(infile) as handle:
            if not is_fasta(infile):
                query = handle.read().lower()

                # Check if query uses the correct alphabet
                if valid_sequence(base, query):
                    if organism == "HIV":
                        print
                    else:
                        print

    else:
        cmd_args = vars(ArgParse().retrieve())
        strain = cmd_args.values()[0]
        region = cmd_args.values()[1]
        start_coord = cmd_args.values()[2]
        end_coord = cmd_args.values()[3]


if __name__ == '__main__':
    seq_locator()
