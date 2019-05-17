"""
Python implementation of HIV Sequence Locator from https://www.hiv.lanl.gov
"""

import argparse
import subprocess
import tempfile
import os
import logging


from Bio import SeqIO


DNA_ALPHABET = "atgc"
AA_ALPHABET = "ARDNCEQGHILKMFPSTWYV"


def parse_args():
    """
    Parses command line arguments
    """

    bases = ["nucl", "prot"]
    ref_organisms = ["K03455.fasta", "M33262.fasta"]

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
    parser_align.add_argument('file',
                              help='Path to the input file.',
                              )
    parser_align.add_argument('base',
                              help='Sequence base type. Allowed bases are: ' + ', '.join(bases),
                              metavar='base',
                              choices=bases
                              )
    parser_align.add_argument('-reference',
                              type=argparse.FileType('w+'),
                              help='FASTA file with the reference genome of HIV strain HXB2 (Accession: K03455) or '
                                   'SIV strain SIVmm239 (Accession: M33262), or any other reference genome.',
                              metavar=''
                              )

    # Create subparser for 'retrieve' mode
    parser_retrieve = subparsers.add_parser('retrieve',
                                            description='Retrieve a sequence in HXB2 '
                                                        'or SIVmm239 from its coordinates',
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                            )
    parser_retrieve.add_argument('-reference',
                                 type=argparse.FileType('r'),
                                 help='FASTA file with reference genomes of HIV strain HXB2 (Accession: K03455) or SIV '
                                      'strain SIVmm239 (Accession: M33262). Options are: ' + ', '.join(ref_organisms),
                                 default="auto",
                                 metavar='',
                                 )
    parser_retrieve.add_argument('-first',
                                 default="start",
                                 help='Starting coordinate of the genomic region'
                                 )
    parser_retrieve.add_argument('-last',
                                 default="end",
                                 help='Ending coordinate of the genomic region'
                                 )
    parser_retrieve.add_argument('-region',
                                 help='List of case sensitive genomic regions. '
                                      'Allowed regions are: ' + ', '.join(regions),
                                 metavar='',
                                 choices=regions
                                 )

    return parser.parse_args()


def is_fasta(infile):
    """
    Checks whether an input file is a FASTA file
    :param infile: the input file
    :return: <true> if the file is a FASTA file; <false> otherwise
    """

    fasta = SeqIO.parse(infile, "fasta")
    return any(fasta)


def make_fasta(handle):
    fasta = []
    seq = ""
    for line in handle:
        if line.startswith('>'):
            header = line
            if len(seq) > 0:
                fasta.append([header, seq])
                seq = ""
                header = line.strip('>\n')
            else:
                seq += line.strip('\n').upper()

    fasta.append([header, seq])
    return fasta


def valid_sequence(base, query):
    """
    Verifies that input sequence uses the correct output
    :param base: the base of the sequence (nucl or prot)
    :param query: input sequence
    :return: <true> if the input sequence uses the correct alphabet; <false> otherwise
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


def align(ref_seq, query):
    """
    Pairwise alignment of input sequence relative to HIV reference genome
    :param: reference sequence
    :param: query sequence
    :return: output of alignment
    """
    with tempfile.NamedTemporaryFile('w+', delete=False) as handle:
        handle.write(ref_seq)
        handle.write('>query sequence\n{}\n'.format(query))

        # Path for the temporary file
        os.path.join(tempfile.gettempdir(), handle.name)

        # Path to find MAFFT
        script_path = os.path.dirname(os.path.abspath(__file__))
        bin_dir = os.path.join(script_path, '/mafft-linux64/bin/')

        if not os.path.isfile(os.path.join(bin_dir, 'mafft-linux64')):
            logging.error("No file exists.")

        raw_output = subprocess.check_output(['mafft', '--quiet', handle.name])

    output = raw_output.decode('utf-8')
    return output


def sequence_align(ref_seq, base, infile):
    """
    Sequence locator for 'align mode'
    :param ref_seq: reference sequence
    :param base: base of the sequence (nucleotide or protein)
    :param infile: input file containing query sequence
    :return: the alignment of the query sequence with the reference genome
    """

    with open(infile) as handle:
        if base == "nucl":
            query = handle.read().lower()
        else:
            query = handle.read().upper()

        # Check if query uses the correct alphabet
        if valid_sequence(base, query):
            print
            # align(ref_seq, query)


def sequence_retrieve(ref_seq, region, start_coord, end_coord):
    """
    Sequence locator for 'retrieve' mode
    :param ref_seq: reference genome sequence
    :param region: genomic region
    :param start_coord: starting coordinate
    :param end_coord: ending coordinate
    :return: return the genomic region defined by the starting and ending coordinates
    """

    return None


def main():
    args = parse_args()
    ref_seq = args.reference.read()

    if args.subcommand == "align":
        base = args.base
        infile = args.file
        sequence_align(ref_seq, base, infile)

    else:
        region = args.region
        start_coord = args.first
        end_coord = args.last
        sequence_retrieve(ref_seq, region, start_coord, end_coord)


if __name__ == '__main__':
    main()
