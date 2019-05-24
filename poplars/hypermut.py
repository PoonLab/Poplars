"""
Python implementation of HYPERMUT 2.0 algorithm
as described in Keele et al. http://www.pnas.org/cgi/content/short/0802203105
Supplementary Materials.
"""

import re
import argparse
from Bio import SeqIO
import scipy.stats as stats


class HyperMut:

    def __init__(self, seq_name=None, num_muts=None, potential_muts=None, ctrl_muts=None,
                 potential_ctrls=None, rate_ratio=None, p_value=None, odds_ratio=None, ctable=None):
        self.seq_name = seq_name
        self.num_muts = num_muts
        self.pot_muts = potential_muts
        self.ctrl_muts = ctrl_muts
        self.potential_ctrls = potential_ctrls
        self.rate_ratio = rate_ratio
        self.p_value = p_value
        self.odds_ratio = odds_ratio
        self.ctable = ctable


mut_pattern = '[AGCT](?=[AG][AGT])'
mut = re.compile(mut_pattern)  # Matches potential mutation sites (GRD)
ctrl_pattern = '[AGCT](?=[CT].|[AG]C)'
ctrl = re.compile(ctrl_pattern)  # Matches control sites (YN|RC)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Classify HIV-1 sequences as being hypermutated '
                    'based on dinucleotide frequencies relative to the consensus.'
    )
    parser.add_argument('fasta', help="<input> path to FASTA file")
    parser.add_argument('--skip', help="<option> number of records to skip")
    return parser.parse_args()


def make_ctable(seq, gees):
    """
    Creates a list of contingency tables for each query sequence
    :param seq: the query sequence
    :param gees: list of positions of G's in the reference sequence
    :return ctable: a contingency table for the query sequence
    """

    # Keele et al. apply Fisher's exact test to a contingency table
    ctable = [
        [0, 0],  # (G->A, G) in GRD context
        [0, 0]  # in other contexts (GYN|GRC is !GRD)
    ]

    # Matches downstream RD motif in query sequence
    ds_motifs = [match.start() for match in mut.finditer(seq)]
    # Matches downstream YN|RC motif in query sequence
    ds_ctrl = [match.start() for match in ctrl.finditer(seq)]

    # Iterate through all positions with G in consensus
    for i in gees:
        nt = seq[i]  # Corresponding base in query sequence

        # G->A mutation in downstream
        row = (0 if nt == 'A' else 1)
        if i in ds_motifs:
            ctable[0][row] += 1

        elif i in ds_ctrl:
            ctable[1][row] += 1

    return ctable


def hypermut(infile):
    # Get reference sequence
    with open(infile) as handle:
        fasta = [list(s) for s in SeqIO.FastaIO.SimpleFastaParser(handle)]

    if args.skip:
        print("skipping first {} records".format(args.skip))
        fasta = fasta[int(args.skip):]

    refseq = fasta[0][1]  # First sequence is the reference sequence

    # Locate GRD motifs in reference sequence
    gees = [i for i, nt in enumerate(refseq) if nt.upper() == 'G']

    results = []

    # Iterate through sequences and identify substitutions from consensus
    for j, seq in enumerate(fasta[1:]):
        header = seq[0]
        ctable= make_ctable(seq[1], gees)

        result = HyperMut()

        # Populate output table using entries from contingency table
        result.seq_name = header
        result.num_muts = ctable[0][0]  # Mutation sites
        result.pot_muts = ctable[0][0] + ctable[0][1]  # Potential mutation sites
        result.ctrl_muts = ctable[1][0]  # Control mutations
        result.potential_ctrls = ctable[1][0] + ctable[1][1]  # Potential controls
        result.rate_ratio = rate_ratio(ctable)  # Rate ratio

        # Calculate one-sided P-value
        odds_ratio, p_value = stats.fisher_exact(ctable, alternative='greater')
        result.p_value = p_value  # Fisher's exact P-value
        result.odds_ratio = odds_ratio
        result.ctable = ctable

        results.append(result)

        # make_data_file(header, seq[1], j + 1, gees, ds_motifs, ds_ctrl)

    return results


def rate_ratio(ctable):
    """
    Calculate rate ratio from contingency table
    @:param ctable: Contingency table
    @:return r_ratio: rate ratio
    """
    muts = float(ctable[0][0])
    pot_muts = float(muts + ctable[0][1])
    ctrl_muts = float(ctable[1][0])
    pot_ctrls = float(ctrl_muts + ctable[1][1])

    try:
        r_ratio = (muts / pot_muts) / (ctrl_muts / pot_ctrls)
    except ZeroDivisionError:
        r_ratio = "undef"

    return r_ratio


def pretty_print(table):
    """
    Print headings in bold
    @:param results table
    """
    print("\033[1m{0:<11} {1:<7} {2:<22} {3:<15} {4:<21} {5:<22} {6:<25} {7:<20}".format
          ("Sequence", "Muts", "Potential Mut Sites", "Control Muts", "Potential Controls",
           "Rate Ratio", "Fisher's Exact P-value", "Odds Ratio\033[0m"))

    # Print values of rows under corresponding headings
    for i in range(len(table)):
        print("{0:<11} {1:<7} {2:<22} {3:<15} {4:<21} {5:<22} {6:<25} {7:<20}".format(*table[i]))

    # Print summary of hypermutated sequences
    print("\n\033[1m\033[4mSummary:\033[0m")
    for i in range(len(table)):
        if table[i][6] <= 0.05:
            print("{} appears to be hypermutated (OR={})".format(table[i][0], table[i][7]))


def make_data_file(header, s, num_seq, gees, ds_motifs, ds_ctrl):
    """
    Writes detailed output of hypermut to a text file
    :param header: the sequence header
    :param s: the sequence
    :param num_seq: number of sequences in the FASTA file (excluding the consensus sequence)
    :param gees: list of position of G's in the consensus sequence
    :param ds_motifs: list of positions that match GRD motif in the query sequence
    :param ds_ctrl: list of positions that match control (YN|RC) motif in the query sequence
    """

    with open("hypermut_output.txt", "a") as output:
        if num_seq == 1:
            output.write("Regex: " + mut_pattern + ", " + ctrl_pattern + "\n")
        output.write("\nSequence Number: " + str(num_seq))
        output.write("\nSequence Name: " + header)
        output.write("\nPos\tMut\n")

        for i in gees:
            nt = s[i]
            if i in ds_motifs:
                if nt == 'A':
                    output.write(str(i + 1) + "\t" + "1\n")
                else:
                    output.write(str(i + 1) + "\t" + "0\n")

        output.write("\nSequence Name: " + header + " control")
        output.write("\nPos\tMut\n")

        for i in gees:
            nt = s[i]
            if i in ds_ctrl:
                output.write("{}\t{}\n".format(
                    str(i + 1),
                    '1' if nt == 'A' else '0'
                ))


if __name__ == '__main__':
    args = parse_args()
    hypermut(args.fasta)
