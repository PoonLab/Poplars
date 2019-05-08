"""
Python implementation of HYPERMUT 2.0 algorithm
as described in Keele et al. http://www.pnas.org/cgi/content/short/0802203105
Supplementary Materials.
"""

import re
import argparse
from Bio import SeqIO
import scipy.stats as stats

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


def hypermut(infile):
    # Get reference sequence
    with open(infile) as handle:
        fasta = [list(s) for s in SeqIO.FastaIO.SimpleFastaParser(handle)]
    if args.skip:
        print ("skipping first {} records".format(args.skip))
        fasta = fasta[int(args.skip):]

    # Open and clear output file
    with open("hypermut_output.txt", "w+") as outfile:
        outfile.truncate(0)

    refseq = fasta[0][1]  # First sequence is the reference sequence

    # Locate GRD motifs in reference sequence
    gees = [i for i, nt in enumerate(refseq) if nt.upper() == 'G']

    results = [[0] * 8 for _ in range(len(fasta) - 1)]

    # Iterate through sequences and identify substitutions from consensus
    for j, seq in enumerate(fasta[1:]):
        header = seq[0]
        s = seq[1]

        # Keele et al. apply Fisher's exact test to a contingency table
        ctable = [
            [0, 0],  # (G->A, G) in GRD context
            [0, 0]  # in other contexts (GYN|GRC is !GRD)
        ]

        # Matches downstream RD motif in query sequence
        ds_motifs = [match.start() for match in mut.finditer(s)]
        # Matches downstream YN|RC motif in query sequence
        ds_ctrl = [match.start() for match in ctrl.finditer(s)]

        # Iterate through all positions with G in consensus
        for i in gees:
            nt = s[i]  # Corresponding base in query sequence

            # G->A mutation in downstream
            row = (0 if nt == 'A' else 1)
            if i in ds_motifs:
                ctable[0][row] += 1

            elif i in ds_ctrl:
                ctable[1][row] += 1

        # Populate output table using entries from contingency table
        results[j][0] = header
        results[j][1] = ctable[0][0]  # Mutation sites
        results[j][2] = ctable[0][0] + ctable[0][1]  # Potential mutation sites
        results[j][3] = ctable[1][0]  # Control mutations
        results[j][4] = ctable[1][0] + ctable[1][1]  # Potential controls
        results[j][5] = rate_ratio(ctable)  # Rate ratio

        # Calculate one-sided P-value
        odds_ratio, p_value = stats.fisher_exact(ctable, alternative='greater')
        results[j][6] = p_value  # Fisher's exact P-value
        results[j][7] = odds_ratio

        make_data_file(header, s, j + 1, gees, ds_motifs, ds_ctrl)

    pretty_print(results)


def rate_ratio(ctable):
    """
    Calculate rate ratio from contingency table
    @param ctable: Contingency table...
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
    # Write hypermut data to a file
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
                    str(i+1), 
                    '1' if nt=='A' else '0'
                ))
                #if nt == 'A':
                #    output.write(str(i + 1) + "\t" + "1\n")
                #else:
                #    output.write(str(i + 1) + "\t" + "0\n")


if __name__ == '__main__':
    args = parse_args()
    hypermut(args.fasta)
