"""
Python implementation of HYPERMUT 2.0 algorithm
as described in Keele et al. http://www.pnas.org/cgi/content/short/0802203105
Supplementary Materials.
"""


import argparse
from poplars.common import convert_fasta
import re

import scipy.stats as stats


class MutationInfo:
    """
    Represents mutation information about how a query sequence compares to a reference sequence
    """

    def __init__(self, seq_name=None, num_muts=None, potential_muts=None, ctrl_muts=None,
                 potential_ctrls=None, rate_ratio=None, p_value=None, odds_ratio=None,
                 ctable=None, mut_sites=None, ctrl_sites=None):
        """
        Generate mutation information for a query sequence
        :param seq_name: the name of the sequence
        :param num_muts: the number of mutated sites, compared to the reference sequence
        :param potential_muts: the number of potential mutation sites (GRD)
        :param ctrl_muts: the number of control sites
        :param potential_ctrls: the number of potential control sites (GYN|GRC)
        :param rate_ratio: the rate ratio for the mutation
        :param p_value: the p-value
        :param odds_ratio: the odds ratio
        :param ctable: a contingency table
        :param mut_sites: <dict> a dictionary where the keys are the locations of potential
                            mutation sites (GRD motifs) and the values are 1 if the sequence
                            is mutated; 0 otherwise
        :param ctrl_sites: <dict> a dictionary where the keys are the locations of potential
                            control sites (GYN|GRC) and the values are 1 if the sequence is
                            mutated; 0 otherwise
        """

        self.seq_name = seq_name
        self.num_muts = num_muts
        self.pot_muts = potential_muts
        self.ctrl_muts = ctrl_muts
        self.potential_ctrls = potential_ctrls
        self.rate_ratio = rate_ratio
        self.p_value = p_value
        self.odds_ratio = odds_ratio
        self.ctable = ctable
        self.mut_sites = mut_sites
        self.ctrl_sites = ctrl_sites


def parse_args():
    parser = argparse.ArgumentParser(
        description='Classify HIV-1 sequences as being hypermutated '
                    'based on dinucleotide frequencies relative to the consensus.'
    )
    parser.add_argument('fasta', help="<input> path to FASTA file")
    parser.add_argument('--skip', help="<option> number of records to skip")
    parser.add_argument('--out', help="<option> write output to the specified file")

    return parser.parse_args()


def make_results(seq, gees):
    """
    Creates a list of contingency tables for each query sequence
    :param seq: the query sequence as a tuple of the form (name, sequence)
    :param gees: list of positions of G's in the reference sequence
    :return ctable: a contingency table for the query sequence
    """

    # FIXME: Regexes do not catch all potential mutation/ control sites
    mut = re.compile('[AGCT](?=[AG][AGT])')  # Matches potential mutation sites (GRD)
    ctrl = re.compile('[AGCT](?=[CT].|[AG]C)')  # Matches control sites (YN|RC)

    # Keele et al. apply Fisher's exact test to a contingency table
    ctable = [
        [0, 0],  # (G->A, G) in GRD context
        [0, 0]  # in other contexts (GYN|GRC is !GRD)
    ]

    s = seq[1]  # query sequence
    # Matches downstream RD motif in query sequence
    ds_motifs = [match.start() for match in mut.finditer(s)]
    # Matches downstream YN|RC motif in query sequence
    ds_ctrl = [match.start() for match in ctrl.finditer(s)]

    mut_sites = {}
    ctrl_sites = {}
    # Iterate through all positions with G in consensus
    for i in gees:
        nt = s[i]  # Corresponding base in query sequence

        # G->A mutation in downstream
        row = (0 if nt == 'A' else 1)
        if i in ds_motifs:
            ctable[0][row] += 1
            # Populate dictionary with locations of mutated sites
            mut_sites[(i + 1)] = 1 if nt == 'A' else 0

        elif i in ds_ctrl:
            ctable[1][row] += 1
            # Populate dictionary with locations of control sites
            ctrl_sites[(i + 1)] = 1 if nt == 'A' else 0

    result = MutationInfo()

    # Set attributes of MutationInfo
    result.seq_name = seq[0]  # header
    result.num_muts = ctable[0][0]  # Mutation sites
    result.pot_muts = ctable[0][0] + ctable[0][1]  # Potential mutation sites
    result.ctrl_muts = ctable[1][0]  # Control mutations
    result.potential_ctrls = ctable[1][0] + ctable[1][1]  # Potential controls
    result.rate_ratio = rate_ratio(ctable)  # Rate ratio

    odds_ratio, p_value = stats.fisher_exact(ctable, alternative='greater')  # Calculate one-sided P-value
    result.p_value = p_value  # Fisher's exact P-value
    result.odds_ratio = odds_ratio
    result.ctable = ctable

    result.mut_sites = mut_sites
    result.ctrl_sites = ctrl_sites

    return result


def hypermut(infile, skip=None):
    """
    Determines if a sequence is hypermutated compared to the reference sequence
    :param infile: the input FASTA file
    :param skip: the number of records to skip
    :return results: a list of MutationInfo Objects
    """

    with open(infile) as handle:
        fasta = convert_fasta(handle)

    if skip:
        print("skipping first {} records".format(skip))
        fasta = fasta[int(skip):]

    # Check that sequences are aligned
    length = len(fasta[0][1])
    for h, s in fasta:
        assert length == len(s), "Sequences are not aligned."

    refseq = fasta[0][1]  # First sequence is the reference sequence
    gees = [i for i, nt in enumerate(refseq) if nt.upper() == 'G']  # Locate GRD motifs in reference sequence
    results = []

    # Iterate through sequences and identify substitutions from consensus
    for j, seq in enumerate(fasta[1:]):
        result = make_results(seq, gees)
        results.append(result)

    return results


def rate_ratio(ctable):
    """
    Calculate rate ratio from contingency table
    :param ctable: Contingency table
    :return r_ratio: rate ratio
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


def pretty_print(results):
    """
    Print results
    :param results: a list of MutationInfo objects
    """
    print("\033[1m{0:<11} {1:<7} {2:<22} {3:<15} {4:<21} {5:<22} {6:<25} {7:<20}".format
          ("Sequence", "Muts", "Potential Mut Sites", "Control Muts", "Potential Controls",
           "Rate Ratio", "Fisher's Exact P-value", "Odds Ratio\033[0m"))

    # Print values of rows under corresponding headings
    for result in results:
        print("{0:<11} {1:<7} {2:<22} {3:<15} {4:<21} {5:<22} {6:<25} {7:<20}"
              .format(result.seq_name, result.num_muts, result.pot_muts,
                      result.ctrl_muts, result.potential_ctrls,
                      result.rate_ratio, result.p_value, result.odds_ratio))

    # Print summary of hypermutated sequences
    print("\n\033[1m\033[4mSummary:\033[0m")
    for result in results:
        if result.p_value <= 0.05:
            print("{} appears to be hypermutated (OR={})".format(result.seq_name, result.odds_ratio))


def make_data_file(file_name, mutation_info_list):
    """
    Writes detailed output of hypermut to a text file
    :param file_name: name of the output file
    :param mutation_info_list: list of MutationInfo objects
    """

    with open(file_name, "w+") as output:
        for mutation_info in mutation_info_list:
            output.write("\nSequence Name: {}".format(mutation_info.seq_name))
            output.write("\nPos\tMut\n")
            for key in mutation_info.mut_sites:
                output.write("{}\t{}\n".format(key, mutation_info.mut_sites[key]))

            output.write("\nSequence Name: {} control".format(mutation_info.seq_name))
            output.write("\nPos\tMut\n")
            for key in mutation_info.ctrl_sites:
                output.write("{}\t{}\n".format(key, mutation_info.ctrl_sites[key]))


def main():
    args = parse_args()
    results = hypermut(args.fasta, args.skip)
    if args.out:
        make_data_file(args.out, results)
    pretty_print(results)


if __name__ == '__main__':
    main()
