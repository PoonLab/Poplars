"""
Python implementation of HYPERMUT 2.0 algorithm
as described in Keele et al. http://www.pnas.org/cgi/content/short/0802203105
Supplementary Materials.
"""

import argparse
import re
import csv
import scipy.stats as stats
from poplars.common import *


def make_results(seq, gees):
    """
    Creates a list of contingency tables for each query sequence
    :param seq: the query sequence as a tuple of the form (name, sequence)
    :param gees: list of positions of G's in the reference sequence
    :return ctable: a contingency table for the query sequence
    """

    mut = re.compile('[AGCT](?=[AG][AGT])')  # Matches potential mutation sites (GRD)
    ctrl = re.compile('[AGCT](?=[CT].|[AG]C)')  # Matches control sites (YN|RC)

    # Keele et al. apply Fisher's exact test to a contingency table
    ctable = [
        [0, 0],  # (G->A, G) in GRD context
        [0, 0]  # in other contexts (GYN|GRC is !GRD)
    ]
    s = seq[1]  # query sequence

    # Matches downstream RD motif in query sequence
    ds_motifs = set([match.start() for match in mut.finditer(s)])
    # Matches downstream YN|RC motif in query sequence
    ds_ctrl = set([match.start() for match in ctrl.finditer(s)])

    mut_sites = set()  # locations of potential mutation sites (GRD motifs)
    ctrl_sites = set()  # locations of potential control sites (GYN|GRC)

    # Iterate through all positions with G in reference sequence
    for i in gees:
        nt = s[i]  # Corresponding base in query sequence
        row = (0 if nt == 'A' else 1)  # G->A mutation in downstream
        if i in ds_motifs:
            ctable[0][row] += 1
            if row:
                mut_sites.update({i+1})
        elif i in ds_ctrl:
            ctable[1][row] += 1
            if row:
                ctrl_sites.update({i+1})

    # Calculate one-sided P-value
    odds_ratio, p_value = stats.fisher_exact(ctable, alternative='greater')

    # prepare output dictionary
    result = {}
    result['seq_name'] = seq[0]  # the name of the sequence
    result['num_muts'] = ctable[0][0]  # the number of mutated sites, compared to the reference sequence
    result['pot_muts'] = ctable[0][0] + ctable[0][1]  # the number of potential mutation sites (GRD)
    result['ctrl_muts'] = ctable[1][0]  # the number of control sites
    result['potential_ctrls'] = ctable[1][0] + ctable[1][1]  # the number of potential control sites (GYN|GRC)
    result['rate_ratio'] = rate_ratio(ctable)  # the rate ratio for the mutation
    result['p_value'] = p_value  # Fisher's exact test P-value
    result['odds_ratio'] = odds_ratio
    result['ctable'] = ctable  # contingency table
    result['mut_sites'] = mut_sites
    result['ctrl_sites'] = ctrl_sites
    return result


def hypermut(infile, cons, skip=None):
    """
    Determines if a sequence is hypermutated compared to the reference sequence
    :param infile: the input FASTA file
    :param cons: <bool> True if the consensus sequence is designated as the reference sequence
                        False if the first sequence is designated as the reference sequence
                 <str> user provides reference sequence
    :param skip:  int, the number of records to skip
    :return results: a list of dict objects
    """

    with open(infile) as handle:
        fasta = convert_fasta(handle)

    if type(cons) is str:
        refseq = cons
        query_seqs = fasta
    else:
        if cons:
            refseq = get_consensus(fasta)
            query_seqs = fasta
        else:
            refseq = fasta[0][1]  # First sequence is the reference sequence
            query_seqs = fasta[1:]

    if skip:
        print("Skipping first {} records...\n".format(skip))
        fasta = fasta[int(skip):]
        query_seqs = fasta

    # Check that sequences are aligned
    length = len(fasta[0][1])
    for h, s in fasta:
        assert length == len(s), "Sequences are not aligned."

    gees = [i for i, nt in enumerate(refseq) if nt.upper() == 'G']  # Locate GRD motifs in reference sequence
    results = []

    # Iterate through sequences and identify substitutions from consensus
    for j, seq in enumerate(query_seqs):
        result = make_results(seq, gees)
        results.append(result)

    return results


def get_consensus(fasta):
    print("Generating consensus sequence...\n")
    refseq = consensus(fasta)
    return refseq


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


def is_hypermutated(result, alpha=0.05):
    return result['p_value'] is not None and result['p_value'] <= alpha


def pretty_print(results):
    """
    Print results
    :param results: a list of MutationInfo objects
    """

    print("{0}\t\t{1}\t\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}\t\t{7}\t\t{8}".format
          ("Sequence", "Muts", "Potential Mut Sites", "Control Muts", "Potential Controls",
           "Rate Ratio", "Fisher's Exact P-value", "Odds Ratio", "Hypermutated"))
    print('-' * 169)

    # Print values of rows under corresponding headings
    for result in results:
        result['is_hypermutated'] = is_hypermutated(result)
        print("{seq_name:8.8}\t\t{num_muts}\t\t\t{pot_muts}\t\t\t\t\t{ctrl_muts}\t\t\t\t\t{potential_ctrls}"
              "\t\t\t\t\t{rate_ratio:2.2}\t\t\t{p_value:6.6}\t\t\t\t\t{odds_ratio:6.6}\t\t{is_hypermutated}"
              .format(**result))

    # Print summary of hypermutated sequences
    print("\nSummary:")
    count = 0
    for result in results:
        if result['p_value'] <= 0.05:
            print("{seq_name} appears to be hypermutated (OR={odds_ratio})".format(**result))
            count += 1

    if count == 0:
        print("No sequences appear to be hypermutated.")



def make_data_file(file_name, results):
    """
    Writes detailed output of hypermut to a text file
    :param file_name: name of the output file
    :param results: list of MutationInfo objects
    """
    fieldnames = ["Sequence", "Muts", "Potential Mut Sites", "Control Muts",
                  "Potential Controls", "Rate Ratio", "Fisher's Exact P-value",
                  "Odds Ratio", "Hypermutated"]
    output = open(file_name, 'w+')
    output.write(','.join(fieldnames) + '\n')  # header row

    for row in results:
        output.write("Results:\n")
        output.write("{seq_name:8.8},{num_muts},{po_muts},{ctrl_muts},"
                     "{potential_ctrls},{rate_ratio:2.2},{p_value:6.6},"
                     "{odds_ratio:6.6},{is_hypermutated}\n".format(**results))

        # Print summary of hypermutated sequences
        output.write("\nSummary:\n")
        for result in results:
            if result.p_value <= 0.05:
                output.write("{seq_name} appears to be hypermutated (OR={odds_ratio})\n".format(**result))
            else:
                output.write("No sequences appear to be hypermutated.\n")
                break

        # Print detailed output
        output.write("\nLocations of matches:\n")
        for result in results:
            output.write("\nSequence Name: {seq_name:8.8}".format(**result))
            output.write("\nPos\tMut\n")
            for key in result.mut_sites:
                output.write("{}\t{}\n".format(key, result['mut_sites'][key]))

            output.write("\nSequence Name: {} control".format(result.seq_name))
            output.write("\nPos\tMut\n")
            for key in result.ctrl_sites:
                output.write("{}\t{}\n".format(key, result['ctrl_sites'][key]))


def main():
    # command line interface
    parser = argparse.ArgumentParser(
        description='Classify HIV-1 sequences as being hypermutated '
                    'based on dinucleotide frequencies relative to the consensus.'
    )
    parser.add_argument('fasta', help='<input> path to FASTA file. If no reference sequence is specified, '
                                      'the program will designate the first sequence as the reference sequence.')
    parser.add_argument('--consensus', action="store_true", default=False,
                        help='<option> the majority-rule consensus sequence is designated as the reference sequence.')
    parser.add_argument('--skip', type=int, help="<option> number of records to skip")
    parser.add_argument('--out', help="<option> write output to the specified file")

    args = parser.parse_args()

    results = hypermut(args.fasta, args.consensus, args.skip)
    if args.out:
        make_data_file(args.out, results)
    pretty_print(results)


if __name__ == '__main__':
    main()
