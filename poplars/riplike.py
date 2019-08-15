# TODO: consistent reference coordinates across outputs

import random
import argparse
import os

from poplars.common import convert_fasta
from poplars.mafft import align


def hamming(bit_aln):
    """
    Convert list of lists into boolean outcomes (difference between query and reference)
    :param bit_aln: aligned sequences converted to bit-strings
    :return: dictionary of boolean lists keyed by reference label
    """
    query = bit_aln.pop('query')

    # Iterate over remaining sequences as references
    results = {}
    for h, s in bit_aln.items():
        result = []
        for i in range(len(query)):
            nt2 = s[i]
            if (query[i] | nt2) & 0B1111:
                result.append(None)
                continue
            result.append(bin(query[i] ^ nt2))
        results.update({h: result})

    return results


# def hamming(fasta):
#     """
#     Convert list of lists into boolean outcomes (difference between query and reference)
#     :param aln_seqs: aligned sequences converted to bit-strings
#     :return: dictionary of boolean lists keyed by reference label
#     """
#     aln = dict(fasta)
#     assert "query" in aln, "Argument <fasta> must contain 'query' entry"
#     query = aln.pop('query')
#     _ = aln.pop('CON_OF_CONS')
#
#     # iterate over remaining sequences as references
#     results = {}
#     for h, s in aln.items():
#         result = []
#         for i, nt1 in enumerate(query):
#             nt2 = s[i]
#             if nt1 == '-' or nt2 == '-':
#                 result.append(None)
#                 continue
#             result.append(int(nt1 != nt2))
#         results.update({h: result})
#
#     return results

def update_alignment(seq, reference):
    """
    Append query sequence <seq> to reference alignment and remove insertions relative to
    global consensus sequence.
    :param seq: the query sequence
    :param reference: the reference sequence
    :return: a list of [header, sequence] lists
    """
    # append query sequence to reference alignment
    fasta = align(seq, reference)

    # eliminate insertions in query relative to references
    try:
        conseq = dict(fasta)['CON_OF_CONS']
    except KeyError:
        print("ERROR: reference alignment in poplars.riplike does not contain CON_OF_CONS entry")
        raise

    skip = [i for i in range(len(conseq)) if conseq[i] == '-']
    fasta2 = []
    for h, s in fasta:
        s2 = [nt for i, nt in enumerate(s) if i not in skip]
        fasta2.append([h, ''.join(s2)])

    return fasta2


def encode(fasta):
    """
    Encodes each nucleotide in a sequence using 4-bits
    :param fatsa: the result of the alignment
    :return: the sequence as a bitstring where each nucleotide is encoded using a 4-bits
    """
    bit_aln = dict(fasta)
    assert "query" in bit_aln, "Argument <fasta> must contain 'query' entry"
    _ = bit_aln.pop('CON_OF_CONS')

    binary_nt = {'A': 0B0001, 'T': 0B0010, 'C': 0B0011, 'G': 0B0100, ' ': 0B0000, '-': 0B1111}

    for h, s in bit_aln.items():
        seq = []
        for nt in s:
            seq.append(binary_nt[nt])
        bit_aln[h] = seq
    return bit_aln


def riplike(seq, reference, window=400, step=5, nrep=100):
    """
    :param seq:  query sequence
    :param reference: the alignment background
    :param window:  width of sliding window in nucleotides
    :param step:  step size of sliding window in nucleotides
    :param nrep:  number of replicates for nonparametric bootstrap sampling
    :return: list of result dictionaries in order of window position
    """

    results = []

    fasta = update_alignment(seq, reference)
    query = dict(fasta)['query']  # aligned query
    seqlen = len(query)
    bit_aln = encode(fasta)
    ham = hamming(bit_aln)

    for centre in range(window // 2, seqlen - (window // 2), step):
        best_p, second_p = 1., 1.  # maximum p-distance
        best_ref, second_ref = None, None
        best_seq = []

        # iterate over reference genomes
        for h, s in ham.items():
            if h == 'query' or h == 'CON_OF_CONS':
                continue

            # slice window segment from reference
            s1 = s[centre - (window // 2): centre + (window // 2)]
            s2 = [x for x in s1 if x is not None]

            # calculate p-distance
            ndiff = sum(s2)
            denom = len(s2)
            if denom == 0:
                # no overlap!  TODO: require minimum overlap?
                continue
            pd = ndiff / denom

            if pd < best_p:
                # query is closer to this reference
                second_p = best_p
                second_ref = best_ref
                best_p = pd
                best_ref = h
                best_seq = s2
            elif pd < second_p:
                # replace second best
                second_p = pd
                second_ref = h

            result = {'centre': centre, 'best_ref': best_ref, 'best_p': best_p,
                      'second_ref': second_ref, 'second_p': None if second_ref is None else second_p}

            quant = None
            if second_ref is not None and nrep > 0:
                # use nonparametric bootstrap to determine significance
                count = 0.
                n = s.count(0B0001)
                sample = random.choices(best_seq, k=n * nrep)
                for rep in range(nrep):
                    boot = sample[rep: rep + n]
                    if sum(boot) / n < second_p:
                        count += 1
                quant = count / nrep

            result.update({'quant': quant})
            results.append(result)

        ### Move result in for loop to bootstrap for every list of booleans values in the dictionary
        ### The value is the list of boolean values representing the number of differences
        ### number of repetitions = number of 1's
        # OR... Sample over only where there are 1's??? ... if this option, how to reconcile random.choices...

        # result = {'centre': centre, 'best_ref': best_ref, 'best_p': best_p,
        #           'second_ref': second_ref, 'second_p': None if second_ref is None else second_p}

        # quant = None
        # if second_ref is not None and nrep > 0:
        #     # use nonparametric bootstrap to determine significance
        #     count = 0.
        #     n = s.count(0B0001)
        #     sample = random.choices(best_seq, k=n*nrep)
        #     for rep in range(nrep):
        #         boot = sample[rep: rep + n]
        #         if sum(boot) / n < second_p:
        #             count += 1
        #     quant = count / nrep
        #
        # result.update({'quant': quant})
        # results.append(result)

    return results


def main():
    parser = argparse.ArgumentParser(
        description='An approximate implementation of the Recombinant Identification '
                    'program by the Los Alamos National Laboratory.'
    )
    parser.add_argument('infile', type=argparse.FileType('r'),
                        help='<input> FASTA file containing sequences to process.')
    parser.add_argument('outfile', type=argparse.FileType('w'),
                        help='<output> file to write CSV results.')
    parser.add_argument('-window', type=int, default=400,
                        help='<optional, int> Window size for p-distances.')
    parser.add_argument('-step', type=int, default=5,
                        help='<optional, int> Window step size.')
    parser.add_argument('-nrep', type=int, default=100,
                        help='<optional, int> Number of bootstrap replicates.')
    parser.add_argument('-custombg', type=argparse.FileType('r'),
                        help='<optional> FASTA file to be used as the alignment background')

    args = parser.parse_args()

    if args.custombg:
        ref_seq = args.custombg
    else:
        # subset of HIV-1 group M subtype references curated by LANL
        seq_path = os.path.dirname(os.path.abspath(__file__))
        ref_seq = os.path.join(seq_path, 'ref_genomes/HIV1_Mgroup.fasta')

    with open(ref_seq) as handle:
        reference = convert_fasta(handle)

    args.outfile.write('qname,pos,rname,pdist,rname2,pdist2,qboot\n')

    fasta = convert_fasta(args.infile)
    for h, s in fasta:
        print(h)  # crude progress monitoring
        results = riplike(s, reference, window=args.window, step=args.step, nrep=args.nrep)
        for result in results:
            args.outfile.write(
                '{},{centre},{best_ref},{best_p},{second_ref},{second_p},{quant}\n'
                .format(h, **result)
            )

    args.outfile.close()


if __name__ == '__main__':
    main()


