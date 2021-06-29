#!/usr/bin/env python3
import random
import re

mixture_dict = {'W': 'AT', 'R': 'AG', 'K': 'GT', 'Y': 'CT', 'S': 'CG',
                'M': 'AC', 'V': 'AGC', 'H': 'ATC', 'D': 'ATG',
                'B': 'TGC', 'N': 'ATGC', '-': 'ATGC'}

ambig_dict = dict(("".join(sorted(v)), k) for k, v in mixture_dict.items())


def convert_fasta(handle):
    """
    Import FASTA-formatted sequences from file stream
    :param handle:  File stream in read mode, or the contents of a file split into lines
    :return: list of lists containing header-sequence pairs
    """
    result = []
    sequence, h = '', ''

    # Verifies files have correct formatting
    for i, line in enumerate(handle):
        if line.startswith('$'):
            continue
        elif line.startswith('>') or line.startswith('#'):
            break
        else:
            print("No header")
            raise NameError

    if hasattr(handle, 'seek'):
        handle.seek(0)

    for line in handle:
        if line.startswith('$'):  # skip header line
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h, sequence])
                sequence = ''  # reset
            h = line.strip('>#\n\r')
        else:
            sequence += line.strip('\n\r').upper()

    result.append([h, sequence])  # handle last entry
    return result


def transpose_fasta(fasta):
    # some checks to make sure the right kind of object is being sent
    if type(fasta) is not list:
        return None
    if type(fasta[0]) is not list or len(fasta[0]) != 2:
        return None

    n_columns = len(fasta[0][1])
    res = []
    for c in range(n_columns):
        res.append([s[c] for h, s in fasta])

    return res


def plurality_consensus(column, alphabet='ACGT', resolve=False):
    """
    Plurality consensus - nucleotide with highest frequency.
    In case of tie, report mixtures.
    """
    freqs = {}
    for char in alphabet:
        freqs.update({char: 0})

    for char in column:
        if char in alphabet:
            freqs[char] += 1
        elif char in mixture_dict:
            # handled ambiguous nucleotides with equal weighting
            resolutions = mixture_dict[char]
            for char2 in resolutions:
                freqs[char2] += 1. / len(resolutions)
        else:
            # unrecognized nucleotide character
            pass

    base = max(freqs, key=lambda n: freqs[n])
    max_count = freqs[base]
    possib = list(filter(lambda n: freqs[n] == max_count, freqs))
    if len(possib) == 1:
        return possib[0]
    elif "-" in possib:
        if resolve:
            possib.remove("-")
            if len(possib) == 0:
                return "-"
            elif len(possib) == 1:
                return possib[0]
            else:
                return ambig_dict["".join(sorted(possib))]
        else:
            # gap character overrides ties
            return "-"
    else:
        if resolve:
            return random.sample(possib, 1)[0]
        else:
            return ambig_dict["".join(sorted(possib))]


def consensus(fasta, alphabet='ACGT', resolve=False):
    """
    Return plurality consensus of alignment.
    """
    consen = []
    columns = transpose_fasta(fasta)

    for column in columns:
        consen.append(plurality_consensus(column, alphabet=alphabet, resolve=resolve))

    newseq = "".join(consen)

    return newseq


def convert_clustal(handle):
    """
    Import Clustal-formatted sequences from file stream
    :param handle:  File stream in read mode, or the contents of a file split into lines
    :return: dictionary containing header-sequence pairs
    """
    result = {'aln': ''}

    pat = re.compile('\w+\s*')
    offset = pat.match(handle[3]).end()

    # Loop over sequences in the Clustal alignment
    for ln in handle[3:]:
        # Skip first three lines
        if len(ln) > 0:
            # Match sequence label and build up sequence
            if ln[0].isalpha():
                line = ln.split()
                header = line[0].strip()
                seq = line[1].strip()
                result.setdefault(header, '')   # Add header to dictionary
                result[header] += seq

            # Add conservation information
            else:
                result['aln'] += ln[offset:].strip('\n\r')

    return result


def resolve_mixtures(seq, replaceN=False):
    """
    Randomly resolve mixtures (ambiguous base calls) to nucleotides
    :param seq:  str, nucleotide sequence
    :param replaceN:  if True, replace N with [ACGT]
    :return:  str, sequence with all mixtures [WRKYSMBDHVN] replaced with
              nucleotides [ACGT]; or None if sequence contains illegal characters
    """
    newseq = ''
    for nt in seq.upper():
        if nt == 'N':
            newseq += random.sample('ACGT', 1)[0] if replaceN else 'N'
        elif nt in 'ACGT-':
            newseq += nt  # unambiguous base call or gap
        elif nt in mixture_dict:
            newseq += random.sample(mixture_dict[nt], 1)[0]
        else:
            print("Unexpected character {} in resolve_mixtures()".format(nt))
            return None
    return newseq
