#TODO: consistent reference coordinates across outputs

import sys
import subprocess
import random
import argparse

from poplars.common import convert_fasta
from poplars.mafft import align



def pdistance(seq1, seq2):
    """
    Calculate p-distance between two aligned sequences
    :param seq1: First sequence
    :param seq2: Second sequence
    :return: <ndiff> is the number of differences, <denom> is the number of valid positions
    """
    denom = 0.  # number of valid columns
    ndiff = 0
    for i, nt1 in enumerate(seq1):
        nt2 = seq2[i]
        if nt1 == '-' or nt2 == '-':
            continue
        denom += 1
        if nt1 != nt2:
            ndiff += 1
    return ndiff, denom


def bootstrap(s1, s2, reps=100):
    """
    Sample positions at random with replacement.
    :param s1:  first sequence
    :param s2:  second sequence (must be aligned to s1)
    :param reps:  number of replicates to generate
    
    :yield: tuples of sequences generated by bootstrap resampling
    """
    seqlen = len(s1)
    assert len(s2)==seqlen, "s1 and s2 must be of same length in bootstrap()"
    
    for rep in range(reps):
        result = []
        bootstrap = [random.randint(0, seqlen-1) for _ in range(seqlen)]        
        b1 = ''.join([s1[i] for i in bootstrap])
        b2 = ''.join([s2[i] for i in bootstrap])
        yield b1, b2


def riplike(seq, window=400, step=5, nrep=100):
    """
    :param seq:  query sequence
    :param outfile:  open file stream in write mode for results
    :param window:  width of sliding window in nucleotides
    :param step:  step size of sliding window in nucleotides
    :param nrep:  number of replicates for nonparametric bootstrap sampling
    """
    
    # subset of HIV-1 group M subtype references curated by LANL
    with open('poplars/ref_genomes/HIV1_Mgroup.fasta') as handle:
        reference = convert_fasta(handle)
    
    results = []
    
    # append query sequence to reference alignment
    fasta = align(seq, reference)
    
    # eliminate insertions in query relative to references
    try:
        conseq = dict(fasta)['CON_OF_CONS']
    except:
        print(fasta)
        raise
    skip = [i for i in range(len(conseq)) if conseq[i] == '-']
    fasta2 = []
    for h, s in fasta:
        s2 = [nt for i, nt in enumerate(s) if i not in skip]
        fasta2.append([h, ''.join(s2)])
    fasta = fasta2
    
    query = dict(fasta)['query']  # aligned query
    seqlen = len(query)
    
    for left in range(0, seqlen-window, step):
        best_p, second_p = 1., 1.  # maximum p-distance
        best_ref, second_ref = None, None
        best_seq = ''
        
        # cut slice from query sequence for this window
        q1 = query[left:(left+window)]
        
        # iterate over reference genomes
        for h, s in fasta:
            if h=='query' or h=='CON_OF_CONS':
                continue 
                
            # slice window segment from reference
            s1 = s[left:(left+window)]
            
            # calculate p-distance
            ndiff, denom = pdistance(s1, q1)
            if denom == 0:
                # no overlap!  TODO: require minimum overlap?
                continue
            pd = ndiff/denom
            
            if pd < best_p:
                # query is closer to this reference
                second_p = best_p
                second_ref = best_ref
                best_p = pd
                best_ref = h
                best_seq = s1
                
            elif pd < second_p:
                # replace second best
                second_p = pd
                second_ref = h
        
        
        if best_ref is None:
            outfile.write('{},{},None,,None,,\n'.format(header, left))
            continue
        
        # use nonparametric bootstrap to determine significance
        if second_ref is not None:
            boot_dist = []
            for bs, bq in bootstrap(best_seq, q1, reps=nrep):
                ndiff, denom = pdistance(bs, bq)
                if denom > 0:
                    boot_dist.append(ndiff/denom)
                    
            # how many are closer than second best?
            quant = list(map(lambda x: x < second_p, boot_dist))
            quant = sum(quant) / float(len(quant))
            
            results.append({
                'left': left, 'best_ref': best_ref, 'best_p': best_p, 
                'second_ref': second_ref, 'second_p': second_p, 'quant': quant
            })
        else:
            # if no valid second best, accept first without bootstrap
            results.append({
                'left': left, 'best_ref': best_ref, 'best_p': best_p, 
                'second_ref': None, 'second_p': None, 'quant': None
            })
        
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

    args = parser.parse_args()
    args.outfile.write('qname,pos,rname,pdist,rname2,pdist2,qboot\n')
    
    fasta = convert_fasta(args.infile)
    for h, s in fasta:
        print(h)  # crude progress monitoring
        res = riplike(s, args.outfile, window=args.window, step=args.step, nrep=args.nrep)
        args.outfile.write(
            '{},{left},{best_ref},{best_p},{second_ref},{second_p},{quant}\n'.format(header, **res)
        )
        
    args.outfile.close()
    

if __name__ == '__main__':
    main()


