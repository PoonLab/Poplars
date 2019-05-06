#TODO: consistent reference coordinates across outputs

import sys
import tempfile
import subprocess
import random
import argparse



def convert_fasta (handle):
    result = []
    sequence = ''
    for line in handle:
        if line.startswith('$'): # skip header line
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h,sequence])
                sequence = ''   # reset
            h = line.strip('>#\n')
        else:
            sequence += line.strip('\n').upper()
            
    result.append([h, sequence])  # handle last entry
    return result


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


def realign(seq, reference):
    # add query sequence to consensus
    with tempfile.NamedTemporaryFile('w', delete=False) as handle:
        handle.write(reference)
        handle.write('>query\n{}\n'.format(seq))
        handle.close()

        output = subprocess.check_output(['mafft', '--quiet', handle.name])
    
    result = output.decode('utf-8')
    return convert_fasta(result.split('\n'))


def bootstrap(s1, s2, reps=100):
    """
    Sample alignment positions at random with replacement.
    Return a list of nested lists (FASTA objects).
    """
    seqlen = len(s1)
    
    for rep in range(reps):
        result = []
        bootstrap = [random.randint(0, seqlen-1) for _ in range(seqlen)]        
        b1 = ''.join([s1[i] for i in bootstrap])
        b2 = ''.join([s2[i] for i in bootstrap])
        yield b1, b2


def riplike(header, seq, outfile, reference, window=400, step=5, nrep=100):
    """
    Assumes query sequence is last
    """
    
    # append query sequence to reference alignment
    fasta = realign(seq, reference)
    
    # eliminate insertions in query relative to references
    conseq = dict(fasta)['CON_OF_CONS']
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
                # no overlap!
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
            
            outfile.write('{},{},{},{},{},{},{}\n'.format(
                header, left, best_ref, best_p, second_ref, second_p, quant)
            )
        else:
            # if no valid second best, accept first without bootstrap
            outfile.write('{},{},{},{},,,\n'.format(
                header,left, best_ref, best_p)
        )


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
    parser.add_argument('-reference', type=argparse.FileType('r'),
                        default='HIV1_Mgroup.fasta',
                        help='<optional> FASTA file with HIV reference genomes.')

    args = parser.parse_args()

    reference = args.reference.read()
    reference = reference.replace('-', '')  # remove all gaps!

    fasta = convert_fasta(args.infile)

    args.outfile.write('qname,pos,rname,pdist,rname2,pdist2,qboot\n')
    for h, s in fasta:
        print(h)  # crude progress monitoring
        riplike(h, s, args.outfile, reference, window=args.window, step=args.step, nrep=args.nrep)


if __name__ == '__main__':
    main()


