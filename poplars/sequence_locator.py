#!/usr/bin/env python3
"""
Python implementation of HIV Sequence Locator from https://www.hiv.lanl.gov
HIV and SIV genomic region coordinates are based on the HXB2 and Mac239
annotation resources from https://www.hiv.lanl.gov/content/sequence/HIV/MAP/annotation.html

Note: The first 256 nucleotides of SIVMM239 correspond to the flanking sequence,
and are included in the complete SIV genome (https://www.ncbi.nlm.nih.gov/nucleotide/M33262)
"""

import re
import sys
import os
import argparse

from poplars.common import convert_fasta, convert_clustal
from poplars.mafft import align

NON_CODING = ["5'LTR", "TAR", "3'LTR"]


class Region(object):
    """
    Stores information about a genome region
    """
    def __init__(self, region_name, genome, ncoords=None, pcoords=None):
        """
        :param region_name: The name of the genomic region
        :param genome: A reference to a genome object
        :param ncoords: A list containing the global start and end indices of the nucleotide region
        :param pcoords: A list containing the global start and end indices of the protein region
        """
        self.region_name = region_name
        self.genome = genome
        self.ncoords = ncoords
        self.pcoords = pcoords
        self.nt_seq = None
        self.aa_seq = None

    def get_coords(self, base):
        if base == 'nucl':
            return self.ncoords
        else:
            return self.pcoords

    def get_sequence(self, base):
        if base == 'nucl':
            return self.nt_seq
        else:
            return self.aa_seq

    def get_sequence_from_genome(self, base):
        seq = ''
        if base == 'nucl':
            seq = self.genome.nt_seq[self.ncoords[0] - 1: self.ncoords[1]]
        else:
            try:
                seq = self.genome.aa_seq[self.region_name]
            except KeyError as e:
                print(e)

        return seq

    def set_coords(self, coord_pair, base):
        if base == 'nucl':
            self.ncoords = coord_pair
        else:
            self.pcoords = coord_pair

    def set_nt_seq_from_genome(self):
        self.nt_seq = self.genome.nt_seq[self.ncoords[0] - 1: self.ncoords[1]]

    def set_aa_seq_from_genome(self):
        self.aa_seq = self.genome.aa_seq[self.region_name]

    def set_sequence(self, base, seq):
        if base == 'nucl':
            self.nt_seq = seq
        else:
            self.aa_seq = seq


class RefRegion(Region):
    """
    Stores information about each region in the reference genome
    """
    def __init__(self, region_name, genome, ncoords=None, pcoords=None):
        super(RefRegion, self).__init__(region_name, genome, ncoords, pcoords)
        self.codon_aln = self.make_codon_aln()

    def make_codon_aln(self):
        """
        Aligns the protein sequence with its associated nucleotide sequence
        """
        if self.nt_seq is not None and self.aa_seq is not None:

            # Account for the frame-shift at position 5572 in vpr (extra T relative to other subtype B)
            if self.region_name == 'Vpr' and self.genome.virus == 'hiv':
                pos = 215   # position of extra T relative to the start of the Vpr CDS (0-based indexing)
                nt_seq = self.nt_seq[:pos] + self.nt_seq[(pos + 1):]
            else:
                nt_seq = self.nt_seq

            codon_aln = []
            codons = [''.join(t) for t in zip(*[iter(nt_seq)] * 3)]
            for aa, codon in zip(self.aa_seq, codons):
                # Check if in-frame stop codon
                if codon == 'TAA' or codon == 'TGA' or codon == 'TAG':
                    codon_aln.append('-*-')
                else:
                    codon_aln.append('-{}-'.format(aa))

            # Account for stop codon at the end of the sequence
            if (len(self.aa_seq) + 1) == len(codons):
                codon_aln.append('-*-')

            self.codon_aln = ''.join(codon_aln)
            return self.codon_aln

    def find_overlap(self, base, q_coords):
        """
        Gets the sequence regions that overlap with the region of interest
        :param base: the base of the sequence (nucleotide or protein)
        :param q_coords: the coordinates of the query region (global coordinates)
        :return overlap: A QueryRegion object that represents the where the query overlaps with the reference region
        """
        ref_coords = self.get_coords(base)  # Coordinates are global coordinates

        # If the ref_coords are in the range of the query region and the q_coords are in the range of the ref region
        if not ((ref_coords[0] < ref_coords[1] < q_coords[0]) or (q_coords[0] < q_coords[1] < ref_coords[0])):

            # If the end of the query region exceeds the end of the reference region
            if q_coords[1] > ref_coords[1]:
                end = ref_coords[1]
            else:
                end = q_coords[1]

            # If the query region starts before the reference region starts
            if q_coords[0] < ref_coords[0]:
                start = ref_coords[0]
            else:
                start = q_coords[0]

            overlap_coords = [start, end]       # Coordinates of the overlap between the query and the reference

            overlap = QueryRegion(self.region_name, self, base, self.genome, q_coords)
            if base == 'nucl':
                overlap.set_coords(overlap_coords, 'nucl')
                overlap.set_nt_seq_from_genome()
                overlap.cds_offset = overlap.set_pos_from_cds()
                overlap.qstart = overlap.set_pos_from_qstart('nucl')

            else:
                overlap.set_coords(overlap_coords, 'prot')
                overlap.set_aa_seq_from_genome()
                overlap.set_pos_from_cds()
                overlap.cds_offset = overlap.set_pos_from_cds()

            # Set the nucleotide and protein regions
            if base == 'nucl':
                overlap.set_coords(overlap_coords, 'nucl')
                overlap.set_nt_seq_from_genome()
                prot_overlap = self.set_protein_equivalents(overlap)
                overlap.qstart = overlap.set_pos_from_qstart('nucl')

                # Set protein coordinates if the overlap is in a coding region
                if overlap.region_name not in NON_CODING:
                    overlap.set_coords(prot_overlap[1], 'prot')
                    overlap.set_sequence('prot', prot_overlap[0])
                    overlap.gstart = overlap.set_pos_from_gstart()
                    overlap.cds_offset = overlap.set_pos_from_cds()
            # else:
            #     nucl_overlap = self.set_nucleotide_equivalents(overlap)
            #     overlap.set_coords(overlap_coords, 'prot')
            #     overlap.set_sequence('prot', self.get_sequence('prot')[start: end])
            #     overlap.set_coords(nucl_overlap[1], 'nucl')
            #     overlap.set_sequence(nucl_overlap[0], 'nucl')

            return overlap
        return None

    def set_protein_equivalents(self, overlap):
        """
        Finds the protein equivalent to the nucleotide sequence
        :param overlap: reference to the region that aligns with the reference region
        :return aa_seq, aa_coords: the sequence and coordinates of the overlapping protein sequence
        """
        aa_seq, aa_coords = '', []
        if self.region_name not in NON_CODING:

            # Slice aligned nucleotide sequence using coordinates of the query region
            if self.codon_aln is not None:
                aligned_prot = self.codon_aln[overlap.query_coords[0]: overlap.query_coords[1]]
                aa_seq = re.sub('[-]', '', aligned_prot)  # Get corresponding protein sequence

                # Set the local amino acid coordinates
                pat = re.compile(aa_seq)
                m = pat.search(self.aa_seq)
                aa_coords = [m.start() + 1, m.end()]

                # Set the pcoords of the overlap region
                overlap.set_coords(aa_coords, 'prot')

        return aa_seq, aa_coords


class QueryRegion(Region):
    """
    Represents information about the query sequence region
    """
    def __init__(self, region_name, ref_region, base, genome, query_coords, ncoords=None, pcoords=None):
        """
        :param region_name: name of the query region
        :param ref-region: reference to the genome region
        :param genome: reference to the genome
        :param query_coords: the local coordinates of the region
        :param ncoords: the global coordinates of the nucleotide sequence
        :param pcoords: the global coordinates of the protein sequence
        """
        super(QueryRegion, self).__init__(region_name, genome, ncoords, pcoords)

        self.ref_region = ref_region
        self.codon_aln = ''
        self.query_coords = query_coords
        self.cds_offset = self.set_pos_from_cds()
        self.gstart = self.set_pos_from_gstart()
        self.qstart = self.set_pos_from_qstart(base)
        self.pstart = self.set_pos_from_pstart()

    def set_pos_from_cds(self):
        """
        Gives the position of a the query sequence relative to the start of the coding sequence
        """
        if self.ncoords and self.ref_region.ncoords:
            query_start = self.ncoords[0]
            query_end = self.ncoords[1]
            ref_start = self.ref_region.ncoords[0]
            ref_end = self.ref_region.ncoords[1]

            if query_start == ref_start:
                start = 1
            else:
                start = query_start - ref_start + 1     # 1-based indexing

            if query_end == ref_end:
                end = query_end - query_start + 1
            else:
                end = start + (query_end - query_start)

            return [start, end]

    def set_pos_from_gstart(self):
        return self.ncoords

    def set_pos_from_qstart(self, base):
        """
        Gives the position of the sequence relative to the start of the region of interest
        :param base: The base of the sequence (nucleotide or protein)
        :return: The position of the reference region relative to the query region
        """
        overlap_coords = self.get_coords(base)
        parent_coords = self.ref_region.get_coords(base)
        overlap_seq = self.get_sequence(base)

        if overlap_coords is not None and parent_coords is not None and overlap_seq is not None:
            start_offset = overlap_coords[0] - parent_coords[0] + 1
            end_offset = start_offset + len(overlap_seq) - 1
            return [start_offset, end_offset]

    def set_pos_from_pstart(self):
        """
        Gives the position of the sequence relative to the start of the protein sequence
        """
        if self.pcoords is not None:
            if 'LTR' in self.region_name:
                return None
            else:
                return self.pcoords

    def set_pcoords_from_ncoords(self):
        """
        Sets protein coordinates relative to the protein start, given the nucleotide coordinates
        """
        if self.cds_offset is not None:
            prot_start = self.cds_offset[0] // 3 + 1
            prot_end = self.cds_offset[1] // 3
            self.pcoords = [prot_start, prot_end]

    def set_ncoords_from_pcoords(self):
        """
        Sets nucleotide coordinates given the protein coordinates relative to the protein start
        :return ncoords: the nucleotide coordinates
        """
        if self.pcoords is not None:
            nucl_start = self.pcoords[0] - 1 * 3
            nucl_end = self.pcoords[1] * 3
            return [nucl_start, nucl_end]


class Genome:
    """
    Stores information about the reference genome
    """
    def __init__(self, virus, nt_coords, nt_seq, aa_seq, reference_sequence, base):
        self.virus = virus
        self.nt_seq = nt_seq
        self.aa_seq = aa_seq  # List of lists
        self.reference_sequence = reference_sequence
        self.ref_base = base
        self.ref_genome_regions = self.make_ref_regions(nt_coords, aa_seq)

    def make_ref_regions(self, nt_coords, aa_seq):
        """
        Reads in the start and end coordinates of the genomic regions and associates the region with its coordinates.
        If no coordinate files are specified, set default nucleotide and protein coordinate files.
        :param nt_coords: File stream containing coordinates for nucleotide regions. Format: region_name,start,end
        :param aa_seq: A list of lists containing the protein sequences
        """
        ref_regions = {}
        # Parse nucleotide region coordinates file
        for nt_line in nt_coords:
            nt_line = nt_line.strip()
            nt_line = nt_line.split(',')
            nucl_coords = [int(nt_line[1]), int(nt_line[2])]
            seq_region = RefRegion(nt_line[0], self, nucl_coords)
            seq_region.set_nt_seq_from_genome()
            ref_regions[nt_line[0]] = seq_region

        # Parse protein coordinates file
        prot_names = []
        for aa_line in aa_seq:
            aa_line = aa_line.strip()
            prot_name = aa_line[aa_line.find('>') + 1: aa_line.find('|')]
            prot_names.append(prot_name)

        # Match protein regions to nucleotide regions
        for region_name in ref_regions:
            for prot_name in prot_names:
                if prot_name.startswith(region_name):
                    # Set protein coordinates and sequences
                    ref_regions[region_name].set_sequence('prot', aa_seq[region_name])
                    ref_regions[region_name].set_coords([1, len(aa_seq[region_name])], 'prot')

        # Align nucleotide sequence with protein sequence
        for name in ref_regions:
            ref_regions[name].make_codon_aln()

        return ref_regions

    def sequence_align(self, query_sequence, outfile=None):
        """
        Aligns the query sequence to the reference genome
        :param query_sequence: The query sequence
        :param outfile: <option> The file stream of the output file in write mode
        """
        if self.ref_base == 'nucl':
            result = align(query_sequence[1], self.reference_sequence[0][1], True)

        else:
            for ref_seq in self.ref_genome_regions:
                result = align(query_sequence[1], ref_seq.get_sequence('prot'))

        if outfile is not None:
            outfile.write("Alignment:\n")
        else:
            print("Alignment:")

        clustal_to_fasta(result, outfile)
        return result

    @staticmethod
    def find_query_match_coords(alignment):
        """
        Gets the indices of regions where the query aligns with the reference sequence
        :param alignment: the aligned query sequence
        :return query_matches: The positions where the query aligns with the reference sequences with no gaps
        """
        pat = re.compile('\w[\w-]*\w')  # Get the query sequence from the alignment
        match = pat.search(alignment)
        query_match_coords = (match.start() + 1, match.end())
        return query_match_coords

    def find_matches(self, base, qmatch_coords, lookup_table):
        """
        Finds the genomic regions where the query sequence aligns with the reference sequence
        :param base: The base of the query sequence
        :param qmatch_coords: tuple containing the coordinates of the query sequence in the alignment
        :param lookup_table: a list that maps coordinates of reference genome to the aligned sequence
        :return query_regions: <dict> with keys as region name and values as the QueryRegion object
        """
        query_regions = {}

        for name in self.ref_genome_regions:
            if name != 'Complete':
                # Convert query coordinates from "alignment space" to "region space"
                converted_coords = (lookup_table[qmatch_coords[0]], lookup_table[qmatch_coords[1]])
                # Find which regions overlap with the query region
                overlap = self.ref_genome_regions[name].find_overlap(base, converted_coords)

                if overlap is not None:
                    query_regions[name] = overlap

        return query_regions

    def retrieve(self, base, region, qstart=1, qend='end'):
        """
        Retrieves a sequence given its coordinates
        :param base: The base of the sequence (nucleotide or protein)
        :param region: The genomic region
        :param qstart: <option> The start coordinate of the query region (given as local coordinate)
        :param qend: <option> The end coordinate of the query region (given as local coordinate)
        :return: The genomic region defined by the starting and ending coordinates
        """
        if region == 'Rev' or 'rev':
            self.retrieve(base, 'Rev(exon1)', qstart, qend)
            self.retrieve(base, 'Rev(exon2)', 8379, 8653)

        elif region == 'Tat' or 'tat':
            self.retrieve(base, 'Tat(exon1)', qstart, 6045)
            self.retrieve(base, 'Tat(exon2)', 8379, 8469)

        elif region == 'LTR3' or region == 'Nef' or region == 'LTR5':
            self.retrieve(base, 'LTR3', qstart, qend)
            self.retrieve(base, 'Nef', qstart, qend)
            self.retrieve(base, 'LTR5', qstart, qend)

        else:
            if region in self.ref_genome_regions:
                ref_region = self.ref_genome_regions[region]
                global_range = ref_region.get_coords(base)
                region_start, region_end = global_range[0], global_range[1]

                if qend == 'end':
                    qend = region_end - region_start

                # Convert global region coordinates to local coordinates
                global_start = region_start + qstart
                global_end = region_end + qend

                # Check if the coordinates are valid
                if global_start < region_start or global_end < region_end:
                    print("Invalid {} coordinates: {} to {}.\nValid range: 1 to {}"
                          .format(region, qstart, qend, (region_end - region_start)))
                    sys.exit(0)

                retrieved_regions = self.find_matches(base, [global_start, global_end])
                query_region = retrieved_regions.pop(region)

                return query_region, retrieved_regions


def make_aln_lookup_table(aln):
    """
    Maps the coordinates of reference genome to the aligned sequence
    :param aln: dictionary of header, sequence pairs
    :param genome: the reference genome
    :return: list where indices correspond to a position in the aligned sequence
              and the value corresponds to the position in the reference sequence
    """
    aln_lookup_table = []
    count = -1
    for pos, nt in enumerate(aln['reference']):
        if nt != '-':   # Account for gaps
            count += 1
        aln_lookup_table.append(count)
    return aln_lookup_table


def clustal_to_fasta(aln, outfile=None):
    """
    Outputs a Clustal-formatted alignment to FASTA to the console by default or to a file
    :param aln: the Clustal formatted alignment as a dictionary of header, sequence pairs
    :param outfile: path to the FASTA output file
    """
    for header in aln:
        if header != 'aln':
            seq = aln[header]

            if outfile:
                with open(outfile) as out_handle:
                    out_handle.write('>{}\n'.format(header))
                    for i in range(0, len(seq), 60):
                        out_handle.write('{}\n'.format(seq[i:i + 60]))

            else:
                print('>{}'.format(header))
                for i in range(0, len(seq), 60):
                    print(seq[i:i+60])


def valid_sequence(base, sequence):
    """
    Verifies that input sequence is valid
    :param base: The base of the sequence (nucl or prot)
    :param sequence: A list of lists containing header and sequence pairs
    :raises ValueError: If the sequence is empty or if it contains invalid characters
    :return: <True> If the input sequence uses the correct alphabet, <False> otherwise
    """
    dna_alphabet = 'ATGC-*XN'
    aa_alphabet = 'ARDNCEQGHILKMFPSTWYV-*X'

    if not sequence:
        print("Invalid sequence: sequence length is 0\n")
        return False

    for h, s in sequence:
        if not s:
            print("Invalid sequence: sequence length is 0\n")
            return False
        elif base == 'nucl':
            if not all(pos in dna_alphabet for pos in s):
                print("Invalid nucleotide sequence:\n{}\n{}\n".format(h, s))
                return False
        else:
            if not all(pos in aa_alphabet for pos in s):
                print("Invalid amino acid sequence:\n{}\n{}\n".format(h, s))
                return False

    return True


def valid_inputs(virus, start_coord, end_coord, region):
    """
    Checks that the start and end coordinates of the genomic region are valid
    :param virus: The reference virus
    :param start_coord: The starting coordinate of the genomic region
    :param end_coord: The ending coordinate of the genomic region
    :param region: The genomic region
    :return: <True> If the coordinates are valid, <False> otherwise
    """

    if type(start_coord) == str:
        print("Invalid start coordinate type: {}".format(type(start_coord)))

    if start_coord < 0:
        print("Invalid start coordinate: {}".format(start_coord))
        return False

    if type(end_coord) != str and type(end_coord) != int:
        print("Invalid end coordinate type: {}".format(type(end_coord)))
        return False

    if type(end_coord) == str:
        if end_coord != "end":
            print("Invalid end coordinate: {}".format(end_coord))
            return False
    else:
        if end_coord <= 0 or start_coord >= end_coord:
            print("Invalid range: {} to {}".format(start_coord, end_coord))
            return False

    if virus == 'hiv' and region == 'Vpx' or virus == 'siv' and region == 'Vpu':
        print("Invalid region: {} in {}".format(region, virus))
        return False

    return True


def get_query(base, query, rev_comp):
    """
    Gets the query sequence and checks that it is valid
    :param base: The base (nucleotide or protein)
    :param query: The query sequence as a string or the file path to the query sequence
    :param rev_comp: Reverse complement flag (False by default)
    :return: A list of lists containing the sequence identifiers and the query sequences
    """

    # If the query is a string
    if not os.path.exists(query):
        query_lines = query.split('\n')
        [line.strip('\n') for line in query_lines]

        # Convert string into a list of lists
        if query_lines[0].startswith('>'):
            header = query_lines[0].strip('>')
            seq = ''.join(query_lines[1:])
        else:
            header = "query"
            seq = ''.join(query_lines)
        query_seq = [[header, seq.upper()]]

    # If the query is a file path
    else:
        with open(query, 'r') as query_handle:
            line = query_handle.readline()
            query_handle.seek(0)  # Reset pointer to beginning

            # Parse query sequence
            if line.startswith('>'):  # Fasta file
                query_seq = convert_fasta(query_handle)

            else:
                query_seq = []
                count = 1
                for line in query_handle:
                    line = line.strip('\n')
                    if len(line) > 0:
                        query_seq.extend([["Sequence{}".format(count), "{}".format(line.upper())]])
                        count += 1

    if not valid_sequence(base, query_seq):
        sys.exit(0)

    # At this point, the sequence is valid
    if base == 'nucl' and rev_comp:
        rc_query = reverse_comp(query_seq[0][1])
        header = query_seq[0][0]
        query_seq = [[header, rc_query]]

    return query_seq


def get_ref_seq(ref_seq, base):
    """
    Converts the reference sequence to a string and checks if the sequence is valid
    :param ref_seq: The path to the reference sequence
    :param base: The base of the sequence
    :return reference_sequence: A list of the form [header, sequence]
    """
    with open(ref_seq, 'r') as ref_handle:
        reference_sequence = convert_fasta(ref_handle)

    if valid_sequence(base, reference_sequence):
        return reference_sequence
    else:
        sys.exit(0)


def reverse_comp(query_sequence):
    """
    Reverses and complements the query sequence
    :param query_sequence: The query sequence
    :return: The reverse complement of the query sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '*': '*', 'N': 'N', '-': '-'}
    rev_comp = "".join(complement.get(nt, nt) for nt in reversed(query_sequence))
    return rev_comp


def output_retrieved_region(region, outfile=None):
    """
    Outputs the retrieved region
    :param region: A list of GenomeRegions where the query sequence aligns with the reference sequence
    :param outfile: The file stream for the output file in write mode
    """

    if outfile is None:
        print("\nRetrieved Region:\t{}".format(region.region_name))
        print("\tNucleotide Sequence:")
        seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
        for line in seq_lines:
            print('\t\t{}'.format(line))

        if region.aa_seq:
            print("\tProtein Sequence:")
            seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
            for line in seq_lines:
                print('\t\t{}\n'.format(line))

        print("\n\tRelative Positions:")

        if region.rel_pos['CDS']:
            print("\t\tNucleotide position relative to CDS start:\t{} --> {}"
                  .format(region.rel_pos['CDS'][0], region.rel_pos['CDS'][1]))
        else:
            print("\t\tNucleotide position relative to CDS start:\tN/A")

        if region.rel_pos['gstart']:
            print("\t\tNucleotide position relative to genome start:\t{} --> {}"
                  .format(region.rel_pos['gstart'][0], region.rel_pos['gstart'][1]))

        if region.rel_pos['pstart']:
            print("\t\tAmino acid position relative to protein start:\t{} --> {}"
                  .format(region.rel_pos['pstart'][0], region.rel_pos['pstart'][1]))

        if region.rel_pos['qstart']:
            print("\t\tPosition relative to query start:\t{} --> {}"
                  .format(region.rel_pos['qstart'][0], region.rel_pos['qstart'][1]))

    else:
        outfile.write("\nRetrieved Region:\t{}".format(region.region_name))
        outfile.write("\tNucleotide Sequence:\n")
        seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
        for line in seq_lines:
            outfile.write('\t\t{}\n'.format(line))

        if region.aa_seq:
            outfile.write("\tProtein Sequence:\n")
            seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
            for line in seq_lines:
                outfile.write('\t\t{}\n'.format(line))

        outfile.write("\n\tRelative Positions:")

        if region.rel_pos['CDS']:
            outfile.write("\t\tNucleotide position relative to CDS start:\t{} --> {}"
                          .format(region.rel_pos['CDS'][0], region.rel_pos['CDS'][1]))
        else:
            outfile.write("\t\tNucleotide position relative to CDS start:\tN/A")

        if region.rel_pos['gstart']:
            outfile.write("\t\tNucleotide position relative to genome start:\t{} --> {}"
                          .format(region.rel_pos['gstart'][0] + 1, region.rel_pos['gstart'][1] + 1))

        if region.rel_pos['pstart']:
            outfile.write("\t\tAmino acid position relative to protein start:\t{} --> {}"
                          .format(region.rel_pos['pstart'][0], region.rel_pos['pstart'][1]))

        if region.rel_pos['qstart']:
            outfile.write("\t\tPosition relative to query start:\t{} --> {}"
                          .format(region.rel_pos['qstart'][0], region.rel_pos['qstart'][1]))


def output_overlap(overlap_regions, outfile=None):
    """
    Outputs the regions that overlap with the query
    :param overlap_regions: A dictionary of GenomeRegions where the keys are the region name and
                                the values are the query sequence aligns with the reference sequence
    :param outfile: The file stream for the output file in write mode
    """

    if outfile is None:

        for key in overlap_regions:
            region = overlap_regions[key]

            print("\nRegion:\t{}".format(region.region_name))
            print("\n\tNucleotide Sequence:")
            seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
            for line in seq_lines:
                print('\t\t{}'.format(line))

            if region.aa_seq is not None:
                print("\n\tProtein Sequence:")
                seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.aa_seq), 60)]
                for line in seq_lines:
                    print('\t\t{}'.format(line))

            print("\n\tRelative Positions:")

            if region.cds_offset:
                print("\t\tNucleotide position relative to CDS start:\t{} --> {}"
                      .format(region.cds_offset[0], region.cds_offset[1]))
            else:
                print("\t\tNucleotide position relative to CDS start: N/A")

            if region.gstart:
                print("\t\tNucleotide position relative to genome start:\t{} --> {}"
                      .format(region.gstart[0], region.gstart[1]))
            if region.pstart:
                print("\t\tAmino acid position relative to protein start:\t{} --> {}"
                      .format(region.pstart[0], region.pstart[1]))
            if region.qstart:
                print("\t\tPosition relative to query start:\t{} --> {}"
                      .format(region.qstart[0], region.qstart[1]))

    else:
        for key in overlap_regions:
            region = overlap_regions[key]
            if region.region_name.startswith('5\'LTR'):
                outfile.write("\t3'LTR\n")

            outfile.write("\nRegion:\t{}".format(region.region_name))
            outfile.write("\n\tNucleotide Sequence:\n")
            seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.nt_seq), 60)]
            for line in seq_lines:
                outfile.write('\t\t{}'.format(line))

            if region.aa_seq is not None:
                outfile.write("\n\tProtein Sequence:\n")
                seq_lines = [region.nt_seq[i:i + 60] for i in range(0, len(region.aa_seq), 60)]
                for line in seq_lines:
                    print('\t\t{}'.format(line))

            outfile.write("\n\tRelative Positions:")

            if region.cds_offset:
                outfile.write("\t\tNucleotide position relative to CDS start:\t{} --> {}"
                              .format(region.cds_offset[0], region.cds_offset[1]))
            else:
                outfile.write("\t\tNucleotide position relative to CDS start: N/A")

            if region.gstart:
                print("\t\tNucleotide position relative to genome start:\t{} --> {}"
                      .format(region.gstart[0] + 1, region.gstart[1] + 1))
            if region.pstart:
                outfile.write("\t\tAmino acid position relative to protein start:\t{} --> {}"
                              .format(region.pstart[0], region.pstart[1]))
            if region.qstart:
                outfile.write("\t\tPosition relative to query start:\t{} --> {}"
                              .format(region.qstart[0], region.qstart[1]))


def handle_args(virus, base):
    """
    Handles the possible execution paths for the program
    :param virus: The reference virus
    :param base: The base of the reference sequence
    :return configs: a list containing the query sequence and paths to the reference nucleotide and protein sequences
    """
    if virus == 'hiv':
        ref_nt = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455.fasta")
        nt_coords = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455_genome_coordinates.csv")
        ref_aa = os.path.join(os.path.dirname(__file__), "ref_genomes/K03455-protein.fasta")

    else:
        ref_nt = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262.fasta")
        nt_coords = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262_genome_coordinates.csv")
        ref_aa = os.path.join(os.path.dirname(__file__), "ref_genomes/M33262-protein.fasta")

    ref_nt_seq = get_ref_seq(ref_nt, 'nucl')
    ref_aa_seq = get_ref_seq(ref_aa, 'prot')

    # Set the reference sequence to be used in the sequence alignment
    if base == 'nucl':
        reference_sequence = ref_nt_seq
    else:
        reference_sequence = ref_aa_seq

    # Convert list of lists to dictionary
    ref_prot_seq = {}
    for name, seq, in ref_aa_seq:
        ref_prot_seq[name.split('|')[0]] = seq      # Remove organism name

    configs = [ref_nt_seq, ref_prot_seq, nt_coords, reference_sequence]

    return configs


def parse_args():
    """
    Parses command line arguments
    """

    parser = argparse.ArgumentParser(
        description='An implementation of the HIV Sequence Locator tool by the Los Alamos National Laboratory.\n\n'
                    'This tool finds the location of the nucleotide or protein query sequence relative to the \n'
                    'HIV or SIV reference genomes; or retrieves a sequence in the HXB2 or SIVmm239 reference \n'
                    'genome from its coordinates.\n\n',
        epilog='For help using a specific sub-command enter:  sequence_locator.py <sub-command> -h',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(title='sub-commands', dest='subcommand')

    # Create sub-parser for 'align' mode
    parser_locate = subparsers.add_parser('locate',
                                          help='find the location of a sequence')
    parser_locate.add_argument('query', metavar='query',
                               help='the query sequence as a string or a FASTA file')
    parser_locate.add_argument('virus', metavar='virus', choices=['hiv', 'siv'],
                               help='the reference virus (choices: hiv, siv)')
    parser_locate.add_argument('base', metavar='base', choices=['nucl', 'prot'],
                               help='sequence base type (choices: \'nucl\' and \'prot\')')
    parser_locate.add_argument('-o', '--out', metavar='FILE', type=argparse.FileType('w'),
                               help='directs the output to a file (default: stdout)')
    parser_locate.add_argument('-rc', '--revcomp', action='store_true',
                               help='aligns the reverse complement of the query with the reference genome')

    # Create sub-parser for 'retrieve' mode
    parser_retrieve = subparsers.add_parser('retrieve',
                                            formatter_class=argparse.RawTextHelpFormatter,
                                            help='retrieve a region by its coordinates')
    parser_retrieve.add_argument('virus', metavar='virus', choices=['hiv', 'siv'],
                                 help='the reference virus (choices: hiv, siv)')
    parser_retrieve.add_argument('base', metavar='base', choices=['nucl', 'prot'],
                                 help='sequence base type (choices: \'nucl\' and \'prot\')')
    parser_retrieve.add_argument('-r', '--region', metavar='REG', default='Complete',
                                 help='list of genomic regions \n'
                                      'Accepted regions are: \n\t '
                                      '5\'LTR         TAR                 Gag-Pol             Gag\t\n\t ' 
                                      'Matrix        Capsid              p2                  Nucleocapsid\t\n\t ' 
                                      'p1            p6                  Pol                 GagPolTF\t\n\t ' 
                                      'Protease      RT                  RNase               Integrase\t\n\t '
                                      'Vif           Vpr                 Tat(with intron)    Tat(exon1)\t\n\t ' 
                                      'Tat(exon2)    Rev(with intron)    Rev(exon1)          Rev(exon2)\t\n\t ' 
                                      'Vpu           Vpx                 Env                 gp160\t\n\t '
                                      'V1            V2                  V3                  V4\t\n\t '
                                      'V5            RRE                 gp120               gp41\t\n\t '
                                      'Nef           3\'LTR               Complete\t\n\t')

    parser_retrieve.add_argument('-s', '--start', metavar='ST', default=1,
                                 help='start coordinate (default: 1)')
    parser_retrieve.add_argument('-e', '--end', metavar='END', default='end',
                                 help='end coordinate, either an integer or \'end\' (deault: end)')
    parser_retrieve.add_argument('-o', '--out', metavar='FILE', type=argparse.FileType('w'),
                                 help='directs the output a file (default: stdout)')

    # Show help text if no arguments or 1 arguments are given
    if len(sys.argv) == 1 or len(sys.argv) == 2:
        print("\033[1mSequence Locator\033[0m")
        parser.print_help()
        print("\n{}".format("-" * 80))
        print("\n\033[1m'locate' sub-command:\033[0m")
        parser_locate.print_help()
        print("\n{}".format("-" * 80))
        print("\n\033[1m'retrieve' sub-command:\033[0m")
        parser_retrieve.print_help()
        sys.exit(2)

    return parser.parse_args()


def main():
    args = parse_args()

    # Ensure proper configuration files are set
    configs = handle_args(args.virus, args.base)
    ref_nt_seq, ref_aa_seq = configs[0][0][1], configs[1]
    nt_coords = configs[2]
    reference_sequence = configs[3]

    with open(nt_coords) as nt_coords_handle:
        # Create genome object based on configuration files
        ref_genome = Genome(args.virus, nt_coords_handle, ref_nt_seq, ref_aa_seq, reference_sequence, args.base)

        if args.subcommand == "locate":
            if args.base == 'nucl':
                query_sequences = get_query(args.base, args.query, args.revcomp)
            else:
                query_sequences = get_query(args.base, args.query, False)

            # Handle multiple query sequences
            for seq in query_sequences:
                alignment = ref_genome.sequence_align(seq, args.out)

                # Find where the query sequence aligns with the reference sequence
                lookup_table = make_aln_lookup_table(alignment)
                query_match_coords = ref_genome.find_query_match_coords(alignment['query'])
                query_regions = ref_genome.find_matches(args.base, query_match_coords, lookup_table)
                output_overlap(query_regions, args.out)

        # Retrieve Mode
        else:
            valid_in = valid_inputs(args.virus, args.start, args.end, args.region)

            if not valid_in:
                sys.exit(0)
            else:
                result = ref_genome.retrieve(args.base, args.region, args.start, args.end)
                query_region = result[0]
                overlap_regions = result[1]

                output_retrieved_region(query_region, args.out)        # Output retrieved region first
                if args.outfile is None:
                    print("\n\nRegions touched by the query sequence:")
                else:
                    args.outfile.write("\n\nRegions touched by the query sequence:\n\n")

                output_overlap(overlap_regions, args.out)


if __name__ == '__main__':
    main()
