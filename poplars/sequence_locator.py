#!/usr/bin/env python3
"""
Python implementation of HIV Sequence Locator from https://www.hiv.lanl.gov
HIV and SIV genomic region coordinates are based on the HXB2 and SMM239
annotation resources from https://www.hiv.lanl.gov/content/sequence/HIV/MAP/annotation.html

Note: The first 256 nucleotides of SMM239 correspond to the flanking sequence,
and are included in the complete SIV genome (https://www.ncbi.nlm.nih.gov/nucleotide/M33262)
"""

import argparse
import os
import re
import sys

from poplars.common import convert_fasta, resolve_mixtures
from poplars.mafft import align

NON_CODING = ["5'LTR", "TAR", "3'LTR"]

CODON_DICT = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
              'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
              'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
              'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
              'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
              'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
              'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
              'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
              'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
              'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
              'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
              'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
              'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
              'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
              'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
              'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}


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
        if base == 'NA':
            return self.ncoords
        else:
            return self.pcoords

    def get_sequence(self, base):
        if base == 'NA':
            return self.nt_seq
        else:
            return self.aa_seq

    def set_coords(self, coord_pair, base):
        if base == 'NA':
            self.ncoords = coord_pair
        else:
            self.pcoords = coord_pair

    def set_nt_seq_from_genome(self):
        self.nt_seq = self.genome.nt_seq[self.ncoords[0] - 1: self.ncoords[1]]

    def set_aa_seq_from_genome(self):
        self.aa_seq = self.genome.aa_seq[self.region_name]

    def set_sequence(self, base, seq):
        if base == 'NA':
            self.nt_seq = seq
        else:
            self.aa_seq = seq


class RefRegion(Region):
    """
    Stores information about each region in the reference genome
    """

    def __init__(self, region_name, genome, ncoords=None, pcoords=None):
        super(RefRegion, self).__init__(region_name, genome, ncoords, pcoords)
        self.codon_aln = None

    def make_codon_aln(self):
        """
        Aligns the protein sequence with its associated nucleotide sequence
        """
        if self.nt_seq is not None and self.aa_seq is not None:

            # Account for the frame-shift at position 5572 in vpr (extra T relative to other subtype B)
            if self.region_name == 'Vpr' and self.genome.virus == 'hiv':
                pos = 215  # position of extra T relative to the start of the Vpr CDS (0-based indexing)
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

    def find_overlap_coords(self, base, global_qcoords):
        """
        Gets the sequence regions that overlap with the region of interest
        :param base: the base of the sequence (nucleotide or protein)
        :param global_qcoords: the global coordinates of the overlap region
        :return overlap: A coordinates indicating where the query overlaps with the reference region
        """
        ref_coords = self.get_coords(base)  # Global coordinates

        # If the ref_coords are in the range of the query region and the q_coords are in the range of the ref region
        if not ((ref_coords[0] < ref_coords[1] < global_qcoords[0]) or
                (global_qcoords[0] < global_qcoords[1] < ref_coords[0])):

            # If the end of the query region exceeds the end of the reference region
            if global_qcoords[1] > ref_coords[1]:
                end = ref_coords[1]
            else:
                end = global_qcoords[1]

            # If the query region starts before the reference region starts
            if global_qcoords[0] < ref_coords[0]:
                start = ref_coords[0]
            else:
                start = global_qcoords[0]

            overlap_coords = [start, end]
            return overlap_coords
        return None

    def local_to_global_coords(self, base, local_coords):
        """
        Convert local coordinates to global coordinates
        :param base: the base of local coordinates
        :param local_coords: the local coordinates
        :return: the global coordinates
        """
        if base == 'NA':
            start = local_coords[0] + self.ncoords[0] - 1
            end = local_coords[1] + self.ncoords[1] - 1
        else:
            start = local_coords[0] + self.pcoords[0] - 1
            end = local_coords[1] + self.pcoords[1] - 1

        return [start, end]


class OverlapRegion(Region):
    """
    Represents information about the overlap sequence region
    """

    def __init__(self, region_name, ref_region, genome, query_coords, ncoords=None, pcoords=None):
        """
        :param region_name: name of the query region
        :param ref-region: reference to the genome region
        :param genome: reference to the genome
        :param query_coords: the local coordinates of the region
        :param ncoords: the global coordinates of the nucleotide sequence
        :param pcoords: the global coordinates of the protein sequence
        """
        super(OverlapRegion, self).__init__(region_name, genome, ncoords, pcoords)

        self.ref_region = ref_region
        self.codon_aln = ''
        self.query_coords = query_coords
        self.aln_coords = None
        self.cds_offset = None
        self.qstart = None

    def set_ncoords_from_pcoords(self):
        """
        Sets nucleotide coordinates given the protein coordinates relative to the protein start
        :return ncoords: the nucleotide coordinates
        """
        if self.pcoords is not None:
            nucl_start = self.pcoords[0] - 1 * 3
            nucl_end = self.pcoords[1] * 3
            self.ncoords = [nucl_start, nucl_end]

    def set_pcoords_from_ncoords(self):
        """
        Sets protein coordinates relative to the protein start, given the nucleotide coordinates
        """
        if self.cds_offset is not None:
            prot_start = self.cds_offset[0] // 3 + 1
            prot_end = self.cds_offset[1] // 3
            self.pcoords = [prot_start, prot_end]

    def find_aln_coords(self, base, lookup_table):
        """
        Gives the position of the overlap relative to the alignment
        :param base: the sequence base type
        :param lookup_table: The mapping of the query and reference sequences to the alignment
        """
        global_q_coords = self.get_coords(base)
        return (lookup_table['reference'].index(global_q_coords[0]),
                lookup_table['reference'].index(global_q_coords[1] - 1))

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
                start = query_start - ref_start + 1  # 1-based indexing

            if query_end == ref_end:
                end = query_end - query_start + 1
            else:
                end = start + (query_end - query_start)

            self.cds_offset = [start, end]

    def set_pos_from_qstart(self, base, lookup_table):
        """
        Gives the position of the sequence relative to the start of the region of interest
        :param base: The base of the sequence (nucleotide or protein)
        :param lookup_table: The mapping of the query and reference sequences to the alignment
        :return: The position of the reference region relative to the query region
        """
        aln_coords = self.find_aln_coords(base, lookup_table)
        self.qstart = [lookup_table['query'][aln_coords[0]], lookup_table['query'][aln_coords[1]] + 1]

    def find_protein_equivalent(self, query_seq):
        """
        Gets the protein equivalent of the query sequence
        :param query_seq: the query sequence
        """
        q_frame = self.cds_offset[0] % 3
        prot_seq = ''

        # Handle partial codons at the beginning
        pos = self.qstart[0] - 1
        if q_frame != 1:  # +0 frame
            prot_seq += 'X'

            # Shift to next codon
            if q_frame == 0:  # +2 frame
                pos += 1
            else:  # + 1 frame
                pos += 2

        # Translate the query sequence
        for i in range(pos, self.qstart[1], 3):
            if self.qstart[1] - i < 3:
                prot_seq += 'X'
            else:
                prot_seq += CODON_DICT[query_seq[i:i + 3]]

        return prot_seq


class Genome:
    """
    Stores information about the reference genome
    """

    def __init__(self, virus, nt_coords, nt_seq, aa_seq, reference_sequence, base):
        self.virus = virus
        self.ref_name = 'HXB2' if virus == 'hiv' else 'SMM239'
        self.nt_seq = nt_seq
        self.aa_seq = aa_seq
        self.reference_sequence = reference_sequence
        self.ref_base = base
        self.ref_genome_regions = self.make_ref_regions(nt_coords, aa_seq)

    def make_ref_regions(self, nt_coords, aa_seq):
        """
        Reads in the start and end coordinates of the genomic regions and associates the region with its coordinates.
        If no coordinate files are specified, set default nucleotide and protein coordinate files.
        :param nt_coords: File stream containing coordinates for nucleotide regions. Format: region_name,start,end
        :param aa_seq: A dictionary of protein sequences
        """
        ref_regions = {}
        # Parse nucleotide region coordinates file
        for nt_line in nt_coords:
            nt_line = nt_line.strip()
            nt_line = nt_line.split(',')
            nucl_coords = [int(nt_line[1]), int(nt_line[2])]
            seq_region = RefRegion(nt_line[0], self, nucl_coords)
            seq_region.set_nt_seq_from_genome()
            seq_region.codon_aln = seq_region.make_codon_aln()
            ref_regions[nt_line[0]] = seq_region

        # Match protein regions to nucleotide regions
        for region_name in ref_regions:
            for prot_name in aa_seq.keys():
                if prot_name.startswith(region_name):
                    # Set protein coordinates and sequences
                    ref_regions[region_name].set_sequence('AA', aa_seq[region_name])
                    ref_regions[region_name].set_coords([1, len(aa_seq[region_name])], 'AA')

        # Align nucleotide sequence with protein sequence
        for name in ref_regions:
            ref_regions[name].make_codon_aln()

        return ref_regions

    def sequence_align(self, query):
        """
        Aligns the query sequence to the reference genome
        :param query: The query object
        """
        result = ''
        if self.ref_base == 'NA':
            result = align(query.query_sequence, self.reference_sequence[0][1], True)
        else:
            for name, region in self.ref_genome_regions.items():
                aa_seq = region.get_sequence('AA')
                if aa_seq is not None:
                    result = align(query.query_sequence, aa_seq, True)

        return result

    def output_alignment(self, query, outfolder):
        """
        Outputs the alignment of the query reference to the reference genome
        :param query: the query object
        :param outfolder: The path to the output folder
        """
        q = query.alignment['query'][query.qcoords[0] - 1: query.qcoords[1]]
        r = query.alignment['reference'][query.qcoords[0] - 1: query.qcoords[1]]
        a = query.alignment['aln'][query.qcoords[0] - 1: query.qcoords[1]]

        if outfolder is None:
            print("\nAlignment of the query sequence to {}".format(self.ref_name))

            for i in range(0, len(q), 60):
                qline = q[i:i + 60]
                rline = r[i:i + 60]
                aline = a[i:i + 60]
                end_pos = query.qcoords[0] - 1 + i + len(qline) - 1
                print('\tQuery\t {}\t{}\n'
                      '\t         {}\n'
                      '\tRef  \t {}\t{}'
                      .format(qline, query.lookup_table['query'][end_pos] + 1, aline, rline,
                              query.lookup_table['reference'][end_pos] + 1))

        else:
            with open(os.path.join(outfolder, '{}.aln'.format(query.name)), 'w') as aln_handle:
                for i in range(0, len(q), 60):
                    qline = q[i:i + 60]
                    rline = r[i:i + 60]
                    aline = a[i:i + 60]
                    end_pos = query.qcoords[0] - 1 + i + len(qline) - 1
                    aln_handle.write('Query\t {}\t{}\n'
                                     '         {}\n'
                                     'Ref  \t {}\t{}\n\n'
                                     .format(qline, query.lookup_table['query'][end_pos] + 1, aline,
                                             rline, query.lookup_table['reference'][end_pos] + 1))


class Query(object):
    """
    Represents the query region
    """
    _count = 0

    def __init__(self, base, genome, outfile=None, query_sequence=None, region_name=None):
        Query._count += 1
        self.query_num = Query._count
        self.base = base
        self.genome = genome

        if query_sequence is not None:
            self.query_sequence = query_sequence[0][1]

            # If the query has no header, create one
            if query_sequence[0][0] == 'query':
                self.name = query_sequence[0][0] + str(Query._count)
            else:
                self.name = query_sequence[0][0]

            self.alignment = self.genome.sequence_align(self)

        self.overlap_regions = None
        self.qregion = region_name
        self.qcoords = self.find_query_match_coords()
        self.lookup_table = self.make_lookup_table()

    def find_query_match_coords(self):
        """
        Gets the indices of regions where the query aligns with the reference sequence
        :return query_matches: The positions where the query aligns with the reference sequences with no gaps
        """
        pat = re.compile('\w[\w-]*\w')  # Get the query sequence from the alignment
        match = pat.search(self.alignment['query'])
        query_match_coords = (match.start() + 1, match.end())
        return query_match_coords

    def make_lookup_table(self):
        """
        Maps the coordinates of reference genome and query sequence to the alignment
        :return: dictionary of of mappings from positions in aligned sequence
                 to positions in the query and reference sequences
        """
        lookup_table = {'reference': [], 'query': []}
        r_count = -1
        q_count = -1

        for r_nt, q_nt in zip(self.alignment['reference'], self.alignment['query']):
            # Account for gaps in the reference sequence
            if r_nt != '-':
                r_count += 1
            lookup_table['reference'].append(r_count)

            # Account for gaps in the query sequence
            if q_nt != '-':
                q_count += 1
            lookup_table['query'].append(q_count)

        return lookup_table

    def find_location(self):
        """
        Finds the location of the query sequence in the genome
        """

        # Convert query coordinates from "alignment space" to "region space"
        # Non-inclusive 0-based indexing
        global_qcoords = (self.lookup_table['reference'][self.qcoords[0]],
                          self.lookup_table['reference'][self.qcoords[1] - 1] + 1)

        overlap_regions = {}
        for name in self.genome.ref_genome_regions:
            if name != 'Complete':
                # Find which regions overlap with the query region
                overlap_coords = self.genome.ref_genome_regions[name].find_overlap_coords(self.base, global_qcoords)

                if overlap_coords is not None:

                    # Exclude overlaps that have insufficient length
                    if overlap_coords[1] - overlap_coords[0] < 3:
                        continue

                    overlap = OverlapRegion(name, self.genome.ref_genome_regions[name], self.genome, overlap_coords)

                    if self.base == 'NA':
                        overlap.set_coords(overlap_coords, self.base)  # Set nt_coords
                        overlap.set_nt_seq_from_genome()
                        overlap.set_pos_from_cds()
                        overlap.set_pos_from_qstart(self.base, self.lookup_table)
                        overlap.aln_coords = overlap.find_aln_coords(self.base, self.lookup_table)

                        # Set protein coordinates if the overlap is in a coding region
                        if overlap.region_name not in NON_CODING:
                            overlap.set_pcoords_from_ncoords()
                            prot_equiv = overlap.find_protein_equivalent(self.query_sequence)
                            overlap.set_sequence('AA', prot_equiv)

                    else:
                        overlap.set_coords(overlap_coords, self.base)  # Set aa_cords
                        overlap.set_aa_seq_from_genome()
                        overlap.set_pos_from_cds()
                    overlap_regions[name] = overlap

        self.overlap_regions = overlap_regions

    def output_overlap(self, outfolder):
        """
        Outputs the regions that overlap with the query
        :param outfolder: The path to the output folder
        """
        if outfolder is None:
            print("\nRegions touched by {}".format(self.name))

            for key in self.overlap_regions:
                region = self.overlap_regions[key]

                print("Region:\t{}".format(region.region_name))

                if region.aa_seq is not None:
                    print("\tProtein Translation of the Query:")
                    seq_lines = [region.aa_seq[i:i + 60] for i in range(0, len(region.aa_seq), 60)]
                    for line in seq_lines:
                        print('\t\t{}\n'.format(line))

                print("\tRelative Positions:")
                if region.cds_offset:
                    print("\t\tNA position relative to CDS start in {}: {},{}"
                          .format(self.genome.ref_name, region.cds_offset[0], region.cds_offset[1]))

                    # Compare query region length and CDS length
                    query_len = region.qstart[1] - region.qstart[0] + 1
                    cds_len = region.cds_offset[1] - region.cds_offset[0] + 1
                    if cds_len < query_len:
                        print('\t\t\tNotice: length of {} portion of query ({}) is greater than its '
                              'length in {} ({}) - possible frameshift'
                              .format(region.region_name, query_len, self.genome.ref_name, cds_len))
                else:
                    print("\t\tNA position relative to CDS start in {}: N/A".format(self.genome.ref_name))

                if region.qstart:
                    print("\t\t{} position relative to query start in {}: {},{}"
                          .format(self.base, self.genome.ref_name, region.qstart[0], region.qstart[1]))
                if region.ncoords:
                    print("\t\tNA position relative to genome start in {}: {},{}"
                          .format(self.genome.ref_name, region.ncoords[0], region.ncoords[1]))
                if region.pcoords:
                    print("\t\tAA position relative to protein start in {}: {},{}\n"
                          .format(self.genome.ref_name, region.pcoords[0], region.pcoords[1]))

        else:
            coords_path = os.path.join(outfolder, 'overlap_coordinates.csv')

            with open(coords_path, 'a+') as out_handle:
                for key in self.overlap_regions:
                    region = self.overlap_regions[key]
                    out_handle.write('{},{},{},'.format(self.name, self.query_sequence, region.region_name))

                    if region.cds_offset:
                        out_handle.write('{}:{},'.format(region.cds_offset[0], region.cds_offset[1]))
                    else:
                        out_handle.write('{},'.format('N/A'))

                    if region.qstart:
                        out_handle.write('{}:{},'.format(region.qstart[0], region.qstart[1]))
                    if region.ncoords:
                        out_handle.write('{}:{},'.format(region.ncoords[0], region.ncoords[1]))
                    if region.pcoords:
                        out_handle.write('{}:{}\n'.format(region.pcoords[0], region.pcoords[1]))

    def retrieve_region(self, region, start, end):
        """
        Retrieves a sequence given its coordinates
        :param region: name of the region to retrieve
        :param start: the start coordinate relative to the start of the region
        :param end: the end coordinate relative to the end of the region
        """
        reg_start, reg_end = 0, 0

        if self.base == 'NA':

            # Set the value of end based on the region name
            if end == 'end':
                if self.qregion.lower() == 'tat':
                    reg_end = 8653  # The last position in the Tat2 exon

                elif self.qregion.lower() == 'rev':
                    reg_end = 8469  # The last position in the Rev2 exon

                else:
                    global_qcoords = self.genome.ref_genome_regions[self.qregion].get_coords(self.base)
                    reg_end = global_qcoords[1] - global_qcoords[0]

            else:
                # Check if the start coordinate falls in the Tat1 exon
                if self.qregion.lower() == 'tat':
                    tat2_coords = self.genome.ref_genome_regions['Tat(exon2)'].get_coords(self.base)

                    if start > len(self.genome.ref_genome_regions['Tat(exon1)'].nt_seq):
                        reg_start = tat2_coords[0]
                        reg_end = end

                # Check if the start coordinate falls in the Tat2 exon
                elif self.qregion.lower() == 'rev':
                    rev2_coords = self.genome.ref_genome_regions['Rev(exon2)'].get_coords(self.base)

                    if start > len(self.genome.ref_genome_regions['Rev(exon1)'].nt_seq):
                        reg_start = rev2_coords[0]
                        reg_end = end

        else:
            reg_start, reg_end = start, end

        # Convert local coordinates to global coordinates
        global_range = self.genome.ref_genome_regions[region].local_to_global(self.base, [reg_start, reg_end])

        overlap_regions = {}
        for name in self.genome.ref_genome_regions:
            if name != 'Complete':
                # Find which regions overlap with the query region
                overlap_coords = self.genome.ref_genome_regions[region].find_overlap_coords(self.base, global_range)

                if overlap_coords is not None:

                    # Exclude overlaps that have insufficient length
                    if overlap_coords[1] - overlap_coords[0] < 3:
                        continue

                    overlap = OverlapRegion(name, self.genome.ref_genome_regions[name], self.genome, overlap_coords)

                    if self.base == 'NA':
                        overlap.set_coords(overlap_coords, self.base)  # Set nt_coords
                        overlap.set_nt_seq_from_genome()
                        overlap.set_pos_from_cds()
                        # overlap.set_pos_from_qstart(self.base)
                        # overlap.aln_coords = overlap.find_aln_coords(self.base, lookup_table)

                        # Set protein coordinates if the overlap is in a coding region
                        if overlap.region_name not in NON_CODING:
                            overlap.set_pcoords_from_ncoords()
                            prot_equiv = overlap.find_protein_equivalent(self.query_sequence)
                            overlap.set_sequence('AA', prot_equiv)

                    else:
                        overlap.set_coords(overlap_coords, self.base)  # Set aa_cords
                        overlap.set_aa_seq_from_genome()
                        overlap.set_pos_from_cds()
                    overlap_regions[name] = overlap

        self.overlap_regions = overlap_regions


def output_retrieved_region(base, region, outfile=None):
    """
    Outputs the retrieved region
    :param base: The sequence base type (NA or AA)
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
        if region.cds_offset:
            print("\t\tNA position relative to CDS start: {} --> {}"
                  .format(region.cds_offset[0], region.cds_offset[1]))
        else:
            print("\t\tNA position relative to CDS start: N/A")

        if region.qstart:
            print("\t\t{} position relative to query start: {} --> {}"
                  .format(base, region.qstart[0], region.qstart[1]))
        if region.ncoords:
            print("\t\tNA position relative to genome start: {} --> {}"
                  .format(region.ncoords[0], region.ncoords[1]))
        if region.pcoords:
            print("\t\tAA position relative to protein start: {} --> {}"
                  .format(region.pcoords[0], region.pcoords[1]))

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
        if region.cds_offset:
            outfile.write("\t\tNA position relative to CDS start: {},{}"
                          .format(region.cds_offset[0], region.cds_offset[1]))
        else:
            outfile.write("\t\tNA position relative to CDS start: N/A")
        if region.qstart:
            outfile.write("\t\t{} position relative to query start: {},{}"
                          .format(base, region.qstart[0], region.qstart[1]))
        if region.ncoords:
            outfile.write("\t\tNA position relative to genome start: {},{}"
                          .format(region.ncoords[0] + 1, region.ncoords[1] + 1))
        if region.pcoords:
            outfile.write("\t\tAA position relative to protein start: {},{}"
                          .format(region.pcoords[0], region.pcoords[1]))


def valid_sequence(base, sequence):
    """
    Verifies that input sequence is valid
    :param base: The base of the sequence (NA or AA)
    :param sequence: A list of lists containing header and sequence pairs
    :param verbatim:  if True, reject any nucleotide sequence with mixtures
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
        elif base == 'NA':
            if not all(pos in dna_alphabet for pos in s):
                print("Invalid nucleotide sequence:\n{}\n{}\n".format(h, s))
                return False
        elif base == 'AA':
            if not all(pos in aa_alphabet for pos in s):
                print("Invalid amino acid sequence:\n{}\n{}\n".format(h, s))
                return False
        else:
            print("Unexpected base argument {} in valid_sequence()".format(base))
            sys.exit()
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


def reverse_comp(query_sequence):
    """
    Reverses and complements the query sequence
    :param query_sequence: The query sequence
    :return: The reverse complement of the query sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '*': '*', 'N': 'N', '-': '-'}
    rev_comp = "".join(complement.get(nt, nt) for nt in reversed(query_sequence))
    return rev_comp


def get_query(base, query, rev_comp, verbatim=False):
    """
    Gets the query sequence and checks that it is valid
    :param base: The base (nucleotide or protein)
    :param query: The query sequence as a string or the file path to the query sequence
    :param rev_comp: Reverse complement flag (False by default)
    :param verbatim: if True, reject any nucleotide sequences with mixtures
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
        if base == 'NA' and not verbatim:
            # attempt to salvage sequence with mixtures
            resolved = []
            for h, s in query_seq:
                rs = resolve_mixtures(s)
                if rs is None:
                    print("Failed to resolve mixtures in {}".format(h))
                    sys.exit()
                resolved.append([h, rs])
            query_seq = resolved
        else:
            sys.exit(0)

    # At this point, the sequence is valid
    if base == 'NA' and rev_comp:
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

    ref_nt_seq = get_ref_seq(ref_nt, 'NA')
    ref_aa_seq = get_ref_seq(ref_aa, 'AA')

    # Set the reference sequence to be used in the sequence alignment
    if base == 'NA':
        reference_sequence = ref_nt_seq
    else:
        reference_sequence = ref_aa_seq

    # Convert list of lists to dictionary
    ref_prot_seq = {}
    for name, seq, in ref_aa_seq:
        ref_prot_seq[name.split('|')[0]] = seq  # Remove organism name

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

    # Create sub-parser for 'locate' mode
    parser_locate = subparsers.add_parser('locate',
                                          help='find the location of a sequence')
    parser_locate.add_argument('virus', metavar='virus', choices=['hiv', 'siv'],
                               help='the reference virus (choices: hiv, siv)')
    parser_locate.add_argument('base', metavar='base', choices=['NA', 'AA'],
                               help='sequence base type (choices: \'NA\' and \'AA\')')
    parser_locate.add_argument('query', metavar='query', nargs='+',
                               help='the query sequence as a string or a FASTA file')
    parser_locate.add_argument('-o', '--out', metavar='DIRECTORY',
                               help='directs the output to the specified directory (default: stdout) '
                                    'and creates an alignment and an output file for each query')
    parser_locate.add_argument('-rc', '--revcomp', action='store_true',
                               help='aligns the reverse complement of the query with the reference genome')
    parser_locate.add_argument('--verbatim', action='store_true',
                               help='no tolerance for ambiguous base calls, i.e., mixtures (R=A/G), '
                                    'exits gracefully')

    # Create sub-parser for 'retrieve' mode
    parser_retrieve = subparsers.add_parser('retrieve',
                                            formatter_class=argparse.RawTextHelpFormatter,
                                            help='retrieve a region by its coordinates')
    parser_retrieve.add_argument('virus', metavar='virus', choices=['hiv', 'siv'],
                                 help='the reference virus (choices: hiv, siv)')
    parser_retrieve.add_argument('base', metavar='base', choices=['NA', 'AA'],
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
                                 help='end coordinate, either an integer or \'end\' (default: end)')
    parser_retrieve.add_argument('-o', '--out', metavar='DIRECTORY',
                                 help='directs the output to the specified directory (default: stdout)'
                                      'and creates an output file for the query')

    # Show help text if no arguments or 1 arguments are given
    if len(sys.argv) == 1 or len(sys.argv) == 2:
        print("\033[1mSequence Locator\033[0m")
        parser.print_help()
        sys.exit(2)

    return parser.parse_args()


def main():
    args = parse_args()

    # Ensure proper configuration files are set
    configs = handle_args(args.virus, args.base)
    ref_nt_seq, ref_aa_seq = configs[0][0][1], configs[1]
    nt_coords = configs[2]
    reference_sequence = configs[3]

    with open(nt_coords, 'r') as nt_coords_handle:
        ref_genome = Genome(args.virus, nt_coords_handle, ref_nt_seq, ref_aa_seq, reference_sequence, args.base)

        # Create output directory if it does not already exist
        if args.out:
            if not os.path.isdir(args.out):
                os.makedirs(args.out)

        # Finds the location of a sequence (locate mode)
        if args.subcommand == "locate":
            # Prepare overlap coordinates output file
            if args.out:
                coords_path = os.path.join(args.out, 'overlap_coordinates.csv')
                with open(coords_path, 'w') as out_handle:
                    out_handle.write('query,query_seq,region,CDS_coords,query_coords,genome_coords,prot_coords\n')

            # Handle multiple query sequences
            for qseq in args.query:
                if args.base == 'NA':
                    query_seq = get_query(args.base, qseq, args.revcomp)
                else:
                    query_seq = get_query(args.base, qseq, False)

                query = Query(args.base, ref_genome, outfile=args.out, query_sequence=query_seq)
                query.find_location()
                query.output_overlap(args.out)
                ref_genome.output_alignment(query, args.out)

        # Retrieves a region by its coordinates (retrieve mode)
        else:
            # Validate input
            valid_in = valid_inputs(args.virus, args.start, args.end, args.region)
            if not valid_in:
                sys.exit(0)

            else:
                query = Query(args.base, ref_genome, outfile=args.out, region_name=args.region)
                output_retrieved_region(query, args.region, args.out)

                if args.outfile is None:
                    print("\n\nRegions touched by the query sequence:")
                else:
                    args.outfile.write("\n\nRegions touched by the query sequence:\n\n")

                query.output_overlap(args.out)


if __name__ == '__main__':
    main()
