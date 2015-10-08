#!/usr/bin/env python3
"""
Read a VCF and a BAM file and phase the variants. The phased VCF is written to
standard output.
"""
"""
 0: ref allele
 1: alt allele
 -: unphasable: no coverage of read that covers at least 2 SNPs
 X: unphasable: there is coverage, but still not phasable (tie)

TODO
* it would be cleaner to not open the input VCF twice
* convert parse_vcf to a class so that we can access VCF header info before
  starting to iterate (sample names)
"""
import os
import logging
import sys
import random
import gzip
import time
import itertools
import platform
from collections import defaultdict, Counter
from contextlib import contextmanager

try:
	from contextlib import ExitStack
except ImportError:
	from contextlib2 import ExitStack  # PY32
from ..vcf import parse_vcf, PhasedVcfWriter, remove_overlapping_variants
from .. import __version__
from ..args import HelpfulArgumentParser as ArgumentParser
from ..core import Read, ReadSet, DPTable, readselection
from ..graph import ComponentFinder
from ..coverage import CovMonitor
from ..bam import MultiBamReader, SampleBamReader, BamIndexingError, SampleNotFoundError, HaplotypeBamWriter


__author__ = "Murray Patterson, Alexander Schönhuth, Tobias Marschall, Marcel Martin"

logger = logging.getLogger(__name__)


def covered_variants(variants, start, bam_read, source_id):
	"""
    Find the variants that are covered by the given bam_read and return a
    core.Read instance that represents those variants. The instance may be
    empty.

    start -- index of the first variant (in the variants list) to check
    """
	core_read = Read(bam_read.qname, bam_read.mapq, source_id)
	i = 0  # index into CIGAR
	j = start  # index into variants
	ref_pos = bam_read.pos  # position relative to reference
	query_pos = 0  # position relative to read

	for cigar_op, length in bam_read.cigar:
		# The mapping of CIGAR operators to numbers is:
		# MIDNSHPX= => 012345678
		if cigar_op in (0, 7, 8):  # we are in a matching region
			# Skip variants that come before this region
			while j < len(variants) and variants[j].position < ref_pos:
				j += 1

			# Iterate over all variants that are in this region
			while j < len(variants) and variants[j].position < ref_pos + length:
				if len(variants[j].reference_allele) == len(variants[j].alternative_allele) == 1:
					# Variant is a SNP
					offset = variants[j].position - ref_pos
					base = bam_read.seq[query_pos + offset]
					allele = None
					if base == variants[j].reference_allele:
						allele = 0
					elif base == variants[j].alternative_allele:
						allele = 1
					if allele is not None:
						# TODO
						# Fix this: we can actually have indel and SNP
						# calls at identical positions. For now, ignore the
						# second variant.
						if variants[j].position in core_read:
							logger.debug("Found two variants at identical positions. Ignoring the second one: %s",
										 variants[j])
						else:
							# Do not use bam_read.qual here as it is extremely slow.
							# If we ever decide to be compatible with older pysam
							# versions, cache bam_read.qual somewhere - do not
							# access it within this loop (3x slower otherwise).
							core_read.add_variant(variants[j].position, allele,
												  bam_read.query_qualities[query_pos + offset])
				elif len(variants[j].reference_allele) == 0:
					assert len(variants[j].alternative_allele) > 0
					# This variant is an insertion. Since we are in a region of
					# matches, the insertion was *not* observed (reference allele).
					qual = 30  # TODO average qualities of "not inserted" bases?
					core_read.add_variant(variants[j].position, allele=0, quality=qual)
				elif len(variants[j].alternative_allele) == 0:
					assert len(variants[j].reference_allele) > 0
					# This variant is a deletion that was not observed.
					# Add it only if the next variant is not located 'within'
					# the deletion.
					deletion_end = variants[j].position + len(variants[j].reference_allele)
					if not (j + 1 < len(variants) and variants[j + 1].position < deletion_end):
						qual = 30  # TODO
						core_read.add_variant(variants[j].position, allele=0, quality=qual)
					else:
						logger.info('Skipped a deletion overlapping another variant at pos. %d', variants[j].position)
						# Also skip all variants that this deletion overlaps
						while j + 1 < len(variants) and variants[j + 1].position < deletion_end:
							j += 1
						# One additional j += 1 is done below
				else:
					assert False, "Strange variant: {}".format(variants[j])
				j += 1
			query_pos += length
			ref_pos += length
		elif cigar_op == 1:  # an insertion
			# Skip variants that come before this region
			while j < len(variants) and variants[j].position < ref_pos:
				j += 1
			if j < len(variants) and variants[j].position == ref_pos and \
							len(variants[j].reference_allele) == 0 and \
							variants[j].alternative_allele == bam_read.seq[query_pos:query_pos + length]:
				qual = 30  # TODO
				assert variants[j].position not in core_read
				core_read.add_variant(variants[j].position, allele=1, quality=qual)
				j += 1
			query_pos += length
		elif cigar_op == 2:  # a deletion
			# Skip variants that come before this region
			while j < len(variants) and variants[j].position < ref_pos:
				j += 1
			# We only check the length of the deletion, not the sequence
			# that gets deleted since we don’t have the reference available.
			# (We could parse the MD tag if it exists.)
			if j < len(variants) and variants[j].position == ref_pos and \
							len(variants[j].alternative_allele) == 0 and \
							len(variants[j].reference_allele) == length:
				qual = 30  # TODO
				deletion_end = variants[j].position + len(variants[j].reference_allele)
				if not (j + 1 < len(variants) and variants[j + 1].position < deletion_end):
					qual = 30  # TODO
					assert variants[j].position not in core_read
					core_read.add_variant(variants[j].position, allele=1, quality=qual)
				else:
					logger.info('Skipped a deletion overlapping another variant at pos. %d', variants[j].position)
					# Also skip all variants that this deletion overlaps
					while j + 1 < len(variants) and variants[j + 1].position < deletion_end:
						j += 1
					# One additional j += 1 is done below
				j += 1
			ref_pos += length
		elif cigar_op == 3:  # a reference skip
			ref_pos += length
		elif cigar_op == 4:  # soft clipping
			query_pos += length
		elif cigar_op == 5 or cigar_op == 6:  # hard clipping or padding
			pass
		else:
			logger.error("Unsupported CIGAR operation: %d", cigar_op)
			sys.exit(1)
	return core_read


class ReadSetReader:
	"""
    Associate VCF variants with BAM reads.
    """

	def __init__(self, paths, mapq_threshold=20):
		self._mapq_threshold = mapq_threshold
		if len(paths) == 1:
			self._reader = SampleBamReader(paths[0])
		else:
			self._reader = MultiBamReader(paths)

	def read(self, chromosome, variants, sample):
		"""
        chromosome -- name of chromosome to work on
        variants -- list of Variant objects (obtained from VCF with parse_vcf)
        sample -- name of sample to work on. If None, read group information is
            ignored and all reads in the file are used.

        Return a ReadSet object.
        """
		# Since variants are identified by position, positions must be unique.
		if __debug__ and variants:
			varposc = Counter(variant.position for variant in variants)
			pos, count = varposc.most_common()[0]
			assert count == 1, "Position {} occurs more than once in variant list.".format(pos)

		# Map read name to a list of Read objects. The list has two entries
		# if it is a paired-end read, one entry if the read is single-end.
		reads = defaultdict(list)

		i = 0  # keep track of position in variants array (which is in order)
		for alignment in self._reader.fetch(reference=chromosome, sample=sample):
			# TODO: handle additional alignments correctly! find out why they are sometimes overlapping/redundant
			if alignment.bam_alignment.flag & 2048 != 0:
				# print('Skipping additional alignment for read ', alignment.bam_alignment.qname)
				continue
			if alignment.bam_alignment.mapq < self._mapq_threshold:
				continue
			if alignment.bam_alignment.is_secondary:
				continue
			if alignment.bam_alignment.is_unmapped:
				continue
			if not alignment.bam_alignment.cigar:
				continue

			# Skip variants that are to the left of this read.
			while i < len(variants) and variants[i].position < alignment.bam_alignment.pos:
				i += 1

			core_read = covered_variants(variants, i, alignment.bam_alignment, alignment.source_id)
			# Only add new read if it covers at least one variant.
			if core_read:
				reads[(alignment.source_id, alignment.bam_alignment.qname)].append(core_read)

		# Prepare resulting set of reads.
		read_set = ReadSet()

		for readlist in reads.values():
			assert 0 < len(readlist) <= 2
			if len(readlist) == 1:
				read_set.add(readlist[0])
			else:
				read_set.add(self._merge_pair(*readlist))
		return read_set

	def _merge_pair(self, read1, read2):
		"""
        Merge the two ends of a paired-end read into a single core.Read. Also
        takes care of self-overlapping read pairs.

        TODO this can be simplified as soon as a variant in a read can be
        modified.
        """
		if read2:
			result = Read(read1.name, read1.mapqs[0], read1.source_id)
			result.add_mapq(read2.mapqs[0])
		else:
			return read1

		i1 = 0
		i2 = 0

		def add1():
			result.add_variant(read1[i1].position, read1[i1].allele, read1[i1].quality)

		def add2():
			result.add_variant(read2[i2].position, read2[i2].allele, read2[i2].quality)

		while i1 < len(read1) or i2 < len(read2):
			if i1 == len(read1):
				add2()
				i2 += 1
				continue
			if i2 == len(read2):
				add1()
				i1 += 1
				continue
			variant1 = read1[i1]
			variant2 = read2[i2]
			if variant2.position < variant1.position:
				add2()
				i2 += 1
			elif variant2.position > variant1.position:
				add1()
				i1 += 1
			else:
				# Variant on self-overlapping read pair
				assert read1[i1].position == read2[i2].position
				# If both alleles agree, merge into single variant and add up qualities
				if read1[i1].allele == read2[i2].allele:
					quality = read1[i1].quality + read2[i2].quality
					result.add_variant(read1[i1].position, read1[i1].allele, quality)
				else:
					# Otherwise, take variant with highest base quality and discard the other.
					if read1[i1].quality >= read2[i2].quality:
						add1()
					else:
						add2()
				i1 += 1
				i2 += 1
		return result

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()

	def close(self):
		self._reader.close()


def find_components(superreads, reads):
	"""
    Return a dict that maps each position to the component it is in. A
    component is identified by the position of its leftmost variant.
    """
	logger.debug('Finding connected components ...')

	# The variant.allele attribute can be either 0 (major allele), 1 (minor allele),
	# or 3 (equal scores). If all_heterozygous is on (default), we can get
	# the combinations 0/1, 1/0 and 3/3 (the latter means: unphased).
	# If all_heterozygous is off, we can also get all other combinations.
	# In both cases, we are interested only in 0/1 and 1/0.
	phased_positions = [v1.position for v1, v2 in zip(*superreads)
						if (v1.allele, v2.allele) in ((0, 1), (1, 0))
	]
	assert phased_positions == sorted(phased_positions)

	# Find connected components.
	# A component is identified by the position of its leftmost variant.
	component_finder = ComponentFinder(phased_positions)
	phased_positions = set(phased_positions)
	for read in reads:
		positions = [variant.position for variant in read if variant.position in phased_positions]
		for position in positions[1:]:
			component_finder.merge(positions[0], position)
	components = {position: component_finder.find(position) for position in phased_positions}
	logger.info('No. of variants considered for phasing: %d', len(superreads[0]))
	logger.info('No. of variants that were phased: %d', len(phased_positions))
	return components


def best_case_blocks(reads):
	"""
    Given a list of core reads, determine the number of phased blocks that
    would result if each variant were actually phased.

    Return the number of connected components and non-singleton components.
    """
	positions = set()
	for read in reads:
		for variant in read:
			positions.add(variant.position)
	component_finder = ComponentFinder(positions)
	for read in reads:
		read_positions = [variant.position for variant in read]
		for position in read_positions[1:]:
			component_finder.merge(read_positions[0], position)
	# A dict that maps each component to the number of SNPs it contains
	component_sizes = defaultdict(int)
	for position in positions:
		component_sizes[component_finder.find(position)] += 1
	non_singletons = [component for component, size in component_sizes.items() if size > 1]
	return len(component_sizes), len(non_singletons)


class StageTimer:
	"""Measure run times of different stages of the program"""

	def __init__(self):
		self._start = dict()
		self._elapsed = defaultdict(float)

	def start(self, stage):
		"""Start measuring elapsed time for a stage"""
		self._start[stage] = time.time()

	def stop(self, stage):
		"""Stop measuring elapsed time for a stage."""
		t = time.time() - self._start[stage]
		self._elapsed[stage] += t
		return t

	def elapsed(self, stage):
		"""
        Return total time spent in a stage, which is the sum of the time spans
        between calls to start() and stop(). If the timer is currently running,
        its current invocation is not counted.
        """
		return self._elapsed[stage]

	def total(self):
		"""Return sum of all times"""
		return sum(self._elapsed.values())

	@contextmanager
	def __call__(self, stage):
		self.start(stage)
		yield
		self.stop(stage)


def ensure_pysam_version():
	from pysam import __version__ as pysam_version
	from distutils.version import LooseVersion

	if LooseVersion(pysam_version) < LooseVersion("0.8.1"):
		sys.exit("WhatsHap requires pysam >= 0.8.1")


def union_sets(index, connection_list, connectivity):
	'''Recursive method which builds the union of the sets in the connection_list
    with index lower than the given index, and which fullfill the given connectivity.
    Result is the connection_list with disjoint sets'''
	#print('Union_sets')
	#print(connection_list)
	# variable to go over the whole list recursively
	i = index - 1
	while (i != -1):
		actual_set = connection_list[index]
		former_set = connection_list[i]

		if len(actual_set.intersection(former_set)) >= connectivity:
			new_union = actual_set.union(former_set)
			connection_list[i] = new_union
			connection_list.pop(index)
			#recursivley check for the sets before the union...
			connection_list = union_sets(i, connection_list, connectivity)
			index -= 1
		i -= 1
	return connection_list


def check_for_connectivity(read_positions, List_of_connections, connectivity):
	connection_found = False
	#print('List_of_connections')
	#print(List_of_connections)
	#
	if len(List_of_connections) > 0:
		#stores the indices with which the read has connected and also the values of the intersection
		connections = []
		for i_c, c in enumerate(List_of_connections):
			#look if connection criteria is fullfilled
			if len(c.intersection(read_positions)) >= connectivity:
				#intersect_value=c.intersection(read_positions)
				union_of_those_sets = c.union(read_positions)
				List_of_connections[i_c] = union_of_those_sets
				#store the indices where in the original list of connection a change has
				# taken place so that later the union set could be called
				connections.append((i_c, union_of_those_sets))
				connection_found = True

		#Could only occure in this setting
		if connection_found:
			#print('DECISION CONNECTION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
			for (ic, union) in connections:
				actual_index_of_element = List_of_connections.index(union)
				union_sets(actual_index_of_element, List_of_connections, connectivity)
		else:
			#print('In FIRST ELSE No Connection')
			#In Order to remove double occurences
			if read_positions not in List_of_connections:
				List_of_connections.append(read_positions)
	else:
		#print('In SECOND ELSE EMpty Connection')
		List_of_connections.append(read_positions)

	return List_of_connections


def analyze_readset(sliced_reads, list_of_bam, connectivity, score):
	print('Working directory')
	print(os.getcwd())
	f = open('Analyzefile.unique_extension', 'w')

	# How many positions are covered by the readset
	f.write('Positions:')
	f.write("\n")
	for pos in sliced_reads.get_positions():
		f.write(str(pos))
		f.write("\t")
	component_finder = ComponentFinder(sliced_reads.get_positions())
	important_positions = set(sliced_reads.get_positions())
	f.write("\n")
	f.write('Length of Readset')
	f.write("\t")
	i = 0
	length_readset = len(sliced_reads)
	f.write(str(length_readset))
	f.write("\n")

	List_of_connections = []
	connected_blocks = {}
	while i != len(sliced_reads):
		positions = [variant.position for variant in sliced_reads[i] if variant.position in important_positions]
		# here the merging of only reads which share one component, could not be updated by the connectivity factor
		for position in positions[1:]:
			component_finder.merge(positions[0], position)
		components = {position: component_finder.find(position) for position in important_positions}

		#search now for components when the connectivity is not 1 :
		read_positions = set(positions)
		#if there is something in the connected blocks, so not the first entry

		List_of_connections = check_for_connectivity(read_positions, List_of_connections, connectivity)
		#print('LIST OF CONNECITONS')
		#print(List_of_connections)
		#
		# #still a boolean if a connection in this read is found
		# decision_connection = False
		# #store the blocks
		# dict_for_intersection_values = {}
		# if len(List_of_connections) > 0:
		# 	#stores the indices with which the read has connected and also the values of the intersection
		# 	connections = []
		# 	for i_c, c in enumerate(List_of_connections):
		# 		print('positions and c')
		# 		print(read_positions)
		# 		print(c)
		# 		#look if connection criteria is fullfilled
		# 		if len(c.intersection(read_positions)) >= connectivity:
		# 			print('IN SECOND IF')
		# 			intersect_value = c.intersection(read_positions)
        #
        #
		# 			#	intersect_val=str(intersect_value)
		# 			#	if intersect_val in dict_for_intersection_values:
		# 			#		dict_for_intersection_values[intersect_val] +=1
		# 			#	else:
		# 			#		dict_for_intersection_values[intersect_val] =1
        #
		# 			connections.append((i_c, intersect_value))
        #
        #
		# 			#Merge the read with the set in the  List
		# 			#print('c.union(read_positions)')
		# 			#print(c.union(read_positions))
		# 			#List_of_connections.pop(i_c)
		# 			#List_of_connections.append(c.union(read_positions))
        #
        #
        #
		# 			decision_connection = True
        #
		# 	#Could only occure in this setting
		# 	if decision_connection:
		# 		print('DECISION CONNECTION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
		# 		for (i_c, value) in connections:
		# 			print('List_of_connection at start')
		# 			print(List_of_connections)
		# 			first_set = List_of_connections[i_c]
		# 			print('FIRST SET')
		# 			print(first_set)
		# 			union_of_those_sets = first_set.union(read_positions)
		# 			print('Union_of sets')
		# 			print(union_of_those_sets)
		# 			List_of_connections[i_c] = union_of_those_sets
		# 			print('List of connection after union')
		# 			print(List_of_connections)
		# 			if str(value) in dict_for_intersection_values:
		# 				print('Founr Some Intersection with more than 1 ')
		# 				dict_for_intersection_values[str(value)] += 1
		# 				print('VALUE for storing in dict')
		# 				print(value)
		# 			#former_list=to_check_for_merge[str(value)]
		# 			#new_list=former_list.add(i_c)
		# 			#to_check_for_merge[str(value)]=new_list
        #
		# 			else:
		# 				print('Foun no intersection')
		# 				dict_for_intersection_values[str(value)] = 1
        #
        #
		# 			#Need to look if the value which are used for combination are used more than one time
		# 			#TODO Other dict_for_intersection_values[value ].... The new is wrotm
		# 			#				for value_for_intersect in dict_for_intersection_values:
		# 			#					if (value_for_intersect > 1):
		# 			#						print('FOUND SOME INTERSECTION WITH MORE THAN 2')
		# 			#					else:
		# 			#						print('FouND NO INTERSECTION ')
        #
        #
        #
        #
        #
		# 	else:
		# 		print('In SECOND ELSE')
		# 		List_of_connections.append(read_positions)
		# 	print('Dictionary')
		# 	print(dict_for_intersection_values)
        #
        #
        #
		# #For the first read in the readset
		# else:
		# 	#print('IN FIRST ELSE')
		# 	List_of_connections.append(read_positions)
        #
        #
        #
        #
        #
        #
        #
        #
        #
        #
        #




		#TODO : HERE something went wrong dont get the same result for connectivity 1 as with the component_finder
        #
        #
		# decision_connection = False
		# #only the first one has to be definitly in the connected_ block,,,, therefore
		# if len(connected_blocks) > 0:
		# 	print('In FIRST IF')
		# 	#boolean if it is connected to some former entry
        #
		# 	#stores the indices of the block, which have later to be merges together....
		# 	connection_indices = []
		# 	#look at every entry of the block an look if it intersects with the newly read positions
		# 	for j in connected_blocks:
        #
		# 		if len(connected_blocks[j].intersection(read_positions)) >= connectivity:
		# 			print('In SECOND If')
		# 			#found some overlap which is larger or equal the connectivity
		# 			decision_connection = True
		# 			#add index of the indey which have to be merged with the read to the list.
		# 			connection_indices.append(j)
		# else:
		# 	print('In FIRST ELSE')
		# 	#append read_position on the blocks if there is no connection possible with the former members of the block
		# 	connected_blocks[len(connected_blocks)] = read_positions
        #
        #
		# if decision_connection:
		# 	print('IN IF WITH BOOLEAN TRUE ')
		# 	#ifthere indexes are stored which have to be merged:
		# 	to_merged_set=read_positions
		# 	for index in connection_indices:
		# 		#merge all positions in the connected blocks together and add it later again at the end
		# 		to_merged_set=to_merged_set.union(connected_blocks[index])
		# 		#TODO VERMUTE HIER GEHT WAS SCHIEF
		# 		del connected_blocks[index]
		# 	connected_blocks[len(connected_blocks)]=to_merged_set
		# else:
		# 	connected_blocks[len(connected_blocks)]=read_positions

		f.write("Read")
		f.write(str(i))
		f.write("\t")
		for variant in sliced_reads[i]:
			pos = variant.position
			qual = variant.quality
			f.write(str(pos))
			f.write("\t")
			f.write(str(qual))
			f.write("\t")

		f.write("\n")
		i = i + 1

	f.write('Connectivity of Blocks')
	f.write("\n")
	f.write(str(connectivity))
	f.write("\n")

	# k=0
	# while k != len(sliced_reads):
	# #at the moment the same informations as in line 427 above
	# 	#positions = [variant.position for variant in sliced_reads[k] if variant.position in important_positions]
	# 	read_positions = set(positions)
	# 	#if there is something in the connected blocks, so not the first entry
	# 	#print('Connected Block length')
	# 	#print(len(connected_blocks))
	# 	decision_connection = False
	# 	if len(connected_blocks) > 1:
	# 		#print('In IF STATEMENT FOR CONNECTED BLOCKS')
	# 		#boolean if it is connected to some former entry
    #
	# 		#stores the indices of the block, which have later to be merges together....
	# 		connection_indices = []
	# 		#look at every entry of the block an look if it intersects with the newly read positions
	# 		for j in connected_blocks:
    #
	# 			if len(connected_blocks[j].intersection(read_positions)) >= connectivity:
	# 				#found some overlap which is larger or equal the connectivity
	# 				decision_connection = True
	# 				#add index of the indey which have to be merged with the read to the list.
	# 				connection_indices.append(j)
	# 	else:
	# 		#append read_position on the blocks if there is no connection possible with the former members of the block
	# 		connected_blocks[len(connected_blocks)] = read_positions
    #
    #
	# 	if decision_connection:
	# 		#ifthere indexes are stored which have to be merged:
	# 		to_merged_set=read_positions
	# 		for index in connection_indices:
	# 			#merge all positions in the connected blocks together and add it later again at the end
	# 			to_merged_set=to_merged_set.union(connected_blocks[index])
	# 			del connected_blocks[index]
	# 		connected_blocks[len(connected_blocks)]=to_merged_set
    #
	# 	k=k+1
	f.write('Connected Blocks by given Connectivity')
	f.write("\n")
	f.write(str(List_of_connections))
	#f.write(str(connected_blocks))
	f.write("\n")

	f.write('Blocks (always with connectivity 1)')
	f.write("\n")
	comset = set(components.values())
	f.write(str(comset))
	f.write("\n")
	f.write('Blocks and their components')
	f.write("\n")
	f.write(str(components))
	f.write("\n")
	f.write("Number of reads in bam files")
	f.write("\n")
	f.write(str(list_of_bam))
	f.write("\n")
	f.write('Which score:')
	f.write("\n")
	f.write(str(score))
	f.write("\n")
	f.close()

	return 0


def run_whatshap(bam, vcf,
				 output=sys.stdout, sample=None, ignore_read_groups=False, indels=True,
				 mapping_quality=20, max_coverage=15, all_heterozygous=True, seed=123, haplotype_bams_prefix=None,
				 connectivity=1, bridge=False, analyze=False, score=0):
	"""
    Run WhatsHap.

    bam -- list of paths to BAM files
    vcf -- path to input VCF
    output -- path to output VCF or sys.stdout
    sample -- name of sample to phase. None means: phase first found sample.
    ignore_read_groups
    mapping_quality -- discard reads below this mapping quality
    max_coverage
    all_heterozygous
    seed -- seed for random numbers
    haplotype_bams_prefix -- name selected reads are stored in a bam file
    connectivity -- number how many overlapping variants define a connection
    bridge -- using bridging
    analyze -- writing additional file with informations
    score -- which scoring method is used

    """
	random.seed(seed)

	class Statistics:
		pass

	stats = Statistics()
	timers = StageTimer()
	timers.start('overall')
	stats.n_homozygous = 0
	stats.n_phased_blocks = 0
	stats.n_best_case_blocks = 0
	stats.n_best_case_nonsingleton_blocks = 0
	stats.n_best_case_blocks_cov = 0
	stats.n_best_case_nonsingleton_blocks_cov = 0
	logger.info("This is WhatsHap %s running under Python %s", __version__, platform.python_version())
	with ExitStack() as stack:
		try:
			bam_reader = stack.enter_context(ReadSetReader(bam, mapq_threshold=mapping_quality))
		except (OSError, BamIndexingError) as e:
			logger.error(e)
			sys.exit(1)
		if output is not sys.stdout:
			output = stack.enter_context(open(output, 'w'))
		command_line = ' '.join(sys.argv[1:])
		vcf_writer = PhasedVcfWriter(command_line=command_line, in_path=vcf, out_file=output)
		vcf_reader = parse_vcf(vcf, sample=sample, indels=indels)
		haplotype_bam_writer = None
		if haplotype_bams_prefix is not None:
			haplotype_bam_writer = HaplotypeBamWriter(bam, haplotype_bams_prefix, sample)
		timers.start('parse_vcf')
		for sample, chromosome, variants in vcf_reader:
			variants = remove_overlapping_variants(variants)
			timers.stop('parse_vcf')
			logger.info('Working on chromosome %s', chromosome)
			logger.info('Read %d variants', len(variants))
			bam_sample = None if ignore_read_groups else sample
			logger.info('Reading the BAM file ...')
			try:
				timers.start('read_bam')
				reads = bam_reader.read(chromosome, variants, bam_sample)
				bam_time = timers.stop('read_bam')
			except SampleNotFoundError:
				logger.error("Sample %r is not among the read groups (RG tags) "
							 "in the BAM header.", bam_sample)
				sys.exit(1)
			logger.info('Read %d reads in %.1f s', len(reads), bam_time)

			with timers('slice'):
				# Sort the variants stored in each read
				# TODO: Check whether this is already ensured by construction
				for read in reads:
					read.sort()
				# Sort reads in read set by position
				reads.sort()
				# defined as 0 at the moment
				score_selection = score
				#TODO: Now time for finding listofbam included into slicing time
				if analyze:
					selected_reads, uninformative_read_count, listofbam = readselection(reads, max_coverage, bridge,
																						analyze, score_selection)
				else:
					selected_reads, uninformative_read_count = readselection(reads, max_coverage, bridge, analyze,
																			 score_selection)
				sliced_reads = reads.subset(selected_reads)

				position_list = reads.get_positions()
				accessible_positions = sliced_reads.get_positions()
				informative_read_count = len(reads) - uninformative_read_count
				unphasable_snps = len(position_list) - len(accessible_positions)
				logger.info('%d variants are covered by at least one read', len(position_list))
				logger.info('Skipped %d reads that only cover one variant', uninformative_read_count)
				if position_list:
					logger.info('%d out of %d variant positions (%.1d%%) do not have a read '
								'connecting them to another variant and are thus unphasable',
								unphasable_snps, len(position_list),
								100. * unphasable_snps / len(position_list)
					)
				if reads:
					logger.info('After read selection: Using %d of %d (%.1f%%) reads that cover two or more variants',
								len(selected_reads), informative_read_count, (
							100. * len(selected_reads) / informative_read_count if informative_read_count > 0 else float(
								'nan'))
					)

			n_best_case_blocks, n_best_case_nonsingleton_blocks = best_case_blocks(reads)
			n_best_case_blocks_cov, n_best_case_nonsingleton_blocks_cov = best_case_blocks(sliced_reads)
			stats.n_best_case_blocks += n_best_case_blocks
			stats.n_best_case_nonsingleton_blocks += n_best_case_nonsingleton_blocks
			stats.n_best_case_blocks_cov += n_best_case_blocks_cov
			stats.n_best_case_nonsingleton_blocks_cov += n_best_case_nonsingleton_blocks_cov
			logger.info('Best-case phasing would result in %d non-singleton phased blocks (%d in total)',
						n_best_case_nonsingleton_blocks, n_best_case_blocks)
			logger.info('... after read selection: %d non-singleton phased blocks (%d in total)',
						n_best_case_nonsingleton_blocks_cov, n_best_case_blocks_cov)
			# in order to not include the time for writing of the file for analysis here the analyze_readset is called
			with timers('analyzing'):
				logger.info('Writing different options in seperate file for later analysis')
				if analyze:
					analyze_readset(sliced_reads, listofbam, connectivity, score)
			with timers('phase'):
				logger.info('Phasing the variants (using %d reads)...', len(sliced_reads))
				# Run the core algorithm: construct DP table ...
				dp_table = DPTable(sliced_reads, all_heterozygous)
				# get the mec score
				mec_score = dp_table.get_optimal_cost()
				# ... and do the backtrace to get the solution
				superreads = dp_table.get_super_reads()

				# output the MEC score of phasing
				logger.info('MEC score of phasing: %d', mec_score)

				n_homozygous = sum(1 for v1, v2 in zip(*superreads)
								   if v1.allele == v2.allele and v1.allele in (0, 1))
				stats.n_homozygous += n_homozygous

			components = find_components(superreads, sliced_reads)
			n_phased_blocks = len(set(components.values()))
			stats.n_phased_blocks += n_phased_blocks
			logger.info('No. of phased blocks: %d', n_phased_blocks)
			if all_heterozygous:
				assert n_homozygous == 0
			else:
				logger.info('No. of heterozygous variants determined to be homozygous: %d', n_homozygous)

			with timers('write_vcf'):
				vcf_writer.write(chromosome, sample, superreads, components)

			if haplotype_bam_writer is not None:
				logger.info('Writing used reads to haplotype-specific BAM files')
				haplotype_bam_writer.write(sliced_reads, dp_table.get_optimal_partitioning(), chromosome)
			logger.info('Chromosome %s finished', chromosome)
			timers.start('parse_vcf')
		timers.stop('parse_vcf')
	logger.info('== SUMMARY ==')
	logger.info('Best-case phasing would result in %d non-singleton phased blocks (%d in total)',
				stats.n_best_case_nonsingleton_blocks, stats.n_best_case_blocks)
	logger.info('... after read selection: %d non-singleton phased blocks (%d in total)',
				stats.n_best_case_nonsingleton_blocks_cov, stats.n_best_case_blocks_cov)
	if all_heterozygous:
		assert stats.n_homozygous == 0
	else:
		logger.info('No. of heterozygous variants determined to be homozygous: %d', stats.n_homozygous)
	timers.stop('overall')
	logger.info('Time spent reading BAM: %6.1f s', timers.elapsed('read_bam'))
	logger.info('Time spent parsing VCF: %6.1f s', timers.elapsed('parse_vcf'))
	logger.info('Time spent slicing:     %6.1f s', timers.elapsed('slice'))
	logger.info('Time spent for analysis:%6.1f s', timers.elapsed('analyzing'))
	logger.info('Time spent phasing:     %6.1f s', timers.elapsed('phase'))
	logger.info('Time spent writing VCF: %6.1f s', timers.elapsed('write_vcf'))
	logger.info('Time spent on rest:     %6.1f s', 2 * timers.elapsed('overall') - timers.total())
	logger.info('Total elapsed time:     %6.1f s', timers.elapsed('overall'))


class NiceFormatter(logging.Formatter):
	"""
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).

    Based on http://stackoverflow.com/a/9218261/715090 .
    """

	def format(self, record):
		if record.levelno != logging.INFO:
			record.msg = '{}: {}'.format(record.levelname, record.msg)
		return super().format(record)


def main():
	ensure_pysam_version()
	logger.setLevel(logging.INFO)
	handler = logging.StreamHandler()
	handler.setFormatter(NiceFormatter())
	logger.addHandler(handler)
	parser = ArgumentParser(prog='whatshap', description=__doc__)
	parser.add_argument('--version', action='version', version=__version__)
	parser.add_argument('--debug', action='store_true', default=False,
						help='Show more verbose output')
	parser.add_argument('-o', '--output', default=sys.stdout,
						help='Output VCF file. If omitted, use standard output.')
	parser.add_argument('--max-coverage', '-H', metavar='MAXCOV', default=15, type=int,
						help='Reduce coverage to at most MAXCOV (default: %(default)s).')
	# maybe it should be possible to change the quality which we selectt for for each bam file individually
	parser.add_argument('--mapping-quality', '--mapq', metavar='QUAL',
						default=20, type=int, help='Minimum mapping quality (default: %(default)s)')
	parser.add_argument('--seed', default=123, type=int, help='Random seed (default: %(default)s)')
	parser.add_argument('--indels', dest='indels', default=False, action='store_true',
						help='Also phase indels (default: do not phase indels)')
	parser.add_argument('--distrust-genotypes', dest='all_heterozygous',
						action='store_false', default=True,
						help='Allow switching variants from hetero- to homozygous in an '
							 'optimal solution (see documentation).')
	parser.add_argument('--ignore-read-groups', default=False, action='store_true',
						help='Ignore read groups in BAM header and assume all reads come '
							 'from the same sample.')
	parser.add_argument('--sample', metavar='SAMPLE', default=None,
						help='Name of a sample to phase. If not given, the first sample in the '
							 'input VCF is phased.')
	parser.add_argument('--haplotype-bams', metavar='PREFIX', dest='haplotype_bams_prefix', default=None,
						help='Write reads that have been used for phasing to haplotype-specific BAM files. '
							 'Creates PREFIX.1.bam and PREFIX.2.bam')
	parser.add_argument('--connectivity', metavar='CONECT', default=1, type=int,
						help='Sets value of connectivity between reads in the'
							 '  unity of overlapping variants (default: %(default)s)')
	parser.add_argument('--bridging', dest='bridge', default=False, action='store_true',
						help='Selecting if in Readselection'
							 ' the bridging is on (default: %(default)s)')
	parser.add_argument('--analyzfile', dest='analyze', default=False, action='store_true',
						help='Write properties of the selected readset,like'
							 'conneted blocks, quality, origin ,... to separate file named Analyzefile.unique_extension '
							 '(default: %(default)s )')
	parser.add_argument('--score', '--score', metavar='SCORE', default=0, type=int,
						help='Select the score you would like to '
							 'use (defualt :%(default)s), score have to be set '
							 '0 = actual best score'
							 '1 = best score for paired_end reads'
							 '2= best score for single ended reads')
	parser.add_argument('vcf', metavar='VCF', help='VCF file')
	parser.add_argument('bam', nargs='+', metavar='BAM', help='BAM file')

	args = parser.parse_args()
	logger.setLevel(logging.DEBUG if args.debug else logging.INFO)
	del args.debug
	run_whatshap(**vars(args))
