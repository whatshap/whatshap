import math
import logging
from collections import defaultdict

from whatshap.coverage import CovMonitor
from .graph import ComponentFinder

logger = logging.getLogger(__name__)


def _construct_indexes(readset):
	''' The parameter readset: is the given ReadSet and returns, all possible variant positions, the vcf_index_ mapping
    and the snp_to_reads_map'''
	positions = readset.get_positions()
	vcf_indices = {position: index for index, position in enumerate(positions)}
	snp_to_reads_map = defaultdict(list)
	for index, read in enumerate(readset):
		for variant in read:
			snp_index = vcf_indices[variant.position]
			snp_to_reads_map[snp_index].append(index)
	return positions, vcf_indices, snp_to_reads_map


cdef priority_type_ptr _update_score_for_reads(priority_type_ptr former_score, ReadSet* readset, index, unordered_set[int]& already_covered_snps):
	'''updatest the score of the read, depending on how many reads are already covered'''
	cdef int first_score = former_score.at(0)
	cdef int second_score = former_score.at(1)
	cdef int quality = former_score.at(2)
	cdef Read* read = readset.get(index)
	cdef unordered_set[int].iterator it
	for i in range(read.getVariantCount()):
		it = already_covered_snps.find(read.getPosition(i))
		if it == already_covered_snps.end():
			first_score -= 1
	cdef priority_type_ptr result = new priority_type()
	result.push_back(first_score)
	result.push_back(second_score)
	result.push_back(quality)
	return result

#TODO not working method here
def _compute_score_depending_on_quality_only (readset,index,vcf_indices):
	'''Method which computes anothe readscore depending only in the quality,mor precisely the
	score of the read is the average quality in the read over all SNP positions covered by the read'''
	read= readset[index]
	#min_quality= 1000
	actual_quality= 0
	for pos in enumerate:
		quality= pos.quality
		actual_quality += quality
	return (actual_quality/ len(read))


cdef priority_type_ptr _compute_score_for_read(ReadSet* readset, int index, vcf_indices):
	'''Compute the score for one read, assuming no other reads have been selected so far
	(after selecting reads, scores can be updated using _update_score_for_reads).
	We use the following scoring scheme: (new - gaps, total - bad, min(qualities)),
	where "new" is the number of SNPs covered by this read and no other selected read (so far),
	"gaps" is the number of SNPs overlapped by (physical) fragment, but not covered by the sequenced part of the read,
	"total" is the total number of SNPs covered by the read, and "min(qualities)" is the minimum
	over all base qualities at SNP positions.
	'''
	cdef Read* read = readset.get(index)
	cdef int min_quality = -1
	cdef int good_score = 0
	cdef int bad_score = 0
	cdef int quality = -1
	cdef int pos = -1
	covered_SNPS = []
	for i in range(read.getVariantCount()):
		quality = read.getVariantQuality(i)
		pos = read.getPosition(i)
		if i == 0:
			min_quality = quality
		else:
			min_quality = min(min_quality, quality)
		SNP_covered = vcf_indices.get(pos)
		if SNP_covered != None:
			covered_SNPS.append(SNP_covered)
			#good_score represents how many SNPs are really covered
			good_score += 1
	#bad_score here means number of SNPs which are not covered...
	if len(covered_SNPS) != (covered_SNPS[len(covered_SNPS) - 1] - covered_SNPS[0] + 1):
		bad_score = (covered_SNPS[len(covered_SNPS) - 1] - covered_SNPS[0] + 1) - len(covered_SNPS)
	cdef priority_type_ptr result = new priority_type()
	#initially new_score is the same as good_score, but new_score is updated later
	result.push_back(good_score - bad_score)
	result.push_back(good_score - bad_score)
	result.push_back(min_quality)
	return result

#TODO work also here
cdef priority_type_ptr _compute_score_including_map_qual(ReadSet* readset, int index, vcf_indices):
	'''Alternative score, maybe for longer reads, just for slicing experiments,and for checking if selection in
	the argument parser works.
	Using only a single score of the number of covered SNPs, but inserted into tupel
	for the moment because the vector structure
	only allows that, has to be redone.
	'''
	cdef Read* read = readset.get(index)
	#Getting the map_quality of a read in form of a vector
	mapquality=read.getMapqs()
	print("MApQuality")
	print(mapquality)
	if (mapquality.front()==60):
		print('FOUND')
	else:
		print('NOT FOUND')
	cdef int min_quality = -1
	cdef int good_score = 0
	cdef int bad_score = 0
	cdef int quality = -1
	cdef int pos = -1
	cdef int score
	covered_SNPS = []
	for i in range(read.getVariantCount()):
		quality = read.getVariantQuality(i)
		pos = read.getPosition(i)
		print('Quality_of_variant')
		print(pos)
		print(quality)
		#if i == 0:
		#	min_quality = quality
		#else:
		#	min_quality = min(min_quality, quality)
		SNP_covered = vcf_indices.get(pos)
		if SNP_covered != None:
			covered_SNPS.append(SNP_covered)
			#_score represents how many SNPs are really covered
			score += 1
	#bad_score here means number of SNPs which are not covered...
	if len(covered_SNPS) != (covered_SNPS[len(covered_SNPS) - 1] - covered_SNPS[0] + 1):
		bad_score = (covered_SNPS[len(covered_SNPS) - 1] - covered_SNPS[0] + 1) - len(covered_SNPS)
	cdef priority_type_ptr result = new priority_type()
	#initially new_score is the same as good_score, but new_score is updated later
	#result.push_back(good_score - bad_score)
	#result.push_back(good_score - bad_score)
	result.push_back(score)
	result.push_back(score)
	result.push_back(score)
	return result

# TODO New Idea for combi situation... Add a factor for paired end, to add them first, if the number of covered positions and
 #the the quality is equal to single end reads









cdef PriorityQueue _construct_priorityqueue(ReadSet* readset, read_indices, vcf_indices,score_selection):
	'''Construct a priority queue containing all given read indicies, each one representing the
	respective read from readset.'''
	cdef PriorityQueue priorityqueue = PriorityQueue()

	for index in read_indices:
		read = readset.get(index)
		if score_selection==0:
			computed_score = _compute_score_for_read(readset, index, vcf_indices)
		else:
			computed_score = _compute_score_including_map_qual(readset, index, vcf_indices)
		#computed_score = _compute_score_for_longer_reads(readset, index, vcf_indices)
		priorityqueue.c_push(computed_score, index)

	return priorityqueue


cdef _slice_read_selection(PriorityQueue pq, coverages, max_cov, ReadSet* readset, vcf_indices, snp_to_reads_map):
	'''Extraction of a set of read indices, where each SNP should be covered at least once, if coverage, or reads are allowing it '''
	# positions of SNPs covered by any read selected so far
	already_covered_snps = set()
	# indices of selected reads
	reads_in_slice = set()
	# indices of reads that cannot be added because doing that would violate coverage constraint
	reads_violating_coverage = set()
	cdef Read* extracted_read = NULL
	cdef queue_entry_type entry
	cdef priority_type_ptr max_score
	cdef int max_item
	cdef int pos
	cdef unordered_set[int] snps_covered_by_this_read
	cdef priority_type_ptr oldscore
	cdef priority_type_ptr newscore
	while not pq.c_is_empty():
		snps_covered_by_this_read.clear()
		entry = pq.c_pop()
		max_score = entry.first
		max_item = entry.second
		extracted_read = readset.get(max_item)
		covers_new_snp = False
		#look if positions covered by this reads are already covered or not
		for i in range(extracted_read.getVariantCount()):
			pos = extracted_read.getPosition(i)
			if pos in already_covered_snps:
				continue
			else:
				covers_new_snp = True
				#stores the positions the read covers
				snps_covered_by_this_read.insert(pos)
		# only add read if it covers at least one new SNP and adding it does not violate coverage constraints
		begin = vcf_indices.get(extracted_read.getPosition(0))
		end = vcf_indices.get(extracted_read.getPosition(extracted_read.getVariantCount() - 1)) + 1
		if coverages.max_coverage_in_range(begin, end) >= max_cov:
			reads_violating_coverage.add(max_item)
		elif covers_new_snp:
			coverages.add_read(begin, end)
			reads_in_slice.add(max_item)
			reads_whose_score_has_to_be_updated = set()

			#again go over the positions in the read and add them to the already_covered_SNP list

			#Only the positions in the read which cover new SNPs are analysed
			for pos in snps_covered_by_this_read:
				already_covered_snps.add(pos)
				reads_whose_score_has_to_be_updated.update(snp_to_reads_map[vcf_indices.get(pos)])

			#find difference between to_decrease_score and selected_reads in order to not to try to decrease score by selected reads
			selected_read_set = set(reads_in_slice)
			decrease_set = reads_whose_score_has_to_be_updated
			d_set = decrease_set.difference(selected_read_set)

			#Catch with None if element is not in the priorityqueue
			for element in d_set:
				oldscore = pq.c_get_score_by_item(element)
				if oldscore != NULL:
					newscore = _update_score_for_reads(oldscore, readset, element, snps_covered_by_this_read)
					pq.c_change_score(element, newscore)
	return reads_in_slice, reads_violating_coverage



def format_read_source_stats(readset, indices):
	"""Creates a string giving information on the source_ids of the reads with the given indices."""
	if len(indices) == 0:
		return 'n/a'
	source_id_counts = defaultdict(int)
	for i in indices:
		source_id_counts[readset[i].source_id] += 1
	present_ids = list(source_id_counts.keys())
	present_ids.sort()
	return ', '.join('{}:{}'.format(source_id, count) for source_id, count in source_id_counts.items())


def readselection(PyReadSet pyreadset, max_cov, bridging,analyze,score_selection):
	'''Return the selected readindices which do not violate the maximal coverage, and additionally usage of a boolean for deciding if
     the bridging is needed or not.'''
	cdef ReadSet* readset = pyreadset.thisptr
	assert readset != NULL

	positions, vcf_indices, snp_to_reads_map = _construct_indexes(pyreadset)

	logger.info('Running read selection for %d reads covering %d variants (bridging %s)', len(pyreadset), len(positions),'ON' if bridging else 'OFF')

	#initialization of Coverage Monitor
	coverages = CovMonitor(len(positions))

	# indices of reads that have been selected
	selected_reads = set()

	# indices of reads that could (potentially) still be selected ,do not consider thes read which only cover one SNP
	undecided_reads = set(i for i, read in enumerate(pyreadset) if len(read) >= 2)

	uninformative_read_count = len(pyreadset) - len(undecided_reads)

	cdef PriorityQueue pq
	cdef Read* read
	loop = 0
	#setting up list for storing how many reads came from which bam file.
	newlist=[]
	while len(undecided_reads) > 0:
		#switch cases : because we maybe need other update rules for other scores.
		#if (score_selection == 0):
		#	pq = _construct_priorityqueue(readset, undecided_reads, vcf_indices,score_selection)
		#else:
		#	pq = _construct_priorityqueue(readset, undecided_reads, vcf_indices,score_selection)
		#
		#		print('FORMER SCORE USED')
		pq = _construct_priorityqueue(readset, undecided_reads, vcf_indices,score_selection)
		reads_in_slice, reads_violating_coverage = _slice_read_selection(pq, coverages, max_cov, readset, vcf_indices, snp_to_reads_map)
		selected_reads.update(reads_in_slice)
		undecided_reads -= reads_in_slice
		undecided_reads -= reads_violating_coverage

		# Create new component finder from reads just selected
		component_finder = ComponentFinder(positions)
		for read_index in reads_in_slice:
			read = readset.get(read_index)
			#only need to build up this list for th e reads in different bam files if we need to write it out
			if analyze:

				readID= read.getSourceID()

				while len(newlist)< (readID+1):
					newlist.append(0)



				score=newlist[readID]
				newscore=score+1
				newlist[readID]=newscore


			for i in range(1, read.getVariantCount()):
				component_finder.merge(read.getPosition(0), read.getPosition(i))
		bridging_reads = set()
		if bridging:
			pq = _construct_priorityqueue(readset, undecided_reads, vcf_indices,score_selection)
			while not pq.is_empty():
				score, read_index = pq.pop()
				read = readset.get(read_index)
				covered_blocks = set()
				for i in range(read.getVariantCount()):
					covered_blocks.add(component_finder.find(read.getPosition(i)))

				# check whether read meets coverage constraints
				begin = vcf_indices.get(read.getPosition(0))
				end = vcf_indices.get(read.getPosition(read.getVariantCount() - 1)) + 1
				if coverages.max_coverage_in_range(begin, end) >= max_cov:
					undecided_reads.remove(read_index)
					continue

				# skip read if it only covers one block
				if len(covered_blocks) < 2:
					continue
				bridging_reads.add(read_index)
				selected_reads.add(read_index)

				coverages.add_read(begin, end)

				undecided_reads.remove(read_index)
				# update component_finder
				for i in range(1, read.getVariantCount()):
					component_finder.merge(read.getPosition(0), read.getPosition(i))
		loop += 1
		logger.info(
			'... iteration %d: selected %d reads (source: %s) to cover positions and %d reads (source: %s) for bridging; %d reads left undecided',
			loop, len(reads_in_slice), format_read_source_stats(pyreadset, reads_in_slice), len(bridging_reads),
			format_read_source_stats(pyreadset, bridging_reads), len(undecided_reads)
		)
	if analyze:
		return selected_reads, uninformative_read_count,newlist
	else:
		return selected_reads, uninformative_read_count
