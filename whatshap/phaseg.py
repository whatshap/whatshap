"""
Phase variants in a Locus file with the WhatsHap algorithm

Read a Locus and reads from GAM file. The optimal partitioning is written to standard output.
"""
import pyfaidx
from xopen import xopen
import stream
import logging
from . import vg_pb2
from collections import Counter
from collections import defaultdict
from .core import ReadSet, Read
from functools import reduce
import operator as op


from contextlib import ExitStack
from .vcf import VcfReader, PhasedVcfWriter
from . import __version__
from .core import ReadSet, readselection, Pedigree, PedigreeDPTable, NumericSampleIds, PhredGenotypeLikelihoods
from .graph import ComponentFinder
from .pedigree import (PedReader, mendelian_conflict, recombination_cost_map,
                       load_genetic_map, uniform_recombination_map, find_recombination)
from .bam import BamIndexingError, SampleNotFoundError
from .timer import StageTimer
from .variants import ReadSetReader, ReadSetError

__author__ = "Shilpa Garg, Tobias Marschall"

logger = logging.getLogger(__name__)

def ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, range(n, n-r, -1))
    denom = reduce(op.mul, range(1, r+1))
    return numer//denom

def vg_reader(locus_file, gam_file):
	"""
	input: reads locus and GAM file output from vg.
	output: sorted readset for core DP.
	assumptions: 
	1. locus file consists of linear ordering of simple bubbles only and hence sorted.
	2. paths in the locus should be covered by atleast one pacbio read.
	2. GAM file is sorted and restricted to locus file.
	3. files consists of all connected components.
	4. add variant only when it identifies the branch uniquely.
	"""
	gam_list = []
	with stream.open(str(gam_file), "rb") as istream:
		for data in istream:
			aln = vg_pb2.Alignment()
			aln.ParseFromString(data)
			gam_list.append(aln)


	locus_list = []
	with stream.open(str(locus_file), "rb") as istream:
		for data in istream:
			aln = vg_pb2.Locus()
			aln.ParseFromString(data)
			locus_list.append(aln)

	# create a dictionary of branches for each locus based on locus file.
	locus_branch_mapping=defaultdict()
	locus_count=0
	alleles_per_pos= defaultdict()
	for l in locus_list:
		per_locus_position=[]
		if len(l.allele) > 1:
			for paths in l.allele:
				tmp=set()
				if len(paths.mapping) > 1: # donot consider homozygous locus, only with multiple paths.
					for mappings in paths.mapping:
						tmp.add(mappings.position.node_id)
					per_locus_position.append(tmp)
			alleles_per_pos[locus_count] = len(per_locus_position)
			locus_branch_mapping[locus_count]=per_locus_position
			locus_count=locus_count+1

	if 3 in locus_branch_mapping: del locus_branch_mapping[3]
	print(locus_branch_mapping)
	# key is the values in locus_branch_mapping and value is triplet(locus, branch, alleles to be flipped)
	# if a simple bubbles consists of complicated path, then there would atleast one node that uniquely determine the branch.
	reverse_mapping= defaultdict(list)
	for k,v in locus_branch_mapping.items():
		if len(v) > 1: # more than one branch
			for i,b in enumerate(v):
				if len(b) > 2:
					for p,j in enumerate(b): # reverse_mapping for every node.
						if p==0 or p==len(b)-1: # first and last element of simple bubble
							reverse_mapping[j]=[[-1,-1,-1]] # handle start and sink node.
							continue
						reverse_mapping[j].append([k,i, len(v)])
	last_index = sorted(reverse_mapping.keys())[-1]
	reverse_mapping[last_index]=[[-1,-1,-1]]


	# extract reads from GAM file associated with the locus and create a sorted readset.
	readset=ReadSet()  
	#TODO: update for phred score (0,10,10....) -- equal to number of alleles.
	for g in gam_list:
		# hard-coded source id, mapping quality and other values.
		read=Read(g.name, 0, 0, 0) # create read for each read alignment
		prev_tmp=[]
		prev_locus= -1
		for i in g.path.mapping: # go over the mapping in a read
			node = i.position.node_id # go over nodes in a mapping
			if node in reverse_mapping and reverse_mapping[node][0][0]!=-1: # handle start and sink node.
				qualities = [10]* reverse_mapping[node][0][2]
				tmp=[]
				node_inf= [tuple(i[0:2]) for i in reverse_mapping[node]]
				for i in node_inf:
					tmp.append(i)
				interset_tmp= list(set(tmp).intersection(set(prev_tmp)))
				if len(tmp)==1 and len(prev_tmp)==0: # simple bubble with only two paths.
					qualities[tmp[0][1]] = 0
					read.add_variant(tmp[0][0], tmp[0][1], qualities) # if any new locus of branch is encountered, enter a variant in read.
					locus = tmp[0][0]
				elif len(prev_tmp) > 0 and len(set(tmp).intersection(set(prev_tmp)))==1: # for complicated bubbles, but with Top-N paths. combination of some nodes uniquely determine branch.
					qualities[interset_tmp[0][1]] = 0
					read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)
					locus= interset_tmp[0][0]
				elif len(prev_tmp) > 0 and len(tmp) ==1: # complicated bubbles, but one node can uniquely identify the branch.
					qualities[tmp[0][1]] = 0
					read.add_variant(tmp[0][0], tmp[0][1], qualities)
					locus = tmp[0][0]

				if prev_locus!=locus:
					prev_tmp = []
				else:
					for i in tmp:
						prev_tmp.append(i)
		readset.add(read)
	for read in readset:
			read.sort()
	readset.sort()
	return readset, alleles_per_pos

def run_phaseg(gam_file, locus_file):
	"""
	Run WhatsHap.

	gam_file -- path to GAM file
	locus_file -- path to input variants
	"""
	timers = StageTimer()
	recombrate=1.26
	all_heterozygous = False
	distrust_genotypes = True
	timers.start('overall')
	#logger.info("This is WhatsHap %s running under Python %s", __version__, platform.python_version())
	with ExitStack() as stack:
		logger.info('Using uniform recombination rate of %g cM/Mb.', recombrate)
		all_reads, alleles_per_pos = vg_reader(locus_file, gam_file)
		selected_reads = readselection(all_reads, 15)
		accessible_positions = sorted(selected_reads.get_positions())
		pedigree = Pedigree(NumericSampleIds())
		# compute the number of alleles at each position.
		alleles_per_accessible_pos =[]
		genotype_likelihoods = []
		for pos in accessible_positions:
			if pos in alleles_per_pos:
				n_alleles = alleles_per_pos[pos]  
				possible_genotypes = n_alleles +  ncr(n_alleles, 2)
				genotype_likelihoods.append(None if all_heterozygous else PhredGenotypeLikelihoods([0]* possible_genotypes))
		# random input of genotypes, since distrust_genotypes is always ON.
		pedigree.add_individual('individual0', [0]* len(accessible_positions), genotype_likelihoods)
		recombination_costs = uniform_recombination_map(recombrate, accessible_positions)
		# Finally, run phasing algorithm
		dp_table = PedigreeDPTable(selected_reads, recombination_costs, pedigree, distrust_genotypes, accessible_positions)
		superreads_list, transmission_vector = dp_table.get_super_reads()
		read_partitions = dp_table.get_optimal_partitioning()
		print(read_partitions)


def add_arguments(parser):
	arg = parser.add_argument
	# Positional arguments
	arg('locus_file', metavar='LOCUS', help='variants in LOCUS file to phase')
	arg('gam_file', metavar='PHASEINPUT',
	    help='read alignments in GAM file ')

def main(args):
	run_phaseg(**vars(args))
