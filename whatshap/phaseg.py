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


from contextlib import ExitStack
from .vcf import VcfReader, PhasedVcfWriter
from . import __version__
from .core import ReadSet, Pedigree, PedigreeDPTable, NumericSampleIds, PhredGenotypeLikelihoods
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
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
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
	"""
	gam_list = []
	with stream.open(locus_file, "rb") as istream:
		for data in istream:
			aln = vg_pb2.Alignment()
			aln.ParseFromString(data)
			gam_list.append(aln)

		
	locus_list = []
	with stream.open(gam_file, "rb") as istream:
		for data in istream:
			aln = vg_pb2.Locus()
			aln.ParseFromString(data)
			locus_list.append(aln)

	  
	# create a dictionary of branches for each locus based on locus file.
	locus_branch_mapping=defaultdict()
	locus_count=0
	alleles_per_pos=[]
	for l in locus_list:
		per_locus_position=[]
		for paths in l.allele:
			tmp=[]
			for mappings in paths.mapping:
				tmp.append(mappings.position.node_id)
			per_locus_position.append(tmp)
		alleles_per_pos[locus_count] = len(per_locus_position)
		locus_branch_mapping[locus_count]=per_locus_position
		locus_count=locus_count+1
		
	# key is the values in locus_branch_mapping and value is triplet(locus, branch, alleles to be flipped)
	reverse_mapping= defaultdict(list)
	for k,v in locus_branch_mapping.items():
		for i,b in enumerate(v):
			reverse_mapping[tuple(b)]=[k,i, len(v)-1]

		
	# extract reads from GAM file associated with the locus and create a sorted readset.
	readset=ReadSet()  

	for g in gam_list:
		# hard-coded source id, mapping quality and other values.
		read=Read(g.name, 0, 0, 0)
		tmp=[]
		for i in g.path.mapping:
			tmp.append(i.position.node_id)
			if tuple(tmp) in reverse_mapping:
				read.add_variant(reverse_mapping[tuple(tmp)][0], reverse_mapping[tuple(tmp)][1], [0]*reverse_mapping[tuple(tmp)][2])
				tmp=[]
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
	logger.info("This is WhatsHap %s running under Python %s", __version__, platform.python_version())
	with ExitStack() as stack:
		logger.info('Using uniform recombination rate of %g cM/Mb.', recombrate)
		all_reads, alleles_per_pos = vg_reader(locus_file, gam_file)
		selected_reads = select_reads(all_reads, 15)
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


def add_arguments(parser):
	arg = parser.add_argument
	# Positional arguments
	arg('locus_file', metavar='LOCUS', help='variants in LOCUS file to phase')
	arg('gam_file', nargs='*', metavar='PHASEINPUT',
	    help='read alignments in GAM file ')

def main(args):
	run_phaseg(**vars(args))
