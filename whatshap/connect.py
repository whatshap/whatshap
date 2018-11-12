"""
Use a reference panel to connect the haplotype blocks

"""
import logging
from .blockparsing import compute_referencepanel, create_blocks, compute_haplotypes, update_haplotypes
from .customcontainer import DefaultOrderedDict
import sys
from cyvcf2 import VCF
import platform
import resource
from .core import scoring_computation
from collections import defaultdict
import collections
import math
from . import __version__
from .timer import StageTimer
from .utils import detect_file_format
logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('variant_file', metavar='BLOCKS', help='VCF file with phased haplotype blocks for one chromosome')
	arg('reference_file', metavar='REFERENCE', help='VCF file that serves as reference panel. This is used to connect haplotype blocks. ')
	arg('output_file', metavar="OUTPUT",help='Output file in VCF format')
	arg('--include-unphased',dest = 'include_unphased', default=False, action='store_true', help = 'Include unphased variant positions for inferring the haplotypes at these sites (Default: Use phased variants only)')
	
def validate(args, parser):
	pass

def included_variants(ref_file, variant_file, include_unphased):
	ref_set = set()
	chromosome_reference_set = set()
	refvcf = VCF(ref_file)
	for variant in refvcf:
		ref_set.add(variant.POS)
		chromosome_reference_set.add(variant.CHROM)
	#for multiple chromosomes, raise an error
	if (len(chromosome_reference_set)>1):
		logger.error("The given reference file is supposed to contain variants from one chromosome only!")
		sys.exit(1)
	#create the list of variants to be considered	
	variant_set = set()
	chromosome_variant_set = set()
	vcf = VCF(variant_file)
	haplo1,haplo2 = "",""	
	for v in vcf:
		#use variant if it is present in the panel, heterozygous and phased
		if (include_unphased==False) and (v.POS in ref_set) and (v.gt_types[0] == 1) and (v.gt_phases[0] == True) and (v.is_snp):	
			if (v.CHROM not in chromosome_reference_set):
				logger.error("Reference and target file must contain variants from the same chromosome!")
				sys.exit(1)
			variant_set.add(v.POS)
			chromosome_variant_set.add(v.CHROM)
		elif include_unphased and (v.POS in ref_set) and (v.gt_types[0] == 1) and (v.is_snp):
			if (v.CHROM not in chromosome_reference_set):
				logger.error("Reference and target file must contain variants from the same chromosome!")
				sys.exit(1)
			variant_set.add(v.POS)
			chromosome_variant_set.add(v.CHROM)
	if (len(chromosome_variant_set)>1):
		logger.error("The given variant file is supposed to contain variants from one chromosome only!")
		sys.exit(1)
	if not variant_set:
		logger.error("No positions found in both the reference and the target file.")
		sys.exit(1)	
	return(variant_set)

def run_connect(variant_file,reference_file,output_file, include_unphased=True):
	timers = StageTimer()
	timers.start('overall')
	logger.info("This is WhatsHap (reference-based phasing) %s running under Python %s", __version__, platform.python_version())
	with timers("variant_list"):
		logger.info("Choosing suitable variant positions.")
		variant_set = included_variants(reference_file, variant_file,include_unphased)

	with timers("ref_panel"):
		fileformat = detect_file_format(reference_file)
		if fileformat == 'VCF':
			logger.info("Computing reference panel. This action may take a few minutes.")
			compute_referencepanel(reference_file,variant_file,variant_set, "ref.txt")
		else:
			logger.error("Wrong file format. Please enter a reference file in VCF format.")
	#reads the phase set information from the VCF and returns a dictionary containing variant positions and its phase set, as well as block ending positions and internal blocks 
	with timers("creating_blocks"):
		logger.info("Reading phased blocks")
		(E_whole, haplo, intE) = create_blocks(variant_file, "blockends.txt", variant_set)
		#writes the original haplotypes from the target into strings
		compute_haplotypes(variant_file, "haplotypes.txt", variant_set)

	#performs the computation of the scoring matrix (dynamic programming) and the backtracing to find the optimal pair of paths through the panel. Paths and corresponding costs are stored in files.
	with timers("dynamic_programming"):
		logger.info("Dynamic programming: Performing haplotype block parsing")
		scoring_computation("haplotypes.txt", "blockends.txt", "ref.txt", 1, "costs.txt", "paths.txt")

	#update the haplotypes accordingly to the paths and resolve any existing switches. Returns two haplotypes that are optimal for the pair of paths.
	with timers("creating_haplotypes"):
		logger.info("Assembling the resulting haplotypes")
		(newhaplo1, newhaplo2) = update_haplotypes("haplotypes.txt", E_whole, "paths.txt", intE)

	with timers("output"):
		logger.info("Writing output file. This action may take a few minutes.")
		new_haplo_dict = DefaultOrderedDict()
		i = 0
		for key in haplo.keys():
			new_haplo_dict[key] = (newhaplo1[i],newhaplo2[i])
			i += 1
		key_dict = DefaultOrderedDict()
		with open(variant_file) as f:
			with open(output_file, 'w') as of:
				counter = 0
				current = 0
				for line in f:				
					counter += 1
					key_found = False
					if line.startswith('#'):
						of.write(line)	
					else:
						for key in list(new_haplo_dict.keys())[current:]:
							if (str(key) in line.split('.')[0]):	
								current = list(new_haplo_dict.keys()).index(key)
								key_found = True
								if ("0|1" in line):
									of.write(line.replace("0|1",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1]))								
								elif ("1|0" in line):
									of.write(line.replace("1|0",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1]))
								break
						if (key_found == False):
							of.write(line)

	logger.info('\n== SUMMARY ==')
	timers.stop("overall")
	if sys.platform == 'linux':
		memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
		logger.info('Maximum memory usage: %.3f GB', memory_kb / 1E6)
	logger.info('Time spent choosing suitable variant positions:                      %6.1f s', timers.elapsed('variant_list'))
	logger.info('Time spent creating reference panel:                      %6.1f s', timers.elapsed('ref_panel'))
	logger.info('Time spent creating haplotype blocks:                      %6.1f s', timers.elapsed('creating_blocks'))
	logger.info('Time spent finding optimal paths through the reference panel:                      %6.1f s', timers.elapsed('dynamic_programming'))
	logger.info('Time spent assembling the connected haplotypes:                      %6.1f s', timers.elapsed('creating_haplotypes'))
	logger.info('Time spent writing output file:                      %6.1f s', timers.elapsed('output'))
	logger.info('Total elapsed time:                      %6.1f s', timers.elapsed('overall'))

def main(args):
	run_connect(**vars(args))	