from blockparser import scoring_computation
from helperfunctions import new_haplos, create_blocks, compute_incomplete_haplos, StrandSeq_improvepaths
import customcontainer
import sys

if (__name__ == '__main__'):
	#argv1: original VCF file
	#argv2: target VCF file
	#argv3: refpanel
	#argv4: mutation parameter
	#argv5: output for the resulting vcf
	#Creates haplotype strings where each position that is missing (in comparison to an original haplotype file) gets an unknown haplotype (1/1 or 0/0 by default)
	(haplo1, haplo2, intE) = compute_incomplete_haplos(sys.argv[1],sys.argv[2], "Strandseq_intE.txt", "Strandseq_haplotypes.txt")
	print("Step 1/4: Haplotypes created")
	#reads the phase set information from the VCF and returns a dictionary containing variant positions and its phase set, as well as block ending positions and internal blocks 
	(E_whole, haplo, intE) = create_blocks(sys.argv[1], "Strandseq_ends.txt")
	print("Step 2/4: Blocks created")
	#performs the computation of the scoring matrix (dynamic programming) and the backtracing to find the optimal pair of paths through the panel. Paths and corresponding costs are stored in files.
	scoring_computation("Strandseq_haplotypes.txt", "Strandseq_intE.txt", sys.argv[3], sys.argv[4], "Strandseq_costs.txt", "Strandseq_paths.txt")
	print("Step 3/4: Paths created")
	#update the haplotypes according to the paths and resolve any existing switches. Returns two haplotypes that are optimal for the pair of paths.
	(newhaplo1, newhaplo2) = StrandSeq_improvepaths("Strandseq_paths.txt", "Strandseq_haplotypes.txt","Strandseq_intE.txt")[1]
	print("Step 4/4: Haplotypes and paths updated. Now writing to output file (can take several minutes).")	
	new_haplo_dict = customcontainer.DefaultOrderedDict()
	i = 0
	for key in haplo.keys():
		new_haplo_dict[key] = (newhaplo1[i],newhaplo2[i])
		i += 1	
	key_dict = customcontainer.DefaultOrderedDict()
	with open(sys.argv[1]) as f:
		with open(sys.argv[5], 'w') as f1:
			for line in f:
				if line.startswith('#'):
					f1.write(line)	
				else:
					for key in new_haplo_dict.keys():
						if (str(key) in line):
							key_dict[key] = line
			for key in key_dict.keys():
				if ("0|1" in key_dict[key]):
					f1.write(key_dict[key].replace("0|1",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1]))
				elif ("1|0" in key_dict[key]):
					f1.write(key_dict[key].replace("1|0",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1]))
				elif ("1/1" in key_dict[key]):
					f1.write(key_dict[key].replace("1/1",new_haplo_dict[key][0]+'|'+new_haplo_dict[key][1]))