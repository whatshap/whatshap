from blockparser import scoring_computation
from helperfunctions import  compute_original, new_haplotypes, compute_switcherrors, create_blocks
import customcontainer
import sys

if (__name__ == '__main__'):
	#argv1: the target file, phased by WhatsHap
	#argv2: the original file (Platinum Genome)
	#blockends.txt: contains the indices where a haplotype block of argv1 ends
	#haplotypes.txt: contains the two haplotype strings extracted from the phasing information in argv1 in the first two lines and the original haplotypes (extracted from argv2) in lines 3 and 4
	#argv3: the refpanel
	#argv4: the  mutation parameter
	
	#reads the phase set information from the VCF and returns a dictionary containing variant positions and its phase set, as well as block ending positions and internal blocks 
	(E_whole, haplo, intE) = create_blocks(sys.argv[1], "blockends.txt")
	print("Step 1/5: Haplotype blocks created")
	#for later evaluation: Extract strings from the given haplotype block sequence and from the original haplotype file and store them in haplotypes.txt
	(originalhaplo1, originalhaplo2) = compute_original(sys.argv[2], sys.argv[1], "haplotypes.txt")[1]
	print("Step 2/5: Original and target haplotype strings created")
	#performs the computation of the scoring matrix (dynamic programming) and the backtracing to find the optimal pair of paths through the panel. Paths and corresponding costs are stored in files.
	scoring_computation("haplotypes.txt", "blockends.txt", sys.argv[3], sys.argv[4], "costs.txt", "paths.txt")
	print("Step 3/5: Paths through the reference panel created")
	#update the haplotypes accordingly to the paths and resolve any existing switches. Returns two haplotypes that are optimal for the pair of paths.
	(newhaplo1, newhaplo2) = new_haplotypes("haplotypes.txt", E_whole, "paths.txt", intE)
	print("Step 4/5: New haplotypes created")
	#comparison to the original haplotypes: Switch error rate within blocks and at block ends is computed.
	(SER, SER_B, number_of_switches,number_of_end_switches, testswitch, switcherrorlist, flipnumber, withinlist) = compute_switcherrors(newhaplo1, originalhaplo1, newhaplo2, originalhaplo2, "blockends.txt")
	print("Step 5/5: Switch error rates computed. Now writing to output file (can take several minutes).")
	print("switch error rate within blocks: ","%.2f" % SER,"%.", " Corresponds to",number_of_switches, "wrongly phased sites out of", len(originalhaplo1))
	print("switch error rate at block ends: ","%.2f" % SER_B, "%.", " Corresponds to",number_of_end_switches, "wrongly phased block ending positions out of",len(E_whole))
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
