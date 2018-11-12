import math
import logging
import vcf
from cyvcf2 import VCF
from collections import defaultdict
from libcpp.vector cimport vector
#from libcpp.set cimport set
#from libcpp.list cimport list
from libcpp.string cimport string
from libcpp.unordered_set cimport unordered_set
from .priorityqueue cimport priority_type, priority_type_ptr, queue_entry_type, PriorityQueue
from core cimport ReadSet
cimport cpp

logger = logging.getLogger(__name__)

cdef extern from "../src/DP_matrix.h":
	float** printlist(vector[int], vector[int], vector[vector[int]],float, set[int] )
	void compute_scoring(char*, char*, char*, float, char*, char*)
	#void compute_score(vector[int] H_A, vector[int] H_B, set[int] E, vector[vector[int]] Ref, string out, float param_mut, float param_switch)


#def list_test(vector[int] H_A, vector[int] H_B, vector[vector[int]] Ref, float param_mut, set[int] E): #, set[int] E, vector[vector[int]] Ref, string out, float param_mut, float param_switch):
#	#compute_score(H_A,H_B,E,Ref,out,param_mut,param_switch)
#	printlist(H_A, H_B, Ref, param_mut,E)



 #distutils: language = c++
#
#import math
#import logging
#import vcf
#from cyvcf2 import VCF
#from collections import defaultdict
#from libcpp.vector cimport vector
#from libcpp.set cimport set
#from libcpp.list cimport list
#from libcpp.string cimport string
#cimport cpp
#logger = logging.getLogger(__name__)
#
#
#cdef extern from "../src/DP_matrix.h":
#	float** printlist(vector[int], vector[int], vector[vector[int]],float, set[int] )
#	void compute_scoring(char*, char*, char*, float, char*, char*)
#	#void compute_score(vector[int] H_A, vector[int] H_B, set[int] E, vector[vector[int]] Ref, string out, float param_mut, float param_switch)
#
#def list_test(vector[int] H_A, vector[int] H_B, vector[vector[int]] Ref, float param_mut, set[int] E): #, set[int] E, vector[vector[int]] Ref, string out, float param_mut, float param_switch):
#	#compute_score(H_A,H_B,E,Ref,out,param_mut,param_switch)
#	printlist(H_A, H_B, Ref, param_mut,E)
#	
#def scoring_computation(haplofile, E_file, panelfile, mutation, out_costfile, out_pathfile):
#	filename_byte_string_haplo = haplofile.encode("UTF-8")
#	filename_byte_string_E = E_file.encode("UTF-8")
#	filename_byte_string_panel = panelfile.encode("UTF-8")
#	filename_byte_string_mutation = mutationfile.encode("UTF-8")
#	filename_byte_string_costs = out_costfile.encode("UTF-8")
#	filename_byte_string_paths = out_pathfile.encode("UTF-8")
#	cdef char* fname_haplo = filename_byte_string_haplo
#	cdef char* fname_E = filename_byte_string_E
#	cdef char* fname_panel = filename_byte_string_panel
#	cdef char* fname_mutation = filename_byte_string_mutation
#	cdef char* fname_costs = filename_byte_string_costs
#	cdef char* fname_paths = filename_byte_string_paths
#	cdef float mutation_parameter = int(mutation)
#	compute_scoring(fname_haplo, fname_E, fname_panel, mutation_parameter, fname_costs, fname_paths)
