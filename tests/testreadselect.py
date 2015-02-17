
from whatshap.readselect import Bestreads
from ..core import Read, ReadSet, DPTable, IndexSet

def test_priority_construction():
    vcf_indices = {}
    vcf_indices[100]= '0'
    vcf_indices[101]='1'
    vcf_indices[103]='2'
    print(vcf_indices)
    #look like in the vcf file
    read1 = ((100, 'G','Entry(MINOR)'),(101,'T','Entry(MINOR)'))
    read2 = ((101, 'A','Entry(MINOR)'),(103,'C','Entry(MINOR)'))
    read3 =((100,'G', 'Entry(MINOR)'),(101, 'A', 'Entry(MINOR)'),(103,'C', 'Entry(MINOR)'))
    readset= set()
    readset.add(read1)
    readset.add(read2)
    readset.add(read3)
    print('readset')
    print(readset)

    br= Bestreads(readset,vcf_indices)
#    br.priorityqueue_construction(pq)



def test_read_selection(self, pq,SNP_dict, max_coverage):


    Bestreads(readset,vcf_indices)
