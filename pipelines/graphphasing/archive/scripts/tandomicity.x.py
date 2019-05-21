import sys
sys.path.append('/local/data/yeast_analysis/scripts/lib/python3.6/site-packages/stream')
import stream
import logging
import vg_pb2
from collections import Counter
from collections import defaultdict
import collections
from collections import OrderedDict, namedtuple
from collections import defaultdict

# tandom and interspersed repeat from both aligned pacbio reads and true_haps
trans_filename = sys.argv[1]
gam_filename = sys.argv[2]
#true_haps_filename = sys.argv[3]
parameter_interspersed = sys.argv[3]
out_filename = sys.argv[4]

out_file = open(out_filename, 'w')


d={}
count=1

bubbles_start = set()
#with stream.open('assembly_graph.P.int.remn2n.X_100.chrXIII.trans' ,"rb") as istream:
with stream.open(str(trans_filename) ,"rb") as istream:
	for data in istream:
		l = vg_pb2.SnarlTraversal()
		l.ParseFromString(data)
		if l.snarl.start.backward == True:
			start_node = l.snarl.end.node_id
		else:
			start_node = l.snarl.start.node_id
		bubbles_start.add(start_node)


multiplicity_bubbles = defaultdict(list)
read_details = defaultdict(list)
#with stream.open('../out.new.gam', "rb") as istream:
with stream.open(str(gam_filename), "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		tmp = []
		for i in range(0,len(g.path.mapping)):
			node = g.path.mapping[i].position.node_id
			if node in bubbles_start:
				multiplicity_bubbles[node].append(g.name)
			if node in tmp:
				out_file.write(str(node)+ '\n')
			tmp.append(node)

#with stream.open(str(true_haps_filename), "rb") as istream:
#	for data in istream:
#		g = vg_pb2.Alignment()
#		g.ParseFromString(data)
#		tmp = []
#		for i in range(0,len(g.path.mapping)):
#			node = g.path.mapping[i].position.node_id
#			if node in bubbles_start:
#				multiplicity_bubbles[node].append(g.name)
#			if node in tmp:
#				out_file.write(str(node)+ '\n')
#			tmp.append(node)

count=0
repeaticity = defaultdict()
repeaticity_read_support = defaultdict()
repeaticity_i = defaultdict()
repeaticity_read_support_i = defaultdict()

for k,v in multiplicity_bubbles.items():
	tmp = defaultdict(set)
	for item,count in collections.Counter(v).items():
		tmp[count].add(item)
	max_val = 0 
	item_val = set()
	l1_max = 0
	# how many reads support multiplicity l1
	for l1,l2 in tmp.items():
		if max_val < len(l2):
			max_val = len(l2)
			item_val = l2
			l1_max = l1
	#if max_val > 5 and l1_max > 1: # minimum read support and minimum multiplicty
	if l1_max > 1:
		out_file.write(str(k)+ '\n')
		#repeaticity[k] = l1_max
		#repeaticity_read_support[k] = item_val
		#print(item_val, max_val, l1_max, k, v)
	if len(set(v)) > int(parameter_interspersed):
		#print(int(parameter_interspersed))
		out_file.write(str(k)+ '\n')
		#cnv = len(set(v))/38.0
		#repeaticity_i[k] = cnv
		#repeaticity_read_support_i[k] = set(v)
		#print('interpersed repeat')
		#true_chrms = set()
		#for i in set(v):
		#	var = var = i.split("_")[2]
		#	true_chrms.add(var)
		#print(set(v), k, cnv, len(true_chrms))
		#if int(round(cnv)) == len(true_chrms):
		#	count=count+1
		

	

