import sys
import stream
import logging
import vg_pb2
from collections import Counter
from collections import defaultdict
import collections
from collections import defaultdict

file_input = sys.argv[1]
#file_out = argv[2]
out = open(file_input + '.gfa', 'w')
#in_file = sys.argv[2]


def reverse_complement(seq):
	seq_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'T', 'c': 'G', 'g': 'C', 't': 'A'}
	return "".join([seq_dict[base] for base in reversed(seq)])

nodes_list = set()
dummy_list = ['0']*1000
orderalignment = defaultdict(list)
orderalignment = defaultdict(lambda: [-1]*10000000, orderalignment)
with stream.open(str(file_input), "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		read_info = g.name.split("_")
		canu_name = read_info[0]

		canu_chunk_num = int(read_info[1].strip('0'))
		orderalignment[canu_name].insert(canu_chunk_num, g)

new_orderalignment = defaultdict(list)
for k,v in orderalignment.items():
        new_orderalignment[k] = [x for x in v if x != -1]

with stream.open(sys.argv[2], 'wb') as ostream:
	count = 0
	for k,v in new_orderalignment.items():
		#new_g = vg_pb2.Alignment()
		for i in range(0,len(v)-1):
			new_g = vg_pb2.Alignment()
			count=count+1
			new_g.name = v[i].name.split("_")[0]+str(count)
			#print(g.path.mapping[0])
			new_g.path.mapping.extend([v[i].path.mapping[-1], v[i+1].path.mapping[0]])
			#new_mapping = new_g.path.mapping.add()
			#new_mapping = g.path.mapping[0]
			#new_g.path.mapping[:] = [g.path.mapping[0]]
			ostream.write(new_g)
			#kprint('{"name": "'+ g.name.split("_")[0] + '", "path": {"mapping": [{"position":{"node_id":'g.path.mapping[0].position)
			#print(g.path.mapping[len(g.path.mapping)-1])
			#print(']}}')
