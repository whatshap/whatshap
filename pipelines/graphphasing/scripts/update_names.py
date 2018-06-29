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
ostream = stream.open(sys.argv[2], 'wb')
with stream.open(str(file_input), "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		read_info = g.name.split("_")
		canu_name = read_info[0]
		new_g = vg_pb2.Alignment()	
		new_g.name = read_info[0]
		new_g.path.mapping.extend([g.path.mapping[0]])
		new_g.query_position = g.query_position
		ostream.write(new_g)
