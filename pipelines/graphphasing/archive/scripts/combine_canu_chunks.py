import sys
import stream
import logging
import vg_pb2
from collections import Counter
from collections import defaultdict
import collections
from collections import defaultdict

file_input = sys.argv[1]
out = open(file_input + '.gfa', 'w')


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

ostream = stream.open(sys.argv[2], 'wb')

for k,v in new_orderalignment.items():
        new_g = vg_pb2.Alignment()
        for i in range(0,len(v)):
                new_g.name = v[i].name.split("_")[0]
                tmp = v[i].name.split("_")[0]
                for j in range(0,len(v[i].path.mapping)):
                        new_g.path.mapping.extend([v[i].path.mapping[j]])
        ostream1 = stream.open(tmp + '.gam', 'wb')
        ostream1.write(new_g)
        ostream1.close()
        ostream.write(new_g)


