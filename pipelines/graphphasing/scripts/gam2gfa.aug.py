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
in_file = sys.argv[2]


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

nodes_done = set()


for k,v in new_orderalignment.items():
        canu_piece_ends=[]
        graph_updated_nodes = []
        for g in v:
                checker0 = '+'
                checker1 = '+'
                if g.path.mapping[0].position.is_reverse == True:
                        checker0 = '-'
                if g.path.mapping[len(g.path.mapping)-1].position.is_reverse == True:
                        checker1 = '-'

                seq1 = g.path.mapping[0].edit[0].sequence
                seq2 = g.path.mapping[len(g.path.mapping)-1].edit[0].sequence
                if seq1 == '':
                        seq1= g.path.mapping[1].edit[0].sequence
                if seq2 == '':
                        seq2 = g.path.mapping[len(g.path.mapping)-2].edit[0].sequence
                if checker0 == '-':
                        seq1 = reverse_complement(g.path.mapping[0].edit[0].sequence)
                        if seq1 == '':
                                seq1 = reverse_complement(g.path.mapping[1].edit[0].sequence)
                if checker1 == '-':
                        seq2 = reverse_complement(g.path.mapping[len(g.path.mapping)-1].edit[0].sequence)
                        if seq2=='':
                                seq2 = reverse_complement(g.path.mapping[len(g.path.mapping)-2].edit[0].sequence)
                canu_piece_ends.append((g.path.mapping[0].position.node_id, checker0))
                canu_piece_ends.append((g.path.mapping[len(g.path.mapping)-1].position.node_id, checker1))
                graph_updated_nodes.append((g.path.mapping[0].position.node_id, seq1))
                graph_updated_nodes.append((g.path.mapping[len(g.path.mapping)-1].position.node_id, seq2))
                nodes_done.add(g.path.mapping[0].position.node_id)
                nodes_done.add(g.path.mapping[len(g.path.mapping)-1].position.node_id)
                print(graph_updated_nodes)
                print(canu_piece_ends)
        for p in range(0,len(canu_piece_ends)-2,2):
                out.write("L" + "\t" + str(canu_piece_ends[p+1][0])+"\t" + str(canu_piece_ends[p+1][1]) + "\t" + str(canu_piece_ends[p+2][0]) + "\t" + str(canu_piece_ends[p+2][1]) +"\t" + "0M" + "\n")

        for s in range(0,len(graph_updated_nodes)):
                out.write("S" + "\t" + str(graph_updated_nodes[s][0]) + "\t" + str(graph_updated_nodes[s][1])+ "\n")

with open(in_file) as fp:
        for line in fp:
                var=line.rstrip().split('\t')
                if var[0] == 'S' and int(var[1]) not in nodes_done:
                        out.write(line)
                if var[0] == 'L':
                        out.write(line)
