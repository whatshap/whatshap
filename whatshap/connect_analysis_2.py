import copy
import cProfile
"""
Connectivity analysis:
Second approach for the connectivity analysis:
Implementation of a graph , where the nodes are the variant positions,the edges respresent the reads,
and we search the connected component.

For the edges the weight represents, how many reads cover this 2 SNPs.

"""

class Node:

    def __init__(self,position):
        self.position = position
        #Initialize the coverage, so how many reads covering this node
        #self.coverage = 0
        #Initilize to which component the node is assigned to
        self.component=None
        #Store the connections to the other node , just store the position of the nodes, this node is connected with
        self.neighbour_nodes=set()
        self.position=position


    def add_connection(self,node):
        '''
        :param node: the node the actual node has a connection to
        :return: nothing
        Add connection
        '''
        #connection is a set  of node and weight

        self.neighbour_nodes.add(node.get_position())


    def get_position(self):
        return self.position

    def get_connections(self):
        #Return set of connected nodes of the actual node
        return self.neighbour_nodes

    def get_component(self):
        return self.component

    def set_component(self,component):
        self.component=component

    def __hash__(self):
        return self.get_position()

class Edge:

    def __init__(self,pos_1,pos_2,i):
        self.min=min(pos_1,pos_2)
        self.max=max(pos_1,pos_2)
        self.index=set()
        self.index.add(i)
        self.weight=1

    def get_min(self):
        return self.min

    def get_max(self):
        return self.max

    def add_edge_for_same_nodes(self,i):
        self.weight+=1
        self.index.add(i)

    def get_weight(self):
        return self.weight

    def get_indices(self):
        return self.index

    def __hash__(self):
        #TODO : NOT SO SURE IF IT IS A PERFECT HASH FUNCTION
        #unique identifier, average value between these two nodes added with the absolute value between them
        return ((self.max +self.min)/2 )+ (self.max - self.min)

class Edges:
    #Edge structure for the graph
    #Stores the nodes,  just by the nodes
    #Stores the positions in a dictionary, where for every node we have a weight assigned to it, how often it is covered


    def __init__(self):
        self.edges_of_nodes={}

   # def add_edge_to_edges(self,node_1,node_2,i):
   #     edge=Edge(node_1.get_position(),node_2.get_position(),i)
   #     self.pro_forma_edges.append(edge)
        #keys_of_nodes=self.edges_of_nodes.keys()
        #if node_1 in keys_of_nodes:
        #    edges_of_node1=self.edges_of_nodes[node_1]
        #    print('Already in the edges of node 1 ')
        #else:
        #    edge=Edge(node_1.get_position(),node_2.get_position(),i)
        #    self.edges_of_nodes[node_1]=[edge]
        #if node_2 in keys_of_nodes:
        #    edges_of_node2=self.edges_of_nodes[node_2]
        #    print('Already in the edges of node 2')
        #else:
        #    edge=Edge(node_1.get_position(),node_2.get_position(),i)
        #    self.edges_of_nodes[node_2]=[edge]


    def add_edge_to_edges(self,node_1,node_2,i):

        if node_1 in self.edges_of_nodes.keys():
            edges_node_1=self.edges_of_nodes[node_1] # list of edges
            #print('Value of the dict edges of nodes')
            #print(edges_node_1)
            run_variable=0
            #go over edges of node_1
            #print('Length of the edges_ node1 %d' %len(edges_node_1))
            for e in edges_node_1:
                #If this edge exists and  for getting the right edge
                if (e.get_min()==node_2.get_position() or e.get_max()==node_2.get_position()) and (e.get_min()==node_1.get_position() or e.get_max()==node_1.get_position()):
                    e.add_edge_for_same_nodes(i)
                else:
                   run_variable+=1
            if run_variable==len(edges_node_1):
                edge=Edge(node_1.get_position(),node_2.get_position(),i)
                self.edges_of_nodes[node_1].append(edge)
        else:
            edge=Edge(node_1.get_position(),node_2.get_position(),i)
            #if (len(self.edges_of_nodes[node_1])==0):
            #    print('Edges of nodes of node_1 is 0 ')
            self.edges_of_nodes[node_1]=[edge] #else:    #    self.edges_of_nodes[node_1].append(edge)
        if node_2 in self.edges_of_nodes.keys():
            edges_node_2=self.edges_of_nodes[node_2]
            run_variable_2=0
            #go over edges of node_2
            for e in edges_node_2:
                #For getting the right edge
                if (e.get_min()==node_2.get_position() or e.get_max()==node_2.get_position()) and (e.get_min()==node_1.get_position() or e.get_max()==node_1.get_position()):
                    e.add_edge_for_same_nodes(i)
                else:
                    run_variable_2+=1
            if  run_variable_2==len(edges_node_2):
                edge=Edge(node_1.get_position(),node_2.get_position(),i)
                self.edges_of_nodes[node_2].append(edge)
        else:
            edge=Edge(node_1.get_position(),node_2.get_position(),i)
            #if (len(self.edges_of_nodes[node_2])==0):
            #    print('Edges of nodes of node 2 is 0')
            self.edges_of_nodes[node_2]=[edge]
            #else:
            #self.edges_of_nodes[node_2].append(edge)

    def get_edges_for_specific_node(self,node):
        #Returns the edges of a specific node
        if node in self.edges_of_nodes.keys():
            return self.edges_of_nodes[node]
        else:
            return None

    def get_all_edges(self):
        return self.edges_of_nodes


class Graph:

    def __init__(self, reads, important_positions):
        #Initialize Graph by the reads and the positions in the readset
        self.edges=Edges()
        self.nodes={}
        self.already_seen=set()
        #Initialize nodes and store in dictionary
        self.components=set()
        for position in important_positions:
            node=Node(position)
            self.nodes[position]= node
        i=0
        #i represents the read index, for assigning edges to read
        while i!=len(reads):
        #Go over reads
            read_positions = [variant.position for variant in reads[i]]
            #Make connections to the positions which are covered
            while len(read_positions)!=0:
                r_p =read_positions.pop()
                #TODO ab hier
                for r_p_2 in read_positions:
                    #Get nodes by position
                    node_of_r_p=self.nodes[r_p]
                    node_of_r_p_2=self.nodes[r_p_2]
                    node_of_r_p.add_connection(node_of_r_p_2)
                    node_of_r_p_2.add_connection(node_of_r_p)
                    self.edges.add_edge_to_edges(node_of_r_p,node_of_r_p_2,i)
            i+=1

    def get_node_by_position(self,position):
        return self.nodes[position]

    def get_edges(self):
        return self.edges

    def get_components(self):
        return self.components

    def get_nodes(self):
        return self.nodes

    def add_components(self,comp):
        #comp is an disjoint element
        self.components.add(comp)



#    def find_components_of_graph(self,factor):
#        nodes=self.get_nodes()
#        #Go over all nodes
#        while len(nodes)!=0:
#            node=nodes.pop()
#            comp=Component(node)
#            #get edges of this node
#            edges_of_this_node=self.get_edges()[node]
#            #Go over edges if we find edge that fullfills the desired connectivity remove node and add node to the component
#            for e in edges_of_this_node:
#                if e.get_weight()>= factor:
#                    if e.get_min()==node.get_position():
#                        node_2_pos=e.get_max()
#                        node_2=self.get_node_by_position(node_2_pos)
#                        comp.add_node(node_2)
#                        #Not sure if pop works
#                        nodes.pop(node_2)
#                    else:
#                        #e.get_max()==node.get_position():
#                        node_2_pos=e.get_min()
#                        node_2=self.get_node_by_position(node_2_pos)
#                        comp.add_node(node_2)
#                        nodes.pop(node_2)
#            comp.covered_node(node)


 #   def add_nodes_and_followers_to_comp(self,nodes,comp,factor):
 #       already_visited_nodes=comp.get_covered_nodes()
 #       #nodes are from the first iteration the nodes of the original component....
#
#        while len(nodes)!=0:
#            actual_node=nodes.pop()
#            edges_of_this_node=self.edges()[actual_node]
#            for e in edges_of_this_node:
#                if e.get_weight()>= factor:

    def find_components_of_graph(self,factor):
        #Solange nicht alle gesehen wurden
        print('Get nodes keys')
        print(len(self.get_nodes().keys()))
        print(len(self.already_seen))
        while (len(self.already_seen) != len(self.get_nodes().keys())):
            print('In While loop')
            nodes=self.get_nodes().keys()
            for node in nodes:
                print('Node in the find_component of graph')
                print(node)
                if node in self.already_seen:
                    continue
                else:
                    print('In else case')
                    new_comp=Component(node)
                    self.already_seen.add(node)
                    self.recursive_method(factor,new_comp,node)
                print('Found here a component')
                print(new_comp)
                self.add_components(new_comp)



    def get_suitable_children_of_node(self,factor,node):
        #Already seen is a set of positions, which are already seen, so where the corresponding nodes are already covered in a comp
        print('In get suitable children method')
        nodes_for_futher_ana=set()
        graph_edges=self.get_edges()
        Nodes=self.get_nodes()
        node_really=Nodes[node]
        edges_of_node=graph_edges.get_edges_for_specific_node(node_really)
        print('Edges of node')
        print(edges_of_node)
        if edges_of_node!=None:
            print('There are edges for this node')
            for e in edges_of_node:
                if e.get_weight()>= factor:
                    #node_1=self.get_node_by_position(e.get_min())
                    #node_2=self.get_node_by_position(e.get_max())

                    #Work on positions not on the nodes
                    if e.get_min() not in self.already_seen:
                        nodes_for_futher_ana.add(e.get_min())
                    if e.get_max() not in self.already_seen:
                        nodes_for_futher_ana.add(e.get_max())
        return nodes_for_futher_ana

    def recursive_method(self,factor,comp,node):
        print('In recursive method with node')
        print(node)
        children=self.get_suitable_children_of_node(factor,node)
        #which means there are nodes to add to the comp
        while len(children)!=0:
            child=children.pop()
            print('child in the while loop')
            print(child)
            if child not in self.already_seen:
                self.already_seen.add(child)
                child_node=self.get_node_by_position(child)
                comp.add_node(child_node)
                self.recursive_method(factor,comp,child)



class Component:

    def __init__(self,node):
        #Node is already the position of the node
        self.nodes=set()
        self.nodes.add(node)
        self.length=1
        self.min=node
        self.max=node
        self.covered_nodes=set()


    def add_node(self,node):
        print('In add node of Component')
        print(node)
        print(node.get_position())
        print('Self nodes before adding ')
        print(self.nodes)
        self.nodes.add(node.get_position())
        print('Self nodes after adding ')
        print(self.nodes)
        self.length+=1
        self.min=min(self.min,node.get_position())
        self.max=max(self.max,node.get_position())
        print('Component after including node')
        print(self)
        print(self.length)
        print(self.min)
        print(self.max)
        print(self.nodes)

#    def covered_node(self,node):
#        self.covered_nodes.add(node)

#    def get_covered_nodes(self):
#        return self.covered_nodes

    def get_nodes_of_component(self):
        return self.nodes

    def get_length(self):
        return self.length

