import copy

"""
Connectivity analysis:
Implementation of a graph, where the nodes are the reads, represented in class element, which store the covered positions
and is identified, by a tuple of the minimal position ,maximal position and the index of the read.

The edges in the graph, include as weight the number of  positions covering both nodes.

"""
class Element:

    def __init__(self,positions,i ):
        self.positions = set(positions)
        self.min = min(positions)
        self.max=max(positions)
        self.index=i
        #include connections, or edges (Element,w) means edge to Element with weight w
        self.connections=[]
        self.component=None

    def add_connection(self,Element,w):
        self.connections.append((Element,w))

    def get_positions(self):
        return self.positions

    def get_min(self):
        return self.min

    def get_max(self):
        return self.max

    def get_index(self):
        return self.index

    def get_connections(self):
        return self.connections

    def get_component(self):
        return  self.component

    def set_component(self,component):
        self.component=component


    #Adding comparison methods for sorting

    #=
    def __eq__(self, other):
        return self.index==other.index

    #x<y
    def __lt__(self, other):
        return (self.min<other.max or (self.min==other.min and self.max<self.max))

    #x<=y
    def __le__(self, other):
        return (self.min<self.other or (self.min==other.min and self.max<=self.max))

    #x!=y
    def __ne__(self, other):
        return (self.min==self.max and self.min==self.max and self.positions==other.positions)

    #x>y
    def __gt__(self, other):
        return (self.min>other.min or (self.min==other.min and self.max>other.max))

    #x>=y
    def __ge__(self, other):
        return (self.min>other.min or(self.min==other.min and self.max >=other.max))


class read_positions_graph:


    def __init__(self, reads):
        #print('In inti positions')
        self.nodes = []
        i=0
        while i != len(reads):
            positions = [variant.position for variant in reads[i]]
            new_node=Element(positions,i)
            self.nodes.append(new_node)
            i+=1
        #Discard here sorting
        #TODO do not know if needed , like in test case also without sorting the nodes are sorted because the reads are already sorted
        #sorted_nodes=sorted(self.nodes)
        #sorted(self.nodes)


        self.components=[]
        self.edges=[]
        #node which connections are added
        for node in self.nodes:
            node_min=node.get_min()
            node_max=node.get_max()
            for sec_node in self.nodes:
                #Discard self loops
                if sec_node==node:
                    continue
                sec_min=sec_node.get_min()
                #TODO If sorted works could make here a cut if sec_node greater then node, or if node_max<sec_node_min ist
                if sec_min<=node_max:
                    sec_max=sec_node.get_max()
                    #cases :
                        #1 : min & max are less than sec_min or bigger than sec_max => no intersection possible - no edges possible - printouts
                        #2 : min < sec_min  but max>=sec_min  - analysis overlap - independent of single or paired end
                        #3 : max>sec_max but min <sec_max - analysis overlap - independent of single or paired end
                        #4 :eiter s_min<min<s_max and and s_max>max>s_min or  min<s_min<max or max>s_max>min
                    #Not implemented in this order
                    if node_min<=sec_min:
                        if node_max>sec_max :
                            weight=len(sec_node.get_positions().intersection(node.get_positions()))
                            if ((sec_node,node,weight)not  in self.edges and weight!=0):
                                node.add_connection(sec_node,weight)
                                sec_node.add_connection(node,weight)
                                self.edges.append((node,sec_node,weight))
                        else:
                            if sec_min<=node_max:
                                weight=len(sec_node.get_positions().intersection(node.get_positions()))
                                if ((sec_node,node,weight)not  in self.edges and weight!=0):
                                    sec_node.add_connection(node,weight)
                                    node.add_connection(sec_node,weight)
                                    self.edges.append((node,sec_node,weight))

                    if sec_min<node_min:
                        if node_max>sec_max:
                            weight=len(sec_node.get_positions().intersection(node.get_positions()))
                            if ((sec_node,node,weight)not  in self.edges and weight!=0):
                                node.add_connection(sec_node,weight)
                                sec_node.add_connection(node,weight)
                                self.edges.append((node,sec_node,weight))

                        else:
                            if sec_max>=node_min:
                                weight=len(sec_node.get_positions().intersection(node.get_positions()))
                                if ((sec_node,node,weight)not  in self.edges and weight!=0):
                                    node.add_connection(sec_node,weight)
                                    sec_node.add_connection(node,weight)
                                    self.edges.append((node,sec_node,weight))



    def get_edges(self):
        return self.edges

    def get_components(self):
        return self.components

    def get_nodes(self):
        return self.nodes

    def add_components(self,comp):
        #comp is an disjoint element
        self.components.append(comp)

class Connect_comp:
    #included_nodes: nodes in the connected components
    # blocks: dictionary of  keys (which are the degrees of the edges) and the values is the number of edges of this degree
    #           in the component
    #stored_for_later: Neighbour nodes, which in the first approach are not added, are stored there to later look again to them
    #Analyze_nodes: Nodes which need to be analyzed in the BFS approach
    #Length : Number of positions in the connected component
    #max = maximal position
    #min = min position


    def __init__(self,node,factor):
        self.positions=node.get_positions() # set
        self.included_nodes =[node]
        self.stored_for_later=[]
        self.analyze_nodes=[]
        self.blocks={}
        neighbours=node.get_connections()

        for (element,w) in neighbours:
            if w>=factor:
                self.analyze_nodes.append(element)
                self.blocks[w]=1
            else:
                self.stored_for_later.append(element)
        #blocks store the number of connection between reads,depending on the degree w
        self.length= len(self.positions)
        self.max=max(self.positions)
        self.min=min(self.positions)


    def get_length(self):
        return self.length

    def get_positions(self):
        return self.positions

    def get_included_nodes(self):
        return self.included_nodes

    def get_analyze_nodes(self):
        return self.analyze_nodes

    def get_blocks(self):
        return self.blocks

    def get_stored_for_later(self):
        return self.stored_for_later

    def get_min(self):
        return self.min

    def get_max(self):
        return self.max

    def update_component(self,node):
        #print('Update component of comp and node')
        #print(self.positions)
        #print(node.get_positions())
        #node added to the actual component
        new_position=node.get_positions()
        self.positions=self.positions.union(new_position)
        self.included_nodes.append(node)
        self.length= len(self.positions)
        self.max=max(self.positions)
        self.min=min(self.positions)

    def expand_component(self,node,factor,not_seen_list):
        #Expand the component by the neighbours of the given node
        neighbours=node.get_connections()
        #TODO need to consider here, that the element positions, which at first iteration was w could have raised - moved to later
        for (element,w) in neighbours:
            #No double occurence of nodes
            if (w>=factor and element in not_seen_list):
                #print('Element added to analyzed')
                #print(element.get_positions())
                if element not in self.analyze_nodes:
                    self.analyze_nodes.append(element)
                if element in self.stored_for_later:
                    self.stored_for_later.remove(element)
                if w not in self.blocks.keys():
                    self.blocks[w]=1
                else:
                    self.blocks[w]+=1
            else:
                if element in not_seen_list and element not in self.analyze_nodes:
                    self.stored_for_later.append(element)

def Find_connected_component_of_this_node(actual_node,graph,factor, not_seen_list):
    #Included some assertions
    new_component=Connect_comp(actual_node,factor)
    assert new_component.get_included_nodes()[0].get_positions()==actual_node.get_positions()
    neighbour_nodes=new_component.get_analyze_nodes()
    assert len(new_component.get_stored_for_later())+ len(new_component.get_analyze_nodes())==len(actual_node.get_connections())
    while ((len(neighbour_nodes)!=0) or (len(new_component.get_stored_for_later())!=0)):
        while len(neighbour_nodes)!=0:
            ana_node=neighbour_nodes.pop()
            new_component.update_component(ana_node)
            not_seen_list.remove(ana_node)
            new_component.expand_component(ana_node,factor,not_seen_list)
            ana_node.set_component(new_component)
        #Need to extend the component by the stored_for_later
        further_analisis_nodes=new_component.get_stored_for_later()
        com_pos=new_component.get_positions()
        if len(further_analisis_nodes) !=0:
            #check if something has changed looking at the component.
            while len(further_analisis_nodes)!=0:
                analyzenode=further_analisis_nodes.pop()
                if analyzenode in not_seen_list:
                    pos_of_node=analyzenode.get_positions()
                    #looking at intersection of read with whole component
                    if (len(pos_of_node.intersection(com_pos))>=factor):
                        new_component.update_component(analyzenode)
                        new_component.expand_component(analyzenode,factor,not_seen_list)
                        not_seen_list.remove(analyzenode)
                        analyzenode.set_component(new_component)
    return not_seen_list,new_component

def find_components_of_graph(graph,factor):
    #Gets the graph and the connection factor as input and computes for the graph the connected components.
    not_seen_list=copy.deepcopy(graph.get_nodes())
    #Till all nodes are assigned to a component
    while len(not_seen_list)!=0:
        #get one node for analysis, is even removed from not seen list
        actual_node=not_seen_list.pop()
        not_seen_list,component=Find_connected_component_of_this_node(actual_node,graph,factor, not_seen_list)
        graph.add_components(component)
    return graph.get_components()



