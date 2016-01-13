"""
Connectivity analysis:
want to implement a graph, where the nodes are the reads, here called element, which store the covered positions
and is identified, by a tuple of the minimal position ,maximal position and tht corresponding index.

The edges in the graph, include as weight the number of  positions covering both nodes.

In order to find the connected components in the graph, use ing either the BFS ir DFS to finde connected components,
but only with higher or even to the given weight for the edges.
During the search the components should be stored ina disjoint data structure.


"""
class Element:
	def __init__(self,positions,i ):
		self.positions = set(positions)
		self.min = min(positions)
		self.max=max(positions)
		self.index=i
		#include connections, or edges (Element,w) means edge to Element with weight w
		self.connections=[]

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



class read_positions_graph:
	"""
	Stores the graph by nodes and edges,
	also supports the DFS for finding connected components depending on the overlapping positions

	"""





	def __init__(self, reads):
		#Need to check if .positions() is the right method
		self.nodes = { x: Element(x.positions(),i) for i,x in enumerate(reads) }
		#Should maybe sort nodes
		#TODO NOT SURe If it works..
		sorted_nodes=sorted(self.nodes, cmp=compare)
		#Stored the components of the graphs.
		self.components=[]
		self.edges=[]
		#node which connections are added
		for node in self.nodes:
			#second node : referenz um edges zu definieren
			node_min=node.get_min()
			node_max=node.get_max()
			for sec_node in self.nodes:
				sec_min=node.get_min()
			#TODO If sorted works could make here a cut if sec_node greater then node, or if node_max<sec_node_min ist
				#if sec_min>node_max:
				#	break
				sec_max=node.get_max()
				#cases :
					#1 : min & max are less than sec_min or bigger than sec_max => no intersection possible - no edges possible - printouts
					#2 : min < sec_min  but max>=sec_min  - analysis overlap - independent of single or paired end
					#3 : max>sec_max but min <sec_max - analysis overlap - independent of single or paired end
					#4 :eiter s_min<min<s_max and and s_max>max>s_min or  min<s_min<max or max>s_max>min
				#Not implemented in this order
				if node_min<sec_min:
					if node_max>sec_max:
						#sec_node positions are included into node positions
						weight=len(sec_node.get_positions().intersection(node.get_positions()))
						node.add_connection(sec_node,weight)
						sec_node.add_connection(node,weight)
						self.edges.append((node,sec_node,weight))
					else:
						if sec_min>node_max:
							print('OUt of range.. no overlap of bothe reads')
						else:
							#max of sec_node  is higher than node max, so  coveraing at least position node_max
							weight=len(sec_node.get_positions().intersection(node.get_positions()))
							node.add_connection(sec_node,weight)
							sec_node.add_connection(node,weight)
							self.edges.append((node,sec_node,weight))
				if sec_min<node_min:
					if node_max>sec_max:
						weight=len(sec_node.get_positions().intersection(node.get_positions()))
						node.add_connection(sec_node,weight)
						sec_node.add_connection(node,weight)
						self.edges.append((node,sec_node,weight))
					else:
						if sec_max<node_min:
							print('Out of range.. no overlap possible')
						else:
							#node completely included into sec node
							#Need to use at all stations the intersection because we could also llok at paired end reads, so the complete including of the read
							#does not mean that alll positions of this read are also covered
							weight=len(sec_node.get_positions().intersection(node.get_positions()))
							node.add_connection(sec_node,weight)
							sec_node.add_connection(node,weight)
							self.edges.append((node,sec_node,weight))
						print('Other Case 2 ')

def compare(item1,item2):
	'''
	Compare function for comparing the node elements
	:param item1 and item2 are  objects of type Element
	:return: if item 1 less ,equal or greater than item 2
	'''
	item1_min=item1.get_min()
	item1_max=item1.get_max()
	item2_min=item2.get_min()
	item2_max=item2.get_max()
	#item 1 has lower position it is less
	if item1_min<item2_min:
		return -1
	#item1 and item 2 have the same  min, so either equal or
	else:
		if item1_min==item2_min:
			if item1_max==item2_max:
				return 0
			else:
				if item1_max<item2_max:
					return -1
				else:
					return 1
		else:
			return 1
