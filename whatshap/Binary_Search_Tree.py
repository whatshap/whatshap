import math
from whatshap._core import PyRead,PyReadSet


#Before starting the algorithm we have to build up the binary search tree.
class Binary_Search_Tree:

    def __init__(self,readset):
        '''Build up a binary search tree, with the attributes: Leaf_list and root_node'''
        #Builds List of possible Leafes
        node_list=self.build_list(readset)
        #resolves double occuring nodes
        leaf_list=self.discover_double_and_sibling(node_list)
        layer_array=[]
        complete_tree=self.building_BST_from_leaf_list(0,len(leaf_list)-1,layer_array,leaf_list)
        self.complete_tree=complete_tree
        self.leaf_list=leaf_list


    def build_list(self,_ana_readset):
        '''Building up the list of the nodes representing the reads...'''
        list_for_reads=[]
         #removing reads which cover less than 2 variants
        indices_of_reads = set(i for i, read in enumerate(_ana_readset) if len(read) >= 2)
        for i in indices_of_reads:
            read_of_index= _ana_readset[i]
            first_pos=read_of_index[0]
            last_pos=read_of_index[len(read_of_index)-1]
            #Initialize the nodes by their value and the index,of the read, sonnect them by sibling
            first_Node=Leaf_node(first_pos.position,i)
            second_Node= Leaf_node(last_pos.position,i)
            first_Node.set_sibling(second_Node)
            second_Node.set_sibling(first_Node)
            list_for_reads.append(first_Node)
            list_for_reads.append(second_Node)
        #sort the list for positions or values so that double occuring are in a row and could be removed
        sorted_list=sorted(list_for_reads,key=lambda node :node.value)
        return sorted_list

    def discover_double_and_sibling(self,sorted_list):
        '''
        Detects double occuring leaf nodes erases later and includes the siblings of the erased nodes to the remaining node,
        :param sorted_list: the list of Leaf nodes out of the given readset
        :return: list of leaf nodes with unique leafs and stored siblings and coverage
        '''
        nodes_to_remove=[]
        #count for keeping track of coverage
        coverage_count=0

        #Go over the leaf nodes
        iter=0
        while iter !=len(sorted_list):
            iter_val= sorted_list[iter].get_value()
            #increase coverage if the new node is a start point of a read
            if (iter_val)<(sorted_list[iter].get_sibling()[0].get_value()):
                coverage_count +=1

            #looking at the next leaf node
            new_var= iter+1
            decrease_coverage_counter= 0
            while (new_var!= len(sorted_list) and iter_val==sorted_list[new_var].get_value() ):

                #If new node start position increase value else remember to decreas it again.
                if (sorted_list[new_var].get_value())<(sorted_list[new_var].get_sibling()[0].get_value()):
                    coverage_count +=1
                else:
                    decrease_coverage_counter+=1
                sorted_list[iter].add_sibling(sorted_list[new_var].get_sibling())
                if new_var not in nodes_to_remove:
                    nodes_to_remove.append(new_var)
                new_var +=1

            #Set coverage of the node
            sorted_list[iter].set_coverage(coverage_count)

            #Also need to set the coverage for the other siblings which have the same coverage
            helping_var=iter
            while helping_var!=new_var:
                sorted_list[helping_var].set_coverage(coverage_count)
                helping_var+=1

            #Reduce counter by the value of the nodes which end in this position
            coverage_count=coverage_count-decrease_coverage_counter

            #decrease score if it is the end of the read
            if (iter_val)>(sorted_list[iter].get_sibling()[0].get_value()):
                coverage_count -=1

            #setting iterable to the point where no doubles occure
            iter=new_var

        for returned_index in range(0,len(nodes_to_remove)):
            to_remove=nodes_to_remove.pop()
            sorted_list.pop(to_remove)
        return sorted_list

    def get_complete_tree(self):
        return self.complete_tree

    def get_leaf_list_of_tree(self):
        return self.leaf_list


    def building_BST_from_leaf_list(self,start,end,arr,node_list):
        '''
        Bottom up method for constructing a balances binary seach tree with the corresponding coverage
        :param start: integer with start point in the node_list
        :param end: integer with end point in the list
        :param arr: array representing the whole tree
        :param node_list: list of leaf nodes
        :return: The Tree structure, where we only get the roo_node as output.

        '''

        if (start>end):
            return 0
        if (start==end):
            #root_node=node_list[start]
            return node_list[start]
        else:
            middel= int((start+end)/2)
            (mini,maxi)=self.coverage_of_range(start,end,node_list)
            root_node=BST_node(mini,maxi)
            arr.append(root_node)
            left_node=self.building_BST_from_leaf_list(start,middel,arr,node_list)
            root_node.set_left_child(left_node)
            left_node.set_parent(root_node)
            left_node.set_is_left_child()
            #print('Set up as left child ')
            #print(left_node.get_coverage())
            #if left_node.isLeaf():
                #print('Left node is Leaf')
                #print(left_node.get_value())
            right_node=self.building_BST_from_leaf_list(middel+1,end,arr,node_list)
            root_node.set_right_child(right_node)
            right_node.set_parent(root_node)
            right_node.set_is_right_child()
            #print('Set up as right child ')
            #print(right_node.get_coverage())
            #if right_node.isLeaf():
            #    print('Right node is Leaf')
            #    print(right_node.get_value())


            return root_node


    def coverage_of_range(self,start,stop,nodelist):
        '''
        Returns the minimum and maximum coverage in a specific range
        '''
        #initializing max and min coverage by -1 because the coverage could not be negative but minimum has to be high
        maximum = -1
        minimum = 50
        #go over the nodelist
        for i in range(start,stop+1):
            if (nodelist[i].isLeaf()):
                coverage=nodelist[i].get_coverage()
                if coverage>maximum:
                    maximum=coverage
                if coverage<minimum:
                    minimum=coverage
            #if inner nodes, we have 2 coverages
            else:
                (min_cov,max_cov)=nodelist[i].get_coverage()
                if min_cov<minimum :
                    minimum=min_cov
                if max_cov>maximum:
                    maximum=max_cov
        return (minimum,maximum)







    #Decides if the read is crucial or not , called from the method in the maxflow.py
    #TODO REad again the criterion
    def is_crucial(self,split_node_coverage, max_coverage,split_balance ):
        '''
        Depends on the split_node and the maximal coverage if the read in this split node is crucial or not
        '''
        d_var=False
        (min_cov,max_cov)=split_node_coverage
        updated_min_coverage= min_cov+split_balance
        updated_max_coverage= max_cov+split_balance
        if updated_min_coverage<= math.floor(max_coverage/2):
            d_var=True
        if updated_max_coverage>max_cov:
            d_var=False
        return d_var



    def synchronize_sibling_with_same_value(self,sibling_list):
        '''
        Assume already sorted sibling list of nodes
        :param sibling_list:
        :return:a updated sibling_list
        '''
        #TODO : NOT EFFICIENT AT ALL But Working
        Leaf_list=self.get_leaf_list_of_tree()
        for actual_node in sibling_list:
            actual_value=actual_node.get_value()
            for j in Leaf_list:
                if j.get_value()==actual_value:
                    actual_node.set_balance(j.get_balance())
                    actual_node.set_parent(j.get_parent())
                    if j.get_is_left_child():
                        actual_node.set_is_left_child()
                    else:
                        actual_node.set_is_right_child()
                    break



    def get__all_nodes(self,node,Set_of_nodes):
        '''
        :param node: Start node of the search, could also be seen as the root
        :return: returns of a node a list of all nodes in this tree
        '''
        n_coverag=node.get_coverage()
        if node.isLeaf():
            if (node.get_value()==50):
                print('Added in get alll nodes')
            Set_of_nodes.add(node)
            return Set_of_nodes
        else:
            Set_of_nodes.add(node)
            r_child=node.get_right_child()
            l_child=node.get_left_child()
            Set_of_nodes.union(self.get__all_nodes(r_child,Set_of_nodes))
            Set_of_nodes.union(self.get__all_nodes(l_child,Set_of_nodes))
            return Set_of_nodes



    def seach_for_split_node(self,start_node,end_node):
        List_of_nodes=set()
        not_found=False
        split_node=None

        while not not_found:
            Need_to_look_if_end_is_right_child=True
            Need_to_look_if_start_is_leaft_child=True
            #TODO have to check if it exists or if it is root.....
            grandparent_start_node=None
            grandparent_end_node=None
            parent_start_node=start_node.get_parent()
            if not (parent_start_node.is_root()):
                grandparent_start_node=parent_start_node.get_parent()
            parent_end_node=end_node.get_parent()
            if (parent_end_node != None):
                grandparent_end_node=parent_end_node.get_parent()
            if (parent_end_node==parent_start_node and (parent_end_node != None and parent_start_node!= None)):
                #print('In first case')
                not_found=True
                split_node=parent_end_node
                List_of_nodes.add(parent_end_node)
                #If the end node is a right child we have to add all nodes left of it also to the list
                if end_node.get_is_right_child():
                    print('ENd node get is right child call all_nodes')
                    Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_right_child)
                #if the start node is a left child we have to add all nodes right of it which do not include the end node
                if start_node.get_is_left_child():
                    print('Start node get is left child call all_nodes')
                    Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_left_child)
            else:
                #Need to be in else, because it the first case occures the fourth is included
                if (grandparent_end_node==grandparent_start_node and (grandparent_end_node != None and grandparent_start_node!= None)):
                    print('In fourth case')
                    not_found=True
                    split_node=grandparent_end_node
                    List_of_nodes.add(grandparent_end_node)
                    if end_node.get_is_right_child():
                        print('End node get is right child after split node found')
                        Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                        List_of_nodes.union(Nodes_of_the_right_child)
                    if start_node.get_is_left_child():
                        print('Start node get is left child after split node found')
                        Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                        List_of_nodes.union(Nodes_of_the_left_child)
                    List_of_nodes.add(parent_end_node)
                    List_of_nodes.add(parent_start_node)


            if (parent_end_node==grandparent_start_node and (parent_end_node != None and grandparent_start_node!= None)):
                #print('In second case')
                not_found=True
                Need_to_look_if_end_is_right_child=False
                split_node=parent_end_node
                List_of_nodes.add(parent_end_node)
                List_of_nodes.add(parent_start_node)
                if start_node.get_is_left_child:
                    Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_left_child)


            if (grandparent_end_node==parent_start_node and (parent_end_node != None and parent_start_node!= None)):
                print('In third case')
                not_found=True
                Need_to_look_if_start_is_leaft_child=False

                split_node=parent_start_node
                #Need_to_look_if_end_is_right_child=False
                List_of_nodes.add(parent_start_node)
                List_of_nodes.add(parent_end_node)
                print('LIST Of NODES')
                print(List_of_nodes)
                new_gelp_list=[l.get_value() for l in List_of_nodes if l.isLeaf()]
                print('Help list')
                print(new_gelp_list)
                help_set=set(new_gelp_list)
                print(50 in help_set)

                if end_node.get_is_right_child():
                    print('in If with end node ')
                    print(end_node.get_value())
                    Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_right_child)

            #Befor resetting them have to look at their siblings according from the parent node

            if (start_node.get_is_left_child() and Need_to_look_if_start_is_leaft_child):
                print('In start node get is left child')
                Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                List_of_nodes.union(Nodes_of_the_left_child)
            if (end_node.get_is_right_child() and Need_to_look_if_end_is_right_child):
                print('In THIS LAST IF ')
                Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                List_of_nodes.union(Nodes_of_the_right_child)

            #reset start and end node

            List_of_nodes.add(start_node)
            List_of_nodes.add(end_node)
            start_node=parent_start_node
            end_node=parent_end_node
            print('BEFOR RETURN')
        return (split_node,List_of_nodes)


class BST_node:

    def __init__(self,minimum,maximum):
        self.balance=0
        self.min_coverage=minimum
        self.max_coverage=maximum
        #will later be overritten except when the node is the root
        self.parent=None
        self.is_left_child=None
        self.is_right_child=None


    def set_left_child(self,left_Node):
        self.left_child=left_Node

    def set_right_child(self,right_Node):
        self.right_child=right_Node

    def set_is_left_child(self):
        self.is_left_child=True
        self.is_right_child=False

    def set_is_right_child(self):
        self.is_left_child=False
        self.is_right_child=True

    def set_balance(self,new_balance):
        self.balance=new_balance

    def set_parent(self,node):
        self.parent= node

    def set_min_coverage(self,value):
        self.min_coverage=value

    def set_max_coverage(self,value):
        self.min_coverage=value


    def get_left_child(self):
        return self.left_child

    def get_right_child(self):
        return self.right_child

    def get_is_left_child(self):
        return self.is_left_child

    def get_is_right_child(self):
        return self.is_right_child

    def get_balance(self):
        return self.balance

    def get_coverage(self):
        return (self.min_coverage,self.max_coverage)

    def get_min_coverage(self):
        return self.min_coverage

    def get_max_coverage(self):
        return self.max_coverage

    def get_parent(self):
        return self.parent


    def isLeaf(self):
        return False

    def is_root(self):
        if (self.parent==None):
            return True
        else:
            return False



class Leaf_node:
    def __init__(self,value, index):
        '''
        :param value: Corresponds to SNP either first or last position in a read
        :param index: Index of the read in the original readset
        :return: A leaf node with at attributes value and index
        '''
        self.value=value
        adding_index=[]
        adding_index.append(index)
        self.index=adding_index
        #set parent to none will later be overwritten
        self.parent=None
        self.balance=0

    def set_sibling(self,sibling):
        #Set sibling node or nodes to the Node
        adding_sibling=[]
        adding_sibling.append(sibling)
        self.sibling= adding_sibling

    def set_parent(self, parent):
        '''Setting the node attributes
        Coverage is exactly the number of sibling the node has
        :param parent: The index of the parent node in the array
        :return
        '''
        self.balance=0
        self.parent=parent

    def set_coverage(self,coverage):
        self.coverage=coverage

    def set_is_left_child(self):
        self.is_left_child=True
        self.is_right_child=False

    def set_is_right_child(self):
        self.is_left_child=False
        self.is_right_child=True

    def set_balance(self,b):
        self.balance=b

    def get_is_left_child(self):
        return self.is_left_child

    def get_is_right_child(self):
        return self.is_right_child

    def get_parent(self):
        return self.parent

    def get_coverage(self):
        return self.coverage

    def get_balance(self):
        return self.balance

    def get_value(self):
        return self.value

    def get_sibling(self):
        #returns list of sibling, either one ore more if the position of the node occured more than once
        return self.sibling

    def add_sibling(self,value_list):
        adding_sibling= self.sibling
        for i in value_list:
            adding_sibling.append(i)
        self.sibling = adding_sibling

    def isLeaf(self):
        return True


    def get_index(self):
        return self.index

    def is_root(self):
        return False

