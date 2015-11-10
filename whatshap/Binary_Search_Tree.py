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
        print('Split node coverage ')
        print(split_node_coverage)
        updated_min_coverage= min_cov+split_balance
        updated_max_coverage= max_cov+split_balance
        print('Updated min coverage %d' %updated_min_coverage)
        print('Updateded max coverage %d' %updated_max_coverage)
        if updated_min_coverage<= math.floor(max_coverage/2):
            print('In updated min coverage less than max')
            d_var=True
        if updated_max_coverage>max_cov:
            print('In max_coverage over max coverage')
            d_var=False
        return d_var
        #start_position = read.getposition()
        #end_positions = read.getposition()
        #TODO Look again if it is really ceil
        #if get_min_cov(start_position, end_positions) <= math.ceil(
        #                get_max_cov(start_position, end_positions, self.max_coverage) / 2):
        #    self.pruned_readset.add(read)
        #return False










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
            #print('In for loop %d' %actual_value)
            #print(actual_node.get_coverage())
            for j in Leaf_list:
                #print('In for loop over leaf list')
                if j.get_value()==actual_value:
                    #print('Found node with same value')
                    actual_node.set_balance(j.get_balance())
                    actual_node.set_parent(j.get_parent())
                    if j.get_is_left_child():
                        #print('IS LEFT CHILD')
                        actual_node.set_is_left_child()
                    else:
                        #print('Is RIGTH CHILD ')
                        actual_node.set_is_right_child()
                    break

            #print('I in sibling list %d' %actual_value)
                #Not sure if this could work  or even if just needed, because these nodes do not exist in the leaf_list
                #sibling_list[i+1].set_sibling(actual_sibling)
        #for sib in sibling_list:
        #    print('In second for ')
        #    print(sib.get_is_left_child())





    def get__all_nodes(self,node,Set_of_nodes):
        '''
        :param node: Start node of the search, could also be seen as the root
        :return: returns of a node a list of all nodes in this tree
        '''
        print('Called get_all_nodes')
        print('With node %d' %node.get_coverage())
        while not node.isLeaf():
            Set_of_nodes.add(node)
           # print('In while appended')
            r_child=node.get_right_child()
            l_child=node.get_left_child()
            Set_of_nodes.union(self.get__all_nodes(r_child,Set_of_nodes))
            Set_of_nodes.add(self.get__all_nodes(l_child,Set_of_nodes))
        Set_of_nodes.add(node)
        #print('Appended node')
        return Set_of_nodes


    def seach_for_split_node(self,start_node,end_node):
        #print('Search for split node ')
        List_of_nodes=set()
        #List_of_nodes.add(start_node)
        #List_of_nodes.add(end_node)
        not_found=False
        split_node=None

        #print('Start node, coverage')
        #print(start_node.get_coverage())
        #print('End node, coverage')
        #print(end_node.get_coverage())

        while not not_found:
            print('In While')
            Need_to_look_if_end_is_right_child=True
            #TODO have to check if it exists or if it is root.....
            grandparent_start_node=None
            grandparent_end_node=None
            parent_start_node=start_node.get_parent()
            print('parent_start_node')
            print(parent_start_node.get_coverage())
            if (parent_start_node != None):
                grandparent_start_node=parent_start_node.get_parent()
            print('grandparent_start_node')
            print(grandparent_start_node.get_coverage())
            parent_end_node=end_node.get_parent()
            print('parent_end_node ')
            print(parent_end_node.get_coverage())
            if (parent_end_node != None):
                grandparent_end_node=parent_end_node.get_parent()
            print('grandparent_end_node')
            print(grandparent_end_node.get_coverage())
            if (parent_end_node==parent_start_node and (parent_end_node != None and parent_start_node!= None)):
                print('In first case')
                not_found=True
                split_node=parent_end_node
                List_of_nodes.add(parent_end_node)
                #If the end node is a right child we have to add all nodes left of it also to the list
                if end_node.get_is_right_child():
                    #e_list=[]
                    Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_right_child)
                #if the start node is a left child we have to add all nodes right of it which do not include the end node
                if start_node.get_is_left_child():
                    #e_list=[]
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
                        Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                        List_of_nodes.union(Nodes_of_the_right_child)
                    if start_node.get_is_left_child():
                        Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                        List_of_nodes.union(Nodes_of_the_left_child)
                    List_of_nodes.add(parent_end_node)
                    List_of_nodes.add(parent_start_node)
                    Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                    Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)



            if (parent_end_node==grandparent_start_node and (parent_end_node != None and grandparent_start_node!= None)):
                print('In second case')
                not_found=True
                Need_to_look_if_end_is_right_child=False
                split_node=parent_end_node
                List_of_nodes.add(parent_end_node)
                List_of_nodes.add(parent_start_node)
                if start_node.get_is_left_child:
                    Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_left_child)


            if (grandparent_end_node==parent_start_node and (parent_end_node != None and parent_start_node!= None)):
                print('In thrid case')
                not_found=True
                split_node=parent_start_node
                List_of_nodes.add(parent_start_node)
                List_of_nodes.add(parent_end_node)
                if end_node.get_is_right_child():
                   Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                   List_of_nodes.union(Nodes_of_the_right_child)

            #Befor resetting them have to look at their siblings accoring from the parent node

            if start_node.get_is_left_child():
                Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                List_of_nodes.union(Nodes_of_the_left_child)
            print('Look if end node is None')
            print(end_node==None)
            print(end_node.isLeaf())
            if (end_node.get_is_right_child() and Need_to_look_if_end_is_right_child):
                print('Called in the end node')
                print(end_node.isLeaf())
                Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                List_of_nodes.union(Nodes_of_the_right_child)


            #reset start and end node

            List_of_nodes.add(start_node)
            List_of_nodes.add(end_node)
            start_node=parent_start_node
            end_node=parent_end_node


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




################################################################################################
#Methods maybe later needed

'''
    def search_for_node_in_subtree(self,end_node,parent_node,Found_node):
        #returns a List_of nodes of the subtree, either the whole subtree or only a part, furthermore a boolean indicating
        #if the node was found. and at last the start point of the search
        #Stores the nodes whose balance should be updated if the read is rejected
        List_of_nodes=[]
        x=parent_node
        while not x.isLeaf():
            r_of_x=x.get_right_child()
            x=r_of_x

        if (x.get_value()==end_node.get_value()):
            List_of_nodes_in_subtree =get__all_nodes(parent_node,List_of_nodes)
            List_of_nodes.append(List_of_nodes_in_subtree)
            Found_node=True
            #In this case the parent_node is the seached split node
            return (List_of_nodes_in_subtree,Found_node,parent_node)

        #value of the rightmost child in the subtree has lower value then the searched node the end_node is not there
        if (x.get_value()<end_node.get_value()):
            List_of_nodes_in_subtree=get__all_nodes(parent_node,List_of_nodes)
            List_of_nodes.append(List_of_nodes_in_subtree)
            Found_node=False
            return (List_of_nodes_in_subtree,Found_node,parent_node)
        #end node has to be in the subtree of the parent node somewhere
        else:
            x_parent=x.get_parent()
            l_of_p_of_x=x_parent.get_left_child()
            #left_child_is aLeaft and has the value we have found the end node
            #Should not do direct comparisoon because l_of_p_of_x could also be an inner node
            if (l_of_p_of_x.isLeaf() and l_of_p_of_x.get_value()==end_node.get_value()):
                List_of_nodes.add(end_node,l_of_p_of_x,x_parent)


#COULD MAYBE USED THERE WHEN WE NEED THE BALANCE
    def get_all_nodes_though_split_list(self,List_to_change,value_of_sibling):
        length_of_list=len(List_to_change)
        while (length_of_list!=0):
            node=List_to_change.pop()
            node.set_balance(node.get_balance()-1)
            if node.isLeaf():
                r_child=node.get_right_child()
                l_child=node.get_left_child()
                #when we have a leaf with an bigger value as the end point the balance there has not to be changed.
                if r_child.isLeaf() and r_child.get_value():
                    continue
                l_child.set_balance(node.get_balance()-1)
                r_child.set_balance(node.get_balance()-1)
            #print('Node :')
            #print(node.isLeaf())
            #print(node.get_coverage())


#TODO does nothing in the moment
    def balance_control(self,split_node_of_read,List_to_change,l_node,value_of_sibling):

        Should update the balance in the tree, because the read starting at l_node and ending on value of sibling is pruned
        out of the original readset
        :param split_node_of_read: Found split node, so we have to check from this point on up an d down
        :param List_to_change: List which stores which inner BST nodes are followed to get to the split node
        :param l_node: start node of the removed read
        :param value_of_sibling:  value of the end node of the read, BUT NOT THE NODE ITSELF
        :return: Nothing a changed tree with changed balance

        #List_of_nodes=self.get_all_nodes_in_the_subtree_without_Those_not_connected_to_sibling(split_node_of_read,value_of_sibling,List_to_change)
        self.get_all_nodes_though_split_list(List_to_change)
        return None


    def get_all_nodes_in_the_subtree_without_Those_not_connected_to_sibling(self,split_node_of_read,value_of_sibling,List_to_change):
        print('One_call_of_the_algorithm')
        balance=split_node_of_read.get_balance()
        split_node_of_read.set_balance(balance-1)
        Right_child_Leaf=False
        Left_child_Leaf=False
       # while not(Right_child_Leaf and Left_child_Leaf):
        #    l_child=split_node_of_read.get_left_child()
         #   r_child=split_node_of_read.get_right_child()
          #  if l_child.isLeaf():
           #     Left_child_Leaf=True
           #     continue
           # if r_child.isLeaf():
           #     Right_child_Leaf=True
           #     continue
            #Right_child_Leaf=True
            #Left_child_Leaf=True


        List_of_nodes_to_change=[]
        for i in List_to_change:
            print('For each node coverage, parent, rchild, and left child')
            print(i.get_coverage())
            print(i.get_parent().get_coverage())
            print(i.get_right_child().get_coverage())
            print(i.get_left_child().get_coverage())
            if i.get_right_child().isLeaf():
                print(i.get_right_child().get_value())
            if i.get_left_child().isLeaf():
                print(i.get_left_child().get_value())
        return None


'''