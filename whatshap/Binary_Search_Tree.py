import math
from whatshap._core import PyRead,PyReadSet


#Before starting the algorithm we have to build up the binary search tree.
class Binary_Search_Tree:



    def __init__(self,readset):
        node_list=self.new_building_of_list(readset)
        l_list=self.set_coverage_of_leaf_list(node_list)
        self.leaf_list=l_list
        layer_list=self.bottom_up_construction_of_BST(l_list)
        inner_nodes=[]
        complete_tree=self.build_balanced_binary_tree(layer_list,inner_nodes)
        self.complete_tree=(complete_tree[len(complete_tree)-1])



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


    #New struct of binary tree
    def new_building_of_list(self,orig_readset):
        list_of_bam=[]
        list_for_leafs=[]
        other_dic={}
        indices_of_reads = set(i for i, read in enumerate(orig_readset) if len(read) >= 2)
        for i in indices_of_reads:
            read=orig_readset[i]
            #read.getID()
            #
            #readID= read.getSourceID()
            #
            #while len(list_of_bam)< (readID+1):
            #    list_of_bam.append(0)
            #
            #    score=list_of_bam[readID]
            #    newscore=score+1
            #    list_of_bam[readID]=newscore

            delimiter_1=(read[0]).position
            delimiter_2=(read[len(read)-1]).position
            #look if the position already exists
            if delimiter_1 in other_dic and delimiter_2 in other_dic:
                f_Node=other_dic[delimiter_1]
                f_Node.add_index(i)
                e_Node=other_dic[delimiter_2]
                e_Node.add_index(i)
                e_Node.add_sibling([f_Node,delimiter_1,i])
                f_Node.add_sibling([e_Node,delimiter_2,i])


            if delimiter_2 in other_dic and delimiter_1 not in other_dic:
                e_Node=other_dic[delimiter_2]
                first_Node=Leaf_node(delimiter_1,i)
                first_Node.set_sibling([e_Node,delimiter_2,i])
                other_dic[delimiter_1]=first_Node
                e_Node.add_sibling([first_Node,delimiter_1,i])
                e_Node.add_index(i)
                list_for_leafs.append(first_Node)


            if delimiter_1 not in other_dic and delimiter_2 not in other_dic:
                first_Node=Leaf_node(delimiter_1,i)
                second_Node= Leaf_node(delimiter_2,i)
                second_Node.set_sibling([first_Node,delimiter_1,i])
                first_Node.set_sibling([second_Node,delimiter_2,i])
                other_dic[delimiter_2]=second_Node
                other_dic[delimiter_1]=first_Node
                list_for_leafs.append(first_Node)
                list_for_leafs.append(second_Node)


            if delimiter_2 not in other_dic and delimiter_1 in other_dic:
                f_Node=other_dic[delimiter_1]
                f_Node.add_index(i)
                second_Node= Leaf_node(delimiter_2,i)
                second_Node.set_sibling([f_Node,delimiter_1,i])
                other_dic[delimiter_2]=second_Node
                list_for_leafs.append(second_Node)
                f_Node.add_sibling([second_Node,delimiter_2,i])
        sorted_list=sorted(list_for_leafs,key=lambda node :node.value)
        self.bam_list=list_of_bam
        return sorted_list


    def set_coverage_of_leaf_list(self,leaf_list):
        coverage_count=0
        iter=0
        while iter !=len(leaf_list):
            val_of_iter=leaf_list[iter].get_value()
            number_ended_reads=0
            siblings_of_node=leaf_list[iter].get_sibling()
            len_siblinling=0
            while len_siblinling<len(siblings_of_node):
                (Node,value,index)=siblings_of_node[len_siblinling]
                #Found beginning of a read
                if (val_of_iter)<value:
                    coverage_count+=1
                else:
                    #Found ending of a read , need to decrease the coverage count in the end
                    number_ended_reads+=1
                len_siblinling+=1
            leaf_list[iter].set_coverage(coverage_count)
            coverage_count=coverage_count-number_ended_reads
            iter+=1
        return leaf_list




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

    def get_bam_list(self):
        return self.bam_list

    def bottom_up_construction_of_BST(self,node_list):
        re_value=1
        number_nodes=len(node_list)
        while re_value>0:
            if math.pow(2,re_value)>number_nodes:
                break
            else:
                re_value+=1
        re_value-=1
        #for_compression is the number of leaves which have to be assigned to an inner node such that the tree get a dlayer of exactly 2^x nodes
        #and then could be build as a normal balanced binary tree
        for_compression=2*int(number_nodes-math.pow(2,re_value))
        i=0
        layer_list=[]
        while for_compression!=0:
            mini=min(node_list[i].get_coverage(),node_list[i+1].get_coverage())
            maxi=max(node_list[i].get_coverage(),node_list[i+1].get_coverage())
            layer_node=BST_node(mini,maxi)
            layer_node.set_left_child(node_list[i])
            layer_node.set_right_child(node_list[i+1])
            node_list[i].set_is_left_child()
            node_list[i+1].set_is_right_child()
            node_list[i].set_parent(layer_node)
            node_list[i+1].set_parent(layer_node)
            layer_list.append(layer_node)
            i+=2
            for_compression-=2
        #adding remaining leafs to the layer list. which then should have exactly 2^re_value elements
        while i!=number_nodes:
            layer_list.append(node_list[i])
            i+=1

        #should give re_value that exceed the number of nodes
        return layer_list

    def build_balanced_binary_tree(self,layer_list,inner_node_list):
        #come to root
        if len(layer_list)!=1:
            next_layer_node_list=[]
            lengt_of_node_list=len(layer_list)
            i=0
            #could assure that because length of the node list is a 2 to power of x case
            while i !=lengt_of_node_list:
                first_node=layer_list[i]
                second_node=layer_list[i+1]
                mini=min(first_node.get_min_coverage(),second_node.get_min_coverage())
                maxi=max(first_node.get_max_coverage(),second_node.get_max_coverage())
                new_BST_node=BST_node(mini,maxi)
                new_BST_node.set_left_child(first_node)
                new_BST_node.set_right_child(second_node)
                first_node.set_parent(new_BST_node)
                second_node.set_parent(new_BST_node)
                first_node.set_is_left_child()
                second_node.set_is_right_child()
                inner_node_list.append(new_BST_node)
                next_layer_node_list.append(new_BST_node)
                i+=2
            return self.build_balanced_binary_tree(next_layer_node_list,inner_node_list)
        else:
            return inner_node_list

#former approach
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
            right_node=self.building_BST_from_leaf_list(middel+1,end,arr,node_list)
            root_node.set_right_child(right_node)
            right_node.set_parent(root_node)
            right_node.set_is_right_child()



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
    #Former Approach
#    def is_crucial(self,split_node_coverage, max_coverage,split_balance ):
#        '''
#        Depends on the split_node and the maximal coverage if the read in this split node is crucial or not
#        '''
#        d_var=False
#        (min_cov,max_cov)=split_node_coverage
#        updated_min_coverage= min_cov+split_balance
#        updated_max_coverage= max_cov+split_balance
#        print("Updatete min coverage %d " %updated_min_coverage)
#        print("Updated max coverage %d " %updated_max_coverage)
#        print("Look if updated min coverage below or equal %d " %(math.floor(max_coverage/2)))
#        if updated_min_coverage<= math.floor(max_coverage/2):
#            d_var=True
#        if updated_max_coverage>max_cov:
#            d_var=False
#
#        return d_var    #
    #adapted version to the idea of really getting thhe coverage of the range query and not of the subtree of the split node
    #thereforegetting the coverage with in the criterion
    def is_crucial_former_approach(self,range_coverage, max_coverage,split_balance ):
        '''
        Depends on the split_node and the maximal coverage if the read in this split node is crucial or not
        '''
        d_var=False
        (min_cov,max_cov)=range_coverage
        updated_min_coverage= min_cov+split_balance
        updated_max_coverage= max_cov+split_balance
        #print("Updatete min coverage %d " %updated_min_coverage)
        #print("Updated max coverage %d " %updated_max_coverage)
        #print("Look if updated min coverage below or equal %d " %(math.floor(max_coverage/2)))
        if updated_min_coverage<= math.floor(max_coverage/2):
            d_var=True
        if updated_max_coverage>max_cov:
            d_var=False

        return d_var
#Idea that for the crucail idea the coverage should not be updatet, because otherwise to many reads will become crucial
    def is_crucial(self,range_coverage, max_coverage,split_balance ):
        '''
        Depends on the split_node and the maximal coverage if the read in this split node is crucial or not
        '''
        d_var=False
        (min_cov,max_cov)=range_coverage
        #updated_min_coverage= min_cov+split_balance
        #updated_max_coverage= max_cov+split_balance
        #print("Updatete min coverage %d " %updated_min_coverage)
        #print("Updated max coverage %d " %updated_max_coverage)
        #print("Look if updated min coverage below or equal %d " %(math.floor(max_coverage/2)))
        if min_cov<= math.floor(max_coverage/2):
            d_var=True
        #if updated_max_coverage>max_cov:
        #    d_var=False

        return d_var






    def get__all_nodes(self,node,Set_of_nodes):
        '''
        :param node: Start node of the search, could also be seen as the root
        :return: returns of a node a list of all nodes in this tree
        '''
        n_coverag=node.get_coverage()
        if node.isLeaf():
            Set_of_nodes.add(node)
            return Set_of_nodes
        else:
            Set_of_nodes.add(node)
            r_child=node.get_right_child()
            l_child=node.get_left_child()
            Set_of_nodes.union(self.get__all_nodes(r_child,Set_of_nodes))
            Set_of_nodes.union(self.get__all_nodes(l_child,Set_of_nodes))
            return Set_of_nodes


#Former
    def seach_for_split_node_former(self,start_node,end_node):
        #for every combination of leaf and end node, so for every range, keep track of the coverage
        mini_coverage=min(start_node.get_min_coverage(),end_node.get_min_coverage())
        maxi_coverage=max(start_node.get_max_coverage(),end_node.get_max_coverage())
        #print('At the start of search split node maxi %d '%maxi_coverage)
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
            parent_end_max_cov=parent_end_node.get_max_coverage()
            parent_end_min_cov=parent_end_node.get_min_coverage()
            parent_start_max_cov=parent_start_node.get_max_coverage()
            parent_start_min_cov=parent_start_node.get_min_coverage()

            if (parent_end_node==parent_start_node and (parent_end_node != None and parent_start_node!= None)):
                #print('In first case')
                #print('Nothing to do here for the coverage in the range ')
                #print('Min coverage %d '%mini_coverage)
                #print('Max coverage %d '%maxi_coverage)

                not_found=True
                split_node=parent_end_node
                List_of_nodes.add(parent_end_node)
                #If the end node is a right child we have to add all nodes left of it also to the list
                if end_node.get_is_right_child():
                    #print('ENd node get is right child call all_nodes')
                    Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_right_child)
                #if the start node is a left child we have to add all nodes right of it which do not include the end node
                if start_node.get_is_left_child():
                    #print('Start node get is left child call all_nodes')
                    Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_left_child)
            else:
                #Need to be in else, because it the first case occures the fourth is included
                if (grandparent_end_node==grandparent_start_node and (grandparent_end_node != None and grandparent_start_node!= None)):
                    #print('In fourth case')
                    not_found=True
                    split_node=grandparent_end_node
                    List_of_nodes.add(grandparent_end_node)

                    if end_node.get_is_right_child():
                        #print('End node get is right child after split node found')
                        Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                        List_of_nodes.union(Nodes_of_the_right_child)
                        #Need to include the coverage pf the parent_end_node subtree
                        mini_coverage=min(parent_end_min_cov,mini_coverage)
                        maxi_coverage=max(parent_end_max_cov,maxi_coverage)

                    if start_node.get_is_left_child():
                        Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                        List_of_nodes.union(Nodes_of_the_left_child)
                        mini_coverage=min(parent_start_min_cov,mini_coverage)
                        maxi_coverage=max(parent_start_max_cov,maxi_coverage)
                    List_of_nodes.add(parent_end_node)
                    List_of_nodes.add(parent_start_node)


            if (parent_end_node==grandparent_start_node and (parent_end_node != None and grandparent_start_node!= None)):
                #print('In second case')
                not_found=True
                Need_to_look_if_end_is_right_child=False
                split_node=parent_end_node
                List_of_nodes.add(parent_end_node)
                List_of_nodes.add(parent_start_node)

                if start_node.get_is_left_child():
                    Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_left_child)
                    mini_coverage=min(parent_start_min_cov,mini_coverage)
                    maxi_coverage=max(parent_start_max_cov,maxi_coverage)


            if (grandparent_end_node==parent_start_node and (parent_end_node != None and parent_start_node!= None)):
                #print('In third case')
                not_found=True
                Need_to_look_if_start_is_leaft_child=False

                split_node=parent_start_node
                #Need_to_look_if_end_is_right_child=False
                List_of_nodes.add(parent_start_node)
                List_of_nodes.add(parent_end_node)

                if end_node.get_is_right_child():
                    #print('in If with end node ')
                    #print(end_node.get_value())
                    Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_right_child)
                    mini_coverage=min(parent_end_min_cov,mini_coverage)
                    maxi_coverage=max(parent_end_max_cov,maxi_coverage)

            #Befor resetting them have to look at their siblings according from the parent node

            if (start_node.get_is_left_child() and Need_to_look_if_start_is_leaft_child):
                #print('In start node get is left child')
                Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                List_of_nodes.union(Nodes_of_the_left_child)
                maxi_coverage=max(maxi_coverage,parent_start_max_cov)
                mini_coverage=min(mini_coverage,parent_start_min_cov)
            if (end_node.get_is_right_child() and Need_to_look_if_end_is_right_child):
                #print('In THIS LAST IF ')
                Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                List_of_nodes.union(Nodes_of_the_right_child)
                maxi_coverage=max(maxi_coverage,parent_end_max_cov)
                mini_coverage=min(mini_coverage,parent_end_min_cov)
            #reset start and end node

            List_of_nodes.add(start_node)
            List_of_nodes.add(end_node)
            start_node=parent_start_node
            end_node=parent_end_node
        coverage_of_range=(mini_coverage,maxi_coverage)
        return (split_node,List_of_nodes,coverage_of_range)


    def seach_for_split_node(self,start_node,end_node):
        #for every combination of leaf and end node, so for every range, keep track of the coverage
        mini_coverage=min(start_node.get_min_coverage(),end_node.get_min_coverage())
        maxi_coverage=max(start_node.get_max_coverage(),end_node.get_max_coverage())
        #print('At the start of search split node maxi %d '%maxi_coverage)
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
            #parent_end_max_cov=parent_end_node.get_max_coverage()
            #parent_end_min_cov=parent_end_node.get_min_coverage()
            #parent_start_max_cov=parent_start_node.get_max_coverage()
            #parent_start_min_cov=parent_start_node.get_min_coverage()
            if start_node.get_is_left_child():
                coverage_start_node_parent_right_child_min=parent_start_node.get_right_child().get_min_coverage()
                coverage_start_node_parent_right_child_max=parent_start_node.get_right_child().get_max_coverage()
            if end_node.get_is_right_child():
                coverage_end_node_parent_left_child_min=parent_end_node.get_left_child().get_min_coverage()
                coverage_end_node_parent_left_child_max=parent_end_node.get_left_child().get_max_coverage()

            if (parent_end_node==parent_start_node and (parent_end_node != None and parent_start_node!= None)):
                #print('In first case')
                #print('Nothing to do here for the coverage in the range ')
                #print('Min coverage %d '%mini_coverage)
                #print('Max coverage %d '%maxi_coverage)

                not_found=True
                split_node=parent_end_node
                List_of_nodes.add(parent_end_node)
                #If the end node is a right child we have to add all nodes left of it also to the list
                if end_node.get_is_right_child():
                    #print('ENd node get is right child call all_nodes')
                    Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_right_child)
                #if the start node is a left child we have to add all nodes right of it which do not include the end node
                if start_node.get_is_left_child():
                    #print('Start node get is left child call all_nodes')
                    Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_left_child)
            else:
                #Need to be in else, because it the first case occures the fourth is included
                if (grandparent_end_node==grandparent_start_node and (grandparent_end_node != None and grandparent_start_node!= None)):
                    #print('In fourth case')
                    not_found=True
                    split_node=grandparent_end_node
                    List_of_nodes.add(grandparent_end_node)

                    if end_node.get_is_right_child():
                        #print('End node get is right child after split node found')
                        Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                        List_of_nodes.union(Nodes_of_the_right_child)
                        #Need to include the coverage pf the parent_end_node subtree
                        mini_coverage=min(coverage_end_node_parent_left_child_min,mini_coverage)
                        maxi_coverage=max(coverage_end_node_parent_left_child_max,maxi_coverage)

                    if start_node.get_is_left_child():
                        Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                        List_of_nodes.union(Nodes_of_the_left_child)
                        mini_coverage=min(coverage_start_node_parent_right_child_min,mini_coverage)
                        maxi_coverage=max(coverage_start_node_parent_right_child_max,maxi_coverage)
                    List_of_nodes.add(parent_end_node)
                    List_of_nodes.add(parent_start_node)


            if (parent_end_node==grandparent_start_node and (parent_end_node != None and grandparent_start_node!= None)):
                #print('In second case')
                not_found=True
                Need_to_look_if_end_is_right_child=False
                split_node=parent_end_node
                List_of_nodes.add(parent_end_node)
                List_of_nodes.add(parent_start_node)

                if start_node.get_is_left_child():
                    Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_left_child)
                    mini_coverage=min(coverage_start_node_parent_right_child_min,mini_coverage)
                    maxi_coverage=max(coverage_start_node_parent_right_child_max,maxi_coverage)


            if (grandparent_end_node==parent_start_node and (parent_end_node != None and parent_start_node!= None)):
                #print('In third case')
                not_found=True
                Need_to_look_if_start_is_leaft_child=False

                split_node=parent_start_node
                #Need_to_look_if_end_is_right_child=False
                List_of_nodes.add(parent_start_node)
                List_of_nodes.add(parent_end_node)

                if end_node.get_is_right_child():
                    #print('in If with end node ')
                    #print(end_node.get_value())
                    Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                    List_of_nodes.union(Nodes_of_the_right_child)
                    mini_coverage=min(coverage_end_node_parent_left_child_min,mini_coverage)
                    maxi_coverage=max(coverage_end_node_parent_left_child_max,maxi_coverage)

            #Befor resetting them have to look at their siblings according from the parent node

            if (start_node.get_is_left_child() and Need_to_look_if_start_is_leaft_child):
                #print('In start node get is left child')
                Nodes_of_the_left_child=self.get__all_nodes(parent_start_node.get_right_child(),List_of_nodes)
                List_of_nodes.union(Nodes_of_the_left_child)
                #coverage_start_node_parent_right_child_min=parent_start_node.get_right_child().get_min_coverage()
                #coverage_start_node_parent_right_child_max=parent_start_node.get_right_child().get_max_coverage()
                maxi_coverage=max(maxi_coverage,coverage_start_node_parent_right_child_max)
                mini_coverage=min(mini_coverage,coverage_start_node_parent_right_child_min)
            if (end_node.get_is_right_child() and Need_to_look_if_end_is_right_child):
                #print('In THIS LAST IF ')
                Nodes_of_the_right_child=self.get__all_nodes(parent_end_node.get_left_child(),List_of_nodes)
                List_of_nodes.union(Nodes_of_the_right_child)
                #coverage_end_node_parent_left_child_min=parent_start_node.get_left_child().get_min_coverage()
                #coverage_end_node_parent_left_child_max=parent_start_node.get_left_child().get_max_coverage()
                maxi_coverage=max(maxi_coverage,coverage_end_node_parent_left_child_max)
                mini_coverage=min(mini_coverage,coverage_end_node_parent_left_child_min)
            #reset start and end node

            List_of_nodes.add(start_node)
            List_of_nodes.add(end_node)
            start_node=parent_start_node
            end_node=parent_end_node
        coverage_of_range=(mini_coverage,maxi_coverage)
        return (split_node,List_of_nodes,coverage_of_range)






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

    def add_index(self,i):
        former_index=self.index
        former_index.append(i)
        self.index=former_index



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

    def get_min_coverage(self):
        return self.coverage

    def get_max_coverage(self):
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
        #for i in value_list:
        #   adding_sibling.append(i)
        adding_sibling.append(value_list)
        self.sibling = adding_sibling

    def isLeaf(self):
        return True


    def get_index(self):
        return self.index

    def __len__(self):
        #Length in this conntext means the number of reads starting or ending in this delimiter
        return len(self.index)

    def is_root(self):
        return False

