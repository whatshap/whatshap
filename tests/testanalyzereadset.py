from whatshap.scripts.whatshap import union_sets, analyze_readset, check_for_connectivity
from phasingutils import string_to_readset

def test_union():
     firstset=set(['a','b','c','d'])
     secondset=set(['c','d','e','f'])
     connectivity_factor=2
     connection_list=[]
     connection_list.append(firstset)
     connection_list.append(secondset)
     new_connection_list=union_sets(1,connection_list,connectivity_factor)
     assert len(new_connection_list[0]) == 1
     finalset=set(['a','b','c','d','e','f'])
     assert(len (finalset.difference(new_connection_list[0][0])) == 0)

def test_union2():
     firstset=set(['a','b','c','d'])
     secondset=set(['c','d','e','f'])
     thirdset=set(['z'])
     fourthset=set(['z','t'])
     connectivity_factor=2
     connection_list=[]
     connection_list.append(firstset)
     connection_list.append(secondset)
     (new_connection_list,index)=union_sets(1,connection_list,connectivity_factor)
     new_connection_list.append(thirdset)
     new_connection_list.append(fourthset)
     (connection_list,index)=union_sets(2,connection_list,connectivity_factor)
     assert len(connection_list) == 3


def test_union3():
     firstset=set(['a','b','c','d'])
     secondset=set(['c','d','e','f'])
     thirdset=set(['z'])
     fourthset=set(['z','t'])
     connectivity_factor=1
     connection_list=[]
     connection_list.append(firstset)
     connection_list.append(secondset)
     (new_connection_list,index)=union_sets(1,connection_list,connectivity_factor)
     new_connection_list.append(thirdset)
     new_connection_list.append(fourthset)
     (connection_list,index)=union_sets(2,connection_list,connectivity_factor)
     assert len(connection_list) == 2

def test_union4():
     firstset=set(['a','b','c','d'])
     secondset=set(['c','d','e','f'])
     thirdset=set(['b','c','e','f','a'])
     fourthset=set(['z','t'])
     connectivity_factor=3
     connection_list=[]
     connection_list.append(firstset)
     connection_list.append(secondset)
     (new_connection_list,index)=union_sets(1,connection_list,connectivity_factor)
     new_connection_list.append(thirdset)
     new_connection_list.append(fourthset)
     (connection_list,index)=union_sets(2,connection_list,connectivity_factor)
     assert len(connection_list) == 2


def test_recursion():
    firstset=set([10,20])
    secondset=set([30,40])
    thirdset=set([50,60])
    fourthset=set([70,50,40])
    connectivity_factor=1
    connection_list=[]
    connection_list.append(firstset)
    connection_list.append(secondset)
    connection_list.append(thirdset)
    connection_list.append(fourthset)
    connection_list=union_sets(3,connection_list,connectivity_factor)
    assert len(connection_list[0]) == 2


def test_recursion2():
    firstset=set([10,20,70])
    secondset=set([30,40])
    thirdset=set([50,60])
    fourthset=set([70,80,90])
    connectivity_factor=1
    connection_list=[]
    connection_list.append(firstset)
    connection_list.append(secondset)
    connection_list.append(thirdset)
    connection_list.append(fourthset)
    connection_list=union_sets(3,connection_list,connectivity_factor)
    assert len(connection_list[0]) == 3



def test_connectivity_1_analysis():
     reads = string_to_readset("""
 	  11
 	   00
 	   111
 	    000
 	    1111
 	   0000
 	     11
 	      00
 	""")
     #connectivity =1 as second factor
     analyze_readset(reads, 1, 1, 0)
     read_positions_r1=[10,20]
     read_positions_r2=[20,30]
     read_positions_r3=[20,30,40]
     read_positions_r4=[30,40,50]
     read_positions_r5=[30,40,50,60]
     read_positions_r6=[20,30,40,50]
     read_positions_r7=[40,50]
     read_positions_r8=[50,60]
     List_of_connections= []
     connectivity=1
     List_of_connections=check_for_connectivity(set(read_positions_r1), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==2
     List_of_connections=check_for_connectivity(set(read_positions_r2), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==3
     List_of_connections=check_for_connectivity(set(read_positions_r3), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==4
     List_of_connections=check_for_connectivity(set(read_positions_r4), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==5
     List_of_connections=check_for_connectivity(set(read_positions_r5), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==6
     List_of_connections=check_for_connectivity(set(read_positions_r6), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==6
     List_of_connections=check_for_connectivity(set(read_positions_r7), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==6
     List_of_connections=check_for_connectivity(set(read_positions_r8), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==6


def test_connectivity_2_analysis():
     reads = string_to_readset("""
 	  11
 	   00
 	   111
 	    000
 	    1111
 	   0000
 	     11
 	      00
 	""")
     #connectivity =1 as second factor
     analyze_readset(reads, 1, 2, 0)
     read_positions_r1=[10,20]
     read_positions_r2=[20,30]
     read_positions_r3=[20,30,40]
     read_positions_r4=[30,40,50]
     read_positions_r5=[30,40,50,60]
     read_positions_r6=[20,30,40,50]
     read_positions_r7=[40,50]
     read_positions_r8=[50,60]
     List_of_connections= []
     connectivity=2
     List_of_connections=check_for_connectivity(set(read_positions_r1), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==2
     List_of_connections=check_for_connectivity(set(read_positions_r2), List_of_connections, connectivity)
     assert len(List_of_connections)==2
     assert len(List_of_connections[0])==2
     List_of_connections=check_for_connectivity(set(read_positions_r3), List_of_connections, connectivity)
     assert len(List_of_connections)==2
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==3
     List_of_connections=check_for_connectivity(set(read_positions_r4), List_of_connections, connectivity)
     assert len(List_of_connections)==2
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==4
     List_of_connections=check_for_connectivity(set(read_positions_r5), List_of_connections, connectivity)
     assert len(List_of_connections)==2
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==5
     List_of_connections=check_for_connectivity(set(read_positions_r6), List_of_connections, connectivity)
     assert len(List_of_connections)==2
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==5
     List_of_connections=check_for_connectivity(set(read_positions_r7), List_of_connections, connectivity)
     assert len(List_of_connections)==2
     assert len(List_of_connections[1])==5
     List_of_connections=check_for_connectivity(set(read_positions_r8), List_of_connections, connectivity)
     assert len(List_of_connections)==2
     assert len(List_of_connections[1])==5

def test_connectivity_3_analysis():
     reads = string_to_readset("""
 	  11
 	   00
 	   111
 	    000
 	    1111
 	   0000
 	     11
 	      00
 	""")
     #connectivity =1 as second factor
     analyze_readset(reads, 1, 3, 0)
     read_positions_r1=[10,20]
     read_positions_r2=[20,30]
     read_positions_r3=[20,30,40]
     read_positions_r4=[30,40,50]
     read_positions_r5=[30,40,50,60]
     read_positions_r6=[20,30,40,50]
     read_positions_r7=[40,50]
     read_positions_r8=[50,60]
     List_of_connections= []
     connectivity=3
     List_of_connections=check_for_connectivity(set(read_positions_r1), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==2
     List_of_connections=check_for_connectivity(set(read_positions_r2), List_of_connections, connectivity)
     assert len(List_of_connections)==2
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==2
     List_of_connections=check_for_connectivity(set(read_positions_r3), List_of_connections, connectivity)
     assert len(List_of_connections)==3
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==2
     assert len(List_of_connections[2])==3
     List_of_connections=check_for_connectivity(set(read_positions_r4), List_of_connections, connectivity)
     assert len(List_of_connections)==4
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==2
     assert len(List_of_connections[2])==3
     assert len(List_of_connections[3])==3
     List_of_connections=check_for_connectivity(set(read_positions_r5), List_of_connections, connectivity)
     assert len(List_of_connections)==4
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==2
     assert len(List_of_connections[2])==3
     assert len(List_of_connections[3])==4
     List_of_connections=check_for_connectivity(set(read_positions_r6), List_of_connections, connectivity)
     assert len(List_of_connections)==3
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==2
     assert len(List_of_connections[2])==5
     List_of_connections=check_for_connectivity(set(read_positions_r7), List_of_connections, connectivity)
     assert len(List_of_connections)==4
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==2
     assert len(List_of_connections[2])==5
     assert len(List_of_connections[3])==2
     List_of_connections=check_for_connectivity(set(read_positions_r8), List_of_connections, connectivity)
     assert len(List_of_connections)==5
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==2
     assert len(List_of_connections[2])==5
     assert len(List_of_connections[3])==2
     assert len(List_of_connections[4])==2



def test_connectivity_paired_end_1_analysis():
     reads = string_to_readset("""
 	  1  1
 	  00
 	  0   1
 	  10  1
 	  1   1
 	    11
 	  0   1
 	  1    1
 	""")
     #connectivity =1 as second factor
     analyze_readset(reads, 1, 1, 0)
     analyze_readset(reads, 1, 2, 0)
     analyze_readset(reads, 1, 3, 0)
     read_positions_r1=[10,40]
     read_positions_r2=[10,20]
     read_positions_r3=[10,50]
     read_positions_r4=[10,20,50]
     read_positions_r5=[10,50]
     read_positions_r6=[30,40]
     read_positions_r7=[10,50]
     read_positions_r8=[10,60]
     List_of_connections= []
     connectivity=1
     List_of_connections=check_for_connectivity(set(read_positions_r1), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==2
     List_of_connections=check_for_connectivity(set(read_positions_r2), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==3
     List_of_connections=check_for_connectivity(set(read_positions_r3), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==4
     List_of_connections=check_for_connectivity(set(read_positions_r4), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==4
     List_of_connections=check_for_connectivity(set(read_positions_r5), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==4
     List_of_connections=check_for_connectivity(set(read_positions_r6), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==5
     List_of_connections=check_for_connectivity(set(read_positions_r7), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==5
     List_of_connections=check_for_connectivity(set(read_positions_r8), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==6



def test_connectivity_paired_end2_analysis():
     reads = string_to_readset("""
 	  1  1
 	  00
 	  0   1
 	  10  1
 	  1   1
 	    11
 	  0   1
 	  1    1
 	""")
     #connectivity =1 as second factor
     read_positions_r1=[10,40]
     read_positions_r2=[10,20]
     read_positions_r3=[10,50]
     read_positions_r4=[10,20,50]
     read_positions_r5=[10,50]
     read_positions_r6=[30,40]
     read_positions_r7=[10,50]
     read_positions_r8=[10,60]
     List_of_connections= []
     connectivity=2
     List_of_connections=check_for_connectivity(set(read_positions_r1), List_of_connections, connectivity)
     assert len(List_of_connections)==1
     assert len(List_of_connections[0])==2
     List_of_connections=check_for_connectivity(set(read_positions_r2), List_of_connections, connectivity)
     assert len(List_of_connections)==2
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==2
     List_of_connections=check_for_connectivity(set(read_positions_r3), List_of_connections, connectivity)
     assert len(List_of_connections)==3
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==2
     assert len(List_of_connections[2])==2
     List_of_connections=check_for_connectivity(set(read_positions_r4), List_of_connections, connectivity)
     assert len(List_of_connections)==2
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==3
     List_of_connections=check_for_connectivity(set(read_positions_r5), List_of_connections, connectivity)
     assert len(List_of_connections)==2
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==3
     List_of_connections=check_for_connectivity(set(read_positions_r6), List_of_connections, connectivity)
     assert len(List_of_connections)==3
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==3
     assert len(List_of_connections[2])==2
     List_of_connections=check_for_connectivity(set(read_positions_r7), List_of_connections, connectivity)
     assert len(List_of_connections)==3
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==3
     assert len(List_of_connections[2])==2
     List_of_connections=check_for_connectivity(set(read_positions_r8), List_of_connections, connectivity)
     assert len(List_of_connections)==4
     assert len(List_of_connections[0])==2
     assert len(List_of_connections[1])==3
     assert len(List_of_connections[2])==2
     assert len(List_of_connections[3])==2

def test_connection_by_deleting_elements():
     reads=string_to_readset("""
     11
      011
     0  0
     """)
     read_positions_r1=[10,20]
     read_positions_r2=[20,30,40]
     read_positions_r3=[10,40]
     read_positions_r4=[10,30,40,50]
     connectivity=2
     List_of_connection=[]
     List_of_connection=check_for_connectivity(set(read_positions_r1),List_of_connection,connectivity)
     assert len(List_of_connection)==1
     assert len(List_of_connection[0])==2
     List_of_connection=check_for_connectivity(set(read_positions_r2),List_of_connection,connectivity)
     assert len(List_of_connection)==2
     assert len(List_of_connection[0])==2
     assert len(List_of_connection[1])==3
     List_of_connection=check_for_connectivity(set(read_positions_r3),List_of_connection,connectivity)
     assert len(List_of_connection)==3
     assert len(List_of_connection[0])==2
     assert len(List_of_connection[1])==3
     assert len(List_of_connection[2])==2
     List_of_connection=check_for_connectivity(set(read_positions_r4),List_of_connection,connectivity)
     assert len(List_of_connection)==1
     assert len(List_of_connection[0])==5
