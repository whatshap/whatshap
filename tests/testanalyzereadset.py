from whatshap.scripts.whatshap import union_sets


def test_union():
    firstset=set(['a','b','c','d'])
    secondset=set(['c','d','e','f'])
    connectivity_factor=2
    connection_list=[]
    connection_list.append(firstset)
    connection_list.append(secondset)
    new_connection_list=union_sets(1,connection_list,connectivity_factor)
    print(new_connection_list)
    print(connection_list)
    assert len(new_connection_list) == 1
    assert(len (set(['a','b','c','d','e','f']).difference(new_connection_list.pop())) == 0)

def test_union2():
    firstset=set(['a','b','c','d'])
    secondset=set(['c','d','e','f'])
    thirdset=set(['z'])
    fourthset=set(['z','t'])
    connectivity_factor=2
    connection_list=[]
    connection_list.append(firstset)
    connection_list.append(secondset)
    new_connection_list=union_sets(1,connection_list,connectivity_factor)
    new_connection_list.append(thirdset)
    new_connection_list.append(fourthset)
    connection_list=union_sets(2,connection_list,connectivity_factor)
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
    new_connection_list=union_sets(1,connection_list,connectivity_factor)
    new_connection_list.append(thirdset)
    new_connection_list.append(fourthset)
    connection_list=union_sets(2,connection_list,connectivity_factor)
    assert len(connection_list) == 2

