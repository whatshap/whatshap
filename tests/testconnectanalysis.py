from phasingutils import string_to_readset
from whatshap.core import readselection
from whatshap.connect_analysis import  Element,read_positions_graph,Connect_comp,find_components_of_graph,Find_connected_component_of_this_node




def test_graph_construction():
    reads = string_to_readset("""
      1111
         111
         1  111
         1     11
        1      11
    """)
    selected_reads, skipped_reads = readselection(reads, max_cov=4, bridging=False, analyze=False, score_selection=0)
    assert selected_reads == {0, 1, 2, 3}, str(selected_reads)
    sliced_reads = reads.subset(selected_reads)
    read_graph= read_positions_graph(sliced_reads)
    edges=read_graph.get_edges()
    nodes=read_graph.get_nodes()
    assert len(nodes)==4
    for node in nodes:
        assert len(node.get_connections())==3
    assert len(edges)==6
    for edge in edges:
        (a,b,i)=edge
        assert i==1

def test_bigger_example():
    #Only single end reads
    reads = string_to_readset("""
      1111
           11111
                1111
       0000
              1111
        00000
               0000
               1111
       11111
    """)
    #Skipped here the selection process
    read_graph= read_positions_graph(reads)
    nodes=read_graph.get_nodes()
    assert len(nodes)==9
    edges=read_graph.get_edges()
    print(len(edges))
    assert len(edges)==17
    weight_1=0
    weight_2=0
    weight_3=0
    weight_4=0
    for (n_1,n_2,weight)in edges:
        if weight==1:
            weight_1+=1
        if weight==2:
            weight_2+=1
        if weight==3:
            weight_3+=1
        if weight==4:
            weight_4+=1
    assert weight_1==3
    assert weight_2==4
    assert weight_3==7
    assert weight_4==3

def test_graph_paired_end():
    reads = string_to_readset("""
      11 11
      0 0000
       11  11
    """)
    read_graph= read_positions_graph(reads)
    nodes=read_graph.get_nodes()
    assert len(nodes)==3
    edges=read_graph.get_edges()
    assert len(edges)==3
    for node in nodes:
        node_min=node.get_min()
        node_max=node.get_max()
        assert node_min<=20
        assert node_max>=50


    #looking if nodes are sorted
    assert nodes[0].get_min()==10
    assert nodes[1].get_min()==10
    assert nodes[2].get_min()==20
    assert nodes[0].get_max()==50
    assert nodes[1].get_max()==60
    assert nodes[2].get_max()==70
    assert nodes[0].get_positions()=={10,20,40,50}
    assert nodes[1].get_positions()=={10,30,40,50,60}
    assert nodes[2].get_positions()=={20,30,60,70}
    assert nodes[0].get_index()==0
    assert nodes[1].get_index()==1
    assert nodes[2].get_index()==2

    #looking if edges are correct
    (node_1,con_node_1,weight_1)=edges[0]
    (node_2,con_node_2,weight_2)=edges[1]
    (node_3,con_node_3,weight_3)=edges[2]
    assert node_1==nodes[0]
    assert con_node_1==nodes[1]
    assert node_2==nodes[0]
    assert con_node_2==nodes[2]
    assert node_3==nodes[1]
    assert con_node_3==nodes[2]
    assert weight_1==3
    assert weight_2==1
    assert weight_3==2
    #assert 0==1


def test_component_construction():
    reads = string_to_readset("""
      11 11
      0 0000
       11  11
    """)
    read_graph= read_positions_graph(reads)
    nodes=read_graph.get_nodes()
    factor=2
    assert len(nodes)==3
    edges=read_graph.get_edges()
    assert len(edges)==3
    #Each node one component
    node_component_0=Connect_comp( nodes[0],factor)
    node_component_1=Connect_comp(nodes[1],factor)
    node_component_2=Connect_comp(nodes[2],factor)

    assert node_component_0.get_length()==4
    assert node_component_1.get_length()==5
    assert node_component_2.get_length()==4

    assert node_component_0.get_positions()=={10,20,40,50}
    assert node_component_1.get_positions()=={10,30,40,50,60}
    assert node_component_2.get_positions()=={20,30,60,70}

    assert node_component_0.get_min()==10
    assert node_component_1.get_min()==10
    assert node_component_2.get_min()==20

    assert node_component_0.get_max()==50
    assert node_component_1.get_max()==60
    assert node_component_2.get_max()==70

    assert len (node_component_0.get_stored_for_later())==1
    assert node_component_0.get_stored_for_later().pop()== nodes[2]
    assert node_component_0.get_blocks()=={3:1}



    assert len(node_component_1.get_stored_for_later())==0
    assert node_component_1.get_blocks()=={2: 1, 3: 1}



    assert len (node_component_2.get_stored_for_later())==1
    assert node_component_2.get_stored_for_later().pop()== nodes[0]
    assert node_component_2.get_blocks()=={2:1}


    assert node_component_0.get_included_nodes().pop()==nodes[0]
    assert node_component_0.get_analyze_nodes()[0]==nodes[1]
    assert len(node_component_0.get_analyze_nodes())==1


    assert node_component_1.get_included_nodes().pop()==nodes[1]
    assert node_component_1.get_analyze_nodes()[0]==nodes[0]
    assert node_component_1.get_analyze_nodes()[1]==nodes[2]
    assert len(node_component_1.get_analyze_nodes())==2


    assert node_component_2.get_included_nodes().pop()==nodes[2]
    assert node_component_2.get_analyze_nodes()[0]==nodes[1]
    assert len(node_component_2.get_analyze_nodes())==1



def test_find_connected_component():
    reads = string_to_readset("""
      11 11
      0 0000
       11  11
    """)
    read_graph= read_positions_graph(reads)
    factor=2
    components=find_components_of_graph(read_graph,factor)
    assert len(components)==1
    only_com=components.pop()
    assert len(only_com.get_included_nodes())==3


def test_accumulation_of_nodes():
    reads = string_to_readset("""
      11111
       0000
        11111
         11  1
          00000
    """)
    read_graph= read_positions_graph(reads)
    factor=3
    components=find_components_of_graph(read_graph,factor)
    assert len(components)==1
    only_component=components.pop()
    assert len(only_component.get_included_nodes())==5
    assert only_component.get_length()==9
    assert only_component.get_max()==90
    assert only_component.get_min()==10
    assert len(only_component.get_stored_for_later())==0
    assert len(only_component.get_analyze_nodes())==0

def test_two_components():
    reads = string_to_readset("""
       00000
         000
      111
         000
           11111
             00000
    """)
    factor=2
    read_graph= read_positions_graph(reads)
    components=find_components_of_graph(read_graph,factor)
    assert len(components)==2
    first_component=components.pop()
    second_component=components.pop()
    assert len(first_component.get_included_nodes())==4
    assert first_component.get_length()==6
    assert first_component.get_max()==60
    assert first_component.get_min()==10
    assert len(first_component.get_stored_for_later())==0
    assert len(first_component.get_analyze_nodes())==0
    assert len(second_component.get_included_nodes())==2
    assert second_component.get_length()==7
    assert second_component.get_max()==120
    assert second_component.get_min()==60
    assert len(second_component.get_stored_for_later())==0
    assert len(second_component.get_analyze_nodes())==0


def test_multiple_reads():
    reads = string_to_readset("""
    00000
    00000
    111111
            000
            11111
            00000
    """)
    factor=2
    read_graph= read_positions_graph(reads)
    components=find_components_of_graph(read_graph,factor)
    assert len(components)==2
    first_component=components.pop()
    second_component=components.pop()
    assert len(first_component.get_included_nodes())==3
    assert first_component.get_length()==6
    assert first_component.get_max()==60
    assert first_component.get_min()==10
    assert len(first_component.get_stored_for_later())==0
    assert len(first_component.get_analyze_nodes())==0
    assert len(second_component.get_included_nodes())==3
    assert second_component.get_length()==5
    assert second_component.get_max()==130
    assert second_component.get_min()==90
    assert len(second_component.get_stored_for_later())==0
    assert len(second_component.get_analyze_nodes())==0

def test_components_of_paired_reads():
    reads = string_to_readset("""
    00    000
       11     111
    """)
    factor=2
    read_graph= read_positions_graph(reads)
    components=find_components_of_graph(read_graph,factor)
    assert len(components)==2
    first_component=components.pop()
    second_component=components.pop()
    assert len(first_component.get_included_nodes())==1
    assert first_component.get_length()==5
    assert first_component.get_max()==90
    assert first_component.get_min()==10
    assert len(first_component.get_stored_for_later())==0
    assert len(first_component.get_analyze_nodes())==0
    assert len(second_component.get_included_nodes())==1
    assert second_component.get_length()==5
    assert second_component.get_max()==130
    assert second_component.get_min()==40
    assert len(second_component.get_stored_for_later())==0
    assert len(second_component.get_analyze_nodes())==0
