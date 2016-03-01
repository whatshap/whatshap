from .phasingutils import string_to_readset
from whatshap.core import readselection
from whatshap.connect_analysis import  Element,read_positions_graph,Connect_comp,find_components_of_graph,Find_connected_component_of_this_node
from whatshap.connect_analysis_2 import Node,Edge,Edges,Component,Graph



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
    node_0=nodes.pop()
    node_1=nodes.pop()
    node_2=nodes.pop()
    assert node_0.get_min()==10
    assert node_1.get_min()==10
    assert node_2.get_min()==20
    assert node_0.get_max()==50
    assert node_1.get_max()==60
    assert node_2.get_max()==70
    assert node_0.get_positions()=={10,20,40,50}
    assert node_1.get_positions()=={10,30,40,50,60}
    assert node_2.get_positions()=={20,30,60,70}
    assert node_0.get_index()==0
    assert node_1.get_index()==1
    assert node_2.get_index()==2

    #looking if edges are correct
    (node_11,con_node_1,weight_1)=edges[0]
    (node_22,con_node_2,weight_2)=edges[1]
    (node_33,con_node_3,weight_3)=edges[2]
    assert node_11==node_0
    assert con_node_1==node_1
    assert node_22==node_0
    assert con_node_2==node_2
    assert node_33==node_1
    assert con_node_3==node_2
    assert weight_1==3
    assert weight_2==1
    assert weight_3==2


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
    node_0=nodes.pop()
    node_1=nodes.pop()
    node_2=nodes.pop()

    node_component_2=Connect_comp(node_0,factor)
    node_component_1=Connect_comp(node_1,factor)
    node_component_0=Connect_comp(node_2,factor)

    assert node_component_0.get_length()==4
    assert node_component_1.get_length()==5
    assert node_component_2.get_length()==4

    assert node_component_2.get_positions()=={10,20,40,50}
    assert node_component_1.get_positions()=={10,30,40,50,60}
    assert node_component_0.get_positions()=={20,30,60,70}

    assert node_component_2.get_min()==10
    assert node_component_1.get_min()==10
    assert node_component_0.get_min()==20

    assert node_component_2.get_max()==50
    assert node_component_1.get_max()==60
    assert node_component_0.get_max()==70

#Look at node_component_0

    assert len (node_component_2.get_stored_for_later())==1
    assert node_component_2.get_stored_for_later().pop()== node_2
    assert node_component_2.get_blocks()=={3:1}
    assert len(node_component_2.get_analyze_nodes())==1
    assert node_component_2.get_analyze_nodes().pop()==node_1
    assert len(node_component_2.get_included_nodes())==1
    assert node_component_2.get_included_nodes().pop()==node_0


    assert len(node_component_1.get_stored_for_later())==0
    assert len(node_component_1.get_analyze_nodes()) ==2
    nodes_of_component_1=node_component_1.get_analyze_nodes()
    assert len(nodes_of_component_1)==2
    assert nodes_of_component_1.pop() == node_0
    assert nodes_of_component_1.pop() == node_2
    assert node_component_1.get_blocks()=={2: 1, 3: 1}
    assert node_component_1.get_included_nodes().pop()==node_1



    assert len (node_component_0.get_stored_for_later())==1
    assert len(node_component_0.get_analyze_nodes())==1
    assert node_component_0.get_analyze_nodes().pop()== node_1
    assert node_component_0.get_stored_for_later().pop()== node_0
    assert node_component_0.get_blocks()=={2:1}
    assert node_component_0.get_included_nodes().pop()==node_2


def test_find_connected_component():
    reads = string_to_readset("""
      11 11
      0 0000
       11  11
    """)
    read_graph= read_positions_graph(reads)
    factor=2
    assert len(read_graph.get_nodes())==3
    assert len(read_graph.get_edges())==3
    components=find_components_of_graph(read_graph,factor)
    assert len(components)==1
    only_com=components.pop()
    assert only_com.get_length()==7
    assert set(only_com.get_positions())==set([10,20,30,40,50,60,70])
    assert only_com.get_blocks() =={2: 1, 3: 1}
    assert only_com.get_min()==10
    assert only_com.get_max()==70
    included_nodes_of_comp=only_com.get_included_nodes()
    #print('included nodes of the component')
    #print(len(included_nodes_of_comp))
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
    assert len(second_component.get_included_nodes())==4
    assert second_component.get_length()==6
    assert second_component.get_max()==60
    assert second_component.get_min()==10
    assert len(second_component.get_stored_for_later())==0
    assert len(second_component.get_analyze_nodes())==0
    assert len(first_component.get_included_nodes())==2
    assert first_component.get_length()==7
    assert first_component.get_max()==120
    assert first_component.get_min()==60
    assert len(first_component.get_stored_for_later())==0
    assert len(first_component.get_analyze_nodes())==0


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
    assert second_component.get_length()==6
    assert second_component.get_max()==60
    assert second_component.get_min()==10
    assert len(second_component.get_stored_for_later())==0
    assert len(second_component.get_analyze_nodes())==0
    assert len(second_component.get_included_nodes())==3
    assert first_component.get_length()==5
    assert first_component.get_max()==130
    assert first_component.get_min()==90
    assert len(first_component.get_stored_for_later())==0
    assert len(first_component.get_analyze_nodes())==0

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
    assert first_component.get_max()==130
    assert first_component.get_min()==40
    assert len(first_component.get_stored_for_later())==0
    assert len(first_component.get_analyze_nodes())==0
    assert len(second_component.get_included_nodes())==1
    assert second_component.get_length()==5
    assert second_component.get_max()==90
    assert second_component.get_min()==10
    assert len(second_component.get_stored_for_later())==0
    assert len(second_component.get_analyze_nodes())==0

def test_finding_of_components_by_multiple_appearances_of_the_same_positions():
    reads = string_to_readset("""
    00
    00
    11
      11
      00
       10
    """)
    factor=1
    read_graph= read_positions_graph(reads)
    assert len(read_graph.get_edges())==6

    assert len(read_graph.get_nodes())==6
    components=find_components_of_graph(read_graph,factor)
    assert len(components)==2

#################Test connect_analysis_2





#TODO check if weight and len(indices in the edges are equal )

def test__connec_ana_2_finding_of_components_by_multiple_appearances_of_the_same_positions():
    reads = string_to_readset("""
    00
    00
    11
      11
      00
       10
    """)
    factor=1
    important_positions = set(reads.get_positions())

    read_graph= Graph(reads,important_positions)

    edges_of_graph=read_graph.get_edges()
    nodes_of_graph=read_graph.get_nodes()
    vals_of_nodes_of_graph=nodes_of_graph.values()
    second_values=list(vals_of_nodes_of_graph)
    node_1=second_values.pop()
    node_2=second_values.pop()
    node_3=second_values.pop()
    node_4=second_values.pop()
    node_5=second_values.pop()
    assert len(nodes_of_graph.keys())==5
    print('nodes_of_graph.keys()')
    print(nodes_of_graph.keys())
    assert len(nodes_of_graph.values())==5
    print('node_1.get_position()')
    print(node_1.get_position())
    assert node_1.get_position()==30
    print('node_2.get_position()')
    print(node_2.get_position())
    assert node_2.get_position()==50
    print('node_3.get_position()')
    print(node_3.get_position())
    assert node_3.get_position()==20
    print('node_4.get_position()')
    print(node_4.get_position())
    assert node_4.get_position()==10
    print('node_5.get_position()')
    print(node_5.get_position())
    assert node_5.get_position()==40
    #TODO : Output, that the nodes in the dictionary are not ordered

    assert node_1.get_connections()=={node_5.get_position()}
    assert node_2.get_connections()=={node_5.get_position()}
    assert node_3.get_connections()=={node_4.get_position()}
    assert node_4.get_connections()=={node_3.get_position()}
    assert node_5.get_connections()=={node_2.get_position(),node_1.get_position()}

    edge_of_node_1=edges_of_graph.get_edges_for_specific_node(node_1)
    edge_of_node_2=edges_of_graph.get_edges_for_specific_node(node_2)
    edge_of_node_3=edges_of_graph.get_edges_for_specific_node(node_3)
    edge_of_node_4=edges_of_graph.get_edges_for_specific_node(node_4)
    edge_of_node_5=edges_of_graph.get_edges_for_specific_node(node_5)


    edge_1=edge_of_node_1.pop() # Edge 30-40
    assert edge_1.get_max()==40
    assert edge_1.get_min()==30
    assert edge_1.get_indices()=={3,4}
    assert edge_1.get_weight()==2

    edge_2=edge_of_node_2.pop() # Edge 50-40
    assert edge_2.get_max()==50
    assert edge_2.get_min()==40
    assert edge_2.get_indices()=={5}
    assert edge_2.get_weight()==1

    edge_3=edge_of_node_3.pop() # Edge 20-10
    assert edge_3.get_max()==20
    assert edge_3.get_min()==10
    assert edge_3.get_indices()=={0,1,2}
    assert edge_3.get_weight()==3

    edge_4=edge_of_node_4.pop() # Edge 10-20
    assert edge_4.get_max()==20
    assert edge_4.get_min()==10
    assert edge_4.get_indices()=={0,1,2}
    assert edge_4.get_weight()==3



    edge_51=edge_of_node_5.pop() # Edge 40-30 & 40-50
    assert edge_51.get_max()==50
    assert edge_51.get_min()==40
    assert edge_51.get_indices()=={5}
    assert edge_51.get_weight()==1

    edge_52=edge_of_node_5.pop() # Edge 40-30 & 40-50
    assert edge_52.get_max()==40
    assert edge_52.get_min()==30
    assert edge_52.get_indices()=={3,4}
    assert edge_52.get_weight()==2

def test_very_small():
    reads = string_to_readset("""
    00
    11
    """)
    factor=1
    important_positions = set(reads.get_positions())
    read_graph= Graph(reads,important_positions)

    edges_of_graph=read_graph.get_edges()
    nodes_of_graph=read_graph.get_nodes()
    vals_of_nodes_of_graph=nodes_of_graph.values()
    second_values=list(vals_of_nodes_of_graph)
    node_1=second_values.pop()
    node_2=second_values.pop()
    assert len(nodes_of_graph.keys())==2
    assert len(nodes_of_graph.values())==2
    assert node_1.get_position()==20
    assert node_2.get_position()==10
    assert node_1.get_connections()=={node_2.get_position()}
    assert node_2.get_connections()=={node_1.get_position()}
    edge_of_node_1=edges_of_graph.get_edges_for_specific_node(node_1)
    edge_of_node_2=edges_of_graph.get_edges_for_specific_node(node_2)

    edge_1=edge_of_node_1.pop()
    edge_2=edge_of_node_2.pop()
    print('Edge_1')
    print(edge_1)
    #print('Edge_2')
    #print(edge_2)
    assert edge_1.get_max()==20
    assert edge_2.get_max()==20
    #Not possible because the __eq__ is not defined
    #assert edge_1==edge_2
    assert edge_1.get_min()==10
    assert edge_2.get_min()==10
    assert edge_1.get_indices()=={0,1}
    assert edge_2.get_indices()=={0,1}
    assert edge_1.get_weight()==2
    assert edge_2.get_weight()==2

def test_small_example_for_component():
    reads = string_to_readset("""
    00
    11
    """)
    factor=1
    important_positions = set(reads.get_positions())
    read_graph= Graph(reads,important_positions)
    Edges=read_graph.get_edges()
    print('Edges')
    print(Edges)
    print(Edges.get_all_edges())

    read_graph.find_components_of_graph(factor)
    components=read_graph.get_components()


    print('Components')
    print(components)
    print('Length of components')
    print(len(components))
    assert len(components)==1
    first_component=components.pop()
    print(first_component.get_nodes_of_component())
    assert len(first_component.get_nodes_of_component())==2
    assert first_component.get_length()==2
    #assert 0==1

def test_bigger_example_of_components():
    reads = string_to_readset("""
    00
    00
    11
      11
      00
       10
    """)
    factor=1
    important_positions = set(reads.get_positions())
    read_graph= Graph(reads,important_positions)
    graph_edges=read_graph.get_edges()
    graph_nodes=read_graph.get_nodes()
    print('graph_edges')
    print(graph_edges.get_all_edges())
    print(len(set(graph_edges.get_all_edges())))
    assert len(graph_edges.get_all_edges())==5
    assert len(graph_nodes.keys())==5

    read_graph.find_components_of_graph(factor)
    components=read_graph.get_components()
    print('components')
    print(components)
    #components=find_components_of_graph(Graph,factor)
    assert len(components)==2
    component_1=components.pop()
    component_2=components.pop()

    #print(component_1.get_length())

    #assert (component_1.get_length())==3
    #assert component_2.get_length()==3

    #assert component_1.get_nodes_of_component()=={10,20}
    #assert component_2.get_nodes_of_component()=={30,40,50}

    #Kann es hier nicht pruefen da dauernd die componentenreihenfolge wechselt.

    assert 0==1