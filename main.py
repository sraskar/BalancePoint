#!/usr/bin/env python
from balence import Swing
import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import write_dot

import random
import logging

import sys,os

logging.getLogger().setLevel(logging.ERROR)

def generate_random_dag(n, p,seed=None):
    random_graph = nx.fast_gnp_random_graph(n, p, directed=True, seed=seed)
    G = nx.DiGraph( [(u, v) for (u, v) in random_graph.edges() if u < v and (u>2 and v>2)])
    # Merge all the leaf
    root=2
    G.add_edges_from([(root,n) for n,d in G.in_degree() if d==0])
    #leaf=G.number_of_nodes()+1
    leaf=max(G.nodes)
    print("#leaf",leaf)
    G.add_edges_from([(n,leaf+1) for n,d in G.out_degree() if d==0])

    assert (nx.is_directed_acyclic_graph(G))

    random.seed(seed)
    for u,v,d in G.edges(data=True):
            d['label'] = random.randint(1,20)

    # G.remove_edge(7,11)
    # G.remove_edge('root',7)
    # G.remove_edge(11,'leaf')
    # G.remove_edge(4,9)
    # G.remove_edge('root',4)
    # G.remove_edge(0,5)
    # G.remove_edge(5,'leaf')
    # G.remove_edge(0,13)
    # G.remove_edge(13,'leaf')
    # G.remove_edge('root',0)

    # G.remove_edge(3,6)
    # G.add_edge('root',6, weight=12)




    # for n in [0,5,13,4,7,11,3]:
    #     G.remove_node(n)

    pos=nx.spring_layout(G)
    labels = nx.get_edge_attributes(G,'label')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
    
    write_dot(G,'66666.dot')
    
    w = { tuple([u,v]): d['label'] for u,v,d in G.edges(data=True)}

    return w



if __name__ == '__main__':

    # Weighted graph to balence
    #Gao fig 8
    w_ref = {
        ('u1', 'u2'): 1,
        ('u2', 'u4'): 3,
        ('u2', 'u3'): 4,
        ('u3', 'u4'): 4,
        ('u4', 'u5'): 1,
        ('u1', 'u5'): 20,
        ('u3', 'u5'): 5
    }

    ##Appl fig 1
    #w_ref = {
    #    ('u1', 'u2'): 1,
    #    ('u1', 'u3'): 1,
    #    ('u2', 'u6'): 90,
    #    ('u2', 'u4'): 1,
    #    ('u4', 'u6'): 80,
    #    ('u6', 'u8'): 1,
    #    ('u3', 'u7'): 20,
    #    ('u3', 'u5'): 1,
    #    ('u5', 'u7'): 1,
    #    ('u7', 'u8'): 1
    #}

    w_ref = generate_random_dag(15,0.1,11) 
    print (w_ref)

    w = dict(w_ref)

    while True:
        s = Swing(w) 
        if not s.non_critical_path:
            print ("Graph Balanced!")
            break
        min_path =  min(s.non_critical_path.items(), key=lambda kv: kv[1])[0]
        set_edge_maybe_on_cp = set(min_path)
        print (min_path) 
        set_edge_on_cp = set()
        for path in s.critical_paths:
            set_edge_on_cp.update(path)

        set_edge_not_on_cp = set_edge_maybe_on_cp - set_edge_on_cp

        from collections import defaultdict
        d_order_trafic = defaultdict(list)
        for edge in set_edge_not_on_cp:
            traffic = s.traffic[edge]
            d_order_trafic[traffic] += [edge]

        max_traffic = max(d_order_trafic)
        max_traffic_edges_not_on_cp = d_order_trafic[max_traffic]

        # Maxium ordering (node at the botom of the graph. Does not seem usefull)
        max_traffic_edge_not_on_cp = max(max_traffic_edges_not_on_cp,key= lambda k: s.ordering(k))
        w[max_traffic_edge_not_on_cp] += 1
        print (max_traffic_edge_not_on_cp, w[max_traffic_edge_not_on_cp])

    #d_opt_buffer = {edge: w[edge] - w_ref[edge] for edge in w}
    d_opt_buffer = {}
    for edge in w:
        buf = w[edge] - w_ref[edge]
        if buf:
            d_opt_buffer[edge] = buf

    gao = Swing(w_ref).opt_edge_buffer


    print("Traffic : ",Swing(w_ref).traffic)
    print ("Us  Buffers : ", sum(d_opt_buffer.values()))
    print ("Gao Buffers : ", sum(gao.values()))

    print ("w", w)
    print ("Us  Bal : ", d_opt_buffer)
    print ("Gao Bal : ", gao)

    print("Calling Simulator : ")
    os.system("./Simulator 3 66666.dot")




    #print (set_edge_not_on_cp)


    #print(s.adjacency_list_topdown)
    #print(f'Max buffer alowed: {max_b}')
    #l_buffer_updated, diff_delay = Swing(w).constrained_edge_buffer(max_b, max_bb_eval)  #adjacency_buffered
    #print(f'list buffer: {l_buffer_updated}')
    #print(f'number of buffer used: {sum(l_buffer_updated.values())}')
    #print(f'max diff delay: {diff_delay}')



