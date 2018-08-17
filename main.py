#!/usr/bin/env python
from balence import Swing
import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import write_dot

import random
import logging

logging.getLogger().setLevel(logging.DEBUG)

def generate_random_dag(n, p,seed=None):
    random_graph = nx.fast_gnp_random_graph(n, p, directed=True, seed=seed)
    G = nx.DiGraph( [(u, v) for (u, v) in random_graph.edges() if u < v])
    # Merge all the leaf
    G.add_edges_from([('root',n) for n,d in G.in_degree() if d==0])
    G.add_edges_from([(n,'leaf') for n,d in G.out_degree() if d==0])
 
    assert (nx.is_directed_acyclic_graph(G))

    random.seed(seed)
    for u,v,d in G.edges(data=True):
            d['weight'] = random.randint(1,20)


    pos=nx.spring_layout(G)
    labels = nx.get_edge_attributes(G,'weight')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
    
    write_dot(G,'66666.dot')
    
    w = { tuple([u,v]): d['weight'] for u,v,d in G.edges(data=True)}

    return w


def increment_edge_weight(w,edge):
    #print(w)
    w[edge] += 1    
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
    w_ref = {
        ('u1', 'u2'): 1,
        ('u1', 'u3'): 1,
        ('u2', 'u6'): 90,
        ('u2', 'u4'): 1,
        ('u4', 'u6'): 80,
        ('u6', 'u8'): 1,
        ('u3', 'u7'): 20,
        ('u3', 'u5'): 1,
        ('u5', 'u7'): 1,
        ('u7', 'u8'): 1
    }

    w_ref = generate_random_dag(15,0.1,11) 
    print (w_ref)
    #del w_ref[(7,11)]
    #del w_ref[('root',7)]
    #w_ref[('root',11)] = 1
    #del w_ref[(4,9)]
    #w_ref[('root',9)] = 1

    #for k,v in w_ref.items():
    #    w_ref[k] = 1

    w = dict(w_ref)

    print (w)

    while True:
        s = Swing(w) 
        if not s.non_critical_path:
            print ("We converge")
            break
        min_path =  min(s.non_critical_path.items(), key=lambda kv: kv[1])[0]
        set_edge_maybe_on_cp = set(min_path)

        set_edge_on_cp = set()
        for path in s.critical_paths:
            set_edge_on_cp.update(path)

        set_edge_not_on_cp = set_edge_maybe_on_cp - set_edge_on_cp
        max_traffic_edge_not_on_cp = max(set_edge_not_on_cp, key= lambda edge: s.traffic[edge]) 

        w[max_traffic_edge_not_on_cp] += 1

    #d_opt_buffer = {edge: w[edge] - w_ref[edge] for edge in w}
    d_opt_buffer = {}
    for edge in w:
        buf = w[edge] - w_ref[edge]
        if buf:
            d_opt_buffer[edge] = buf

    gao = Swing(w_ref).opt_edge_buffer
    print ("Us", sum(d_opt_buffer.values()))
    print ("Gao", sum(gao.values()))

    print ("w", w)
    print ("Us", d_opt_buffer)
    print ("Gao", gao)
    #print (set_edge_not_on_cp)

    # for non critical paths 
    #     #find path with min weight 
    #         find edge with max_traffic on those paths
    #         increment edge weight until path weight == critical path weight 


    #print(s.path)
    #w[('u1','u2')] = 2

    #s2 =Swing(w)
    # optarion

    #print(s.adjacency_list_topdown)
    #print(f'Max buffer alowed: {max_b}')
    #l_buffer_updated, diff_delay = Swing(w).constrained_edge_buffer(max_b, max_bb_eval)  #adjacency_buffered
    #print(f'list buffer: {l_buffer_updated}')
    #print(f'number of buffer used: {sum(l_buffer_updated.values())}')
    #print(f'max diff delay: {diff_delay}')



