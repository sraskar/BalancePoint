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
    #print("#leaf",leaf)
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
    
    write_dot(G,'original.dot')
    
    w = { tuple([u,v]): d['label'] for u,v,d in G.edges(data=True)}

    return w

def draw_graph(w,tech):
    G=nx.DiGraph()
    #print("Drawing : ",w)
    #print( u for k in w.items())
    for k,v in w.items():
        #print(k[0],k[1],":",v)
        G.add_edge(k[0],k[1],label=v)

    #print(u)
    write_dot(G,tech+'.dot')

    return

def max_traffic (w_ref,B):


    w = dict(w_ref)

    while True:
        s = Swing(w) 
        if not s.non_critical_path:
            print ("Graph Balanced!")
            break
        min_path =  min(s.non_critical_path.items(), key=lambda kv: kv[1])[0]
        set_edge_maybe_on_cp = set(min_path)
        #print (min_path) 
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
        #print (max_traffic_edge_not_on_cp, w[max_traffic_edge_not_on_cp])

    
    draw_graph(w,"max_traffic")
    #d_opt_buffer = {edge: w[edge] - w_ref[edge] for edge in w}
    d_opt_buffer = {}
    for edge in w:
        buf = w[edge] - w_ref[edge]
        if buf:
            d_opt_buffer[edge] = buf

    gao = Swing(w_ref).opt_edge_buffer


    gao_graph={}
    for edge in w_ref:
        if( edge not in gao):
            gao_graph[edge] =  w_ref[edge]
        else:
            gao_graph[edge] = gao[edge] + w_ref[edge]

    draw_graph(gao_graph,"gao_linear")
    print("Gao Bal : ",gao_graph)
    # print("Traffic : ",Swing(w_ref).traffic)
    print ("Us  Buffers : ", sum(d_opt_buffer.values()))
    print ("Gao Buffers : ", sum(gao.values()))

    # print ("w", w)
    print ("Us  Bal : ", d_opt_buffer)
    print ("Gao Bal : ", gao)

    # print("Calling Simulator : ")

    os.system("./Simulator 10 original.dot")
    os.system("./Simulator 10 gao_linear.dot")
    os.system("./Simulator 10 max_traffic.dot")
    
    #print("Gao      Cycles : ",os.system("./Simulator 3 gao_linear.dot"))
    #print("Max Traf Cycles : ",os.system("./Simulator 3 max_traffic.dot"))
    return
    

def brute_force(w_ref,B):
    import itertools as it
    
    s=Swing(w_ref)

    ''' Begin Brute Force Method '''
    
    E = len(s.l_edges)
    #B = 4   # uncomment to overide. 

    individual_buffers=[]
    individual_buffers = (list(i for i in it.combinations_with_replacement(range(B+1),E) if sum(i) == B))
    #print(individual_buffers)
    
    # for bb in range(1,B+1):
    #     #print(bb)
    #     ind_buffers = (list(i for i in it.combinations_with_replacement(range(bb+1),E) if sum(i) == bb))
    #    individual_buffers.extend(ind_buffers)
    
    # for ib in individual_buffers:
    #     print (f"individual_buffers: {ib}")
    #     print (f"placement on edges: {list(it.permutations(ib,E))}")

    ### Python noob
    # l = []
    # for ib in individual_buffers:
    #     for config in it.permutations(ib,E):
    #         l.append(config)
    ### Python hard to read
    # l = [config for ib in individual_buffers for config in it.permutations(ib,E)]
    # print (l)

    # Function programming approach
    # from functools import partial
    # from itertools import chain

    # f = partial(it.permutations,r=E)
    # l = map(f, individual_buffers)
    # print (list(chain.from_iterable(l)))
    # print (l)

    ### Python noob 2
    l2 = []
    for ib in individual_buffers:
        l2.extend(it.permutations(ib,E))

    
    #print("Original Graph : ",w_ref)
    import subprocess 
    cur_map={}
    i=0
    # Mapping edges to numbers 
    for v in s.l_edges:
        cur_map[i]=v
        i=i+1

    cur_w={}
    min_w={}
    counter=1
    min_cycle = -1

    for item in l2:
        cur_w={}

        c=0
        for ee in item:
            cur_w[cur_map[c]] = w_ref[cur_map[c]] + ee
            c = c + 1
        
        counter = counter + 1
        
        # Simulating the current configuration.
        draw_graph(cur_w,"brute_"+str(counter-1))
        command = ['./Simulator', '3' , 'brute_'+str(counter-1)+'.dot']
        p = subprocess.Popen(command, stdout=subprocess.PIPE)
        text = p.stdout.read().decode()
        retcode = p.wait()
        c_count, p_count = text.split(" ",1)
        
        if min_cycle == -1:
            min_cycle = c_count

        if c_count < min_cycle:
            min_cycle = c_count
            min_w = cur_w

        #print( min_cycle , item , cur_w )

    print("Min Cycles   : ",min_cycle)
    print("Min Config   : ",min_w)
    print("Combinations : ",len(l2))

    os.system("rm brute*")

    return 


if __name__ == '__main__':

    # Flags for selective execution 
    verify = 1
    gen_rand = 0
    m_traffic = 0
    B = 10

    # Weighted graph to balence
    #Gao fig 8
    # w_ref = {
    #     ('u1', 'u2'): 1,
    #     ('u2', 'u4'): 3,
    #     ('u2', 'u3'): 4,
    #     ('u3', 'u4'): 4,
    #     ('u4', 'u5'): 1,
    #     ('u1', 'u5'): 20,
    #     ('u3', 'u5'): 5
    # }

    # w_ref = {
    #     ('2', '3'): 1,
    #     ('2', '4'): 1,
    #     ('3', '5'): 1,
    #     ('4', '5'): 1
    # }

    ##Appl fig 1
    # w_ref = {
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
    # }

    w_ref = {
   ('2', '3'): 1,
   ('2', '4'): 1,
   ('3', '7'): 90,
   ('3', '5'): 1,
   ('5', '7'): 80,
   ('7', '9'): 1,
   ('4', '8'): 20,
   ('4', '6'): 1,
   ('6', '8'): 1,
   ('8', '9'): 1
    }

    if (gen_rand == 1):
        w_ref = generate_random_dag(15,0.1,11) 
        print ("Original : ",w_ref)
    if (verify == 1 ):
        brute_force(w_ref,B)
    if (m_traffic == 1 ):
        max_traffic (w_ref,B)

    #print (set_edge_not_on_cp)


    #print(s.adjacency_list_topdown)
    #print(f'Max buffer alowed: {max_b}')
    #l_buffer_updated, diff_delay = Swing(w).constrained_edge_buffer(max_b, max_bb_eval)  #adjacency_buffered
    #print(f'list buffer: {l_buffer_updated}')
    #print(f'number of buffer used: {sum(l_buffer_updated.values())}')
    #print(f'max diff delay: {diff_delay}')



