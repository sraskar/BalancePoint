#!/usr/bin/env python
from collections import defaultdict
from itertools import chain
from functools import lru_cache, partial

# Python Typing
from typing import List, Dict, Tuple, Set, NewType
Node = NewType('Node', str)
Edge = Tuple[Node,Node]

import numpy as np


def lazy_property(f):
    return property(lru_cache()(f))

from itertools import tee
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

class B(object):
    def __init__(self, weigh):
        self.weigh = weigh # Top-> bottom by default. Only once

    @lazy_property
    def adjacency_topdown(self) -> Dict[Node, List[Node] ]:
        d = defaultdict(list)
        for i,o in self.weigh:
            d[i].append(o)
        return d

    @lazy_property
    def adjacency_bottomup(self) -> Dict[Node, List[Node] ] :
        d = defaultdict(list)
        for k,lv in self.adjacency_topdown.items():
            for v in lv:
                d[v].append(k)
        return d

    @lazy_property
    def l_node(self) -> Set[Node]:
        return self.adjacency_topdown.keys() |  self.adjacency_bottomup.keys()

    @lazy_property
    def order(self) -> Dict[Node,int]:
        return {k:i for i,k in enumerate(self.l_node)}

    @lazy_property
    def leaf(self) -> Node:
        '''Assume only one leaf for now'''
        return (self.l_node - self.adjacency_topdown.keys()).pop()


    @lazy_property
    def delta_degree(self) -> Dict[Node, int]:
        in_ =  self.adjacency_topdown
        out_ = self.adjacency_bottomup
        return { n:(len(in_[n]) - len(out_[n])) for n in self.l_node}


    @lru_cache(None)
    def path_node(self, cur) -> List[Node]:
        '''Compute all the node form cur to the root of the DAG.
            Assume a unique root'''
        d =  self.adjacency_bottomup

        if d[cur]:
            it =  chain.from_iterable(self.path_node(p) for p in d[cur])
        else:
            it = [ [] ]

        return [ k + [cur] for k in it ]

    @lazy_property
    def path(self) -> Dict[Tuple[Edge], int]:
        '''
        Compute for all the path ( of the form (e1->e2),(e2->e3) )
        and the sum of th weigh  from the leaf to the read
        '''
        d = defaultdict(list)

        for p in self.path_node(self.leaf):
                path_edge = tuple(pairwise(p))
                d[path_edge] = sum(map(self.weigh.get,path_edge))

        return d

    @lazy_property
    def critical_weigh(self) -> int:
        return max(self.path.values())

    @lazy_property
    def non_critical_path(self) -> Dict[Tuple[Edge], int]:
        c = self.critical_weigh
        return {p:w for p,w in self.path.items() if w != c}

    #                        _                        
    # |  o ._   _   _. ._   |_) ._ _   _  ._ _. ._ _  
    # |_ | | | (/_ (_| |    |   | (_) (_| | (_| | | | 
    #                                  _|             
    #
    # https://doi.org/10.1016/0743-7315(89)90041-5
    #
    # Maximize \sum_i u_i (outdegree(i) - indegree(i))
    # Subject to
    #       u_i - u_j \leq -w_{ij} for all (i,j) \subseteq E
    #       u_n - u_1 = w_{st}
    #       u_i >= 0 
     
    @lazy_property
    def lp_constrain_matrix(self):

        od = self.order
        n_node, n_edge = map(len, (self.l_node,self.weigh))

        A = np.zeros( (n_edge, n_node), float)
        for idx, (i,o) in enumerate(self.weigh):
                A[idx,od[i] ] = 1
                A[idx,od[o] ] = -1

        # Expant with positif contrains
        A_pos = np.identity(n_node) * -1
        return np.concatenate((A, A_pos), axis=0)
    
    @lazy_property
    def lp_constrain_vector(self):
        b = np.array([-k for k in self.weigh.values()],dtype=float)
    
        # Expend with posotif contrains
        b_pos = np.zeros(len(self.l_node))
        return np.concatenate( (b,b_pos) )

    @lazy_property
    def lp_objective_vector(self):
        # Float is required by CVXopt
        return np.array([-1*self.delta_degree[i] for i in self.order],dtype=float) # -1 <=> Minizise
    
    @lazy_property
    def lp_firing_buffered(self):
        '''This function use GAO algorithm. It will generated the optimal buffer 
           need to minimze the graph
        '''
        
        from cvxopt import matrix, solvers
        # We add the constrain that each objectif need to me postiif
        A,b = map(matrix, (self.lp_constrain_matrix,self.lp_constrain_vector))
        c = matrix(self.lp_objective_vector)

        #Integer solution use GLPK and turn off verbose output
        solvers.options['glpk'] = {'msg_lev' : 'GLP_MSG_OFF'}
        sol=solvers.lp(c,A,b, solver='glpk')

        sol['x_int'] = [int(round(i)) for i in sol['x'] ]
        assert all(i == f for i,f in zip(sol['x_int'],sol['x']))
        return dict(zip(self.order,sol['x_int'])) # Assume ordered dict

    @lazy_property
    def opt_adjacency_buffer(self) -> Dict[Tuple[str,str],int]:
        b = self.lp_firing_buffered
        w = self.weigh
        return { (i,o): (b[o] - b[i]) - w for (i,o),w in w.items() if (b[o] - b[i]) != w}

     #                   
     #|\/| o ._  |\/|  _.    
     #|  | | | | |  | (_| >< 
     #                                 
     # Miminize  max(u_n - u_1)
     # Subject to
     #       u_n - u_1 \geq 0
     #       u_i unrestricted
    def bb(self, x, A, C, max_b) -> int:
        '''This function compute the function to minimize.
            f(x) = max(C-sum(A*x))  (1)     
            Subject to:
                sum(x) > max_b      (2)
                C-sum(A*x) > 0      (3)

        ''' 

        dim = x.get_n()
        b = [x.get_coord(i) for i in range(dim)]

        delta = C - np.sum(np.multiply(A,b),axis=1)

        x.set_bb_output(0, max(delta))     # (1)
        x.set_bb_output(1, sum(b) - max_b) # (2)

        for i,v in enumerate(delta,2):
            x.set_bb_output(i,-v)          #(3)
 
        return 1

    def constrained_adjacency_buffer(self,max_b) -> Tuple[Dict[Tuple[Edge],int], int]:
        '''
        We will optimize min(max(abs(C-p))) where p are the total weigh of the graph
        under the constrain the sum(e) < max_b
        '''
        # Initialization
        nc_path = self.non_critical_path 
        opt_adj_buf = self.opt_adjacency_buffer 

        # Edge matrix and constrain vector
        A = np.array( [ [ e in pairs for e in opt_adj_buf ] for pairs, w in nc_path.items() ] , dtype=int) 
        C = self.critical_weigh - np.array(list(nc_path.values()),dtype=int)

        bb = partial(self.bb,max_b=max_b,A=A,C=C)

        # Nomad parameter  
        lb = [0] * len(opt_adj_buf) 
        ub = list(opt_adj_buf.values())
        x0 = lb

        # Each path will give us a constrain (dela >0)
        params = [f'BB_OUTPUT_TYPE OBJ EB {" ".join(["EB"] * len(nc_path))}',
                  'MAX_BB_EVAL 100',   # Parameter
                  'LOWER_BOUND * 0',   # Buffer are postif
                  'BB_INPUT_TYPE * I', # All integer
                  'DISPLAY_DEGREE 0']

        import PyNomad
        x_return, f_return, h_return, nb_evals, nb_iters, stopflag = PyNomad.optimize(bb,x0,lb,ub,params)

        # Map back to the correct name
        l_buffer_updated = {k:int(v) for k,v in zip(opt_adj_buf,x_return)}
        return l_buffer_updated, f_return


#Gao fig 7
w = {('u1','u2'): 1,
     ('u2','u4'): 3,
     ('u2','u3'): 4,
     ('u3','u4'): 4,
     ('u4','u5'): 1,
     ('u1','u5'): 20,
     ('u3','u5'): 5}

#Appl fig 1
w = {('u1','u2'): 1,
     ('u1','u3'): 1,
     ('u2','u6'): 90,
     ('u2','u4'): 1,
     ('u4','u6'): 80,
     ('u6','u8'): 1,
     ('u3','u7'): 20,
     ('u3','u5'): 1,
     ('u5','u7'): 1,
     ('u7','u8'): 1}

max_b = 100

print (f'Max buffer alowed: {max_b}')
l_buffer_updated, diff_delay = B(w).constrained_adjacency_buffer(max_b) #adjacency_buffered
print (f'list buffer: {l_buffer_updated}')
print (f'number of buffer used: {sum(l_buffer_updated.values())}')
print (f'max diff delay: {diff_delay}')

