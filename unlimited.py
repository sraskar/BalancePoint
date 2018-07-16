#!/usr/bin/env python
from functools import lru_cache
def lazy_property(f):
    return property(lru_cache()(f))

import numpy as np
from collections import defaultdict

class B(object):
    def __init__(self, weigh):
        self.weigh = weigh # Top-> bottom by default. Only once

    @lazy_property
    def adjacency_topdown(self):
        d = defaultdict(list)
        for i,o in self.weigh:
            d[i].append(o)
        return d

    @lazy_property
    def adjacency_bottomup(self):
        from collections import defaultdict
        d = defaultdict(list)
        for k,lv in self.adjacency_topdown.items():
            for v in lv:
                d[v].append(k)
        return d

    @lazy_property
    def l_node(self):
        return self.adjacency_topdown.keys() |  self.adjacency_bottomup.keys()

    @lazy_property
    def order(self):
        return {k:i for i,k in enumerate(self.l_node)}

    @lazy_property
    def delta_degree(self):
        in_ =  self.adjacency_topdown
        out_ = self.adjacency_bottomup
        return { n:(len(in_[n]) - len(out_[n])) for n in self.l_node}

    @lazy_property
    def constrain_matrix(self):

        od = self.order
        n_node, n_edge = map(len, (self.l_node,self.weigh))

        c = np.zeros( (n_edge, n_node), float)
        for idx, (i,o) in enumerate(self.weigh):
                c[idx,od[i] ] = 1
                c[idx,od[o] ] = -1

        return c
    
    @lazy_property
    def constrain_vector(self):
        return np.array([-k for k in self.weigh.values()],dtype=float)

    @lazy_property
    def objective_vector(self):
        # Float is required by CVXopt
        return np.array([-1*self.delta_degree[i] for i in self.order],dtype=float) # -1 <=> Minizise

    
    def add_pos_contrain(self,A,b):

        # Number of variable
        n = A.shape[1]

        # Add the contrain in the A matrix
        A_pos = np.identity(n) * -1
        A_tot = np.concatenate((A, A_pos), axis=0)

        # Add that they should positive
        b_pos = np.zeros(n)
        b_tot = np.concatenate( (b,b_pos) )

        return A_tot, b_tot

    @lazy_property
    def firing_buffered(self):

        from cvxopt import matrix, solvers
        A,b = map(matrix,self.add_pos_contrain(self.constrain_matrix,self.constrain_vector))
        c = matrix(self.objective_vector)

        #Integer solution use GLPK and turn off verbose output
        solvers.options['glpk'] = {'msg_lev' : 'GLP_MSG_OFF'}
        sol=solvers.lp(c,A,b, solver='glpk')

        sol['x_int'] = [int(round(i)) for i in sol['x'] ]
        assert all(i == f for i,f in zip(sol['x_int'],sol['x']))
        return dict(zip(self.order,sol['x_int'])) # Assume ordered dict

    @lazy_property
    def adjacency_buffer(self):
        b = self.firing_buffered
        w = self.weigh
        return { (i,o): (b[o] - b[i]) - w for (i,o),w in w.items() }

w = {('u1','u2'): 1,
     ('u2','u4'): 3,
     ('u2','u3'): 4,
     ('u3','u4'): 4,
     ('u4','u5'): 1,
     ('u1','u5'): 20,
     ('u3','u5'): 5}

d_b= B(w).adjacency_buffer
from operator import itemgetter
print (dict(filter(itemgetter(1), d_b.items())))
