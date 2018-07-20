#!/usr/bin/env python

import numpy as np
from collections import defaultdict
from itertools import chain
from operator import itemgetter

from functools import lru_cache
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
    def adjacency_topdown(self):
        d = defaultdict(list)
        for i,o in self.weigh:
            d[i].append(o)
        return d

    @lazy_property
    def adjacency_bottomup(self):
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
    def leaf(self):
        '''Assume only one leaf for now'''
        return (self.l_node - self.adjacency_topdown.keys()).pop()

    #  __
    # /__  _.  _
    # \_| (_| (_)
    #

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
        '''Extend the matrix A,b with positif constrain
        '''
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
        '''This function use GAO algorithm. It will generated the optimal buffer 
           need to minimze the graph
        '''
        
        from cvxopt import matrix, solvers
        # We add the constrain that each objectif need to me postiif
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

    #                        
    # |\ |  _  ._ _   _.  _| 
    # | \| (_) | | | (_| (_| 
    #                        

    @lru_cache(None)
    def path(self, cur):
        '''Compute all the path form cur to the root of the DAG. 
            Assume a unique root'''
        d =  self.adjacency_bottomup
        
        if d[cur]:
            it =  chain.from_iterable(self.path(p) for p in d[cur])
        else:
            it = [ [] ]

        return [ k + [cur] for k in it ]

    @lazy_property
    def weigted_path(self):
        '''
        Compute for all the path ( of the form (e1->e2),(e2->e3) )
        and the sum of th weigh  from the leaf to the read 
        '''
        l_edge= [ list(pairwise(p)) for p in self.path(self.leaf) ]
        return sorted([ (sum(map(self.weigh.get, e)), e) for e in l_edge ],key=itemgetter(0))


    def bb(self, x, A, C, max_b):
        '''This function compute the function to minimize.
            f(x) = max(abs(C-sum(A*x)))      
            where sum(x) > max_b
        ''' 

        #Maybe add the constrain that each path should be >0

        dim = x.get_n()
        b = [x.get_coord(i) for i in range(dim)]

        m = np.multiply(A,b)
        s = np.sum(m,axis=1)
        f = max(abs(C - s))
        x.set_bb_output(0, f)

        c = sum(b) - max_b
        x.set_bb_output(1,c)
        return 1

    def nomad(self,max_b):
        '''
        We will optimize min(max(abs(C-p))) where p are the total weigh of the graph
        under the constrain the sum(e) < max_b
        '''
        import PyNomad
        
        # Initialization
        wp = self.weigted_path
        critical_path = self.weigted_path[-1][0]

        optimal_buffer = {k:v for k,v in self.adjacency_buffer.items() if v}
 
        # In order to reduce the complexity of the algorithm,
        # We will optimized only the edges that Gao's algotihm touched.
        get_unik = lambda l_edge: [e for e in l_edge if e in optimal_buffer]
        wp_filter = list(filter(itemgetter(1),((w, get_unik(p)) for w,p in wp )))
        l_edge_unique = set(chain.from_iterable(map(itemgetter(1), wp_filter)))

        # We construct the matrix of edges that we want to optimzation
        A = np.array([  [e in l_edge for e in l_edge_unique] for _, l_edge in wp_filter ], dtype=int)
        # This is the constrain for each edge: They cannot grow bigger than the critical path
        C = np.array([critical_path-w for w,_ in wp_filter],dtype=int)

        # Param for nomad

        # Function
        from functools import partial
        bb = partial(self.bb,max_b=max_b,A=A,C=C)
        
        # Bound and initial value
        lb = [0] * len(l_edge_unique) 
        ub = [optimal_buffer[e] for e in l_edge_unique]
        x0 = lb

        params = ['BB_OUTPUT_TYPE OBJ EB',
                  'MAX_BB_EVAL 100',
                  'LOWER_BOUND * 0',
                  'BB_INPUT_TYPE * I',
                  'DISPLAY_DEGREE 0']

        x_return, f_return, h_return, nb_evals, nb_iters, stopflag = PyNomad.optimize(bb,x0,lb,ub,params)

        l_buffer_updated = {k:int(v) for k,v in zip(l_edge_unique,x_return)}
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

max_b = 94 

print (f'Max buffer alowed: {max_b}')
l_buffer_updated, diff_delay = B(w).nomad(max_b) #adjacency_buffered
print (f'list buffer: {l_buffer_updated}')
print (f'number of buffer used: {sum(l_buffer_updated.values())}')
print (f'max diff delay: {diff_delay}')

#from operator import itemgetter
#print (dict(filter(itemgetter(1), d_b.items())))
