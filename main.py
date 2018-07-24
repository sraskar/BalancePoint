#!/usr/bin/env python
from balence import Swing
if __name__ == '__main__':

    # Weighted graph to balence
    #Gao fig 7
    w = {
        ('u1', 'u2'): 1,
        ('u2', 'u4'): 3,
        ('u2', 'u3'): 4,
        ('u3', 'u4'): 4,
        ('u4', 'u5'): 1,
        ('u1', 'u5'): 20,
        ('u3', 'u5'): 5
    }

    #Appl fig 1
    w = {
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

    max_b = 92

    print(f'Max buffer alowed: {max_b}')
    l_buffer_updated, diff_delay = Swing(w).constrained_edge_buffer(max_b)  #adjacency_buffered
    print(f'list buffer: {l_buffer_updated}')
    print(f'number of buffer used: {sum(l_buffer_updated.values())}')
    print(f'max diff delay: {diff_delay}')



