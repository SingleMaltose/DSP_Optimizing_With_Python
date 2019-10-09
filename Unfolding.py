# -*- coding: utf-8 -*-
"""
Created on Sun Sep 29 11:03:06 2019

This is the straight-forward implementation of algorithms in Chaptr 5, "Unfolding"
of "VLSI Digital Signal Processing Systems".

Unfolding creates a new graph that equals to the original graph with multiple 
iterations, based on the unfolding factor J.It is also called loop unfolding.

For a unforling graph, every latch is J-slowed.

Unfolding can shorten Sample_Cycle to achieve Loop_Bound_Cycle, when retiming or 
pipline is futile. 

Function "unfolding" implementated the unfolding algorithm, which transformed the
original 2-dimensional graph into 4-dimensional graph. The first 2 dimensions both
equal to J, and the last 2 dimensions is same as the original graph, since the 
unfolding graph is a J-time-copy of original graph.

For example, assuming the original graph is G and the unfolding graph is G',
the G'(i,j,m,n) means the path from the m-th node of i-th copy of G to the n-th
node of j-th copy of G.

The first 2 dimension is removable for calculating Loop_Bound or Critic_Path.
@author: Singlemaltmaltose
"""

import numpy as np
import math

class Graph(object):
    
    def __init__(self, graph_matrix, node_calc_time):
        
        
        # graph_matrix是图graph的路径构成的矩阵M, M(i,j)为从第i个节点为尾，指向第j
        # 个节点的, 不经过其它节点的箭头的路径对应的权重。若不存在i到j的箭头，则M(i,j)=None
        self.graph_matrix = graph_matrix
        
        # node_calc_time是图graph的节点的计算时间，形式是一维np.array。长度等于graph_matrix
        # 的长度
        self.node_calc_time = node_calc_time
        
    
    def print_graph(self):
        print(self.graph_matrix)
    
    def if_looped_graph(self):
        # 判定graph_matrix是否包含环路
        g_mat = self.graph_matrix
        w = g_mat.shape[1]
        g_mat_reference = np.zeros([w,w])
        while np.not_equal(g_mat_reference,g_mat).any():
            g_mat_reference=np.copy(g_mat)
            for i in range(w):  
                if (g_mat[:,i] == None).all():
                    g_mat[i,:] = None
        if (g_mat==None).all():
            return g_mat, False
        else:
            return g_mat, True
    
    def unfolding(self,J):
        
        g_mat = np.copy(self.graph_matrix)
        dim = g_mat.shape[1]
        Unfolding_graph = float('inf')*np.ones([J,J,dim,dim])
        for i in range(J):
            for j in range(dim):
                for k in range(dim):
                    if g_mat[j,k]!= None:
                        w = g_mat[j,k]
                        Unfolding_graph[i,(i+w)%J,j,k] = math.floor((i+w)/J)
        
        print(Unfolding_graph)
        return Unfolding_graph

        
    
        
    
