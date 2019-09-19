# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 16:07:29 2019

This is the straight-forward implementation of finding shortest path in graph
based on Bellman_Ford algorithm and Floyd_Warshall algorithm

The B-F algorithm is O(n^3) complexity for one-to-all node pair's shortest path
computation, and for all-to-all node pair it requires O(n^4) complexity

The F-W algorithm is O(n^3) complexity for all-to-all node pair's shortest path
computation, while consume n times of storage.  

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
    
    def SPF_with_Bellman_Ford(self,U):
        # Bellman_Ford算法用于计算从U点到图中任意一点的最短路径，U应该不大于w
        g_mat = np.copy(self.graph_matrix)
        w = g_mat.shape[1]   
        
        #将输入路径矩阵中的None转化为inf
        for i in range(w):
            for j in range(w):
                if g_mat[i,j] == None:
                    g_mat[i,j] = float('inf')

        # 构建r矩阵
        r = np.zeros([w,w-1])
        # Bellman_Ford算法
        for k in range(w):
            if k != U:
                r[k,0] = g_mat[U,k]
        for k in range(w-2):
            for v in range(w):
                r[v,k+1] = r[v,k]
                for z in range(w):
                    if r[v,k+1] > r[z,k] + g_mat[z,v]:
                        r[v,k+1] = r[z,k] + g_mat[z,v]
        for v in range(w):
            for k in range(w):
                if r[v,w-2] > r[k,w-2] + g_mat[k,v]:
                    print("包含负回路")
                    return
        print("不包含负回路")
        return r
    
                    
    def SPF_with_Floyd_Warshall(self,graph_matrix=np.array([])):
        # Floyd_Warshall算法用于计算图中所有点对的最短路径, 该算法共构建w+1
        # 个矩阵，并由第w+1个矩阵给出最短路径。对于第w+1个矩阵R，R[j,k]即第
        # j个点到第k个点的最短路径
            
        if graph_matrix.size != 0:
            g_mat = np.copy(graph_matrix)
        else:
            g_mat = np.copy(self.graph_matrix)
        w = g_mat.shape[1]   
        
        #将输入路径矩阵中的None转化为inf
        for i in range(w):
            for j in range(w):
                if g_mat[i,j] == None:
                    g_mat[i,j] = float('inf')      
        R = []
        for i in range(w+1):
            R.append(np.zeros([w,w]))
        for V in range(w):
            for U in range(w):
                R[0][U,V] = g_mat[U,V]
        for k in range(w):
            for V in range(w):
                for U in range(w):
                    R[k+1][U,V] = R[k][U,V]
                    if R[k+1][U,V] > R[k][U,k] + R[k][k,V]:
                        R[k+1][U,V] = R[k][U,k] + R[k][k,V]
        
        for k in range(w):
            for U in range(w):
                if R[k][U,U] < 0:
                    print("包含负回路")
                    return
        print("不包含负回路")
        return R
    
    def retiming_with_minimized_cycle(self,target_cycle):
        g_mat = np.copy(self.graph_matrix)
        node = np.copy(self.node_calc_time)
        
        t_max = np.max(node)
        node_number = node.shape[0]
        M = t_max*node_number
        g_mat_new = np.copy(self.graph_matrix)
        for i in range(g_mat.shape[0]):
            for j in range(g_mat.shape[1]):
                if g_mat_new[i,j] != None:
                    g_mat_new[i,j] = M*g_mat_new[i,j] - node[i]
        
        R = self.SPF_with_Floyd_Warshall(graph_matrix = g_mat_new)
        
        S_UV = R[-1]
        W_UV = np.zeros([g_mat.shape[0],g_mat.shape[1]])
        D_UV = np.diag(node)
        for i in range(g_mat.shape[0]):
            for j in range(g_mat.shape[1]):
                if i!=j:
                    W_UV[i,j] = math.ceil(S_UV[i,j]/M)
                    D_UV[i,j] = M*W_UV[i,j] - S_UV[i,j] + node[j]
            
        inequation_matrix = float('inf')*np.ones([g_mat.shape[0]+1,g_mat.shape[1]+1])
        for i in range(g_mat.shape[0]):
            for j in range(g_mat.shape[1]):
                if g_mat[i,j]!=None:
                    inequation_matrix[j,i] = g_mat[i,j]
        
        for i in range(D_UV.shape[0]):
            for j in range(D_UV.shape[1]):
                if D_UV[i,j] > target_cycle:
                    if W_UV[i,j]-1 <= inequation_matrix[j,i]:
                        inequation_matrix[j,i] = W_UV[i,j] - 1
        
        inequation_matrix[-1,:] = 0       
        Retiming_result = self.SPF_with_Floyd_Warshall(graph_matrix = inequation_matrix)
        return Retiming_result[-1][-1,:-1]
        
        
        
                
                