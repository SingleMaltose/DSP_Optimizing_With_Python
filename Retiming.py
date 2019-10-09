# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 16:07:29 2019

This is the straight-forward implementation of algorithms in Chaptr 4, "Retiming"
of "VLSI Digital Signal Processing Systems".

Retiming aims to rearrange latches to shorten critical path or reduce latch numbers.
This module achieved critical path optimization with target cycle, which is based on 
Bellman_Ford algorithm and Floyd_Warshall algorithm

The B-F algorithm is O(n^3) complexity for one-to-all node pair's shortest path
computation, and for all-to-all node pair it requires O(n^4) complexity

The F-W algorithm is O(n^3) complexity for all-to-all node pair's shortest path
computation, while consume n times of storage.  

Amendment1: 2019/10/09

1.Create a new function "calculation of S_W_D" to eperate calculation of S,W,D from 
"retiming_with_minimized_cycle" function.

2.Create a new function "calc_critic_path" and "build_critic_path"to calculate 
critic_path and its calculating time(latency or cycle).

@author: Singlemaltmaltose
"""
import numpy as np
import math
import copy

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
                    print("包含负回路,该优化目标不可实现")
                    return [np.zeros([w,w])]
        print("不包含负回路")
        return R
    
    def calculation_of_S_W_D(self):
        # 构建DSP-graph的S矩阵，W矩阵，D矩阵
        # S矩阵为原始graph转换而成的新graph，G'的任意两个节点的最短路径矩阵
        # W矩阵为任意两个节点任意路径上寄存器最少的数目
        # D矩阵为任意两个节点以对应W为权重的，所有路径中最长的计算时间
        g_mat = np.copy(self.graph_matrix)
        node = np.copy(self.node_calc_time)
        
        t_max = np.max(node)
        node_number = node.shape[0]
        M = t_max*node_number
        # 构建G'
        g_mat_new = np.copy(self.graph_matrix)
        for i in range(g_mat.shape[0]):
            for j in range(g_mat.shape[1]):
                if g_mat_new[i,j] != None:
                    g_mat_new[i,j] = M*g_mat_new[i,j] - node[i]
        
        # 使用最短路径问题求解方法求解G'的最短路径矩阵
        R = self.SPF_with_Floyd_Warshall(graph_matrix = g_mat_new)
        
        self.S_UV = R[-1]
        print(self.S_UV)
        # 构建W，D
        self.W_UV = np.zeros([g_mat.shape[0],g_mat.shape[1]])
        self.D_UV = np.diag(node)+np.zeros([node.shape[0],node.shape[0]])
        for i in range(g_mat.shape[0]):
            for j in range(g_mat.shape[1]):
                if i!=j:
                    if self.S_UV[i,j] != float('inf'):
                        self.W_UV[i,j] = math.ceil(self.S_UV[i,j]/M)
                        self.D_UV[i,j] = M*self.W_UV[i,j] - self.S_UV[i,j] + node[j]
                    else:
                        self.W_UV[i,j] = float('inf')
                        self.D_UV[i,j] = float('inf')

    
    def retiming_with_minimized_cycle(self,target_cycle):
        # 对原DSP-graph进行重定时，通过重新安排寄存器的位置，在不影响电路功能的情况下，
        # 实现最小化关键路径，从而降低时钟周期
        g_mat = np.copy(self.graph_matrix)
        inequation_matrix = float('inf')*np.ones([g_mat.shape[0]+1,g_mat.shape[1]+1])
        for i in range(g_mat.shape[0]):
            for j in range(g_mat.shape[1]):
                if g_mat[i,j]!=None:
                    inequation_matrix[j,i] = g_mat[i,j]
        
        # 重定时的电路满足以下约束：
        # 可行性约束：重定时后的任意边延时不为负数，即r(U)-r(V)<=w(e)
        # 关键路径约束：对于计算时间超过目标周期的路径，使其重定时后的路径延迟大于1，
        # 即r(U)-r(V)<=W(U,V)-1
        for i in range(self.D_UV.shape[0]):
            for j in range(self.D_UV.shape[1]):
                if self.D_UV[i,j] > target_cycle:
                    if self.W_UV[i,j]-1 <= inequation_matrix[j,i]:
                        inequation_matrix[j,i] = self.W_UV[i,j] - 1
        
        inequation_matrix[-1,:] = 0
        # 根据约束条件得出的不等式，构建图，并使用最短路径法进行求解
        Retiming_result = self.SPF_with_Floyd_Warshall(graph_matrix = inequation_matrix)
        

        
        return Retiming_result[-1][-1,:-1]
        
    def calc_critic_path(self):
        # 求解DSP-graph的关键路径长度，并标出关键路径经由的具体单元。
        # 求解方法：1.根据W矩阵确定所有路径延迟为0的路径，并记录其头尾节点
        # 2.基于D矩阵，计算路径延迟为0的路径的计算时间，并取其中的最大值，即为关键路径延迟
        # 3.对于路径延迟等于关键路径延迟的路径，确定为关键路径，并记录这些关键路径的头尾节点
        # 4.根据关键路径的头尾节点对，遍历其对应的所有0延迟路径，并计算该路径计算时间
        # 5.对于计算时间等于关键路径延迟的路径，即为关键路径，记录进self并返回
        
        zero_path_cord = np.where(self.W_UV == 0)
        zero_path_cord_pair = list(zip(list(zero_path_cord[0]),list(zero_path_cord[1])))
        zero_path_calc_time = self.D_UV[zero_path_cord]        
        path_to_calc_time_dict = dict(zip(zero_path_cord_pair,zero_path_calc_time))
        self.critic_path_calc_time = max(zero_path_calc_time)
        #print(path_to_calc_time_dict)
        self.critic_path_cord = []
        for key in path_to_calc_time_dict.keys():
            if path_to_calc_time_dict[key] == self.critic_path_calc_time:
                self.critic_path_cord.append(key)
        for cord in self.critic_path_cord:
            self.build_critic_path(cord)
        return self.critic_path_calc_time, self.critic_path_cord, self.critic_path_list
    
        
    def build_critic_path(self,cord_pair):
        # 对应求解关键路径的第4、5步，使用递归的方法遍历路径延时为0，计算时间为关键路
        # 径延迟的路径首尾对，并记录对应的路径具体节点到self.critic_path_list中
        self.zero_path.append(cord_pair[0])
        if cord_pair[0] == cord_pair[1]:
            if sum([self.node_calc_time[x] for x in self.zero_path]) == self.critic_path_calc_time:
                critic_path = copy.deepcopy(self.zero_path)
                self.critic_path_list.append(critic_path)
            #print(self.critic_path_list)
            #print(self.zero_path)
        for index,i in enumerate(self.graph_matrix[cord_pair[0],:]):
            if i == 0:
                self.build_critic_path(cord_pair=[index,cord_pair[1]])
            if index == self.graph_matrix.shape[0]-1:
                self.zero_path.pop()
        
        
        
                
                
