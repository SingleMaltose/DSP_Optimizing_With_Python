# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:18:53 2019

@author: x1carbon
"""
import numpy as np

class graph(object):
    def __init__(self,node_calc_time, graph_matrix, graph_delay_matrix):
        self.node_calc_time = node_calc_time
        
        # graph_matrix是以计算单元为节点的图对应的矩阵M, M(i,j)为从第i个计算单元为尾
        # 指向第j个计算单元的箭头的延迟数，若不存在i到j的箭头，则M(i,j)=-1
        self.graph_matrix = graph_matrix
        
        # graph_delay_matrix是以延迟单元为节点的图对应的矩阵graph_delay_matrix, 
        # graph_delay_matrix(i,j)为从第i个延迟单元为尾，指向第j个延迟单元的箭头，
        # 且中间不包含其他延迟单元的路径中，最大的总计算时长，
        # 若不存在i到j且中间无其他延迟单元的箭头，则graph_delay_matrix(i,j)=-1
        self.graph_delay_matrix = graph_delay_matrix
    
    def print_graph(self):
        print(self.graph_matrix)
        
    def if_looped_graph(self):
        g_mat = self.graph_matrix
        w = g_mat.shape[1]
        for i in range(w):
            if np.max(g_mat[:,i])<0:
                g_mat[i,:] = -1
        if np.max(g_mat)<0:
            return g_mat,False
        else:
            return g_mat,True
        
    def calc_loop_bound_LPM(self):
        graph_delay_matrix = self.graph_delay_matrix
        w = graph_delay_matrix.shape[1]
        L = [graph_delay_matrix]
        for i in range(1,w):
            L_temp = np.zeros([w,w],dtype=int)
            for j in range(w):
                for k in range(w):
                    L_temp[j][k] = -1
                    for m in range(w):
                        if L[0][j][m] != -1 and L[i-1][m][k] != -1:
                            if L[0][j][m] + L[i-1][m][k] > L_temp[j][k]:
                                L_temp[j][k] = L[0][j][m] + L[i-1][m][k]
            L.append(L_temp)
        bound_candidates=[]
        for i in range(w):
            bound_candidates += list(np.diag(L[i])/(i+1))
        loop_bound = max(bound_candidates)
        return loop_bound
    
    def calc_loop_bound_MCM(self):
        graph_delay_matrix = self.graph_delay_matrix
        w = graph_delay_matrix.shape[1]
        f0 = 100000 * np.ones([w,1],dtype=int)
        f0[0] = 0
        f = [f0]
        for i in range(w):
            f_temp = np.zeros([w,1],dtype=int)
            for j in range(w):
                f_temp[j] = int(min([f[i][k]-graph_delay_matrix[k,j] \
                                    for k in range(w) \
                                    if graph_delay_matrix[k,j] >=0]))
                if f_temp[j] > 10000:
                    f_temp[j] = 100000;
            f.append(f_temp)
        bound_candidates=np.concatenate([(f[w]-f[i])/(w-i) for i in range(w)],axis=1)
        loop_bound = -min(bound_candidates[range(w),np.argmax(bound_candidates,axis=1)])
        return loop_bound
        
    

                
        
        