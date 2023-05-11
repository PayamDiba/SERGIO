import numpy as np
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import networkx as nx
import random
import scipy
import itertools

from collections import deque
import itertools
import warnings
import scipy
import argparse
from functools import lru_cache
from multiprocessing import Pool
import sys
import logging
from numba import jit
import os
#from correlation_cavity2 import correlation_parallel,correct_correlation
logging.basicConfig(level=logging.WARNING)
@jit(nopython=True)
def cavity_grad_numba(P_loc,inter,T,theta,K,J0):
    '''
    Dynamic programming is run using iterative calls. See documentation of cavity_single_numba
    '''

    pos = inter.copy()
    neg = inter.copy()
    pos[pos<0]=0#apply theta function
    neg[neg>0]=0
    m = np.zeros((K+1,K+1))

    h_tilde =np.arange(np.sum(neg),np.sum(pos)+1)#the list of possible values of the local field we get
    for ind in range(len(h_tilde)):
        m[K,ind]=1/4/T/(np.cosh((h_tilde[ind]-theta)*J0/2/T))**2#fill the bottom row of the matrix first 
    offset = np.sum(neg)#this is the offset to map \tilde{h} to index of m matrix
    for l,low,top in list(zip(np.arange(1,K),np.cumsum(neg)[:-1],np.cumsum(pos)[:-1]+1))[::-1]:
        #print(l)
        for h in np.arange(low-offset,top-offset):
            m[l,h] = P_loc[l]*m[l+1,h+inter[l]]+(1-P_loc[l])*m[l+1,h]

        #print(m)
    return P_loc[0]*m[1,inter[0]-offset]+(1-P_loc[0])*m[1,-offset]

@jit(nopython=True)
def cavity_single_numba(P_loc,inter,T,theta,K,J0):
    '''
    Dynamic programming is run using iterative calls. It creates a matrix for this task. The row i of the matrix corresponds to taking the average over the first i-1 nodes.
    The columns of the matrix are mappped according to the K+1 values of the partial sum of the interaction in ascending order. For example if inter = [+1,-1,+1,+1] then the
    the column index  maps to [-1,0,1,2,3].
    i\j  -1,0,1,2,3

    1
    2
    3
    4
    5  
    '''

    pos = inter.copy()
    neg = inter.copy()
    pos[pos<0]=0#apply theta function
    neg[neg>0]=0
    m = np.zeros((K+1,K+1))

    h_tilde =np.arange(np.sum(neg),np.sum(pos)+1)#the list of possible values of the local field we get
    for ind in range(len(h_tilde)):
        m[K,ind]=1/2*(1+np.tanh((h_tilde[ind]-theta)*J0/2/T))#fill the bottom row of the matrix first 
    offset = np.sum(neg)#this is the offset to map \tilde{h} to index of m matrix
    for l,low,top in list(zip(np.arange(1,K),np.cumsum(neg)[:-1],np.cumsum(pos)[:-1]+1))[::-1]:
        #print(l)
        for h in np.arange(low-offset,top-offset):
            m[l,h] = P_loc[l]*m[l+1,h+inter[l]]+(1-P_loc[l])*m[l+1,h]

        #print(top)
    return P_loc[0]*m[1,inter[0]-offset]+(1-P_loc[0])*m[1,-offset]
def optimise_caller(J_transpose,T,theta,N,precision=1e-4,J0 = 'auto'):
    '''J here is a 1d vector with N**2 elements. scipy optimise works with 1d vectors'''
    js = [np.arange(N)[el!=0] for el in J_transpose]# list of list, structure is [el[i]] where el[i]
    # is the list of  predecessors of gene i ( the index)
    interaction = [el[el!=0] for el in J_transpose]  # list of list, structure is [el[i]]
    # where el[i] is the list of  predecessors of gene i (interaction strength with sign)
    Ks = np.array([len(neigh) for neigh in js])  # in degree of each gene
    max_outdegree = max(Ks)
    max_recursions = int((max_outdegree + 1) * (max_outdegree + 2) / 2)
    '''
    if max_recursions > sys.getrecursionlimit():
        print("Warning! maximum degree larger than default recursion limit, I 'll update recursion limit to",
              max_recursions)
        sys.setrecursionlimit(max_recursions)
    '''
    if J0 =='auto':
        avg_degree = np.mean(Ks)
        J0 = 1/ np.sqrt(avg_degree)
    P_thr,grad =  cavity(np.random.rand(N), js, T, interaction, N, Ks, theta,J0,precision)
    return P_thr,grad
def cavity(P, js, T, interaction, N, Ks, theta,J0, precision=1e-4, max_iter=50):
    """
    This runs the dynamical cavity without recursive calls but using iterative calls. It creates instead a matrix. This works only if couplings are in the form  \pm J.
    If couplings are in a different form, use it cavity_general
     It computes the node activation probability for a  directed network.
    :param P_init: list of floats of length N
    :param T: float
    :param js: list of list, structure is [el[i]] where el[i] is the list of  predecessors of gene i ( the index)
    :param interaction:  list of list, structure is [el[i] for i in range(N)]
            where el[i] is the list of  predecessors of gene i (interaction strength with sign)
    :param theta: float (in units of 1/sqrt(<K>))
    :param max_iter: int
    :param precision: float
    :return: P_new it is a  list of dimensions N which contains the probability of active state for each gene.
    ----NOTES------
    In order to help storing, couplings are taken to be +-1, at the end the local field is rescaled by 1/sqrt(<|J_{ij}|>)
    Even though code runs for any directed network, results  are exact for fully asymmetric networks only.
    """

    if T == 0:
        return cavity_zero_T(P, js, interaction, N, Ks, theta)
    avg_degree = np.mean(Ks)
    for count in range(max_iter):
        P_new = np.zeros(N)
        for i in range(N):
            j = np.array(js[i],dtype = int)
            bias = 0
            K = len(j)
            if len(j) ==0:
                P_new[i]=0.5
            else:
                #print(count,j)
                inter = np.array(interaction[i],dtype = int)
                P_new[i]=cavity_single_numba(P[j],np.array(inter),T,theta,K,J0)
        if max(np.abs(np.array(P) - np.array(P_new))) < precision:
            P = P_new
            #print('finishing after', count, 'iterations')
            break
        if count == max_iter:
            print("Maximum number of repetition reached, but target  precision has not been reached. Precision reached is "+str(max(np.abs(np.array(P) - np.array(P_new)))))

        P = np.array(P_new)
    P = np.array(P)
    return P,[]
    #now do the gradient
    #grad = {'row':np.zeros(len(J.data)),'col':np.zeros(len(J.data)),'data':np.zeros(len(J.data))}#initialise the gradient
    grad = defaultdict(list) #gradient matrix matches the shape of  J 
    for i in range(N):
        bias = 0
        K = Ks[i]
        if K ==0:
            grad[i]=[]
        if K == 1:
            grad[i]= [1/4/T/(np.cosh((interaction[i][0])*J0/2/T))**2]
        else:
            j = deque(js[i])#use deque to rotate over index and loop over all neighbours. It is a trick to give
            #the set of nodes neighbouring i excluding js[i][ind]
            inter = deque(interaction[i])
            for ind in range(K):
                neigh = np.array(list(itertools.islice(j, 1, len(j))))
                coup = np.array(list(itertools.islice(inter, 1, len(inter))))
                a = cavity_grad_numba(P[neigh],coup,T,theta-inter[0],K-1,J0)
                grad[i]+=[a]
                j.rotate(-1)
                inter.rotate(-1)
    return P,grad


def binarise_matrix(A,thr = 0.1):
    A[A>thr]=1
    A[A<-thr]=-1
    A[np.abs(A)<1]=0
    return A
def draw_net(G,offset = 0.05,**kwargs):
    c =np.array([ c['weight'] for a,b,c in list(G.edges(data=True))])
    edge_color=np.where(c>0,'green','red')
    #nx.draw_circular(G,with_labels=True,edge_color=edge_color,alpha = 0.7,arrowsize = 15,**kwargs)
    nodePos = nx.circular_layout(G)
    nx.draw_networkx_nodes(G,pos=nodePos, label=True)
    nx.draw_networkx_labels(G,pos = nodePos)
    edges = list(G.edges)
    bi_edges = [(a,b) for a,b in edges  if (b,a)in edges]
    non_bi_edges =list(set(edges)-set(bi_edges))
    a,b,c  =zip(*[ (a,b,c['weight']) for a,b,c in list(G.edges(data=True))])
    dic = {(start,stop):col for start,stop,col in zip(a,b,c)}
    #nx.draw_networkx_edges(G,pos = nodePos,edgelist=edges[weight>0],edge_color= 'green',arrowsize = 15,)
    non_bi_weight = np.array([dic[endpoints] for endpoints in non_bi_edges ])
    bi_weight = np.array([dic[endpoints] for endpoints in bi_edges ])
    nx.draw_networkx_edges(G,pos = nodePos,edgelist=np.array(non_bi_edges)[non_bi_weight>0],edge_color= 'green',arrowsize = 20,)
    nx.draw_networkx_edges(G,pos = nodePos,edgelist=np.array(non_bi_edges)[non_bi_weight<0],edge_color= 'red',arrowsize = 10,arrowstyle='-[',alpha = 0.6)
    #draw bi-directional links parallel one another such that they do not overlap
    unique_bi_edges = []
    for start,stop in bi_edges:
        if (stop, start) not in unique_bi_edges:
            unique_bi_edges+=[(start,stop)]
    offset = 0.05
    new_nodePos={}
    for start,end in unique_bi_edges:
        new_nodePos[start] = nodePos[start]-[0,offset]
        new_nodePos[end] = nodePos[end]-[0,offset]
    nx.draw_networkx_edges(G,pos = new_nodePos,edgelist=np.array(bi_edges)[bi_weight>0],edge_color= 'green',arrowsize = 20,)
    new_nodePos={}
    for start,end in unique_bi_edges:
        new_nodePos[start] = nodePos[start]+[0,offset]
        new_nodePos[end] = nodePos[end]+[0,offset]

    nx.draw_networkx_edges(G,pos = new_nodePos,edgelist=np.array(bi_edges)[bi_weight<0],edge_color= 'red',arrowsize = 10,arrowstyle='-[',alpha = 0.6,width = 1.5)


def belief_propagation_pert(P_i,P_perturb,P_perturb_crispa, T,  theta,lam,beta, J0=1,precision=1e-4, max_iter=100, using_pert = True):
    """
    Belief propagation using pertubation. This implementation uses all the factor nodes. This implementation has shown to present instabilities in convergence. It is advisable to use the belief_propagation_pert_short.
    :param P_i: list of floats of length N. They are the target Probabilities
    :param P_w: 3d array of shape (3,N,N)
    :param T: float
    :param theta: float (in units of 1/sqrt(<K>))
    :param max_iter: int
    :param precision: float
    :return: P_w_new it is a 3d array of shape (3,N,N) which contains the probability of positive,negative,or zero weight for every link.
    ----NOTES------
    In order to help storing, couplings are taken to be +-1
    Even though code runs for any directed network, results  are exact for fully asymmetric networks only.
    """
    N = len(P_i)
    '''
    max_recursions = int((N + 1) * (N + 2) / 2)
    if max_recursions > sys.getrecursionlimit():
        print("Warning! maximum degree larger than default recursion limit, I 'll update recursion limit to",
              max_recursions)
        sys.setrecursionlimit(max_recursions)
    '''
    def update_belief_pert(P_w_i,P_i,i,j,w,theta,J0=1):
        '''
        P_w_i is P_w[:,i,:]. Because for every i, only links pointing to it matters
        '''
        bias = 0
        @lru_cache(maxsize=None)
        def recursion(bias, l):
            '''at the end of recursions, bias is the local field'''

            if (l == N):
                
                bias = (bias - theta) *J0
                mismatch = (0.5 + 0.5 *np.tanh((bias-theta) / 2 / T)-P_i[i])**2#difference betwen theoretical and target prob.
                return np.exp(-beta*mismatch)
            elif (l == i):
                #avoid self interaction
                return recursion(bias,l+1)
            elif (l == j):
                return recursion(bias+w,l+1)*P_i[l]+recursion(bias,l+1)*(1-P_i[l])
            activing = P_i[l] * recursion(bias + 1, l + 1)*P_w_i[0,l]  #  node l  active with prob. P_i[l] and link w_[il] positive
            inhibiting = P_i[l] * recursion(bias - 1, l + 1)*P_w_i[2,l]  #  node l  active with prob. P_i[l] and link w_[il] negative
            inactive =  (1-  P_i[l]+P_i[l]*P_w_i[1,l])*recursion(bias , l + 1)# node l sends no contribution, either because it is active, or because the link w_[il] is zero
            return activing+inhibiting+inactive

        N = len(P_i)
        result = recursion(bias, 0)
        recursion.cache_clear()
        return result

    if using_pert:
        M = len(P_perturb)+len(P_perturb_crispa)+1#number of experimental conditions
        input_data = np.concatenate((np.expand_dims(P_i,0),P_perturb,P_perturb_crispa),axis = 0)#iterator on the different experimental conditions
    else:
        raise ValueError('Cavity iteration not implemented yet')
        M = 1
        input_data = [P_i]
    #initialise the P_w
    #rho= np.random.rand(3,M,N)
    P_w = np.zeros((3,N,N))#P_w[0] is P(w_ij=1),P_w[1] is P(w_ij=0),P_w[2] is P(w_ij=-1)
    store = np.empty(N,dtype = object)
    error = []#np.zeros((M,N))
    for i in range(N):
        rho_nu = np.random.rand(3,M,N)#the rho for the cases w =1,0,-1
        rho_nu = rho_nu/rho_nu.sum(axis = 0)#normalise rho_nu
        rho_nu_old = rho_nu.copy()
        P_w_nu = rho_nu.copy()        
        P_w_nu[np.array([0,2])]=P_w_nu[np.array([0,2])]*np.exp(-lam)#penality term for the w = \pm 1     
        P_w_nu = P_w_nu/P_w_nu.sum(axis = 0)#normalise P_w_nu
        error = []
        for count in range(max_iter):         
            for nu,P_i_nu in enumerate(input_data):
                if (nu == i+1)|(nu == i+1+N):#(nu == i+1) corresponds to gene i removed in crispi, (nu == i+1+N) gene i in crispa 
                    rho_nu[:,nu,:]= 1/3#if the condition corresponds to crisp perturbation of i, data do not provide any info on the link, so we use equiprobability
                    continue
                for j in list(set(range(N))-{i}):
                    #this loop updates the rho_nu
                    for w_index in [0,1,2]:
                        w = np.array([1,0,-1])[w_index]#W_{ij}
                        rho_nu[w_index,nu,j] = update_belief_pert(P_w_nu[:,nu,:],P_i_nu,i,j,w,theta,J0)
                        #print(rho_nu[w_index,nu,j])
                    #if rho_nu[:,nu,j].sum(axis = 0) == 0:
                    #    raise ValueError('sum rho_nu = 0 at i,j,nu',i,j,nu)
            rho_nu[:,:,i] =np.array([[0,1,0]]*M).T#no self interaction for j=i
            rho_nu[rho_nu>1]=1#avoid  numeric explosion
            rho_nu[rho_nu<0]=0
            rho_nu = rho_nu/rho_nu.sum(axis = 0)#normalise rho
            loc_error = np.sum(np.abs(rho_nu - rho_nu_old),axis = 0)
            error += [loc_error]#sum over the w={1,0,-1}
            rho_nu_old = rho_nu.copy()
  
            if np.max(loc_error) < precision:
                logging.info('finishing after '+str(count)+' iterations')
                
                break
            if count == max_iter-1:
                logging.warning("Maximum number of repetition reached, but target  precision has not been reached for node "+str(i) +'\n with error'+str(np.max(loc_error)))

            for nu in range(M):
                #this loop updates the P_nu once the loop over the rho_nu is over
                rho_cav = np.concatenate((rho_nu[:,:nu,:],rho_nu[:,nu+1:,:]),axis= 1)#exclude condition nu, i.e. remove the data associated to condition nu
                P_w_nu[:,nu,:] = np.prod(rho_cav,axis = 1)#multiplication is performed on the rho_nu for all mu neq nu
                
            P_w_nu[[0,2]]=P_w_nu[[0,2]]*np.exp(-lam)#penality term for the w = \pm 1                            
            P_w_nu = P_w_nu/P_w_nu.sum(axis = 0)#normalise P_w_nu

        store[i] = error
        #END of cavity, now do the marginal inside the loop over i
        for j in range(N):
            P_w[:,i,j] = np.prod(rho_nu[:,:,j],axis = 1)#multiplication is performed on the rho_nu for all mu
            
    P_w[[0,2]]=P_w[[0,2]]*np.exp(-lam)#penality term for the w = \pm 1
    P_w = P_w/P_w.sum(axis = 0)#normalise P_w_nu

    '''
    if np.max(error) < precision:
        print('finishing after', count, 'iterations')
        break
    if count == max_iter-1:
        print("Maximum numberprint(l) of repetition reached, but target  precision has not been reached. ")
    '''
    return P_w,store
@jit(nopython=True)
def kl_divergence(p, q):
    '''p,q are arrays containing the probability of gene activation for all nodes.'''
    eps = 1e-12 #to cure for pathological cases where p or q are 0 or 1 
    return p * np.log((p + eps) / (q + eps)) + (1 - p) * np.log((1 - p + eps) / (1 - q + eps))
@jit(nopython=True)
def kl_symmetric(p, q):
    '''p,q are arrays containing the probability of gene activation for all nodes.'''
    eps = 1e-12 #to cure for pathological cases where p or q are 0 or 1 
    return (kl_divergence(p,q)+kl_divergence(q,p))/2

@jit(nopython=True)
def l2_norm(p,q):
    return (p-q)**2
    
@jit(nopython=True)
def update_belief_iterative(P_w_i,P_i,i,j,w,T,beta,theta,metric=kl_divergence,J0=1):
    '''
     Compute average of mismatch using iterative calls. If using with @jit, it is much faster than recursive calls, otherwise it is slightly slower.
    P_w_i is P_w[:,i,:]. Because for every i, only links pointing to it matters.
    Returns the scalar value
    '''
    N  = len(P_i)
    bias = 0
    #neg and pos below are preappended with 0, making them two vectors of lenght N+1.This helps to match the dimensions of matrix m
    if w<0:
        neg = np.concatenate((np.array([0]),-np.ones(j),np.array([w]),-np.ones(N-j-1)))#tracks negative links
        pos = np.concatenate((np.array([0]),np.ones(j),np.array([0]),np.ones(N-j-1)))#put 0 in correspondance of negative links
    else:
        neg = np.concatenate((np.array([0]),-np.ones(j),np.array([0]),-np.ones(N-j-1)))
        pos = np.concatenate((np.array([0]),np.ones(j),np.array([w]),np.ones(N-j-1)))
    neg[i+1]= 0 #exclude self interaction
    pos[i+1] = 0
    neg = neg.astype(np.int32)
    pos = pos.astype(np.int32)

    h_tilde =np.arange(np.sum(neg),np.sum(pos)+1) #possible values excluding nodes i,j
    m = np.zeros((N+1,len(h_tilde)))#one link ( i.e. j->i) is frozen
    for ind in range(len(h_tilde)):
        mismatch = metric(1/2*(1+np.tanh((h_tilde[ind]*J0-theta)/2/T)),P_i[i])
        m[-1,ind]=np.exp(-beta*mismatch)
    offset = np.sum(neg)#this is the offset to map \tilde{h} to index of m matrix
    for l,low,top in list(zip(np.arange(0,N),np.cumsum(neg),np.cumsum(pos)+1))[::-1]:
        for h in np.arange(low-offset,top-offset):#h here is the column index of m, not the h_tilde itself 
            
            if l == i:
                #avoid self interaction
                m[i,h]=m[i+1,h]
            elif (l == j):

                m[j,h]=  P_i[l] * m[l+1,h+w]+(1-P_i[l])* m[l+1,h]
            else:
                activating = P_i[l] * m[l+1,h+1]*P_w_i[0,l]  #  node l  active with prob. P_i[l] and link w_[il] positive
                inhibiting = P_i[l] * m[l+1,h-1]*P_w_i[2,l]  #  node l  active with prob. P_i[l] and link w_[il] negative
                inactive =  (1-  P_i[l]+P_i[l]*P_w_i[1,l])*m[l+1,h]# node l sends no contribution, either because it is active, or because the link w_[il] is zero
                m[l,h] = activating+ inhibiting+inactive
    return m[0,-offset]
def update_belief_pert(P_w_i,P_i,i,j,w,T,beta,theta,J0=1):
    '''
    Compute average of mismatch using recursive calls, whenever possible use the iterative option instead.
    P_w_i is P_w[:,i,:]. Because for every i, only links pointing to it matters
    '''
    bias = 0
    @lru_cache(maxsize=None)
    def recursion(bias, l):
        '''at the end of recursions, bias is the local field'''

        if (l == N):
            
            bias = (bias - theta) *J0
            mismatch = metric(1/2*(1+np.tanh((h_tilde[ind]*J0-theta)/2/T)),P_i[i])
            #mismatch = (0.5 + 0.5 *np.tanh((bias-theta) / 2 / T)-P_i[i])**2#difference betwen theoretical and target prob.
            return np.exp(-beta*mismatch)
        elif (l == i):
            #avoid self interaction
            return recursion(bias,l+1)
        elif (l == j):
            return recursion(bias+w,l+1)*P_i[l]+recursion(bias,l+1)*(1-P_i[l])
        activing = P_i[l] * recursion(bias + 1, l + 1)*P_w_i[0,l]  #  node l  active with prob. P_i[l] and link w_[il] positive
        inhibiting = P_i[l] * recursion(bias - 1, l + 1)*P_w_i[2,l]  #  node l  active with prob. P_i[l] and link w_[il] negative
        inactive =  (1-  P_i[l]+P_i[l]*P_w_i[1,l])*recursion(bias , l + 1)# node l sends no contribution, either because it is active, or because the link w_[il] is zero
        return activing+inhibiting+inactive

    N = len(P_i)
    result = recursion(bias, 0)
    recursion.cache_clear()
    return result

def belief_propagation_pert_short(P_i,P_perturb,P_perturb_crispa, T,  theta,lam,beta,distance_name = 'kl-divergence', J0=1,precision=1e-4, max_iter=100, using_pert = True):
    """
    Belief propagation using pertubation. This implementation uses only the unperturbed, crispri, and crispra condition .
    :param P_i: list of floats of length N. They are the target Probabilities
    :param P_w: 3d array of shape (3,N,N)
    :param T: float
    :param theta: float (in units of 1/sqrt(<K>))
    :param max_iter: int
    :param precision: float
    :return: P_w_new it is a 3d array of shape (3,N,N) which contains the probability of positive,negative,or zero weight for every link.
    ----NOTES------
    In order to help storing, couplings are taken to be +-1
    Even though code runs for any directed network, results  are exact for fully asymmetric networks only.
    """
    N = len(P_i)
    '''
    max_recursions = int((N + 1) * (N + 2) / 2)
    if max_recursions > sys.getrecursionlimit():
        print("Warning! maximum degree larger than default recursion limit, I 'll update recursion limit to",
              max_recursions)
        sys.setrecursionlimit(max_recursions)
    '''


    if distance_name.lower() in ['l2-norm','l2','l2_norm']:
        metric = l2_norm
    elif distance_name.lower() in ['kl-divergence','kl','kl_divergence']:
        metric = kl_divergence
    elif distance_name.lower() in ['kl-symmetric','kl_symmetric']:
        metric = kl_symmetric
    else:
        raise ValueError("Don't understand what metric you want")



    if using_pert:
        M = len(P_perturb)+len(P_perturb_crispa)+1#number of experimental conditions
        input_data = np.concatenate((np.expand_dims(P_i,0),P_perturb,P_perturb_crispa),axis = 0)#iterator on the different experimental conditions
    else:
        raise ValueError('Cavity iteration not implemented yet')
        M = 1
        input_data = [P_i]
    #initialise the P_w
    #rho= np.random.rand(3,M,N)
    P_w = np.zeros((3,N,N))#P_w[0] is P(w_ij=1),P_w[1] is P(w_ij=0),P_w[2] is P(w_ij=-1)
    store = np.empty(N,dtype = object)
    error = np.zeros((M,N))
    for i in range(N):
        #data type rho_nu has dimensions (3,3,N). It is structured as follows:
        # 1st dimension maps to link weight, i.e. indices 0,1,2 corrsponds to J = 1,0,-1
        # 2nd dimension maps to the perturbation experiment considered, i.e. indices  0,1,2 corresponds to  wild type, crispra, crispri
        # 3rd dimensions maps node j
        rho_nu = np.random.rand(3,3,N)
        rho_nu = rho_nu/rho_nu.sum(axis = 0)#normalise rho_nu
        rho_nu_old = rho_nu.copy()
        P_w_nu = rho_nu.copy()        
        P_w_nu[np.array([0,2])]=P_w_nu[np.array([0,2])]*np.exp(-lam)#penality term for the w = \pm 1     
        P_w_nu = P_w_nu/P_w_nu.sum(axis = 0)#normalise P_w_nu
        for count in range(max_iter):         
            for j in list(set(range(N))-{i}):
                for nu,P_i_nu in enumerate(input_data[[0,j+1,j+1+N]]):#take only the conditions corresponding to wild type, crispra on gene j, crispri on gene j
                #this loop updates the rho_nu
                    for w_index,w in enumerate([1,0,-1]):
                        rho_nu[w_index,nu,j] = update_belief_iterative(P_w_nu[:,nu,:],P_i_nu,i,j,w,T,beta,theta,metric,J0)

            rho_nu[:,:,i] =np.array([[0,1,0]]*3).T#no self interaction for j=i
            rho_nu = rho_nu/rho_nu.sum(axis = 0)#normalise rho

            error = np.sum(np.abs(rho_nu - rho_nu_old),axis = 0)#sum over the w={1,0,-1}
            rho_nu_old = rho_nu.copy()

            store[i] = np.append(store[i],error)
            if np.max(error) < precision:
                logging.info('finishing after '+ str(count) +' iterations')
                
                break
            if count == max_iter-1:
                logging.warning("Maximum number of repetition reached, but target  precision has not been reached for node "+str(i) +'\n with error'+str(np.max(error)))

            for nu in range(3):
                #this loop updates the P_nu once the loop over the rho_nu is over
                rho_cav = np.concatenate((rho_nu[:,:nu,:],rho_nu[:,nu+1:,:]),axis= 1)#exclude condition nu, i.e. remove the data associated to condition nu
                P_w_nu[:,nu,:] = np.prod(rho_cav,axis = 1)#multiplication is performed on the rho_nu for all mu neq nu
                
            P_w_nu[[0,2]]=P_w_nu[[0,2]]*np.exp(-lam)#penality term for the w = \pm 1                            
            P_w_nu = P_w_nu/P_w_nu.sum(axis = 0)#normalise P_w_nu


        #END of cavity, now do the marginal inside the loop over i
        for j in range(N):
            P_w[:,i,j] = np.prod(rho_nu[:,:,j],axis = 1)#multiplication is performed on the rho_nu for all mu
            
    P_w[[0,2]]=P_w[[0,2]]*np.exp(-lam)#penality term for the w = \pm 1
    P_w = P_w/P_w.sum(axis = 0)#normalise P_w_nu

    '''
    if np.max(error) < precision:
        print('finishing after', count, 'iterations')
        break
    if count == max_iter-1:
        print("Maximum numberprint(l) of repetition reached, but target  precision has not been reached. ")
    '''
    return P_w,store
