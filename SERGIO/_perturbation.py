#perturbation
import numpy as np
from SERGIO.GRN import grn_from_human,grn_from_networkx
from SERGIO._sergio import sergio
from SERGIO.MR import mrProfile
import networkx as nx
from collections import defaultdict
import pandas as pd
from SERGIO.GRN._grn import GRN
from SERGIO.GRN._components import Gene, SingleInteraction
from collections import Counter
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import copy
import random
import scperturb
import logging
import pickle
import os
#logging.basicConfig(level=logging.INFO)
class perturbation(object):
    '''
    Use method self.perturbation_all to simulate perturbations.
    The outcome is stored in self.wt and self.store
    '''
    def __init__(self) -> None:
        pass
    def create(self,N):
        #create GRN object
        G = self.create_random_graph(N)
        self.grn0 = grn_from_networkx(G)
        self.G = self.grn0.to_networkx()
        logging.info('Number of nodes is '+str(len(self.G)))
        #grn0 = copy.deepcopy(self.grn)#it is the grn that I pass to perturbation experiments
    @staticmethod
    def _wild_type(grn,mr_profs,nCells):
        '''Return a numpy array of shape nGenes,nCells'''
        sim = sergio(grn)
        sim.simulate(nCells = nCells, noise_s = 1, safety_iter = 150, scale_iter = 10)
        return  sim.getSimExpr().values
    @staticmethod
    def _single_perturbation(grn,mr_profs,target_gene,basal_prod,nCells,cutting = True):
        '''
        Function Description:

        This function perturbs a given gene in a gene regulatory network (GRN) and returns
        the simulated dynamics of the network in response to the perturbation. The perturbation involves
        changing the production rate of the target gene to a given basal level (positive for CRISPRa,
        negative or 0 for CRISPRi), and optionally, cutting incoming edges to the target gene.
        The function returns the simulated gene expression dynamics of the perturbed network.

        Parameters:

        grn: an instance of the GRN class representing the GRN to be perturbed
        mr_profs: an instance of the mrProfile class representing the master regulator profiles of the system
        target_gene: the name of the target gene to be perturbed
        basal_prod: the basal production rate of the target gene after perturbation (positive for CRISPRa, negative or 0 for CRISPRi)
        cutting: a boolean variable indicating whether incoming edges to the target gene should be cut (True) or not (False)
        in response to the perturbation. Default is True.
        Returns:

        sim: an instance of the sergio class representing the simulated gene expression dynamics of the perturbed network.
        
        '''
        new_grn = copy.deepcopy(grn)
        new_grn.attr_['mrs'].add(target_gene)
        
        if cutting:
            #identifies the links pointing to gene, and removes them from the net
            links_2_remove = list(filter(lambda x:x.split('-')[1]==target_gene,new_grn.attr_['interactions'].keys()))
            [new_grn.attr_['interactions'].pop(link) for link in links_2_remove]#remove interaction from GRN
        
        crispr = copy.deepcopy(mr_profs)
        crispr.profile[target_gene]=basal_prod*np.ones(1)
        #clear up previous stationary point trajectory
        for g in new_grn.attr_['genes'].values():
            g.sim_conc_ = defaultdict(list) 
        new_grn.init(crispr, update_half_resp = False)#we set update_half_resp = False because the half_response has already been set in the wt.
        #If we set = True, the perturbation does not propagate to regulated genes 

        sim = sergio(new_grn)
        sim.simulate(nCells = nCells, noise_s = 1, safety_iter = 150, scale_iter = 10)
        return sim
    def perturbation_all(self,nCells):
        nGenes = len(self.grn0.attr_['genes'])
        self.store = np.zeros(shape = (2,nGenes,nGenes,nCells))
        #dimension 0 describe the type of perturbation experiment, 0:Crispri, 1: crispra
        #dimension 1 describes the index target gene
        #dimension 2 describe the index of genes whose expression is measured
        #dimension 3 describes the index of cell
        
        '''initialise the simulator'''
        grn = copy.deepcopy(self.grn0)#it is the grn that I pass to perturbation experiments
        mrs = grn.get_mrs()
        mr_profs = mrProfile(MR_names = mrs, n_types = 1)
        mr_profs.build_rnd(range_dict={'L': [1, 2.5], 'H': [3.5, 5]})
        grn.init(mr_profs, update_half_resp = True)
        '''now simulate wild type'''
        self.wt = self._wild_type(grn=grn,mr_profs=mr_profs,nCells=nCells)

        for i,target_gene in enumerate(grn.attr_['genes'].keys()):
            #simulate Crispri
            self.store[0,i]= self._single_perturbation(grn = grn,mr_profs=mr_profs,target_gene=target_gene,basal_prod=0.2,nCells=nCells,cutting=True).getSimExpr().values
            #simulate Crispra
            self.store[1,i]= self._single_perturbation(grn = grn,mr_profs=mr_profs,target_gene=target_gene,basal_prod=7,nCells=nCells,cutting=True).getSimExpr().values

        
    @staticmethod
    def draw_net(G,offset = 0.05,node_size = 1000,**kwargs):
        c =np.array([ c['weight'] for a,b,c in list(G.edges(data=True))])
        edge_color=np.where(c>0,'green','red')
        #nx.draw_circular(G,with_labels=True,edge_color=edge_color,alpha = 0.7,arrowsize = 15,**kwargs)
        nodePos = nx.circular_layout(G)
        nx.draw_networkx_nodes(G,pos=nodePos, label=True,node_size=node_size,node_color='none',edgecolors='#1f78b4')
        nx.draw_networkx_labels(G,pos = nodePos)
        edges = list(G.edges)
        bi_edges = [(a,b) for a,b in edges  if (b,a)in edges]
        non_bi_edges =list(set(edges)-set(bi_edges))
        a,b,c  =zip(*[ (a,b,c['weight']) for a,b,c in list(G.edges(data=True))])
        dic = {(start,stop):col for start,stop,col in zip(a,b,c)}
        #nx.draw_networkx_edges(G,pos = nodePos,edgelist=edges[weight>0],edge_color= 'green',arrowsize = 15,)
        non_bi_weight = np.array([dic[endpoints] for endpoints in non_bi_edges ])
        bi_weight = np.array([dic[endpoints] for endpoints in bi_edges ])
        nx.draw_networkx_edges(G,pos = nodePos,edgelist=np.array(non_bi_edges)[non_bi_weight>0],edge_color= 'green',arrowsize = 20,node_size=node_size)
        nx.draw_networkx_edges(G,pos = nodePos,edgelist=np.array(non_bi_edges)[non_bi_weight<0],edge_color= 'red',arrowsize = 10,node_size=node_size,arrowstyle='-[',alpha = 0.6)
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

        nx.draw_networkx_edges(G,pos = new_nodePos,edgelist=np.array(bi_edges)[bi_weight<0],edge_color= 'red',arrowsize = 10,arrowstyle='-[',alpha = 0.6,width = 1.5,node_size=node_size)
    @staticmethod
    def create_random_graph(N):
        '''N:number of nodes
        Note that the number of nodes in the graph will likely be smaller than N
        Weights sampled from 2 uniform distributions
        '''
        weight_pos = 3
        weight_neg = -weight_pos
        band=1
        frac_pos = 0.5 #fraction of positive links over total
        G = nx.generators.random_k_out_graph(n = N,k = 4,alpha = 0.9)
        largest = max(nx.strongly_connected_components(G), key=len)
        G = G.subgraph(largest).copy()#filter graph to maximum  SCC
        G = nx.DiGraph(G)#remove multi links
        M = len(G.edges())#n. of links
        w1 = np.random.uniform(weight_neg-band,weight_neg+band,int(M*frac_pos))#sample weights from 2 uniform distributions
        w2 = np.random.uniform(weight_pos-band,weight_pos+band,M-len(w1))
        w = np.concatenate([w1, w2])
        random.shuffle(w)#avoid having correlation in the position
        J = nx.adjacency_matrix(G)
        J.data = w
        G = nx.from_scipy_sparse_matrix(J,create_using=nx.DiGraph())
        return G
    def save(self,filename):
        folder = 'data'
        if not os.path.exists(folder):
            logging.debug('creating folder '+folder)
            os.makedirs(folder)
        nGenes = len(self.grn0.attr_['genes'])
        nCells = self.wt.shape[-1]
        filename = 'pert_'+str(nGenes)+'_'+str(nCells)
        with open(folder+'/'+filename+'.pkl','wb') as file:
            pickle.dump(self.__dict__,file)
     
    def load(self,filename):
        """try load filename.pkl"""
        folder = 'data'
        with open(folder+'/'+filename+'.pkl','rb') as file:
            dataPickle = pickle.load(file)
        self.__dict__ = dataPickle
        return self