#perturbation
import numpy as np
from SERGIO.GRN import grn_from_human,grn_from_networkx
from functools import partial
from SERGIO._sergio import sergio
from SERGIO.MR import mrProfile
import networkx as nx
from collections import defaultdict
import pandas as pd
from SERGIO.GRN._grn import GRN
from SERGIO.GRN._components import Gene, SingleInteraction
from collections import Counter
import itertools
import time
import copy
import random
import logging
from tqdm import tqdm
import pickle
import os
if nx.__version__>'3.':
    nx.from_scipy_sparse_matrix = nx.from_scipy_sparse_array
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class Perturbation(sergio):
    '''
    Use method self.perturbation_all to simulate perturbations.
    The outcome is stored in self.wt, self.crispra_data, and self.crispri_data
    '''
    def __init__(self,relative_path='SERGIO') -> None:
        '''
        Load the simulation as a sergio object that has been already run
        '''
        self.data_path = os.path.join(relative_path,'data')#path where simulation are saved/loaded
    def initialise_parameters(self,grn,basal_prod_crispri = 0.2,basal_prod_crispra = 7):
        '''
        basal_prod_crispri, basal_prod_crispra is the basal production of the nodes that are perturbed with crispri and crispra respectively
        '''
        self.grn_ = grn
        self.gNames_ = self._get_gene_names(grn)
        G = self.grn_.to_networkx()
        #finds the nodes with out-degree >0
        n,d = zip(*dict(G.out_degree()).items())
        self.nodes_2perturb = np.array(n)[np.array(d)>0]#it is the array of nodes which regulates at least a node 
        logger.info('Number of nodes will be perturbed is '+str(len(self.nodes_2perturb)))
        self.basal_prod_crispri = basal_prod_crispri
        self.basal_prod_crispra = basal_prod_crispra
    @staticmethod
    def create(N):
        '''
        Static method to create a Gene Regulatory Network (GRN) object.

        Parameters:
        - N (int): Number of nodes in the generated random graph.

        Returns:
        - grn_ (GRN): Generated Gene Regulatory Network object.

        Description:
        This static method creates a Gene Regulatory Network (GRN) object.
        It generates a random directed graph with a specified number of nodes using the
        `create_random_graph` function from the perturbation module.
        Then, it converts the generated graph to a GRN object using the `grn_from_networkx` function.
        The resulting GRN object is returned.
        '''

        G = perturbation.create_random_graph(N)
        grn_ = grn_from_networkx(G)
        return grn_
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
    def perturbation_all(self,nCells,multiprocess = True):
        '''creates the self.crispri_data,self.crispra_data,self.wt:
        Self.crispri_data,self.crispra_data are multidimensional np.arrays with: 
        - dimension 0 describes the index target gene
        - dimension 1 describe the index of genes whose expression is measured
        - dimension 2 describes the index of cell
        self.wt is a 2d np.array with
        - dimension 1 describe the index of genes whose expression is measured
        - dimension 2 describes the index of cell
 
        '''
        if multiprocess:
            from pathos.multiprocessing import ProcessingPool as Pool
            #from multiprocessing.pool import Pool

        grn = copy.deepcopy(self.grn_)#it is the grn that I pass to perturbation experiments
        nGenes = len(grn.attr_['genes'])
        self.nCells_ = [nCells]#it is only one cell type, for more, do [nCells]*#celltypes
        
        '''initialise the simulator'''
        mrs = grn.get_mrs()
        mr_profs = mrProfile(MR_names = mrs, n_types = 1)
        mr_profs.build_rnd(range_dict={'L': [1, 2.5], 'H': [3.5, 5]})
        grn.init(mr_profs, update_half_resp = True)
        '''now simulate wild type'''
        start_time = time.time()
        self.wt = self._wild_type(grn=grn,mr_profs=mr_profs,nCells=nCells)

        print("--- %s seconds ---" % (time.time() - start_time))
        logger.info('Time of execution for wt is  %s seconds'% (time.time() - start_time))

        def wrapper_pert_funct(target_gene,basal_prod):
            return self._single_perturbation(grn = grn,mr_profs=mr_profs,target_gene=target_gene,basal_prod=basal_prod,nCells=nCells,cutting=True).getSimExpr().values            

        if multiprocess:
            with Pool() as pool:
                self.crispri_data = np.array(list(pool.map(partial(wrapper_pert_funct,basal_prod=self.basal_prod_crispri),self.nodes_2perturb)))
                self.crispra_data = np.array(list(pool.map(partial(wrapper_pert_funct,basal_prod=self.basal_prod_crispra),self.nodes_2perturb)))
        else:
            self.crispri_data = np.zeros(shape = (len(self.nodes_2perturb),nGenes,nCells))
            self.crispra_data = np.zeros(shape = (len(self.nodes_2perturb),nGenes,nCells))
            for i,target_gene in tqdm(enumerate(self.nodes_2perturb)):
                #simulate Crispri
                self.crispri_data[i]= self._single_perturbation(grn = grn,mr_profs=mr_profs,target_gene=target_gene,basal_prod=self.basal_prod_crispri,nCells=nCells,cutting=True).getSimExpr().values
                #simulate Crispra
                self.crispra_data[i]= self._single_perturbation(grn = grn,mr_profs=mr_profs,target_gene=target_gene,basal_prod=self.basal_prod_crispra,nCells=nCells,cutting=True).getSimExpr().values
    def technical_noise(self,crispri = None, crispra = None, wt = None,multiprocess = True):
        '''
        Apply technical noise to the simulated expression data from perturbation.


        Returns:
        - wt (list): List of expression data with technical noise applied for the wild type.
        - crispri (list): List of expression data with technical noise applied for CRISPRi experiments.
        - crispra (list): List of expression data with technical noise applied for CRISPRa experiments.

        Description:
        This method applies technical noise to the expression data, including outlier genes, library size effect,
        and dropouts. The noise is added separately for each experiment type (wild type, CRISPRi, CRISPRa).
        If `multiprocess` is set to True, multiprocessing is utilized to speed up the computation.
        '''
        if (wt is None ) or (crispri is None) or (cirspra is None):
            wt = self.wt
            crispri = self.crispri_data
            crispra = self.crispra_data
        wt = self._apply_technical_noise(wt)
        if multiprocess:
            from pathos.multiprocessing import ProcessingPool as Pool
            with Pool() as pool:
                crispri = pool.map(self._apply_technical_noise,crispri_data)
                crispra = pool.map(self._apply_technical_noise,crispri_data)                
             
        else:
            crispri = list(map(self._apply_technical_noise,self.crispri_data))
            crispra = list(map(self._apply_technical_noise,self.crispri_data))                
            
        return wt,crispri,crispra
    def _apply_technical_noise(self,expr):
        '''
        Apply technical noise to the given expression data.

        Parameters:
        - expr (array-like): Input expression data.

        Returns:
        - count_matrix (array-like): Expression data with technical noise applied.

        Description:
        This method adds technical noise to the input expression data by simulating outlier genes,
        library size effect, and dropouts. It returns the expression data with noise applied.
        '''
        """
        Add outlier genes
        """
        expr_O = self.outlier_effect(expr, outlier_prob = 0.01, mean = 0.8, scale = 1)

        """
        Add Library Size Effect
        """
        libFactor, expr_O_L = self.lib_size_effect(expr_O, mean = 4.6, scale = 0.4)

        """
        Add Dropouts
        """
        binary_ind = self.dropout_indicator(expr_O_L, shape = 6.5, percentile = 82)
        expr_O_L_D = np.multiply(binary_ind, expr_O_L)

        """
        Convert to UMI count
        """
        count_matrix = self.convert_to_UMIcounts(expr_O_L_D)
        return count_matrix
    
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
        '''
        Generate a random directed graph with a specified number of nodes.

        Parameters:
        - N (int): Number of nodes in the graph.

        Returns:
        - G (networkx.DiGraph): Random directed graph.

        Description:
        This function generates a random directed graph with a specified number of nodes.
        The graph is generated using the `random_k_out_graph` function from the NetworkX library,
        which creates a directed graph where each node has a fixed out-degree.
        The graph is then filtered to retain only the largest strongly connected component.
        The weights of the edges are sampled from two uniform distributions,
        with a specified fraction of positive links.
        The final graph is returned as a NetworkX DiGraph object.

        Notes:
        - The number of nodes in the graph may be smaller than the specified N.
        - This function may generate graphs with fewer nodes or edges than expected
        due to the filtering process to retain only the largest strongly connected component.
        - The weights of the edges are sampled from two uniform distributions
        to introduce randomness and variability in the graph structure.
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
    def save(self):
        folder = self.data_path
        if not os.path.exists(folder):
            logger.debug('creating folder '+folder)
            os.makedirs(folder)
        nGenes = len(self.nodes_2perturb)
        nCells = self.wt.shape[-1]
        filename = 'pert_'+str(nGenes)+'_'+str(nCells)
        with open(self.data_path+'/'+filename+'.pkl','wb') as file:
            pickle.dump(self.__dict__,file)
     
    def load(self,nGenes, nCells):
        """try load filename.pkl"""
        filename = 'pert_'+str(nGenes)+'_'+str(nCells)
        with open(self.data_path+'/'+filename+'.pkl','rb') as file:
            dataPickle = pickle.load(file)
        self.__dict__.update(dataPickle)#add loaded dictionary to methods of the class
        return self