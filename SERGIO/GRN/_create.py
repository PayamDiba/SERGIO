import warnings
from ._components import Gene, SingleInteraction
from ._grn import GRN
import pandas as pd
from ._utils import grnParam, parameterize_grn
import networkx as nx
import numpy as np
import zipfile
import random


def grn_from_v1(path, n = 2, decay = 0.8):
    """
    if set n = None, it's read from file
    """
    ret = GRN()
    with open(path,'r') as f:
        for r in f.readlines():
            r = r.split(',')
            nInt = int(float(r[1]))
            tar = r[0]
            regs = r[2:2+nInt]
            ks = r[2+nInt : 2+2*nInt]
            ks = [float(i) for i in ks]
            if not n:
                ns = r[2+2*nInt : 2+3*nInt]
                ns = [float(i) for i in ns]
            else:
                ns = [n] * nInt

            for currR,currK,currN in zip(regs,ks,ns):
                reg = Gene(name = currR, decay = decay)
                tar = Gene(name = tar, decay = decay)
                currInter = SingleInteraction(reg = [reg], tar = tar, k = currK, h = None, n = currN)
                ret.add_interaction(currInter)
    ret.setMRs()
    return ret


def grn_from_file(path, parameterize = False, k_act = [1,5], k_rep = [-5,-1], n = 2, decay = 0.8):
    names = ['index','reg','coop','tar','k','n','h','reg_decay','tar_decay']
    net = pd.read_csv(path, header = 0, index_col = 0, names = names)
    # it will hanadle it even if file does not contain k,n or h

    if parameterize:
        param = grnParam(k_act, k_rep, n, decay)
        net = parameterize_grn(net, param)


    ret = GRN()
    for _,r in net.iterrows():
        if r.coop:
            raise ValueError('Cooperative interactions are not implemented yet') # TODO: implement multiple interactions
        reg = Gene(name = str(r.reg), decay = r.reg_decay)
        tar = Gene(name = str(r.tar), decay = r.tar_decay)
        currInter = SingleInteraction(reg = [reg], tar = tar, k = r.k, h = r.h, n = r.n)
        ret.add_interaction(currInter)

    ret.setMRs()
    return ret

def grn_from_human(nGenes = 400,k_act = 1,k_rep = -1,n = 2, decay = 0.8):
    def download_network():
        '''It checks if data exists, if not it downloads from https://regnetworkweb.org/download/RegulatoryDirections.zip'''
        import os
        import ssl
        import urllib.request

        url = "https://regnetworkweb.org/download/RegulatoryDirections.zip"

        if not os.path.exists(file_path):
            print("File does not exist. Downloading...")
            #Disable SSL verification, othwerwise it gives error
            context = ssl.create_default_context()
            context.check_hostname = False
            context.verify_mode = ssl.CERT_NONE

            with urllib.request.urlopen(url, context=context) as u, open(file_path, "wb") as f:
                file_data = u.read()
                f.write(file_data)

            print("Download complete!")
        else:
            print("File already exists.")

    file_path = "SERGIO/GRN/ref/RegulatoryDirections.zip"
    file_name = "new_kegg.human.reg.direction.txt"

    download_network()

    """
    read edges
    loads the data into a pandas dataframe and clean up
    """
    
    with zipfile.ZipFile(file_path, "r") as zip_file:
        with zip_file.open(file_name) as file:
            tb = pd.read_csv(file, sep=" ",usecols=[0,1,2,3,4],names = ['#TF', 'ID', 'Target', 'ID.1','sign'],skiprows=1)
    tb.replace({'-->':1,'--|':-1},inplace = True)
    tb.drop(labels= tb.index[~tb.sign.isin([k_act,k_rep])],axis = 'index', inplace=True)
    
    g = nx.DiGraph()
    g.add_weighted_edges_from(tb[['ID','ID.1','sign']].values)
    if nGenes>len(g.nodes()):
        nGenes = len(g.nodes())
        warnings.warn('You want to sample '+str(nGenes)+' genes, but the full network contains '+str(len(g.nodes()+'. Sampling all genes available')))

    """
    sample genes
    Uses a method that starts from a node and follows links, it aims to return a connected network
    """

    subset = _sample_network(G = g,nGenes=nGenes)
    subgraph = nx.subgraph(g,subset).copy()
    grn = grn_from_networkx(g = subgraph,parametrize=False)
    return grn

def grn_from_Ecoli(nGenes = 400, k_act = [1,5], k_rep = [-5,-1], n = 2, decay = 0.8):
    """
    read genes
    """
    df = pd.read_csv('SERGIO/GRN/ref/EcoliNet.v1.txt', sep = '\t', names = ['reg','tar','score'], header = None)
    edges = df[['reg','tar']].drop_duplicates()

    """
    remove cycles
    """
    g = nx.DiGraph()
    g.add_edges_from(edges.values)

    try:
        a = nx.find_cycle(g)
        while len(a) > 0:
            g.remove_edges_from(a)
            try:
                a = nx.find_cycle(g)
            except:
                a = []
    except:
        pass

    ind = 0
    modified_edges = []
    for i in list(g.edges()):
        modified_edges.append([ind, i[0], 0, i[1]])
        ind += 1

    modified_edges = pd.DataFrame(modified_edges)
    modified_edges.columns = ['index','reg','coop','tar']

    """
    sample genes
    """
    genes = modified_edges.reg.unique().flatten().tolist() + modified_edges.tar.unique().flatten().tolist()
    genes = np.unique(genes)
    genes = np.random.choice(genes, nGenes, replace = False)

    edges = modified_edges.loc[(modified_edges.reg.isin(genes)) & (modified_edges.tar.isin(genes))]

    """
    parameterize
    """
    param = grnParam(k_act, k_rep, n, decay)
    net = parameterize_grn(edges, param)

    ret = GRN()
    for _,r in net.iterrows():
        if r.coop:
            raise ValueError('Cooperative interactions are not implemented yet') # TODO: implement multiple interactions
        reg = Gene(name = str(r.reg), decay = r.reg_decay)
        tar = Gene(name = str(r.tar), decay = r.tar_decay)
        currInter = SingleInteraction(reg = [reg], tar = tar, k = r.k, h = r.h, n = r.n)
        ret.add_interaction(currInter)

    ret.setMRs()
    return ret

def _sample_network(G,nGenes):
    '''nGenes: size of the sample

    Sample a connected subset of nodes from a directed network using a random walk.

    The function takes as input a directed network, and samples a connected subset of nodes from the network using a random walk starting from a random node. The subset is guaranteed to be connected, and is sampled with replacement until it reaches the desired size.


    Arguments:
        G    nx.DiGraph()
        nGenes      The number of nodes to sample from the network.

    Output:
        A list of node IDs representing the sampled subset of the network.


    Notes:
        - The network file should contain one edge per line, with the source and target node IDs separated by a space.
        - The script uses a random walk starting from a random node to sample the subset, and may take some time to complete for large networks.
        - The sampling process is repeated with replacement until the desired number of nodes is reached, so the resulting subset may contain duplicate nodes.
        - The sampling process is biased towards nodes with a high degree, so the resulting subset may not be representative of the network as a whole.
    '''

    # Initialize the subset
    subset = set()

    # Loop until the subset has nGenes nodes
    while len(subset) < nGenes:
        # Select a random start node from a weakly connected component not already in the subset
        start_nodes = set(G.nodes()) - subset
        start_node = None
        while not start_node:
            if not start_nodes:
                break
            node = random.choice(list(start_nodes))
            for comp in nx.weakly_connected_components(G.subgraph(start_nodes)):
                if node in comp:
                    start_node = node
                    break
            if not start_node:
                start_nodes.remove(node)

        # If no more weakly connected components are available, break the loop
        if not start_node:
            break

        # Initialize the subset with the start node and its descendants
        descendants = set(nx.descendants(G, start_node))
        subset.update(descendants.union({start_node}))

        # Loop until the subset has nGenes nodes or the weakly connected component is fully sampled
        while len(subset) < nGenes and len(descendants) < len(G):
            # Find all nodes adjacent to the subset in the weakly connected component
            adjacent_nodes = set()
            for node in subset:
                adjacent_nodes.update(set(G.successors(node)))
            adjacent_nodes -= subset

            # If there are no more adjacent nodes, break the loop and move on to the next weakly connected component
            if not adjacent_nodes:
                break

            # Select a random adjacent node to add to the subset
            new_node = random.choice(list(adjacent_nodes))
            subset.add(new_node)
            descendants.add(new_node)

            # Remove all edges between the subset and its adjacent nodes in the weakly connected component
            for node in adjacent_nodes:
                if node in descendants:
                    G.remove_edge(node, new_node)

    return subset
def grn_from_networkx(g,k_act = [1,5], k_rep = [-5,-1],parametrize = False,hill_coef = 2., decay = 0.8):
    '''
    Function Name: grn_from_networkx

    Description:
    This function converts a networkx.DiGraph() object into a Gene Regulatory Network (GRN) object of SERGIO. The interaction strengths between genes are loaded from the weights of the links in the graph. If the parametrize parameter is set to True or the input graph is not weighted, the function assigns random weights to links. The weights are sampled uniformly from the domain given by the k_act and k_rep parameters, with a 75% probability of selecting k_act and a 25% probability of selecting k_rep.

    Parameters:

    g: A networkx.DiGraph() object to be converted to a GRN object.
    k_act (default=[1,5]): A list containing the lower and upper bounds of the domain from which the function samples activation strengths.
    k_rep (default=[-5,-1]): A list containing the lower and upper bounds of the domain from which the function samples repression strengths.
    parametrize (default=False): If True, the function assigns random weights to links.
    hill_coef (default=2.): The Hill coefficient used to calculate the strength of the interaction.
    decay (default=0.8): The decay rate of the genes in the GRN.
    Returns:

    ret: A GRN() object containing the converted GRN.
    '''
    def _remove_cycles(G):
        #removing cycles
        while True:
            try:
                cycle = nx.find_cycle(G)
                G.remove_edges_from([cycle[0]])
            except nx.exception.NetworkXNoCycle:
                break
        return G

    ret = GRN()
    G = g.copy()# I am making a copy because the input graph is modified ( in _remove_cycled, if parametrize==True, or if g is not weighted)
    '''
    Removes cycles
    '''
    G = _remove_cycles(G)
    if parametrize or not nx.is_weighted(G):
        n_edges = len(G.edges())
        cond = np.random.uniform(size = n_edges) < 0.75 #sample True with prob 75% 
        kRange = np.where(cond, np.array([k_act]*n_edges).T,np.array([k_rep]*n_edges).T)#sample k_act to True and k_rep
        k = np.random.uniform(low = kRange[0], high = kRange[1])
        reg,tar = zip(*list(g.edges()))
        G = nx.DiGraph()
        G.add_weighted_edges_from(zip(reg,tar,k))#create a new networkx.DiGraph() with weighted edges
        
    for reg,tar,weight in G.edges(data = True):
        k = weight['weight']# extract numerical value
        reg = Gene(name = str(reg), decay = decay)
        tar = Gene(name = str(tar), decay = decay)
        currInter = SingleInteraction(reg = [reg], tar = tar, k = k, h = None, n = hill_coef)
        ret.add_interaction(currInter)

    ret.setMRs()
    return ret