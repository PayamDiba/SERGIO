from ._components import Gene, SingleInteraction
from ._grn import GRN
import pandas as pd
from ._utils import grnParam, parameterize_grn
import networkx as nx
import numpy as np

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

def grn_from_human(nGenes = 400, k_act = [1,5], k_rep = [-5,-1], n = 2, decay = 0.8):
    """
    read genes
    """
    df = pd.read_csv('SERGIO/GRN/ref/human.node', sep = '\t', names = ['id','symb','unk'], header = None)
    genes = []
    for _,r in df.iterrows():
        curr = [str(r.id), str(r.symb)]
        if str(r.unk) != 'nan' and str(r.unk) != str(r.id):
            curr[-1] = str(r.unk)

        if 'miR' not in curr[-1] and 'hsa' not in curr[-1]:
            genes.append(curr)
    genes = pd.DataFrame(np.array(genes).reshape(-1,2))
    genes.columns = ['id','symb']

    """
    read edges
    """
    df = pd.read_csv('SERGIO/GRN/ref/human.source', sep = '\t', header = None, dtype='string', names = ['reg_symb', 'reg_id', 'tar_symb', 'tar_id'])
    edges = df.loc[(df['reg_id'].isin(genes.id)) & (df['tar_id'].isin(genes.id))][["reg_id","tar_id"]].drop_duplicates()

    """
    remove cycles
    """
    g = nx.DiGraph()
    g.add_edges_from(edges.values)


    a = nx.find_cycle(g)
    while len(a) > 0:
        g.remove_edges_from(a)
        try:
            a = nx.find_cycle(g)
        except:
            a = []

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
