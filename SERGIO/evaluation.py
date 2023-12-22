import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import pandas as pd
def draw_net(G,offset = 0.05,node_size = 3000,arrowsize = 20,as_matrix = False):
    '''Option as_matrix displays the adjacency matrix, where element a_ij refers to link i->j'''
    if as_matrix:
        import matplotlib
        import seaborn as sns
        cmap = matplotlib.colors.ListedColormap(["red", "whitesmoke", "green"])
        sns.heatmap(pd.DataFrame(nx.adjacency_matrix(G).todense(),index = G.nodes(),columns=G.nodes()),ax=plt.gca(),cmap = cmap)
    else:
        c =np.array([ c['weight'] for a,b,c in list(G.edges(data=True))])
        edge_color=np.where(c>0,'green','red')
        #nx.draw_circular(G,with_labels=True,edge_color=edge_color,alpha = 0.7,arrowsize = 15,**kwargs)
        nodePos = nx.kamada_kawai_layout(G,dist = dict(nx.shortest_path_length(G, weight = None)))
        nx.draw_networkx_nodes(G,pos=nodePos, label=True,node_size=node_size,node_shape='o',node_color='w',edgecolors='k',linewidths=2)
        nx.draw_networkx_labels(G,pos = nodePos,font_size=14)
        edges = list(G.edges)
        bi_edges = [(a,b) for a,b in edges  if (b,a)in edges]
        non_bi_edges =list(set(edges)-set(bi_edges))
        a,b,c  =zip(*[ (a,b,c['weight']) for a,b,c in list(G.edges(data=True))])
        dic = {(start,stop):col for start,stop,col in zip(a,b,c)}
        #nx.draw_networkx_edges(G,pos = nodePos,edgelist=edges[weight>0],edge_color= 'green',arrowsize = 15,)
        non_bi_weight = np.array([dic[endpoints] for endpoints in non_bi_edges ])
        bi_weight = np.array([dic[endpoints] for endpoints in bi_edges ])
        nx.draw_networkx_edges(G,pos = nodePos,edgelist=np.array(non_bi_edges)[non_bi_weight>0],edge_color= 'green',arrowsize = arrowsize*1.5,node_size=node_size)
        nx.draw_networkx_edges(G,pos = nodePos,edgelist=np.array(non_bi_edges)[non_bi_weight<0],edge_color= 'red',arrowsize = arrowsize,arrowstyle='-[',alpha = 0.6,node_size=node_size)
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
        nx.draw_networkx_edges(G,pos = new_nodePos,edgelist=np.array(bi_edges)[bi_weight>0],edge_color= 'green',arrowsize = arrowsize*1.5,node_size=node_size)
        new_nodePos={}
        for start,end in unique_bi_edges:
            new_nodePos[start] = nodePos[start]+[0,offset]
            new_nodePos[end] = nodePos[end]+[0,offset]

        nx.draw_networkx_edges(G,pos = new_nodePos,edgelist=np.array(bi_edges)[bi_weight<0],edge_color= 'red',arrowsize = arrowsize,arrowstyle='-[',alpha = 0.6,width = 1.5,node_size=node_size)


def compare_network_to_ref(J_ref,J_pred):
    def statistics(J_ref,J_pred):
        prec = (np.count_nonzero(J_pred==J_ref)-N)/(N*(N-1))
        #print('precision of reconstruction:',prec)
        tp = np.count_nonzero((J_pred==J_ref)[J_ref!=0])/np.count_nonzero(J_ref!=0)# true positive
        #print('true positive rate',tp)
        
        tn = (np.count_nonzero((J_pred==J_ref)[J_ref==0])-N)/(np.count_nonzero(J_ref==0)-N)# true positive
        #print('true negative rate',tn)
        
        return prec,tp,tn
    def create_random_matrix(number_links ):
        
        indx = random.sample(range(N*(N-1)),k=number_links)
        A=np.arange(N**2)
        indx_0 = list(set(A)-set(indx))
        A[indx]=np.where(np.random.rand(len(indx))>0.5,1,-1)
        A[indx_0]=0
        return A.reshape((N,N))
    N = J_ref.shape[0]
    val = [np.array(list(zip(*[statistics(J_ref,create_random_matrix(np.count_nonzero(J_pred)))   for _ in range(5000)]))).mean(axis = 1)]
    val = np.concatenate((val,[statistics(J_ref,J_pred)]))
    result = pd.DataFrame({'precision':val[:,0],'true_positive':val[:,1],'true_negative':val[:,2]},index=['random','model'])
    
    return result
def network(P_w,quantile=0.75,graphical = True):
    """
    Builds a directed network from a probability matrix P_w.

    Parameters:
    - P_w (np.array): a 3xNxN matrix where N is the number of nodes. The first slice of P_w represents the probability of positive links, the second slice of P_w represents the probability of neutral links, and the third slice of P_w represents the probability of negative links.
    - quantile (float, optional): the quantile used to extract the top-quantile links from P_w (default=0.75).
    - graphical (bool, optional): whether to display a histogram of P_w values and the selected quantile (default=True).

    Returns:
    A directed networkx graph.

    Note:
    The extracted links are weighted with a positive weight for positive links, 0 weight (no link) for neutral links, and a negative weight for negative links. The function uses the seaborn library to display the histogram.
    """

    #extract links that are more likely from P_w 
    # look at positive links
    def sample_links_from_P_w(A,title):
        '''A is the matrix P_w, i.e. the prob of links, for either active or inhibitor links'''
        q = np.quantile(A.flatten(),quantile)#take the top 1% of links by default
        if  graphical:
            
            plt.figure()
            sns.histplot(A.flatten())
            plt.axvline(q,ls = '--',c = 'k')
            plt.xlabel('$P(w)$')
            plt.title(title)

        yv,xv = np.meshgrid(nodeset,nodeset)
        return list(zip(xv[A>q],yv[A>q]))
    nodeset=np.arange(0,len(P_w[0]))
    pos_links = sample_links_from_P_w(P_w[0],'positive links')
    neg_links = sample_links_from_P_w(P_w[2],'negative links')

    G = nx.DiGraph()
    G.add_nodes_from(nodeset)
    G.add_weighted_edges_from(list(zip(*([*zip(*pos_links)]+[[1]*len(pos_links)]))))
    G.add_weighted_edges_from(list(zip(*([*zip(*neg_links)]+[[-1]*len(neg_links)])))) 
    return G
def crispr_raw_network_gen(P_crispra,P_crispri,perc):
    '''It is a simple estimation of the network from crispr experiments.
    Note:
    The formula np.log(u*(1-d)/(d*(1-u))) below is derived through a mean field assumption'''
    N = np.shape(P_crispri)[1]
    u= P_crispra
    d = P_crispri
    A = (np.log(u*(1-d)/(d*(1-u))))
    np.fill_diagonal(A,0)

    non_diag = np.where(~np.eye(N,dtype=bool))
    #plt.hist(A[a],bins = 50)
    x = A[non_diag]
    low, up =np.percentile(x,[perc*100,100*(1-perc)])#threshold for low and up links
    J_reconstruct = np.where(A>up,1,np.where(A<low,-1,0))
    return J_reconstruct
def PR_ROC(P_w,J,graphical=True):
    '''Returns:
     auPR,auROC. Considers all element of J but the diagonal ( which are set to 0 by construction, so they are not part of the prediction
     '''
    from sklearn.preprocessing import label_binarize
    from sklearn import metrics
    from sklearn.metrics import precision_recall_curve, roc_curve
    N = J.shape[1]
    row,col = np.concatenate((np.triu_indices(N,k = 1),np.tril_indices(N,k = -1)),axis = 1)#all indices but the diagonal
    y_true = J.toarray()[col,row]
    one_hot_enc = label_binarize(y_true,classes = [1,0,-1])
    precision = dict()
    recall = dict()
    fpr = dict()
    tpr = dict()
    auROC = dict()
    auPR = dict()
    if graphical:
        plt.figure(figsize = (4,3))
    for w_index,w in enumerate([1,0,-1]):
        prob_pred = P_w[w_index,row,col]
        precision[w],recall[w],_ = precision_recall_curve(one_hot_enc[:,w_index],prob_pred)
        auPR[w] = metrics.auc( recall[w],precision[w])        

        if graphical:
            plt.plot(recall[w], precision[w], lw=2, label='J =  {}'.format(w))
            print('Random model AUPR for class '+str(w)+' is: '+str(np.mean(one_hot_enc[:,w_index])))
    if graphical:
        plt.xlabel("recall")
        plt.ylabel("precision")
        plt.legend(loc="best")
        plt.title("precision vs. recall curve")
        plt.text(0.35,0.2,'  AUPRC\nJ= 1: '+str(round(auPR[1],2))+'\nJ= 0: '+str(round(auPR[0],2))+'\nJ=-1: '+str(round(auPR[-1],2))+'\n')

        plt.figure(figsize = (4,3))

    for w_index,w in enumerate([1,0,-1]):
        prob_pred = P_w[w_index,row,col]
        if graphical:
            fpr[w], tpr[w], _ = roc_curve(one_hot_enc[:,w_index],prob_pred)
            plt.plot(fpr[w], tpr[w], lw=2, label='J= {}'.format(w))
        auROC[w] = metrics.roc_auc_score(one_hot_enc[:,w_index], prob_pred)
    if graphical:
        plt.xlabel("false positive rate")
        plt.ylabel("true positive rate")
        plt.legend(loc="best",ncol = 2)
        plt.title("ROC curve")
        plt.plot([0,1],[0,1],ls = '--',c = 'k')
        plt.tight_layout()
        plt.text(0.8,0.2,'  AUROC\nJ= 1: '+str(round(auROC[1],2))+'\nJ= 0: '+str(round(auROC[0],2))+'\nJ=-1: '+str(round(auROC[-1],2))+'\n')
        plt.tight_layout()
    return  auPR,auROC