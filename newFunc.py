import pandas as pd 
import networkx as nx 
import numpy as np
import os 
import scipy
from sknetwork.clustering import Louvain, modularity, bimodularity
import sknetwork as skn
from itertools import count
import igraph as ig
import leidenalg as la

import sklearn.metrics

from scipy.stats import bootstrap




def Optimazed_Build_G_with_Treshold(matrix1,Treshold):
    G=nx.algorithms.bipartite.matrix.from_biadjacency_matrix(scipy.sparse.csr_matrix(matrix1))
    to_remove = [(a,b) for a, b, attrs in G.edges(data=True) if attrs["weight"] < Treshold]
    G.remove_edges_from(to_remove)
    
    return G








def GetResol(resolution,matrix1,namefile,Algo,T):  
    
    optimiser = la.Optimiser()
    totSiluette = []
    ClSil = []
    iteration = 100
    #T = 2.0
    N_non_spurius_C=[]

    rng = np.random.default_rng()
    dfSelection = []


    for res in resolution:
        print("resolution :" , res)
        
        #Build graph with Treshold T
        B = Optimazed_Build_G_with_Treshold(matrix1,T)
        
        #RelabelNode
        sample2_name = list(matrix1.index)
        sample1_name= list(matrix1.columns)
        names = sample2_name+ sample1_name
        mapping = dict(zip(B, names))
        B = nx.relabel_nodes(B, mapping)
        
        if(Algo == "Leiden"):
            h = ig.Graph.from_networkx(B)

            p_01, p_0, p_1 = la.CPMVertexPartition.Bipartite(h,resolution_parameter_01=res ,weights='weight',types='bipartite')
            diff = optimiser.optimise_partition_multiplex([p_01, p_0, p_1],layer_weights=[1, -1, -1],n_iterations=-1)
                
                #MemebershipTot
            d = {'node': h.vs['_nx_name'],'clusterOverall' : p_01.membership} 
            df = pd.DataFrame(data= d).set_index('node')
        
         
        if(Algo == "LeidenNormale"):
            h = ig.Graph.from_networkx(B)
            part = la.find_partition(h, la.ModularityVertexPartition)
                #MemebershipTot
            d = {'node': h.vs['_nx_name'],'clusterOverall' : part.membership} 
            df = pd.DataFrame(data= d).set_index('node')
        
        
        
        if(Algo == "Luvain"):
            
        #Clustering Total
            biadjacency =  nx.algorithms.bipartite.matrix.biadjacency_matrix(B, row_order = list(matrix1.index), column_order=list(matrix1.columns))
                
            louvain = Louvain(resolution= res,random_state = 42)
            louvain.fit(biadjacency,force_bipartite=True)
            labels_row = louvain.labels_row_
            labels_col = louvain.labels_col_
            lbz = list(labels_row)+list(labels_col)
            d = {'node': list(B.nodes()),'clusterOverall' : lbz}
            df = pd.DataFrame(data= d).set_index('node')
    
    
    
        labels_tot = df['clusterOverall'].to_numpy() 
        n_cell_tot = len(df)
        columnsz = df.columns[2:]
        numM = np.zeros((n_cell_tot,n_cell_tot))
        denM = np.zeros((n_cell_tot,n_cell_tot))
        
        for ite in range(iteration):
            
            #Sampling 80% of nodes
            H = B.subgraph(list(pd.DataFrame(B.nodes()).sample(frac=0.8)[0])) 
            if(Algo == "Luvain"):  
                top_nodes = {n for n, d in H.nodes(data=True) if d["bipartite"] == 0}
                bottom_nodes = set(H) - top_nodes
                biadjacency_temp =  nx.algorithms.bipartite.matrix.biadjacency_matrix(H, row_order = list(top_nodes),    column_order=list(bottom_nodes))
            
            #Clustering on sample matrix 
                louvain = Louvain(resolution= res,random_state = 42)
                louvain.fit(biadjacency_temp,force_bipartite=True)
        
                labels_row = louvain.labels_row_
                labels_col = louvain.labels_col_
                lbz = list(labels_row)+list(labels_col)

                d_temp = {'node': list(H.nodes()),'clusterIt'+str(ite) : lbz}
                
                
            if(Algo == "Leiden"):
                h = ig.Graph.from_networkx(H)

                p_01, p_0, p_1 = la.CPMVertexPartition.Bipartite(h,resolution_parameter_01=res ,weights='weight',types='bipartite')
                diff = optimiser.optimise_partition_multiplex([p_01, p_0, p_1],layer_weights=[1, -1, -1],n_iterations=-1)
                
                #MemebershipTot
                d_temp = {'node': h.vs['_nx_name'],'clusterIt'+str(ite) : p_01.membership}
                 
            if(Algo == "LeidenNormale"):
                h = ig.Graph.from_networkx(H)
                part = la.find_partition(h, la.ModularityVertexPartition)
                #MemebershipTot
                d_temp = {'node': h.vs['_nx_name'],'clusterIt'+str(ite) : part.membership} 
                
        
        
            df_temp = pd.DataFrame(data= d_temp).set_index('node')
            df = pd.concat([df,df_temp],axis = 1)
        
        
            a = df['clusterIt'+str(ite)].to_numpy()
            b = df['clusterIt'+str(ite)].to_numpy()
            j_temp = a[:, None] == b
            j_temp  = np.where(j_temp==True, 1, 0)



            a = df['clusterIt'+str(ite)].to_frame()
            a[a>=0] = 1
            a = a.fillna(0)
            a = a['clusterIt'+str(ite)].to_numpy()
            b = a 
            ll_temp=a[:, None] + b
            ll_temp[ll_temp <= 1]= 0
            ll_temp[ll_temp ==2.] = 1
        
            numM = numM + j_temp
            denM = denM + ll_temp
        
        
        N_ofC_tot = len(df["clusterOverall"].value_counts())
        N_ofC = len(df["clusterOverall"].value_counts()[df["clusterOverall"].value_counts() >1])
        D1 = np.divide(numM, denM)
        D = 1-D1 
    	
        pd.DataFrame(D).to_csv("ReM/"+Algo+"/residualmatrix-"+namefile+"-"+str(res)+".csv")
        df.to_csv("ReM/"+Algo+"/cluster-"+namefile+"-"+str(res)+".csv")
        
    return 0
    


def getResol1(namefile,resolution):
    
    rng = np.random.default_rng()
    dfSelection = []

    for res in resolution:
    
    
        D = pd.read_csv("ReM/Luvain/residualmatrix-"+namefile+"-"+str(res)+".csv",index_col=[0]).to_numpy()
        df =pd.read_csv("ReM/Luvain/cluster-"+namefile+"-"+str(res)+".csv")[["node","clusterOverall"]].set_index("node")
    
        countV = df["clusterOverall"].value_counts() 
        countV = countV[countV ==1] 
        listrename = list(countV.index)
        df["clusterOverall"] = df["clusterOverall"].replace(listrename,"unasigned")
        labels_tot = list(df["clusterOverall"])
        df["Silux"] = sklearn.metrics.silhouette_samples(D, labels_tot)
        df1 = df.groupby(["clusterOverall"])["Silux"].median().to_frame()
        df1 = df1.drop(["unasigned"])
        N_ofC_tot = len(df["clusterOverall"].value_counts())
        N_ofC = len(df["clusterOverall"].value_counts()[df["clusterOverall"].value_counts() >1])
        data = (df1.to_numpy().flatten(),)  
        medianz = np.median(df1.to_numpy().flatten())
        boot = bootstrap(data, np.median, confidence_level=0.95,random_state=rng,n_resamples=25000)
        lower_b = boot.confidence_interval[0]
        upper_b = boot.confidence_interval[1]
        oUtofCI = list(df1[~df1['Silux'].between(lower_b,upper_b)]["Silux"])
        NCininterval = len(df1) - len(oUtofCI)
        NMI = 0
        dfSelection.append([res,N_ofC,medianz,lower_b,N_ofC_tot,NMI,upper_b,oUtofCI,NCininterval])
    
        
    dfg = pd.DataFrame(dfSelection)
    dfg = dfg.loc[dfg[2] >max(dfg[3]) ]
    dfg = dfg.loc[dfg[8] ==max(dfg[8]) ] 
    residual = dfg[0].values
    
    return residual,dfSelection



