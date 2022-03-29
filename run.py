#Before make shure yuo run PvalueM_HGT.R



import pandas as pd 
import networkx as nx 
import scipy
from sklearn.metrics.cluster import normalized_mutual_info_score
from sknetwork.clustering import Louvain, modularity, bimodularity
from sknetwork.linalg import normalize
import sklearn
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
from newFunc import * 


import scanpy as sc



adataB = sc.read_h5ad('DataPerPy/SeuratObj/BT346.h5ad')
adataS = sc.read_h5ad('DataPerPy/SeuratObj/BT400.h5ad')

print(len(list(adataB.obs.index)))
print(len(list(adataS.obs.index)))

adataB = adataB[adataB.obs["celltypesign"] != "unasigned"]
adataS = adataS[adataS.obs["celltypesign"] != "unasigned"]

l1 = list(adataB.obs.index)
l2 = list(adataS.obs.index)
l1 = list(pd.Series(l1).str.replace("-","."))
matrix1 = pd.read_csv("matrixPvalue/BT400-BT346-dimz-50.csv",index_col = [0])

matrix1 = matrix1.loc[l2]
matrix1 = matrix1.loc[:, matrix1.columns.isin(l1)]


resolution = 0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,3.0,4.0 #[2.0,3.0,4.0,5.0]
GetResol(resolution,matrix1,"BT400"+"-"+"BT346","Louvain",2.0)









GetResol(resolution,matrix1,"BT400"+"-"+"BT346","Louvain",2.0)


