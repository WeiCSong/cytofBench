from sklearn.neighbors import kneighbors_graph, BallTree
from sklearn.utils.graph import graph_laplacian
import scipy.sparse
import numpy as np
import time

def add_value(cor, key, val):
    if key not in cor:
        cor[key] = 0.0
    cor[key] += val 
            
def dict_2_cs(cor):
    data, indices, indptr = [], [], [0]
    ptr_curr = 0
    for key, item in cor.items():
        indices += item.keys()
        data += item.values()
        ptr_curr += len(item)
        indptr.append(ptr_curr)
        
    return data, indices, indptr

def knn_neighbors(data, lpts, n_neighbor):
    X_semi = np.concatenate([lpts, data], axis = 0)
    btree = BallTree(X_semi)
    nbr_id = btree.query(X_semi, k=n_neighbor + 1, return_distance=False)
    return nbr_id
    
def compute_BTLu(y_lpt, nbr_id, n_sample, n_type, n_lpt):
    
    cor_Lu = {i:{} for i in range(n_sample)}
    cor_BT = {i:{} for i in range(n_type)}

    for i in range(n_lpt):
        for j in nbr_id[i, :]:
            if j >= n_lpt:
                add_value(cor_BT[y_lpt[i]], j - n_lpt, -1)
                add_value(cor_Lu[j - n_lpt], j - n_lpt, 1)

    for i in range(n_sample):
        for j in nbr_id[n_lpt + i, 1:]:
            add_value(cor_Lu[i], i, 1)
            if j < n_lpt:
                add_value(cor_BT[y_lpt[j]], i, -1)
                
            elif j >= n_lpt:
                add_value(cor_Lu[j - n_lpt], i, -1)
                add_value(cor_Lu[i], j - n_lpt, -1)
                add_value(cor_Lu[j - n_lpt], j - n_lpt, 1)

    BT = scipy.sparse.csc_matrix(dict_2_cs(cor_BT), shape=(n_sample, n_type), dtype = np.float)
    Lu = scipy.sparse.csc_matrix(dict_2_cs(cor_Lu), shape=(n_sample, n_sample), dtype = np.float)
    return BT, Lu

def rm_classify(data, lpts, y_lpt, n_neighbor):
    ## data: (n_sample, n_feature) ndarray. unlabeled data
    ## lpts: (n_lpt, n_feature) ndarray. labeled data
    ## y_lpts: (n_lpt, ) ndarray, labels for the labeled data
    ## n_neighbors: number of neighbors
    ## labels: order of labels
    
    
    n_lpt = lpts.shape[0]
    n_type = len(set(y_lpt))
    n_sample = data.shape[0]
    
    idx2lpt =sorted(set(y_lpt))
    lpt2idx = {l:idx for idx, l in enumerate(idx2lpt)}
    
    y_lpy_tmp = np.array([lpt2idx[l] for l in y_lpt])     
    
    nbr_id = knn_neighbors(data, lpts, n_neighbor)

    BT, Lu = compute_BTLu(y_lpy_tmp, nbr_id, n_sample, n_type, n_lpt)
    
    lp = np.zeros((n_sample, n_type)) 
    for i in range(n_type):
        result = scipy.sparse.linalg.bicg(Lu, BT[:, i].toarray(), tol = 1e-6)
        lp[:, i] = -result[0]
        
    
    y_pred = np.argmax(lp, axis = 1)
    y_pred = np.array([idx2lpt[l] for l in y_pred])
    return lp, y_pred


