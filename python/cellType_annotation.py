import pandas as pd
import numpy as np

from sklearn import preprocessing
#from sklearn.mixture import GMM
from sklearn.mixture import GaussianMixture as GMM
from scipy.stats import norm
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as sch

from collections import Counter



def get_label(table):
    return table.axes[0], table.axes[1]


def sort_feature(data):
    pds = np.asarray(data)
    Y = pdist(pds)  # ~ N^2 / 2
    ind = sch.leaves_list(sch.linkage(Y))  # N-1 
    return data[ind, :]


def compute_marker_model(df, table, thres):
    ## compute 2-gaussian mixture model for each marker
    ## output mean, variance, and weight for each marker
    mk_model = {}
    for mk in table.axes[1]:
        gmm = GMM(n_components=2, n_init=10)
        tmp = df[mk].as_matrix()
        gmm.fit(tmp[tmp > thres, np.newaxis])
        index = np.argsort(gmm.means_, axis = 0)    
        mk_model[mk] = (gmm.means_[index], gmm.covariances_[index], gmm.weights_[index])
        
    return mk_model

def appr_fun(x, xroot, slope):
    ## heuristic for axxproimating the decision boundary
    y = np.exp( (x - xroot) * slope)
    y = y / (1.0 + y)
    return np.squeeze(y)

def compute_paras(mk_model, mk):
    ## compute the critical point of the decision bounary
    ## critical point is defined by P(theta = 1|x) = P(theta = 0|x) 
    ## given a 1-D two-mixture model
    ##
    ## inputs:
    ## mk_model: precomputed 2-gassian mixture model for all markers
    ## mk: one selected marker
    ##
    ## outputs:
    ## xroot: location of the critical point
    ## slope: slope of the decision boundary at the critical point
    
    mus, sigmas, ws = mk_model[mk]
    a = (-0.5 / sigmas[0] + 0.5 /sigmas[1])
    b = mus[0] / sigmas[0] - mus[1] / sigmas[1]
    c = 0.5 * (-mus[0] **2 / sigmas[0] + mus[1] ** 2 / sigmas[1]) + np.log(ws[0] / ws[1]) + 0.5 * np.log(sigmas[1] / sigmas[0])
    xroot = (-b - np.sqrt(b ** 2 - 4.0 * a * c) ) / (2.0 * a)
    slope = 0.5 * (xroot - mus[0]) / sigmas[0] - 0.5 * (xroot - mus[1]) / sigmas[1]
    return xroot, slope

def score_fun(x, mk_model, mk):
    ## compute the score function
    ## the score should roughly proportional to P(theta = + | x)

    xroot, slope = compute_paras(mk_model, mk)
    return appr_fun(x, xroot, slope)



def get_score_mat(X, weights, table, y_cluster, mk_model):
    ## compute the cluster x annotation score matrix
    ## X: sample x feature data matrix
    ## table: a list of (ind, sign)
    if y_cluster:
        clusters = np.unique(y_cluster)
        mu = np.vstack([np.mean(X[y_cluster==cl, :], axis = 0) for cl in clusters])
        score = get_score_mat_main(mu, weights, table, mk_model, thres)
    else:
        score = get_score_mat_main(X, weights, table, mk_model)
    return score

def get_score_mat_main(X, weights, table, mk_model):
    ## compute the cluster x annotation score matrix
    ## X: sample x feature data matrix
    ## table: a list of (ind, sign)
    score = np.zeros((X.shape[0], len(table.axes[0])))
    for i, ct in enumerate(table.index):
        #score_tmp = np.zeros(X.shape[0])
        #count = 0.0
        score_tmp = np.ones(X.shape[0])
        count = 1.0
        for j, mk in enumerate(table.columns):
            if table.loc[ct, mk] > 0:
                score_tmp =  np.min([score_tmp, score_fun(X[:, j], mk_model, mk)], axis = 0)
                #score_tmp +=  score_fun(X[:, j], mk_model, mk)
                #count += 1.0
            elif  table.loc[ct, mk] < 0:
                score_tmp =  np.min([score_tmp, 1.0 - score_fun(X[:, j], mk_model, mk)], axis = 0)
                #score_tmp +=  1.0 - score_fun(X[:, j], mk_model, mk)
                #count += 1.0
        score[:, i]= score_tmp / count
    return score



def get_unique_index(X, score, table, thres):
    
    ct_score = np.abs(table.as_matrix()).sum(axis = 1)
    ct_index = np.zeros((X.shape[0], len(table.index) + 1))
    
    for i, ct in enumerate(table.index):
        ct_index[:, i] = (score[:, i] > thres) * 1

    for i in np.where(ct_index.sum(axis = 1) > 1)[0]:
        mk_tmp = np.where(ct_index[i, :] == 1)[0]
        score_tmp = ct_score[mk_tmp]
        ct_index[i, :] = 0
        if np.sum(score_tmp == score_tmp.max()) == 1:
            ct_index[i, mk_tmp[np.argmax(score_tmp)]] = 1    

    ct_index[:, -1] = (score[:, -1] > thres) * 1
    return ct_index

def get_landmarks(X, score, ct_index, idx2ct, obj, thres):

    res_c = {}
    for idx, ct in enumerate(idx2ct):
        if np.any(ct_index[:, idx] == 1):
            item = X[ct_index[:, idx] == 1, :]
        else:
            item = X[score[:, idx] > thres, :]

        if item.shape[0] > 60:
            res = obj.cluster(item, k=30, directed=False, prune=False, min_cluster_size=10, jaccard=True,
                    primary_metric='euclidean', n_jobs=-1, q_tol=1e-3)
            res_c[ct] = return_center(item, res[0])
        elif item.shape[0]> 0:
            res_c[ct] = item.mean(axis = 0)[np.newaxis, :]

    return res_c


def soft_max(x):
    # x: sample x feature matrix
    # w: 1 x weights
    return 1.0 / (1.0 + np.exp(-x))


def return_center(X, labels):
    clusters = np.unique(labels)
    centers = []
    for i in clusters:
        centers.append(np.mean(X[labels == i, :], axis = 0))
    return np.vstack(centers)

def select_centers(res_center, X, table, mk_model):
    # X: sample x feature data matrix
    # table: a list of (ind, sign)
    
    res_score = {}
    res_select = {} 
                
    for ct, item in res_center.items():
        if item.shape[0] > 1:
            score = np.ones(item.shape[0])
            for j, mk in enumerate(table.axes[1]):
                if table.loc[ct, mk]== 1:
                    score = np.min([score, score_fun(item[:, j], mk_model, mk)], axis = 0)
                if table.loc[ct, mk]== -1:
                    score = np.min([score, 1.0 - score_fun(item[:, j], mk_model, mk)], axis = 0)
                    
            res_score[ct] = score
            res_select[ct] = res_center[ct][np.argmax(score), :][np.newaxis, :]
        else:
            res_select[ct] = res_center[ct].copy()
    return res_score, res_select

def output_feature_matrix(res, label=None):
    ## turn a dictionary data structure into a matrix
    ## input:
    ## res: a dictionary that map a key (cluster) to an (n_sample, n_features) array
    ## label: a list of keys. output will be arragned as the order of keys in this list
    ## output:
    ## feature_mat: a (n_sample, n_features) matrix
    ## ct_mat: a (n_sample,) that each element is a key (cluster)
    
    ct_mat = []
    feature_mat = []
    if not label:
        label = res.keys()
        
    for key in label:
        if key in res:
            item = res[key]
            if item.shape[0] > 0:
                feature_mat.append(item)
                ct_mat += [key] * item.shape[0]
                
    feature_mat = np.vstack(feature_mat)
    return feature_mat, ct_mat
