import numpy as np


def jaccard_dis(x, c):
    '''
    Calculates Jaccard distance between two numpy arrays using matrix algebra.
    Note, input arrays should have identical n sizes (.shape[1]).

    Parameters
    ----------
    x  : numpy.array
        A numpy array.
    c : numpy.array
        A second numpy array to calculate the distance from/to.

    Returns
    ----------
    score : lists of numpy.float64
        An array of Jaccard distance of shape 
        x.m, c.m (x.shape[0], c.shape[0]).
    '''

    tot = x.shape[1]

    m11 = np.dot(x, c.T)
    jd = 1 - (m11/tot)
    
    return jd


class HierarchicalClustering():
	pass

class PartitionClustering():
    '''
    ***Under Construction***
    K-means clustering
    '''
    
    def kmeans(data, n_clusters):
    
        # Determine dimmensions of input matrix
        m, n = data.shape[0], data.shape[1]

        # Initailize matrix of random centroids
        c = np.random.randint(0, 2, (n_clusters, n)).astype('float64')

        # Calculate distances to centroids
        jd = jaccard_dis(data, c)

        # Determine cluster labels
        p_lab = np.argmax(jd, axis=1) + 1
        c_lab = np.arange(0,3)
    
        return p_lab, c_lab

    pass
