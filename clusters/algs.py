from csv import reader
import random
import numpy as np
from scipy.stats import mode


def ligand_reader(file_name):
    '''
    Reads .csv file of the following format and imports individual ligands
    into the Ligand class.

    Format:
    | Ligand ID | Score | SMILES | OnBits |
    |     0     |   10  |  OC#N  | 53,624 |

    Parameters
    ----------
    ligand_id : string
        Path to .csv file.

    Returns
    ----------
    ligand_list : list of ligand class instances
        List of instances of the ligand class from .csv file.
    '''
    ligand_list = []

    with open(file_name) as file:
        next(file)
        for line in reader(file):
            ligand_id = int(line[0])
            score = float(line[1])
            smiles = line[2]
            onbits = [int(i) for i in line[3].split(',')]
            ligand_list.append(Ligand(ligand_id, score, smiles, onbits))

    return ligand_list


class Ligand():
    '''
    Class for ligand data containing id's, docking scores, molecular id's,
    and ECFP onbits.
    
    This class also contains the function to transform onbits from sting 
    format to np.array format.
    
    Parameters
    ----------
    ligand_id : int
        Numerical indentification of ligands. 
    score : float
        Autodock Vina score.
    smiles : string
        A chemical species string in simplified molecular-input 
        line-entry system (SMILES) format.
    onbits : list of ints
        Extended-Connectivity Fingerprints (ECFP) characterization.

    Attributes
    ----------
    ligand_id : int
        Numerical indentification of ligands. 
    score : float
        Autodock Vina score.
    smiles : string
        A chemical species string in simplified molecular-input 
        line-entry system (SMILES) format.
    onbits : list of ints
        Extended-Connectivity Fingerprints (ECFP) characterization.
    onbits_array : np.array
        Long form ECFP binary array of shape [n, 1024].

    References
    ----------
    1. Weininger, David. "SMILES, a chemical language and information system. 
    1. Introduction to methodology and encoding rules." Journal of chemical 
    information and computer sciences 28.1 (1988): 31-36.
    
    2. Rogers, David, and Mathew Hahn. "Extended-connectivity fingerprints."
    Journal of chemical information and modeling 50.5 (2010): 742-754.

    Examples
    --------
    >>> Ligand.ligand_id 
    8 
    >>> Ligand.score
    -2.6
    >>> Ligand.smiles
    'O=C(C)C'
    >>> Ligand.onbits
    [4, 6, 33, 397, 650, 807, 893, 1017]
    >>> ligand_list = ligand_reader('ligand_information.csv')
    >>> [ligand.list2array() for ligand in ligand_list]
    >>> Ligand.onbits_array[0, 0:10]
    array([0., 0., 0., 0., 1., 0., 1., 0., 0., 0.])
    '''
    def __init__(self, ligand_id, score, smiles, onbits):
        self.ligand_id = ligand_id
        self.score = score
        self.smiles = smiles
        self.onbits = onbits
        self.onbits_array = []
    
    def list2array(self):
        '''
        Convert onbits list to onbits binary array.
        '''
        onbits_array = np.zeros(1024)
        onbits_array[self.onbits] = 1
        self.onbits_array = onbits_array


def jaccard_dis(x, c):
    '''
    Calculates Jaccard/Tanimoto distance between two numpy arrays using 
    matrix algebra. Note, input arrays should have identical n sizes 
    (.shape[1]).

    Parameters
    ----------
    x  : numpy.array
        A numpy array.
    c : numpy.array
        A second numpy array to calculate the distance from/to.

    Returns
    ----------
    score : lists of numpy.float64
        An array of Jaccard/Tanimoto distance of shape 
        x.m, c.m (x.shape[0], c.shape[0]).
    '''

    n = x.shape[1]

    m11 = np.dot(x, c.T)
    m00 = np.dot(np.abs(x - 1), np.abs(c.T - 1))
    jd = 1 - (m11/(n-m00))

    return jd


class HierarchicalClustering():
    '''
    Single-linkage hierarchical clustering.
    Cluster members are determined through Jaccard/Tanimoto distance.
    
    Parameters
    ----------    
    x : np.array
        Values to be clustered. Values in array must be either zero or one.
    thresh : int
        Threshold value for determining cluster membership.
        
    Attributes
    ----------
    x : np.array
        Values to be clustered. Values in array must be either zero or one.
    thresh : int
        Threshold value for determining cluster membership.
    dis : np.array
        Jaccard/Tanimoto distance matrix.
    min_ind_list : list
        List of step-wise clustering.
    bl_array : np.array
        Array of dendogram branch lengths.
    dendo : list
        List of dendogram hierarchy.
    p_lab : np.array
        Categorical cluster labels for cluster membership for data points.

    References
    ----------
    1. Gower, John C., and Gavin JS Ross. 
    "Minimum spanning trees and single linkage cluster analysis." 
    Journal of the Royal Statistical Society: 
    Series C (Applied Statistics) 18.1 (1969): 54-64.

    Examples
    --------
    >>> SLC = algs.HierarchicalClustering(test, 0.42)
    >>> p_lab = SLC.cluster()
    >>> p_lab[0:10]
    array([0, 1, 1, 1, 2, 0, 0, 0, 1, 0])
    '''
    
    def __init__(self, x, thresh):
        self.x = x
        self.thresh = thresh
        self.dis = np.array([])
        self.min_ind_list = []
        self.bl_array = np.array([])
        self.dendo = []
        self.p_lab = np.array([])
       
    def _single_linkage(self):
        '''
        Generate lists of branch indices and lengths 
        through single-linkage clustering.
        A helper function for the cluster function in the PartionClustering
        class. Function should not be called directly.

        Returns
        ----------
        min_ind_list : list
            List of step-wise clustering.
        bl_array : np.array
            Array of dendogram branch lengths.
        '''

        # Determine dimmensions of input matrix
        m, n = self.dis.shape[0], self.dis.shape[1]

        # Initialize lists of branch indices and lengths
        bl_list = [0]
        min_ind_list = []

        for i in range(0, m-1):
            # Determine the index of the minimum value 
            # on the right side of the diagonal
            ind = np.where(self.dis == self.dis.min())
            min_ind = np.array([ind[0][0], ind[1][0]])

            # Calculate the branch length
            bl = np.min(self.dis)/2

            # Store values
            min_ind_list.append(min_ind)
            bl_list.append(bl)

            # Merge clusters
            up = np.min(self.dis[min_ind], axis=0)

            # Update distance matrix
            self.dis[min_ind[0], :] = up
            self.dis[:, min_ind[0]] = up
            self.dis = np.delete(self.dis, min_ind[1], 0)
            self.dis = np.delete(self.dis, min_ind[1], 1)
            np.fill_diagonal(self.dis, np.inf)

        # Convert branch length list to np.array
        bl_array = np.array(bl_list)
        
        return min_ind_list, bl_array
    
    
    def _hierarch(self):
        '''
        Generate hierachical clustering structure from 
        list of branch indices and lengths.
        A helper function for the cluster function in the PartionClustering
        class. Function should not be called directly.

        Returns
        ----------
        dendo : list
            List of dendogram hierarchy.
        '''
        # Generate first list of labels, where each point is it's own cluster
        label_list = [[i] for i in range(0, len(self.min_ind_list)+1)]
        # Inialize list of lists of labels
        dendo = [label_list.copy()]

        # Generate list of hierarchical structures.
        for j,k in self.min_ind_list:
            node_list = label_list.copy()

            # Update labels given new cluster assignments
            label_list[j] = label_list[j] + node_list[k]

            # Remove labels that have been merged
            label_list.pop(k)
            dendo.append(label_list.copy())

        return dendo


    def _label_gen(self):
        '''
        Generates cluster labels.
        A helper function for the cluster function in the PartionClustering
        class. Function should not be called directly.

        Returns
        ----------
        p_lab: np.array
            Categorical cluster labels for cluster membership for data points.
        '''
        # Determine clustes given threshold
        cluster_ind = np.where(self.bl_array <= self.thresh)[0][-1]
        label_list = self.dendo[cluster_ind] 
        
        # Initialize labels array
        p_lab = np.zeros(len(self.dendo))

        for n, cluster in enumerate(label_list):
            for p in cluster:
                p_lab[p] = n

        return p_lab
    
    
    def cluster(self):
        '''
        Single-linkage hierarchical clustering.
        
        Returns
        ----------
        self.p_lab: np.array
            Categorical cluster labels for cluster membership for data points.
        '''
        # Handeling of negative thresh input values
        if self.thresh < 0: self.thresh = 0

        # Calculate distance matrix
        self.dis = jaccard_dis(self.x, self.x)
        np.fill_diagonal(self.dis, np.inf)

        # Generate lists of branch indices and lengths
        self.min_ind_list, self.bl_array = self._single_linkage()

        # Generate hierachical clustering structure
        self.dendo = self._hierarch()
        
        # Determine cluster labels
        self.p_lab = self._label_gen()

        return self.p_lab

class PartitionClustering():
    '''
    k-modes clustering. Centroids are initialized though k-means++.
    Cluster members are determined through Jaccard/Tanimoto distance.
    
    Parameters
    ----------    
    x : np.array
        Values to be clustered. Values in array must be either zero or one.
    n_clusters : int
        Number of clusters.

    Attributes
    ----------
    x : np.array
        Values to be clustered. Values in array must be either zero or one.
    n_clusters : int
        Number of clusters.
    c : np.array
        Cluster centroids.
    c_lab : np.array
        Categorical cluster labels for cluster membership for centroids.
    p_lab : np.array
        Categorical cluster labels for cluster membership for data points.

    References
    ----------
    1. Arthur, David, and Sergei Vassilvitskii. k-means++: 
    The advantages of careful seeding. Stanford, 2006.

    Examples
    --------
    >>> kmodes = PartitionClustering(ligand_array, n_clusters=5)
    >>> c, c_lab, p_lab = kmodes.cluster()
    >>> c_lab
    array([0, 1, 2, 3, 4])
    >>> p_lab[0:10]
    array([0, 1, 1, 2, 1, 2, 2, 1, 1, 1])
    '''

    def __init__(self, x, n_clusters):
        self.x = x
        self.n_clusters = n_clusters
        self.c = np.array([])
        self.c_lab = np.array([])
        self.p_lab = np.array([])

    
    def _kmeansplusplus(self):
        '''
        Centroid initialization algorithm for k-means type clustering as
        described in "k-means++: The Advantages of Careful Seeding".
        A helper function for the cluster function in the PartionClustering
        class. Function should not be called directly.
        
        Returns
        ----------
        self.c: np.array
            An array of initial centroid as determined by the k-means++.
        
        '''
        # Select a data point at random to be the first centroid
        rand_int = random.randrange(len(self.x))
        self.c = self.x[rand_int: rand_int + 1]

        # Determine dimmensions of input matrix
        m, n = self.x.shape[0], self.x.shape[1]

        # Determine other centroids through a 
        # distance-based probablity distribution 
        for i in range(0, self.n_clusters-1):

            # Calculate distances to centroids
            jd = jaccard_dis(self.x, self.c)

            # Calculate the minimum distance to a centoid for all data points
            if i > 0:
                jd = np.min(jd, axis=1)

            # Calculate weighted distance probability and standardize values
            prob = (jd**2)/(np.sum(jd**2))
            c_ind = np.random.choice(np.arange(0,m), 1, p=prob.ravel())
            self.c = np.concatenate((self.c, self.x[c_ind[0]:c_ind[0]+1]))

        return self.c


    def _assign_labels(self):
        '''
        Calculates cluster membership for each point according to 
        Jaccard/Tanimoto distance.
        A helper function for the cluster function in the PartionClustering
        class. Function should not be called directly.
        
        Returns
        ----------
        self.p_lab: np.array
            Categorical cluster labels for cluster membership for data points.
        '''

        # Calculate distances to centroids
        jd = jaccard_dis(self.x, self.c)

        # Determine cluster labels
        self.p_lab = np.argmin(jd, axis=1)
        n_labels = len(np.unique(self.p_lab))

        # If a cluster contatains no points, 
        # randomly move points to the empty cluster
        if not n_labels == self.n_clusters:
            unq = np.unique(p_lab)
            no_lab = np.setdiff1d(c_lab, unq) # clusters without points
            for i in no_lab:
                self.p_lab[random.randrange(len(self.p_lab))] = i

        return self.p_lab


    def _calculate_centroids(self):
        '''
        Calculates centroids accorinding to the mode of cluster members.
        A helper function for the cluster function in the PartionClustering
        class. Function should not be called directly.
        
        Returns
        ----------
        self.c : np.array
            Cluster centroids.
        '''

        # Calculate new centoids
        for k in range(0, self.n_clusters):
            ind = np.concatenate(np.argwhere(self.p_lab == k))
            self.c[k,:] = mode(self.x[ind], axis=0)[0]

        return self.c


    def cluster(self):
        '''
        k-modes clustering function.
        
        Returns
        ----------
        c : np.array
            Cluster centroids.
        c_lab : np.array
            Categorical cluster labels for cluster membership for centroids.
        p_lab : np.array
            Categorical cluster labels for cluster membership for data points.
        '''
        # Determine dimmensions of input matrix
        m, n = self.x.shape[0], self.x.shape[1]

        # Initailize centroids matrix with kmeans++
        self.c = self._kmeansplusplus()
        self.c_lab = np.arange(0, self.n_clusters)

        iteration_count = 0
        c_old = np.zeros([self.n_clusters, n])

        while not np.array_equal(c_old, self.c):
            # Variable setting
            iteration_count += 1
            c_old = np.copy(self.c)

            ## Assign labels
            self.p_lab = self._assign_labels()

            # Calculate centroids
            self.c = self._calculate_centroids()

        return self.c, self.c_lab, self.p_lab
