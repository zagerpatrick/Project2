from csv import reader
import numpy as np


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
    ligand_id_ : int
        Numerical indentification of ligands. 
    score_ : float
        Autodock Vina score.
    smiles_ : string
        A chemical species string in simplified molecular-input 
        line-entry system (SMILES) format.
    onbits_ : list of ints
        Extended-Connectivity Fingerprints (ECFP) characterization.
    onbits_array : np.array
        Long form ECFP binary array of shape [n, 1024].
        
    Notes
    -----
    ...

    References
    ----------
    1. Weininger, David. "SMILES, a chemical language and information system. 
    1. Introduction to methodology and encoding rules." Journal of chemical 
    information and computer sciences 28.1 (1988): 31-36.
    
    2. Rogers, David, and Mathew Hahn. "Extended-connectivity fingerprints."
    Journal of chemical information and modeling 50.5 (2010): 742-754.


    Examples
    --------
    >>> Ligand.ligand_id_ 
    8 
    >>> Ligand.score_
    -2.6
    >>> Ligand.smiles_
    'O=C(C)C'
    >>> Ligand.onbits_
    [4, 6, 33, 397, 650, 807, 893, 1017]
    >>> Ligand.onbits_array_[0, 0:10]
    array([0., 0., 0., 0., 1., 0., 1., 0., 0., 0.])
    '''
    def __init__(self, ligand_id, score, smiles, onbits):
        self.ligand_id_ = ligand_id
        self.score_ = score
        self.smiles_ = smiles
        self.onbits_ = onbits
        self.onbits_array_ = []
    
    def list2array(self):
        '''
        Convert onbits list to onbits binary array.
        '''
        onbits_array = np.zeros(1024)
        onbits_array[self.onbits_] = 1
        self.onbits_array_ = onbits_array


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
    ***Under Construction***j
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
