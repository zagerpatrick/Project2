# Project 2 - Clustering and Drug Discovery
## Due 02/17/2021

[![HW2](https://github.com/zagerpatrick/Project2/actions/workflows/test.yml/badge.svg)](https://github.com/zagerpatrick/Project2/actions/workflows/test.yml)

## API

    ligand_reader(file_name):

        Reads .csv file of the following format and imports individual ligands
        into the Ligand class.
        Format:
        | Ligand ID | Score | SMILES | OnBits |
        |     0     |   10  |  OC#N  | 53,624 |
        Parameters
        ----------
        file_name : string
            Path to .csv file.
        Returns
        ----------
        ligand_list : list of ligand class instances
            List of instances of the ligand class from .csv file.

    
    Ligand():

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
    
    
        list2array(self)
    
            Convert onbits list to onbits binary array.
    
    
    jaccard_dis(x, c):

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


    silhouette_score(x, x_lab):
    
	    Calculates silhouette score.
	    Parameters
	    ----------
	    x  : numpy.array
	        Cluster values.
	    x_lab : numpy.array
	        A numpy array of cluster labels.
	    Returns
	    ----------
	    si : float
	        Per point silhouette score.


    rand_index(x, y):

        Calculates the Rand Index between to arrays of cluster labels.
        Parameters
        ----------
        x  : numpy.array
            A numpy array of cluster labels.
        y : numpy.array
            A numpy array of cluster labels.
        Returns
        ----------
        r_ind : float
            Rand Index value.


	HierarchicalClustering():

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


        _single_linkage(self):

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
        
        
        _hierarch(self):

            Generate hierachical clustering structure from 
            list of branch indices and lengths.
            A helper function for the cluster function in the PartionClustering
            class. Function should not be called directly.
            Returns
            ----------
            dendo : list
                List of dendogram hierarchy.


        _label_gen(self):

            Generates cluster labels.
            A helper function for the cluster function in the PartionClustering
            class. Function should not be called directly.
            Returns
            ----------
            p_lab: np.array
                Categorical cluster labels for cluster membership for data points.


        cluster(self):

            Single-linkage hierarchical clustering.
            Returns
            ----------
            self.p_lab: np.array
                Categorical cluster labels for cluster membership for data points.
        
        
    class PartitionClustering():

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


        _kmeansplusplus(self):

            Centroid initialization algorithm for k-means type clustering as
            described in "k-means++: The Advantages of Careful Seeding".
            A helper function for the cluster function in the PartionClustering
            class. Function should not be called directly.
            Returns
            ----------
            self.c: np.array
                An array of initial centroid as determined by the k-means++.


        _assign_labels(self):

            Calculates cluster membership for each point according to 
            Jaccard/Tanimoto distance.
            A helper function for the cluster function in the PartionClustering
            class. Function should not be called directly.

            Returns
            ----------
            self.p_lab: np.array
                Categorical cluster labels for cluster membership for data points.


        _calculate_centroids(self):

            Calculates centroids accorinding to the mode of cluster members.
            A helper function for the cluster function in the PartionClustering
            class. Function should not be called directly.
            Returns
            ----------
            self.c : np.array
                Cluster centroids.


        cluster(self):

            k-modes clustering function.
            Returns
            ----------
            c : np.array
                Cluster centroids.
            c_lab : np.array
                Categorical cluster labels for cluster membership for centroids.
            p_lab : np.array
                Categorical cluster labels for cluster membership for data points.
