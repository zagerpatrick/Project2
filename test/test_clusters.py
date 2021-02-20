import pytest
import numpy as np
import clusters.algs as algs

def test_hierarchical():
	# Setup
    thresh = 0.42
    
    x_input_hc = np.array([[0., 0., 1., 0., 0.],
                        [1., 0., 0., 0., 0.],
                        [0., 1., 0., 1., 1.],
                        [0., 1., 0., 0., 0.],
                        [0., 0., 0., 1., 0.]])
    
    desired_p_lab_hc = np.array([0, 1, 2, 2, 2])

    # Exercise
    HC = algs.HierarchicalClustering(x_input_hc, thresh)
    p_lab_hc = HC.cluster()

    # Verify
    np.testing.assert_array_equal(p_lab_hc, desired_p_lab_hc)


def test_partitioning():
    # Setup
    n_clusters = 1
    
    x_input_km = np.array([[0., 0., 1., 0., 0.],
                        [1., 0., 0., 0., 0.],
                        [0., 1., 0., 1., 1.],
                        [0., 1., 0., 0., 0.],
                        [0., 0., 0., 1., 0.]])
    
    desired_p_lab_km = np.array([0, 0, 0, 0, 0])

    # Exercise
    kmodes = algs.PartitionClustering(x_input_km, n_clusters)
    c, c_lab, p_lab_km = kmodes.cluster()
    

    # Verify
    np.testing.assert_array_equal(p_lab_km, desired_p_lab_km)