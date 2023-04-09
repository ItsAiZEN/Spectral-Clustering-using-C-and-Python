"""
    This project implements the spectral clustering algorithm.
    Most of the code is written in C in order to achieve better performance.
    Both in C and in Python parts of the code can be run, in C there are the following functions:
    1. wam - weighted adjacency matrix
    2. ddg - degree diagonal matrix
    3. gl - graph laplacian
    4. jacobi - jacobi algorithm

    In Python there are the following functions:
    1. spkmeans - spectral k-means algorithm
    2. all the functions from C

    Files and inputs are given in the following format:
    1. C code - the first argument is the goal (wam, ddg, gl, jacobi), the second argument is the file name
    2. Python code - the first argument is the number of clusters (optional), second is the goal (wam, ddg, gl, jacobi, spk), third is the file name

    The input file is a csv file with the following format:
    1. each line represents a vector, each value is separated by a comma
    2. file ends with a new line

    Output for jacobi goal is a matrix with the following format:
    1. first line is the eigenvalues
    2. the rest of the lines are the eigenvectors (as columns)

    Output for spk goal has the following format:
    1. first line is the observation indices of the chosen vectors (k eigenvectors as columns chosen by the algorithm) for initial centroids
    2. the rest of the lines are the centroids

    all the other goals output a matrix with the following format:
    1. each line represents a vector, each value is separated by a comma

    The code was written in ANSI-C and Python 3.
"""

import math
import sys
import pandas as pd
import numpy as np
import mykmeanssp

# to import mykmeanssp, have spkmeans.py, mykmeanssp.c, spkmeansmodule.c, setup.py in the same file,and run in terminal:
# python setup.py build_ext --inplace
# if it doesn't work, try:
# python3 setup.py build_ext --inplace

# TODO: make the code more elegant
# TODO: final tests and documentation

# seed set for uniformity, can be removed if needed
np.random.seed(0)

if len(sys.argv) == 3:
    goal = sys.argv[1]
    file = open(sys.argv[2], "r")
    # in this case k will be denoted from eigengap_heuristic
elif len(sys.argv) == 4:
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file = open(sys.argv[3], "r")
else:
    print("An Error Has Occurred")
    exit()

# read the csv file into a pandas dataframe
vectors1 = pd.read_csv(file, header=None)
vector_dimension = vectors1.shape[1]
num_of_vectors = vectors1.shape[0]
file.close()

# convert the pandas dataframe into a list of lists
vectors = vectors1.values.tolist()

# checks the validity of the input
if len(sys.argv) == 4:
    if k >= num_of_vectors or k <= 1:
        print("An Error Has Occurred")
        exit()


def euclidean_distance(vector1, vector2):
    """calculates the euclidean distance between two vectors"""
    total_sum = 0
    for i in range(len(vector1)):
        total_sum += (vector1[i] - vector2[i]) ** 2
    return math.sqrt(total_sum)


def kmeans_pp(u_matrix1, k_dimension):
    """performs kmeans++ algorithm"""
    # make u_matrix a pandas dataframe
    u_matrix = pd.DataFrame(u_matrix1)
    # add a left column of indices to the dataframe
    u_matrix.insert(0, "index", range(0, len(u_matrix)))
    original_vec_list = np.array(u_matrix)
    keys = original_vec_list[:, 0]
    centroids = []
    chosen_index = np.random.choice(keys)
    chosen_vector = original_vec_list[np.where(original_vec_list[:, 0] == chosen_index)[0][0], :]
    centroids.append(chosen_vector[1:].tolist())
    choices_indices = []
    choices_indices.append(chosen_index)
    for i in range(1, k_dimension):
        distances = []
        distances_sum = 0
        for vector in original_vec_list:
            min_dist = float("inf")
            for centroid in centroids:
                min_dist = min(min_dist, euclidean_distance(centroid, vector[1:]))
            distances.append(min_dist)
        for distance in distances:
            distances_sum += distance
        for j in range(len(distances)):
            distances[j] = distances[j] / distances_sum
        chosen_index = np.random.choice(keys, p=distances)
        chosen_vector = original_vec_list[np.where(original_vec_list[:, 0] == chosen_index)[0][0], :]
        centroids.append(chosen_vector[1:].tolist())
        choices_indices.append(chosen_index)
    for i in range(len(choices_indices) - 1):
        print(str(int(choices_indices[i])) + ",", end="")
    print(int(choices_indices[len(choices_indices) - 1]))
    final_vec_list = []
    for i in range(original_vec_list.shape[0]):
        final_vec_list_inner = []
        for j in range(original_vec_list.shape[1] - 1):
            final_vec_list_inner.append(original_vec_list[i, j + 1])
        final_vec_list.append(final_vec_list_inner)
    mykmeanssp.spk(k_dimension, original_vec_list.shape[0], final_vec_list, centroids)


if goal == "wam":
    wam = mykmeanssp.wam(vectors, num_of_vectors, vector_dimension)
    for i in range(num_of_vectors):
        for j in range(num_of_vectors):
            print('%.4f' % wam[i][j], end="")
            if j != num_of_vectors - 1:
                print(",", end="")
        print()

elif goal == "ddg":
    wam = mykmeanssp.wam(vectors, num_of_vectors, vector_dimension)
    ddg = mykmeanssp.ddg(wam, num_of_vectors)
    for i in range(num_of_vectors):
        for j in range(num_of_vectors):
            print('%.4f' % ddg[i][j], end="")
            if j != num_of_vectors - 1:
                print(",", end="")
        print()

elif goal == "gl":
    wam = mykmeanssp.wam(vectors, num_of_vectors, vector_dimension)
    ddg = mykmeanssp.ddg(wam, num_of_vectors)
    gl = mykmeanssp.gl(wam, ddg, num_of_vectors)
    for i in range(num_of_vectors):
        for j in range(num_of_vectors):
            print('%.4f' % gl[i][j], end="")
            if j != num_of_vectors - 1:
                print(",", end="")
        print()

elif goal == "jacobi":
    jacobi = mykmeanssp.jacobi(vectors, num_of_vectors)
    for i in range(num_of_vectors + 1):
        for j in range(num_of_vectors):
            print('%.4f' % jacobi[i][j], end="")
            if j != num_of_vectors - 1:
                print(",", end="")
        print()


elif goal == "spk":
    wam = mykmeanssp.wam(vectors, num_of_vectors, vector_dimension)
    ddg = mykmeanssp.ddg(wam, num_of_vectors)
    gl = mykmeanssp.gl(wam, ddg, num_of_vectors)
    jacobi = mykmeanssp.jacobi(gl, num_of_vectors)
    if len(sys.argv) == 3:
        k = mykmeanssp.eigengap_heuristic(jacobi, num_of_vectors)
    u_matrix = mykmeanssp.calculateUmatrix(jacobi, num_of_vectors, k)

    kmeans_pp(u_matrix, k)

else:
    print("An Error Has Occurred")
    exit()
