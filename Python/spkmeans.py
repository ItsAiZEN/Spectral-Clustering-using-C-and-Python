import math
import sys
import pandas as pd
import numpy as np
import mykmeanssp

# TODO: either check memory of eigengap, calculateUmatrix and kmeanspp or move them to the c api.
# TODO: make the code run on NOVA
# TODO: make the code more elegant
# TODO: final tests and documentation


np.random.seed(0)

if len(sys.argv) == 3:
    goal = sys.argv[1]
    file = open(sys.argv[2], "r")
    # k denoted from eigengap
elif len(sys.argv) == 4:
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file = open(sys.argv[3], "r")
else:
    print("An Error Has Occurred")
    exit()

vectors1 = pd.read_csv(file, header=None)
vector_dimension = vectors1.shape[1]
N = vectors1.shape[0]

# transpose the matrix
vectors = vectors1.values.tolist()

if len(sys.argv) == 4:
    if k >= N or k <= 1:
        print("An Error Has Occurred")
        exit()


def euclidean_distance(vector1, vector2):
    total_sum = 0
    for i in range(len(vector1)):
        total_sum += (vector1[i] - vector2[i]) ** 2
    return math.sqrt(total_sum)


def kmeans_pp(u_matrix1, k_dimension):
    # make u_matrix1 a pandas dataframe
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
    wam = mykmeanssp.wam(vectors, N, vector_dimension)
    size = len(wam)
    # print the wam matrix, seperate by commas and output with 4 decimal places
    for i in range(size):
        for j in range(size):
            print('%.4f' % wam[i][j], end="")
            if j != N - 1:
                print(",", end="")
        print()

elif goal == "ddg":
    wam = mykmeanssp.wam(vectors, N, vector_dimension)
    ddg = mykmeanssp.ddg(wam, len(wam))
    size = len(ddg)
    # print the ddg matrix, seperate by commas and output with 4 decimal places
    for i in range(size):
        for j in range(size):
            print('%.4f' % ddg[i][j], end="")
            if j != N - 1:
                print(",", end="")
        print()

elif goal == "gl":
    wam = mykmeanssp.wam(vectors, N, vector_dimension)
    ddg = mykmeanssp.ddg(wam, len(wam))
    gl = mykmeanssp.gl(wam, ddg, len(wam))
    size = len(gl)
    # print the gl matrix, seperate by commas and output with 4 decimal places
    for i in range(size):
        for j in range(size):
            print('%.4f' % gl[i][j], end="")
            if j != N - 1:
                print(",", end="")
        print()

elif goal == "jacobi":
    jacobi = mykmeanssp.jacobi(vectors, N)
    size = len(jacobi) - 1
    # print the jacobi matrix, seperate by commas and output with 4 decimal places
    for i in range(size + 1):
        for j in range(size):
            print('%.4f' % jacobi[i][j], end="")
            if j != N - 1:
                print(",", end="")
        print()


elif goal == "spk":
    wam = mykmeanssp.wam(vectors, N, vector_dimension)
    ddg = mykmeanssp.ddg(wam, len(wam))
    gl = mykmeanssp.gl(wam, ddg, len(wam))
    jacobi = mykmeanssp.jacobi(gl, len(gl))
    if len(sys.argv) == 3:
        k = mykmeanssp.eigengap_heuristic(jacobi, len(jacobi) - 1)
    u_matrix = mykmeanssp.calculateUmatrix(jacobi, len(jacobi) - 1, k)

    kmeans_pp(u_matrix, k)

file.close()
