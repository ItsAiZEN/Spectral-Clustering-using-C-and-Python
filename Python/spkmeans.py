import math
import sys
import pandas as pd
import numpy as np
import mykmeanssp

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

vectors = pd.read_csv(file, header=None)
vector_dimension = vectors.shape[1]
N = vectors.shape[0]

if len(sys.argv) == 4:
    if k >= N or k <= 1:
        print("Invalid number of clusters!")
        exit()


def euclidean_distance(vector1, vector2):
    total_sum = 0
    for i in range(len(vector1)):
        total_sum += (vector1[i] - vector2[i]) ** 2
    return math.sqrt(total_sum)


if goal == "wam":
    wam = mykmeanssp.wam(vectors, N, vector_dimension)
    size = len(wam)
    # print the wam matrix, seperate by commas and output with 4 decimal places
    for i in range(size):
        for j in range(size):
            print('%.4f' % wam[i][j], end="")
            if j != vector_dimension - 1:
                print(",", end="")
        print()

elif goal == "ddg":
    ddg = mykmeanssp.ddg(mykmeanssp(vectors, N, vector_dimension), len(vectors))
    size = len(ddg)
    # print the ddg matrix, seperate by commas and output with 4 decimal places
    for i in range(size):
        for j in range(size):
            print('%.4f' % ddg[i][j], end="")
            if j != vector_dimension - 1:
                print(",", end="")
        print()

elif goal == "gl":
    wam = mykmeanssp.wam(vectors, N, vector_dimension)
    ddg = mykmeanssp.ddg(mykmeanssp(vectors, N, vector_dimension), len(vectors))
    gl = mykmeanssp.gl(wam, ddg, len(wam))
    size = len(gl)
    # print the gl matrix, seperate by commas and output with 4 decimal places
    for i in range(size):
        for j in range(size):
            print('%.4f' % gl[i][j], end="")
            if j != vector_dimension - 1:
                print(",", end="")
        print()

elif goal == "jacobi":
    # NEXT TIME CONTINUE HERE


def kmeans_pp(vectors, k, iter, epsilon):
    original_vec_list = np.array(vectors)
    keys = [vector for vector in vectors[0]]
    centroids = []
    chosen_index = np.random.choice(keys)
    chosen_vector = original_vec_list[np.where(original_vec_list[:, 0] == chosen_index)[0][0], :]
    centroids.append(chosen_vector[1:].tolist())
    choices_indices = []
    choices_indices.append(chosen_index)
    for i in range(1, k):
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
    mykmeanssp.fit(k, iter, original_vec_list.shape[1] - 1, original_vec_list.shape[0], final_vec_list, epsilon,
                   centroids)


kmeans_pp(vectors, k, iter, epsilon)
