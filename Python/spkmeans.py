import math
import sys
import pandas as pd
import numpy as np
import mykmeanssp

np.random.seed(0)

if len(sys.argv) == 6:
    file1 = open(sys.argv[4], "r")
    file2 = open(sys.argv[5], "r")
    iter = int(sys.argv[2])
    epsilon = float(sys.argv[3])

elif len(sys.argv) == 5:
    file1 = open(sys.argv[3], "r")
    file2 = open(sys.argv[4], "r")
    iter = 300
    epsilon = float(sys.argv[2])
else:
    print("An Error Has Occurred")
    exit()

k = int(sys.argv[1])

f1 = pd.read_csv(file1, header=None)
f2 = pd.read_csv(file2, header=None)
vectors = pd.merge(f1, f2, on=0)
vector_dimension = f1.shape[1]
vectors.sort_values(by=[vectors.columns[0]], inplace=True)

N = vectors.shape[0]

if k >= N or k <= 1:
    print("Invalid number of clusters!")
    exit()

if iter <= 1 or iter >= 1000:
    print("Invalid maximum iteration!")
    exit()


def euclidean_distance(vector1, vector2):
    total_sum = 0
    for i in range(len(vector1)):
        total_sum += (vector1[i] - vector2[i]) ** 2
    return math.sqrt(total_sum)


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
