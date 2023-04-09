## This project implements the spectral clustering algorithm.
### Most of the code is written in C in order to achieve better performance, The C code is imported into Python using a C API module
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
