#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double euclidean_distance(double *vector1, double *vector2, int max_d);

double **wam(double **vectors, int num_of_vectors, int vector_dimension);

double **ddg(double **wam_matrix, int num_of_vectors);

double **gl(double **ddg_matrix, double **wam_matrix, int num_of_vectors);

int *find_largest_absolute_value_coordinates(double **matrix, int num_of_vectors);

double **calculate_rotation_matrix(double **mat, int num_of_vectors, int *coordinates);

int check_convergence(double **matrix, double **previous_matrix, int num_of_vectors, double eps);

double **jacobi(double **matrix, int num_of_vectors);

double **multiply_matrices(double **matrix1, double **matrix2, int num_of_vectors);

int sign(double x);

int sortFunc(const void *a, const void *b);

int eigengap_heuristic(double **jacobi_matrix, int num_of_vectors);

double **calculateUmatrix(double **jacobi_matrix, int num_of_vectors, int k);

void kmeanspp(int num_of_clusters, int num_of_iterations, int vector_dimension, int count,
              double vector_list[][vector_dimension], double eps, double init_centroids[][vector_dimension]);

void print_error();