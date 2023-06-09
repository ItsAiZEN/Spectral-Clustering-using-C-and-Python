#include "spkmeans.h"


int main(int argc, char **argv) {
    double **returned_matrix;
    double **returned_jacobi_matrix;
    double **vector_list;
    double **wam_returned;
    double **ddg_returned;
    char *goal;
    char *file_name;
    char c;
    int vector_dimension, vector_count, i, j, scan;
    FILE *file;

    vector_dimension = 0;
    vector_count = 0;
    scan = 0;
    scan = scan + 1;

    /* checks if the number of arguments is correct */
    if (argc != 3) {
        print_error();
    }
    goal = argv[1];
    file_name = argv[2];

    /* read the file and count the number of vectors and their dimension */
    file = fopen(file_name, "r");
    if (file == NULL) {
        print_error();
    }
    while ((c = fgetc(file)) != EOF) {
        if (c == '\n') {
            vector_count++;
        }
    }
    fseek(file, 0, SEEK_SET);
    while ((c = fgetc(file)) != EOF) {
        if (c == ',') {
            vector_dimension++;
        }
        if (c == '\n') {
            break;
        }
    }
    vector_dimension++;
    vector_list = (double **) malloc(vector_count * sizeof(double *));
    if (vector_list == NULL) {
        print_error();
    }
    for (i = 0; i < vector_count; i++) {
        vector_list[i] = (double *) malloc(vector_dimension * sizeof(double));
        if (vector_list[i] == NULL) {
            print_error();
        }
    }
    /* read the file and store the vectors in a matrix */
    fseek(file, 0, SEEK_SET);
    for (i = 0; i < vector_count; i++) {
        for (j = 0; j < vector_dimension; j++) {
            c = fgetc(file);
            if (c != ',' && c != '\n') {
                /* if the character is not a comma or a new line, it is a number, therefore we need to go back one character */
                fseek(file, -1, SEEK_CUR);
            }
            scan = fscanf(file, "%lf", &vector_list[i][j]);
        }
    }
    fclose(file);

    if (strcmp(goal, "wam") == 0) {
        returned_matrix = wam(vector_list, vector_count, vector_dimension);
        for (i = 0; i < vector_count; i++) {
            for (j = 0; j < vector_count; j++) {
                if (j == vector_count - 1)
                    printf("%.4f\n", returned_matrix[i][j]);
                else
                    printf("%.4f%c", returned_matrix[i][j], ',');
            }
        }
        for (i = 0; i < vector_count; i++) {
            free(returned_matrix[i]);
        }
        free(returned_matrix);
    }
    else if (strcmp(goal, "ddg") == 0) {
        wam_returned = wam(vector_list, vector_count, vector_dimension);
        returned_matrix = ddg(wam_returned, vector_count);
        for (i = 0; i < vector_count; i++) {
            free(wam_returned[i]);
        }
        free(wam_returned);
        for (i = 0; i < vector_count; i++) {
            for (j = 0; j < vector_count; j++) {
                if (j == vector_count - 1)
                    printf("%.4f\n", returned_matrix[i][j]);
                else
                    printf("%.4f%c", returned_matrix[i][j], ',');
            }
        }
        for (i = 0; i < vector_count; i++) {
            free(returned_matrix[i]);
        }
        free(returned_matrix);
    }
    else if (strcmp(goal, "gl") == 0) {
        wam_returned = wam(vector_list, vector_count, vector_dimension);
        ddg_returned = ddg(wam_returned, vector_count);
        returned_matrix = gl(ddg_returned, wam_returned, vector_count);
        for (i = 0; i < vector_count; i++) {
            free(wam_returned[i]);
        }
        free(wam_returned);
        for (i = 0; i < vector_count; i++) {
            free(ddg_returned[i]);
        }
        free(ddg_returned);
        for (i = 0; i < vector_count; i++) {
            for (j = 0; j < vector_count; j++) {
                if (j == vector_count - 1)
                    printf("%.4f\n", returned_matrix[i][j]);
                else
                    printf("%.4f%c", returned_matrix[i][j], ',');
            }
        }
        for (i = 0; i < vector_count; i++) {
            free(returned_matrix[i]);
        }
        free(returned_matrix);
    }
    else if (strcmp(goal, "jacobi") == 0) {
        returned_jacobi_matrix = jacobi(vector_list, vector_count);
        for (i = 0; i < vector_count + 1; i++) {
            for (j = 0; j < vector_count; j++) {
                if (j == vector_count - 1)
                    printf("%.4f\n", returned_jacobi_matrix[i][j]);
                else
                    printf("%.4f%c", returned_jacobi_matrix[i][j], ',');
            }
        }
        for (i = 0; i < vector_count + 1; i++) {
            free(returned_jacobi_matrix[i]);
        }
        free(returned_jacobi_matrix);
    }
    else {
        print_error();
    }


    /*free memory*/
    for (i = 0; i < vector_count; i++) {
        free(vector_list[i]);
    }
    free(vector_list);

    return 0;
}

void print_error() {
    printf("An error has occurred\n");
    exit(1);
}

/* returns the euclidean distance between two vectors */
double euclidean_distance(double *vector1, double *vector2, int max_d) {
    double total_sum;
    int i;
    total_sum = 0;
    for (i = 0; i < max_d; i++) {
        total_sum += pow(vector1[i] - vector2[i], 2);
    }
    return sqrt(total_sum);
}

/* wam = weighted adjacency matrix */
double **wam(double **vectors, int num_of_vectors, int vector_dimension) {
    double **wam_matrix;
    int i, j;
    wam_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (wam_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        wam_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (wam_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        for (j = 0; j < num_of_vectors; j++) {
            if (i == j) {
                wam_matrix[i][j] = 0;
            } else {
                wam_matrix[i][j] = exp(-pow(euclidean_distance(vectors[i], vectors[j], vector_dimension), 2) / 2);
            }
        }
    }
    return wam_matrix;
}

/*ddg = degree diagonal matrix (diagonal matrix with the sum of each row as the diagonal elements)*/
double **ddg(double **wam_matrix, int num_of_vectors) {
    double **ddg_matrix;
    double sum;
    int i, j;
    ddg_matrix = (double **) calloc(num_of_vectors, sizeof(double *));
    if (ddg_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        ddg_matrix[i] = (double *) calloc(num_of_vectors, sizeof(double));
        if (ddg_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        sum = 0;
        for (j = 0; j < num_of_vectors; j++) {
            sum += wam_matrix[i][j];
        }
        ddg_matrix[i][i] = sum;
    }
    return ddg_matrix;
}

/* gl = graph laplacian (ddg - wam) */
double **gl(double **ddg_matrix, double **wam_matrix, int num_of_vectors) {
    double **gl_matrix;
    int i, j;
    gl_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (gl_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        gl_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (gl_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        for (j = 0; j < num_of_vectors; j++) {
            gl_matrix[i][j] = ddg_matrix[i][j] - wam_matrix[i][j];
        }
    }
    return gl_matrix;
}

/*finds Aij for jacobi (off diagonal element with the largest absolute value)*/
int *find_largest_absolute_value_coordinates(double **matrix, int num_of_vectors) {
    double max;
    int *coordinates;
    int i, j;
    max = 0;
    coordinates = (int *) malloc(2 * sizeof(int));
    if (coordinates == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        for (j = 0; j < num_of_vectors; j++) {
            if (i >= j) {
                continue;
            }
            if (fabs(matrix[i][j]) > max) {
                max = fabs(matrix[i][j]);
                coordinates[0] = i;
                coordinates[1] = j;
            }
        }
    }
    return coordinates;
}

/*returns the sign of a number (positive or negative, for 0 returns 1)*/
int sign(double x) {
    if (x >= 0) {
        return 1;
    }
    return -1;
}

/*calculates P matrix*/
double **calculate_rotation_matrix(double **mat, int num_of_vectors, int *coordinates) {
    double **rotation_matrix;
    double aii, ajj, aij, t, c, s, theta;
    int i, j;
    rotation_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (rotation_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        rotation_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (rotation_matrix[i] == NULL) {
            print_error();
        }
    }
    aii = mat[coordinates[0]][coordinates[0]];
    ajj = mat[coordinates[1]][coordinates[1]];
    aij = mat[coordinates[0]][coordinates[1]];
    theta = (ajj - aii) / (2 * aij);
    t = sign(theta) / (fabs(theta) + sqrt(1 + pow(theta, 2)));
    c = 1 / sqrt(1 + pow(t, 2));
    s = c * t;
    for (i = 0; i < num_of_vectors; i++) {
        for (j = 0; j < num_of_vectors; j++) {
            if (i == j) {
                rotation_matrix[i][j] = 1;
            } else {
                rotation_matrix[i][j] = 0;
            }
        }
    }
    rotation_matrix[coordinates[0]][coordinates[0]] = c;
    rotation_matrix[coordinates[1]][coordinates[1]] = c;
    rotation_matrix[coordinates[0]][coordinates[1]] = s;
    rotation_matrix[coordinates[1]][coordinates[0]] = -s;

    return rotation_matrix;
}

/* checks if the difference between the sum of the squares of the off diagonal elements of the current matrix
     and the previous matrix is smaller than eps*/
int check_convergence(double **matrix, double **previous_matrix, int num_of_vectors, double eps) {
    double sum, previous_sum;
    int i, j;
    previous_sum = 0;
    sum = 0;
    for (i = 0; i < num_of_vectors; i++) {
        for (j = 0; j < num_of_vectors; j++) {
            if (i == j) {
                continue;
            }
            previous_sum += pow(previous_matrix[i][j], 2);
            sum += pow(matrix[i][j], 2);
        }
    }
    if (previous_sum - sum <= eps) {
        return 1;
    }
    return 0;
}

/* multiplies two square matrices */
double **multiply_matrices(double **matrix1, double **matrix2, int num_of_vectors) {
    double **result_matrix;
    int i, j, k;

    result_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (result_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        result_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (result_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        for (j = 0; j < num_of_vectors; j++) {
            result_matrix[i][j] = 0;
            for (k = 0; k < num_of_vectors; k++) {
                result_matrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return result_matrix;
}

/* calculates the eigenvalues and eigenvectors of a matrix using the jacobi method*/
double **jacobi(double **original_matrix, int num_of_vectors) {
    double **matrix; /* contains eigenvalues on diagonal at the end of the iteration (A')*/
    double **jacobi_matrix; /* contains eigenvalues in first row and eigenvectors as columns beneath (returned matrix)*/
    double **rot_mat; /* P for each iteration*/
    double **final_matrix; /* the matrix that contains the eigenvectors as columns at the end of the iteration (multiplication of all P matrices)*/
    double **previous_matrix; /* "matrix" from the previous iteration (A)*/
    double **mul_matrices;
    double eps, s, c;
    int *coordinates;
    int num_of_iterations, i, j, k;
    matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        for (j = 0; j < num_of_vectors; j++) {
            matrix[i][j] = original_matrix[i][j];
        }
    }
    previous_matrix = (double **) malloc(num_of_vectors * sizeof(double *));
    if (previous_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        previous_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (previous_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        for (j = 0; j < num_of_vectors; j++) {
            previous_matrix[i][j] = matrix[i][j];
        }
    }

    eps = 0.00001;
    num_of_iterations = 0;
    final_matrix = (double **) calloc(num_of_vectors, sizeof(double *));
    if (final_matrix == NULL) {
        print_error();
    }
    for (k = 0; k < num_of_vectors; k++) {
        final_matrix[k] = (double *) calloc(num_of_vectors, sizeof(double));
        if (final_matrix[k] == NULL) {
            print_error();
        }
    }
    for (k = 0; k < num_of_vectors; k++) {
        final_matrix[k][k] = 1;
    }
    while (num_of_iterations < 100) { /* 100 is set as the maximum number of iterations*/
        num_of_iterations++;
        coordinates = find_largest_absolute_value_coordinates(matrix, num_of_vectors);
        i = coordinates[0];
        j = coordinates[1];
        rot_mat = calculate_rotation_matrix(matrix, num_of_vectors, coordinates);
        s = rot_mat[i][j];
        c = rot_mat[i][i];
        for (k = 0; k < num_of_vectors; k++) {
            /* calculates the new matrix (A') using P and A instead of multiplying using P and A */
            matrix[i][k] = c * previous_matrix[i][k] - s * previous_matrix[j][k];
            matrix[j][k] = s * previous_matrix[i][k] + c * previous_matrix[j][k];
            matrix[k][i] = matrix[i][k];
            matrix[k][j] = matrix[j][k];
        }
        matrix[i][i] = pow(c, 2) * previous_matrix[i][i] - 2 * c * s * previous_matrix[i][j] +
                       pow(s, 2) * previous_matrix[j][j];
        matrix[j][j] = pow(s, 2) * previous_matrix[i][i] + 2 * c * s * previous_matrix[i][j] +
                       pow(c, 2) * previous_matrix[j][j];
        matrix[i][j] = (pow(c, 2) - pow(s, 2)) * previous_matrix[i][j] +
                       c * s * (previous_matrix[i][i] - previous_matrix[j][j]);
        matrix[j][i] = matrix[i][j];
        /* add current P to final matrix calculation*/
        mul_matrices = multiply_matrices(final_matrix, rot_mat, num_of_vectors);
        for (i = 0; i < num_of_vectors; i++) {
            for (j = 0; j < num_of_vectors; j++) {
                final_matrix[i][j] = mul_matrices[i][j];
            }
        }
        for (i = 0; i < num_of_vectors; i++) {
            free(mul_matrices[i]);
        }
        free(mul_matrices);
        for (i = 0; i < num_of_vectors; i++) {
            free(rot_mat[i]);
        }
        free(rot_mat);
        free(coordinates);
        if (check_convergence(matrix, previous_matrix, num_of_vectors, eps)) { /* check convergence*/
            break;
        }
        /* need to reset prev matrix for next iteration*/
        for (i = 0; i < num_of_vectors; i++) {
            for (j = 0; j < num_of_vectors; j++) {
                previous_matrix[i][j] = matrix[i][j];
            }
        }
    }
    jacobi_matrix = (double **) malloc((num_of_vectors + 1) * sizeof(double *));
    if (jacobi_matrix == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors + 1; i++) {
        jacobi_matrix[i] = (double *) malloc(num_of_vectors * sizeof(double));
        if (jacobi_matrix[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) { /* set the eigenvalues as the first row */
        jacobi_matrix[0][i] = matrix[i][i];
    }
    for (i = 1; i < num_of_vectors + 1; i++) {
        for (j = 0; j < num_of_vectors; j++) {
            jacobi_matrix[i][j] = final_matrix[i - 1][j];
        }
    }
    /*free memory*/
    for (i = 0; i < num_of_vectors; i++) {
        free(matrix[i]);
    }
    free(matrix);
    for (i = 0; i < num_of_vectors; i++) {
        free(final_matrix[i]);
    }
    free(final_matrix);
    for (i = 0; i < num_of_vectors; i++) {
        free(previous_matrix[i]);
    }
    free(previous_matrix);
    /* check if any of the eigenvalues between -0.0000 to -0.0001, if so multiply by -1 */
    for (i = 0; i < num_of_vectors; i++) {
        if (jacobi_matrix[0][i] < 0 && jacobi_matrix[0][i] > -0.0001) {
            jacobi_matrix[0][i] = 0;
            for (j = 1; j < num_of_vectors + 1; j++) {
                jacobi_matrix[j][i] = jacobi_matrix[j][i] * -1;
            }
        }
    }

    return jacobi_matrix;
}

/* implements a comparator function for qsort*/
int sortFunc(const void *a, const void *b) {
    double val;
    val = (*(double *) a - *(double *) b);
    if (val > 0) {
        return 1;
    }
    if (val < 0) {
        return -1;
    }
    return 0;
}

/* calculates k using the eigengap heuristic (if k is not provided as an argument) */
int eigengap_heuristic(double **jacobi_matrix, int num_of_vectors) {
    double *eigenvalues;
    double k;
    int i, index;
    index = 0;
    k = 0;
    eigenvalues = (double *) malloc(num_of_vectors * sizeof(double));
    if (eigenvalues == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        eigenvalues[i] = jacobi_matrix[0][i];
    }
    qsort(eigenvalues, num_of_vectors, sizeof(double), sortFunc);
    for (i = 0; i < num_of_vectors / 2; i++) {
        if (fabs(eigenvalues[i + 1] - eigenvalues[i]) > k) {
            k = fabs(eigenvalues[i + 1] - eigenvalues[i]);
            index = i;
        }
    }
    /* free memory*/
    free(eigenvalues);
    return index + 1;
}

/* calculates the U matrix, returns a matrix with k eigenvectors (with smallest eigenvalues) as its columns */
double **calculateUmatrix(double **jacobi_matrix, int num_of_vectors, int k) {
    double **U;
    double *eigenvalues;
    int i, j, l;

    eigenvalues = (double *) malloc(num_of_vectors * sizeof(double));
    if (eigenvalues == NULL) {
        print_error();
    }
    U = (double **) malloc(num_of_vectors * sizeof(double *));
    if (U == NULL) {
        print_error();
    }
    for (i = 0; i < num_of_vectors; i++) {
        U[i] = (double *) malloc(k * sizeof(double));
        if (U[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_vectors; i++) {
        eigenvalues[i] = jacobi_matrix[0][i];
    }
    qsort(eigenvalues, num_of_vectors, sizeof(double), sortFunc);
    for (i = 0; i < k; i++) {
        for (j = 0; j < num_of_vectors; j++) {
            if (eigenvalues[i] == jacobi_matrix[0][j]) {
                jacobi_matrix[0][j] = 2147483647; /* used to avoid duplicates of the same eigenvector for eigenvalues with multiple occurrences*/
                for (l = 0; l < num_of_vectors; l++) {
                    U[l][i] = jacobi_matrix[l + 1][j];
                }
            }
        }
    }
    /* free memory*/
    free(eigenvalues);
    return U;
}

/* implementation of kmeans++ */
void kmeanspp(int num_of_clusters, int num_of_iterations, int vector_dimension, int count,
              double **vector_list, double eps, double **init_centroids) {

    double ***clusters;
    double **centroids;
    double *temp_centroid;
    int *cluster_sizes_copy;
    int *cluster_sizes;
    double max_distance, min_distance1, min_distance2, distance1, distance2, sum;
    int min_index1, min_index2, i, j, m, n, l, k, p, q;

    max_distance = -1;
    centroids = calloc(num_of_clusters, sizeof(double *));
    if (centroids == NULL) {
        print_error();
    }
    cluster_sizes_copy = calloc(num_of_clusters, sizeof(int));
    if (cluster_sizes_copy == NULL) {
        print_error();
    }
    temp_centroid = calloc(vector_dimension, sizeof(double));
    if (temp_centroid == NULL) {
        print_error();
    }
    cluster_sizes = calloc(num_of_clusters, sizeof(int));
    if (cluster_sizes == NULL) {
        print_error();
    }
    clusters = calloc(num_of_clusters, sizeof(double **));
    if (clusters == NULL) {
        print_error();
    }

    for (i = 0; i < num_of_clusters; i++) {
        centroids[i] = calloc(vector_dimension, sizeof(double));
        if (centroids[i] == NULL) {
            print_error();
        }
    }
    for (i = 0; i < num_of_clusters; i++) {
        for (j = 0; j < vector_dimension; j++) {
            centroids[i][j] = init_centroids[i][j];
        }
    }
    for (i = 0; i < num_of_iterations; i++) {

        for (j = 0; j < count; j++) {
            min_distance1 = euclidean_distance(vector_list[j], centroids[0], vector_dimension);
            min_index1 = 0;
            for (k = 1; k < num_of_clusters; k++) {
                distance1 = euclidean_distance(vector_list[j], centroids[k], vector_dimension);
                if (distance1 < min_distance1) {
                    min_distance1 = distance1;
                    min_index1 = k;
                }
            }
            cluster_sizes[min_index1] += 1;
        }
        for (j = 0; j < num_of_clusters; j++) {
            cluster_sizes_copy[j] = cluster_sizes[j];
        }
        for (j = 0; j < num_of_clusters; j++) {
            clusters[j] = calloc(cluster_sizes[j], sizeof(double *));
            if (clusters[j] == NULL) {
                print_error();
            }
            for (m = 0; m < cluster_sizes[j]; m++) {
                clusters[j][m] = calloc(vector_dimension, sizeof(double));
                if (clusters[j][m] == NULL) {
                    print_error();
                }
            }
        }
        for (j = 0; j < count; j++) {
            min_distance2 = euclidean_distance(vector_list[j], centroids[0], vector_dimension);
            min_index2 = 0;
            for (p = 1; p < num_of_clusters; p++) {
                distance2 = euclidean_distance(vector_list[j], centroids[p], vector_dimension);
                if (distance2 < min_distance2) {
                    min_distance2 = distance2;
                    min_index2 = p;
                }
            }
            for (q = 0; q < vector_dimension; q++) {
                clusters[min_index2][cluster_sizes[min_index2] - 1][q] = vector_list[j][q];
            }
            cluster_sizes[min_index2] -= 1;
        }
        max_distance = -1;
        for (j = 0; j < num_of_clusters; j++) {
            for (l = 0; l < vector_dimension; l++) {
                temp_centroid[l] = centroids[j][l];
            }
            for (m = 0; m < vector_dimension; m++) {
                sum = 0;
                for (n = 0; n < cluster_sizes_copy[j]; n++) {
                    sum += clusters[j][n][m];
                }
                centroids[j][m] = sum / cluster_sizes_copy[j];

            }
            if (euclidean_distance(temp_centroid, centroids[j], vector_dimension) > max_distance) {
                max_distance = euclidean_distance(temp_centroid, centroids[j], vector_dimension);
            }
        }
        for (i = 0; i < num_of_clusters; i++) {
            for (j = 0; j < cluster_sizes_copy[i]; j++) {
                free(clusters[i][j]);
            }
            free(clusters[i]);
        }
        if (max_distance <= eps) {
            break;
        }
    }
    for (i = 0; i < num_of_clusters; i++) {
        for (j = 0; j < vector_dimension; j++) {
            if (j == vector_dimension - 1)
                printf("%.4f\n", centroids[i][j]);
            else
                printf("%.4f%c", centroids[i][j], ',');
        }
    }

    /* free memory */
    for (i = 0; i < num_of_clusters; i++) {
        free(centroids[i]);
    }
    free(centroids);
    free(temp_centroid);
    free(cluster_sizes);
    free(clusters);
    free(cluster_sizes_copy);
}


