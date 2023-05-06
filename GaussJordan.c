#include <stdio.h>

#define MAX_N 10 // maximum size of the matrix

void print_matrix(double A[][MAX_N+1], int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void gaussian_elimination(double A[][MAX_N+1], int n, int m) {
    int num_pivots = (n < m) ? n : m; // number of pivots to find
    for (int i = 0; i < num_pivots; i++) {
        // find pivot row
        int max_row = i;
        for (int j = i+1; j < n; j++) {
            if (A[j][i] > A[max_row][i]) {
                max_row = j;
            }
        }
        // swap rows if necessary
        if (max_row != i) {
            for (int k = i; k < m+1; k++) {
                double temp = A[i][k];
                A[i][k] = A[max_row][k];
                A[max_row][k] = temp;
            }
        }
        // perform elimination
        double pivot = A[i][i];
        for (int j = i; j < m+1; j++) {
            A[i][j] /= pivot;
        }
        for (int j = 0; j < n; j++) {
            if (j != i) {
                double factor = A[j][i];
                for (int k = i; k < m+1; k++) {
                    A[j][k] -= factor * A[i][k];
                }
            }
        }
    }
}

int main() {
    int n, m;
    double A[MAX_N][MAX_N+1];

    printf("Enter the number of rows of the matrix (up to %d): ", MAX_N);
    scanf("%d", &n);

    printf("Enter the number of columns of the matrix (up to %d): ", MAX_N);
    scanf("%d", &m);

    printf("Enter the coefficients of the matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m+1; j++) {
            scanf("%lf", &A[i][j]);
        }
    }

    printf("Original matrix:\n");
    print_matrix(A, n, m);

    gaussian_elimination(A, n, m);

    printf("Matrix in reduced row echelon form:\n");
    print_matrix(A, n, m);

    printf("Solution:\n");
    for (int i = 0; i < n; i++) {
        printf("x%d = %f\n", i+1, A[i][m]);
    }

    return 0;
}
