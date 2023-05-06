#include <stdio.h>

#define MAX_N 10 // maximum size of the matrix

void print_matrix(double A[][MAX_N*2], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n*2; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void gaussian_elimination(double A[][MAX_N*2], int n) {
    for (int i = 0; i < n; i++) {
        // find pivot row
        int max_row = i;
        for (int j = i+1; j < n; j++) {
            if (A[j][i] > A[max_row][i]) {
                max_row = j;
            }
        }
        // swap rows if necessary
        if (max_row != i) {
            for (int k = 0; k < n*2; k++) {
                double temp = A[i][k];
                A[i][k] = A[max_row][k];
                A[max_row][k] = temp;
            }
        }
        // perform elimination
        double pivot = A[i][i];
        for (int j = i; j < n*2; j++) {
            A[i][j] /= pivot;
        }
        for (int j = 0; j < n; j++) {
            if (j != i) {
                double factor = A[j][i];
                for (int k = i; k < n*2; k++) {
                    A[j][k] -= factor * A[i][k];
                }
            }
        }
    }
}

int main() {
    int n;
    double A[MAX_N][MAX_N*2];

    printf("Enter the size of the matrix (up to %d): ", MAX_N);
    scanf("%d", &n);

    printf("Enter the coefficients of the matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%lf", &A[i][j]);
        }
        for (int j = n; j < n*2; j++) {
            A[i][j] = (i == j-n) ? 1.0 : 0.0;
        }
    }

    printf("Original matrix:\n");
    print_matrix(A, n);

    gaussian_elimination(A, n);

    printf("Inverse matrix:\n");
    print_matrix(A, n);


    return 0;
}
