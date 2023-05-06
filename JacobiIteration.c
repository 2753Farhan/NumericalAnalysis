#include <stdio.h>
#include <math.h>

#define MAX_N 10 // maximum size of the matrix
#define MAX_ITER 100 // maximum number of iterations
#define TOL 1e-6 // tolerance for convergence

void print_matrix(double A[][MAX_N+1], int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m+1; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void jacobi_iteration(double A[][MAX_N+1], int n, double x[]) {
    double x_new[MAX_N];
    int iter = 0;
    double diff = TOL+1;
    while (iter < MAX_ITER && diff > TOL) {
        // compute new solution
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += A[i][j] * x[j];
                }
            }
            x_new[i] = (A[i][n] - sum) / A[i][i];
        }
        // compute difference between old and new solution
        diff = 0.0;
        for (int i = 0; i < n; i++) {
            diff += fabs(x_new[i] - x[i]);
            x[i] = x_new[i];
        }
        iter++;
    }
    printf("Converged in %d iterations with error %f\n", iter, diff);
}

int main() {
    int n;
    double A[MAX_N][MAX_N+1];
    double x[MAX_N];

    printf("Enter the number of equations (up to %d): ", MAX_N);
    scanf("%d", &n);

    printf("Enter the coefficients of the matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n+1; j++) {
            scanf("%lf", &A[i][j]);
        }
    }

    printf("Enter the initial solution:\n");
    for (int i = 0; i < n; i++) {
        scanf("%lf", &x[i]);
    }

    printf("Original matrix:\n");
    print_matrix(A, n, n);

    jacobi_iteration(A, n, x);

    printf("Solution:\n");
    for (int i = 0; i < n; i++) {
        printf("x%d = %f\n", i+1, x[i]);
    }

    return 0;
}
