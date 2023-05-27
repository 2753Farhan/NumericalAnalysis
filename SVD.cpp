#include <bits/stdc++.h>
using namespace std;
#define M 10

void jacobi_eigenvalue(int N, double A[M][M], double V[M][M], double d[])
{
    int i, j, k, nrot;
    double tresh, theta, tau, t, sm, s, c, p, g, h, tmp;
    double b[N], z[N];

    // Initialize V to the identity matrix
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            V[i][j] = 0.0;
        }
        V[i][i] = 1.0;
    }

    // Initialize b and d to the diagonal of A
    for (i = 0; i < N; i++)
    {
        b[i] = d[i] = A[i][i];
        z[i] = 0.0;
    }

    nrot = 0;
    for (i = 0; i < 50; i++)
    {
        sm = 0.0;
        for (j = 0; j < N - 1; j++)
        {
            for (k = j + 1; k < N; k++)
            {
                sm += fabs(A[j][k]);
            }
        }
        if (sm == 0.0)
        {
            // A is diagonal, so we're done
            return;
        }
        if (i < 3)
        {
            tresh = 0.2 * sm / (N * N);
        }
        else
        {
            tresh = 0.0;
        }
        for (j = 0; j < N - 1; j++)
        {
            for (k = j + 1; k < N; k++)
            {
                g = 100.0 * fabs(A[j][k]);
                if (i > 3 && (fabs(d[j]) + g == fabs(d[j]))
                        && (fabs(d[k]) + g == fabs(d[k])))
                {
                    A[j][k] = 0.0;
                }
                else if (fabs(A[j][k]) > tresh)
                {
                    h = d[k] - d[j];
                    if (fabs(h) + g == fabs(h))
                    {
                        t = A[j][k] / h;
                    }
                    else
                    {
                        theta = 0.5 * h / A[j][k];
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta));
                        if (theta < 0.0)
                        {
                            t = -t;
                        }
                    }
                    c = 1.0 / sqrt(1 + t*t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * A[j][k];
                    z[j] -= h;
                    z[k] += h;
                    d[j] -= h;
                    d[k] += h;
                    A[j][k] = 0.0;
                    for (int l = 0; l < j; l++)
                    {
                        tmp = A[l][j];
                        A[l][j] = c*tmp - s*A[l][k];
                        A[l][k] = s*tmp + c*A[l][k];
                    }
                    for (int l = j + 1; l < k; l++)
                    {
                        tmp = A[j][l];
                        A[j][l                    ] = c*tmp - s*A[l][k];
                        A[l][k] = s*tmp + c*A[l][k];
                    }
                    for (int l = k + 1; l < N; l++)
                    {
                        tmp = A[j][l];
                        A[j][l] = c*tmp - s*A[k][l];
                        A[k][l] = s*tmp + c*A[k][l];
                    }
                    for (int l = 0; l < N; l++)
                    {
                        tmp = V[l][j];
                        V[l][j] = c*tmp - s*V[l][k];
                        V[l][k] = s*tmp + c*V[l][k];
                    }
                    nrot++;
                }
            }
        }
        for (j = 0; j < N; j++)
        {
            b[j] += z[j];
            d[j] = b[j];
            z[j] = 0.0;
        }
    }

    printf("Too many iterations in Jacobi eigenvalue algorithm.");
}

void sort_eigenvectors(double eigenvalues[], double eigenvectors[][M], int n)
{
    // Create an array of pairs, where each pair contains the eigenvalue and its corresponding eigenvector
    pair<double, vector<double>> pairs[n];
    for (int i = 0; i < n; i++)
    {
        vector<double> temp;
        for(int j=0; j<n; j++)
        {
            temp.push_back(eigenvectors[j][i]);
        }
        pairs[i] = make_pair(eigenvalues[i], temp);
    }

    // Sort the array of pairs in descending order of the eigenvalues
    sort(pairs, pairs + n, greater<pair<double, vector<double>>>());

    // Create a copy of the eigenvector matrix
    double sorted_eigenvectors[n][M];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            sorted_eigenvectors[j][i] = pairs[i].second[j];
        }
    }

    // Copy the sorted eigenvectors into the original eigenvector matrix
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            eigenvectors[j][i] = sorted_eigenvectors[j][i];
        }
    }
    sort(eigenvalues,eigenvalues+n,greater<double>());
}



int main()
{
    int n,m, i, j, p, q, flag;
    double mat[M][M], A[M][M], B[M][M];
    double temp[M][M], theta, zero=1e-5, Max, pi=3.141592654;

    cout << "Enter the number of rows and columns of the matrix: ";
    cin >> n >> m;

    double matrix[n][m];
    cout << "Enter the elements of the matrix:" << endl;
    for(int i = 0; i <n; i++)
    {
        for(int j = 0; j <m; j++)
        {
            cin >> matrix[i][j];
        }
    }

    double transpose[m][n];
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j <n; j++)
        {
            transpose[i][j] = matrix[j][i];
        }
    }


    cout << "The transpose of the matrix is:" << endl;
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j <n; j++)
        {
            cout << transpose[i][j] << " ";
        }
        cout << endl;
    }

    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < m; j++)
        {
            A[i][j] = 0;
            for(int k = 0; k <n; k++)
            {
                A[i][j] += matrix[k][j] * transpose[i][k];
            }
        }
    }


    cout << "The product of the transpose and the matrix is:" << endl;
    for(int i = 0; i <m; i++)
    {
        for(int j = 0; j <m; j++)
        {
            cout <<  A[i][j] << " ";
        }
        cout << endl;
    }

    double V[M][M], d[M],Vt[M][M];
    jacobi_eigenvalue(m,A, V, d);


    sort_eigenvectors(d, V, m);

    printf("Eigenvalues:\n");
    for (int i = 0; i < m; i++)
    {
        printf("%f\n", d[i]);
    }

    printf("So V is:\n");
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            printf("%f ", V[i][j]);
        }
        printf("\n");
    }
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j <m; j++)
        {
            Vt[i][j] = V[j][i];
        }
    }
    printf("\n\nMatrix Vt :\n............................\n");
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            printf("%f ", Vt[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    double w[M][M];
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
            if(i==j)
            {
                w[i][j]=sqrt(d[i]);
            }
            else
            {
                w[i][j]=0;
            }
        }
    }

    printf("\n\nMatrix w :\n............................\n");

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
            cout << w[i][j] << " ";
        }
        cout << "\n";
    }
    printf("\n");


    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            B[i][j] = 0;
            for(int k = 0; k <m; k++)
            {
                B[i][j] += matrix[i][k] * transpose[k][j];
            }
        }
    }
    cout << "The product of the matrix and its transpose is:" << endl;
    for(int i = 0; i <n; i++)
    {
        for(int j = 0; j <n; j++)
        {
            cout <<  B[i][j] << " ";
        }
        cout << endl;
    }


    double U[M][M],d1[M];
    jacobi_eigenvalue(n,B, U, d1);

    sort_eigenvectors(d1, U,n );

    printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%f\n", d1[i]);
    }

    printf("\n\nMatrix U :\n............................\n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%f ", U[i][j]);
        }
        printf("\n");
    }

    return 0;
}
