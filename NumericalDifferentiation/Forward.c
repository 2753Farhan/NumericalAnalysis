#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define F(x)  ((x)*(x) + 1)

int main() {
  int i;
  double a = 0.0, b = 2.0, h, k;
  int n = 10;
  double x[n + 2];
  printf("i    x[i]      f(x[i])     f'(x[i])\n\n");
  x[0] = a; x[n] = b;
  for (i = 1; i < n; i++) {
    h = (b - a)/n;
    x[i] = x[i-1] + h;
  }
  for (i = 0; i <= n; i++) {
    k = (F(x[i + 1]) - F(x[i]))/h;
    if(i != n) printf("%d   %lf   %lf   %lf\n", i, x[i], F(x[i]), k);
    else {
        printf("%d  %lf   %lf   ", i, x[i - 1] + h, F(x[i]));
        printf("--------");
    }
 }
 return 0;
}
