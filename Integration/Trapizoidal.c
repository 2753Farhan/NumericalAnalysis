#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define F(x)  ((x)*(x) + 1)

int main() {
  int i;
  double a = 0.0, b = 2.0, h, k;
  int n = 10;
  double x[n + 2];
  x[0] = a; x[n] = b;
  h = (b - a)/n;
  for (i = 1; i < n; i++) {
    x[i] = x[i-1] + h;
  }
  double y[n+2], sum = 0;
  for (i = 0; i <= n; i++) {
    y[i] = F(x[i]);
    sum += y[i];
  }
  double sum_x = y[0] + y[n];
  double sum_i = sum - sum_x;

  double ans = h*( sum_x + 2*sum_i )/2;

  printf("%lf\n", ans);
  return 0;
}
