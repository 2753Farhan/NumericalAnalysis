#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define F(x)  ((x)*(x) + 1)

int main() {
  int i;
  double a = 0.0, b = 2.0, h, k;
  int n = 10;
  if(n % 2 != 0) n++;
  double x[n + 2];
  x[0] = a; x[n] = b;
  h = (b - a)/n;
  for (i = 1; i < n; i++) {
    x[i] = x[i-1] + h;
  }
  double sum_odd = 0, sum_even = 0;
  double y[n+2], sum = 0;
  for (i = 0; i <= n; i++) {
    y[i] = F(x[i]);
    if(i > 0 && i < n) {
        if(i % 2 == 0) sum_even += y[i];
        else sum_odd += y[i];
    }
  }
  double sum_x = y[0] + y[n];

  double ans = h*( sum_x + 4*sum_odd + 2*sum_even)/3;

  printf("%lf\n", ans);
  return 0;
}
