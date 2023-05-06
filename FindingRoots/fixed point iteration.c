#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAXN 200
#define F(x)  ((x)*(x)*(x) + 4*(x)*(x) - 10)
#define g(x)  sqrt(10 - (x)*(x)*(x))/2  //F(x) ke g(x) e convert kora hoise
//for example : x^3 + 4x^2 - 10 = 0; so x^2 = (10 - x^3)/4; x = sqrt((10 - x^3)/4)

int main (void) {
  int i = 1;
  double tol = 1.0e-5;
  double x0 = 1.5;
  printf("-  -  -  -  -  -  -  -  -  -  -  - - - - - - \n");
  printf("i     x0         x1         g(x0)       f(x1)\n");
  printf("-  -  -  -  -  -  -  -  -  -  -  - - -  - - -\n");

  while(i < MAXN) {
    printf("%d   %lf   ",i, x0);
    double x1 = g(x0);
    printf("%lf   %lf     %lf\n", x1, g(x0), F(x1));
    if (fabs(F(x1)) <= tol ) {
      printf("Root = %lf, Itr = %d\n", g(x0), i);
      return 1;
    }
    x0 = x1;
    i++;
  }
  printf("Itr Overflow ...\n");
  return 0;
}
