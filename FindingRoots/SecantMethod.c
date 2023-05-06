//root NR.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAXN 200
#define F(x)  ( (x)*(x)*(x) + 4*(x)*(x) - 10 )

int main (void) {
  int i;
  double x1=1.5, x2=2, x3, tol = 10.0e-6;
  double fx1,fx2,fx3;

  fx1 = F(x1);   fx2 = F(x2);
printf("i       x0                x1                x2              F(x0)            F(x1)           F(x2)\n");
  for (i=1; i< MAXN; i++) {
    x3 = x2 - fx2*(x2 - x1)/(fx2 - fx1);
    fx3 = F(x3);
    printf("%d      %lf         %lf         %lf           %lf         %lf       %lf\n",i,x1,x2,x3,fx1,fx2,fx3);

    if (fabs(fx3) <= tol ) {
      printf("Root = %lf, Itr = %d\n", x3,i );
      return EXIT_SUCCESS;
    }

    x1 = x2;
    x2 = x3;
    fx1 = fx2;
    fx2 = fx3;
  }

  printf("Itr Overflow ...\n");
  return EXIT_FAILURE;
}
