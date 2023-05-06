
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAXN 200
#define F(x)  ( (x)*(x)*(x) + 4*(x)*(x) - 10 )

int main (void)
{
    int i;
    double a=1.25, b=1.5, tol = 1.0e-5, c;
    double fa,fb,fc;

    fa = F(a);
    fb = F(b);

    if (fa*fb >= 0 )
    {
        printf("No Root ...\n");
        return EXIT_FAILURE;
    }
    printf("i       a              b            c           F(a)          F(b)         F(c)\n");
    for (i=1; i< MAXN; i++ )
    {

        c= (a*fb - b*fa)/(fb - fa);
        fc = F(c);
        printf("%d      %lf     %lf     %lf     %lf     %lf     %lf\n",i,a,b,c,F(a),F(b),F(c));

        if (fabs(fc) <= tol )
        {
            printf("\nAbsolute Root = %lf, Itr = %d\n", c,i );
            //return EXIT_SUCCESS;
            break;
        }

        if (fa*fc < 0 )
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    //printf("Itr Overflow ...\n");
    //return EXIT_FAILURE;
}


