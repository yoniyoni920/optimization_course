#define SIZE 9
#include <math.h>
#include <stdio.h>



double numeric(double (*fp)(double), double x, double epsilon) {
    double h = epsilon;
    return (fp(x + h) - fp(x - h)) / (2.0 * h);


}
double opt_inv_quad(double (*fp)(double), double low, double high, double h, double epsilon) {
    double x0=0, x1=0, x2=5;
    double y0, y1, y2;
    double x_new=0;
    int found_bracket = 0;


    for (double x = low + h; x < high; x += h) {
        double left = fp(x - h);
        double mid = fp(x);
        double right = fp(x + h);

        if (mid < left && mid < right) {
            x0 = x - h;
            x1 = x;
            x2 = x + h;
            found_bracket = 1;
            break;
        }
    }

    if (!found_bracket) {
        printf("no minimum in the range [%lf, %lf]", low, high);
        return -1.0; // ערך שגיאה
    }


    while (fabs(x_new - x2) > epsilon) {

        y0 = numeric(fp, x0,epsilon);
        y1 = numeric(fp, x1,epsilon);
        y2 = numeric(fp, x2,epsilon);


        if (fabs(y0 - y1) < epsilon || fabs(y0 - y2) < epsilon || fabs(y1 - y2) < epsilon) {
            break;
        }


        double m0 = (y1 * y2) / ((y0 - y1) * (y0 - y2)) * x0;
        double m1 = (y0 * y2) / ((y1 - y0) * (y1 - y2)) * x1;
        double m2 = (y0 * y1) / ((y2 - y0) * (y2 - y1)) * x2;

        x_new = m0 + m1 + m2;

        x0 = x1;
        x1 = x2;
        x2 = x_new;
    }

    return x2;
}


double obj_f(double x)
    {
        double xp, sum;

        sum = x*x*x*x;
        sum = sum - 4.0*x*x;
        sum = sum - 0.2;

        return sum;
    } // obj_f


int main()
    {
        double low, high, xstar, epsilon, h;

        printf("\nEnter low, high, h, epsilon\n");
        scanf("%lf %lf %lf %lf", &low, &high, &h, &epsilon);

        xstar =  opt_inv_quad(obj_f,  low, high, h, epsilon);

        printf("\nxstar = %lf\n", xstar);

    } // main

