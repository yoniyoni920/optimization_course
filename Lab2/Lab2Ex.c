//
// Created by Yoni on 18/11/2025.
//
#define SIZE 9
#include <math.h>
#include <stdio.h>
double power_of_x(double x, int power) ;
double f(double x);
/*
coeffs - are the parameters of the polynom
x      - the power of the polimon
n      - length of the array
*/
double compute_polynom(double coeffs[], double x, int n) {
    double solution = 0;
    for (int i = 0; i <= n; i++){
        solution = solution + coeffs[n-i-1] * power_of_x(x,i);
    }
    return solution;
}
double power_of_x(double x, int power) {
    if (power == 0)
        return 1;
    double powered_x = x;
    for (int i = 1; i < power; i++)
        powered_x = powered_x * x;
    return powered_x;
}
double scaled_compute_polynom(double coeffs[], double x, int n) {
    double cmax = 0.0;
    for (int i = 0; i <= n && i < SIZE; i++) {
        if (cmax < fabs(coeffs[i]))
            cmax = fabs(coeffs[i]);
    }
    if (cmax > 1.0) {
        double Pn_cmax_x0 = coeffs[0] / cmax;
        for (int i = 1; i < n; i++)
            Pn_cmax_x0 = Pn_cmax_x0 * x + coeffs[i] / cmax;
        return Pn_cmax_x0 * cmax;
    } else {

        return compute_polynom(coeffs, x, n);
    }
}
double inv_quad(double (*fp)(double), double x0, double x1, double x2, double epsilon) {
    double f0, f1, f2, nominator_a, denominator_a, a, b, c;

    while (fabs(x2 - x1) >= epsilon) {
        // Evaluate the function at the current points
        f0 = fp(x0);
        f1 = fp(x1);
        f2 = fp(x2);

        /* We calculated by hand the formulas to find the value for each coefficient */

        // Compute coefficient 'a'
        nominator_a = (x2 - x0) * (f1 - f0) - (x1 - x0) * (f2 - f0);
        denominator_a = (f2 - f0) * (f1 - f0) * (f2 - f1);
        a = nominator_a / denominator_a;

        // Compute coefficient 'b'
        b = (x1 - x0 - a * (pow(f1, 2) - pow(f0, 2))) / (f1 - f0);

        // Compute coefficient 'c'
        c = x0 - pow(f0, 2) * a - f0 * b;

        // Shift the points for the next iteration
        x0 = x1;
        x1 = x2;
        x2 = c;
    }
}

int main() {
    double coeffs[9] = {10.1, 0.0, -20.2, 0.0, 30.3, 0.0, -60.6};
    double arr_x[3];

    x = sqrt(2.0);

    val2 = scaled_compute_polynom(coeffs, x, 7);

    printf("val1 = %lg, val2 = %lg\n", val1, val2);

    return 0;
}

int main()
{
    double x0, x1, x2, xstar, epsilon;

    printf("\nEnter x0, x1, x2, epsilon\n");
    scanf("%lf %lf %lf %lf", &x0, &x1, &x2, &epsilon);

    xstar =  inv_quad(f, x0, x1, x2, epsilon);

    printf("\nxstar = %lf\n", xstar);


} // main
double f(double x)
{
   double xp, sum;

  sum = x*x*x*x;
  sum = sum - 2.0*x*x*x;
  sum = sum + 3.0*x*x;
  sum = sum - 4.0*x;
  sum = sum - 4.0;

  return sum;
} // f

