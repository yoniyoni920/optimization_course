#define SIZE 9
#include <math.h>
#include <stdio.h>
double power_of_x(double x, int power) ;
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

int main() {
    double coeffs[9] = {10.1, 0.0, -20.2, 0.0, 30.3, 0.0, -60.6};
    double x, val1, val2;

    x = sqrt(2.0);

    val1 = compute_polynom(coeffs, x, 7);
    val2 = scaled_compute_polynom(coeffs, x, 7);

    printf("val1 = %lg, val2 = %lg\n", val1, val2);

    return 0;
}
