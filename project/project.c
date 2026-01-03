#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NMAX 100



typedef double (*FUN_PTR)(double[]); 
typedef void (*GRAD_FUN_PTR)(double grad[], double x[]); 
typedef int VECTOR_CONVERGENCE_TEST(double arr[], int n, double epsilon);


FUN_PTR objective_function;
double grad_vector[NMAX];
double xnm1[NMAX];
double gtemp[NMAX];
int vector_n;
double diff[NMAX];



double falpha(double alpha);
double golden(double (*fp)(double), double a, double b, double eps);
void steepest(double xn[], double x0[], int n, FUN_PTR f, GRAD_FUN_PTR grad, double epsilon, VECTOR_CONVERGENCE_TEST v);


double f(double x[])
{
    int i, j;

    return (-sin(x[0]+3*x[1] + 15*x[2]) - sin(3*x[0]+3*x[1] + 9*x[2]) - sin(x[0]+x[1] + 18*x[2]) );



} // f
void copy_vector(double dest[], double source[], int n)
{
    int i;
    for(i=0; i < n; i++)
        dest[i] = source [i];
}

int vector_convergence_test(double arr[], int n, double epsilon)
{
    int i;

    for(i=0; i < n; i++)
        if (fabs(arr[i]) > epsilon)
            return 0;
    return 1;
}

double approx_partial_derivative(double (*obj_f)(double x[]),
     int i, double x[])
{
    double temp1,temp2, xi_orig, result, h;
    double eps_const = 1048576.0;

    xi_orig = x[i];
    h = x[i]/eps_const;

    x[i] =  xi_orig + h;

    temp1 = (*obj_f)(x);

    x[i] =  xi_orig - h;

    temp2 = (*obj_f)(x);

    result =  (temp1 - temp2)/(2*h);

    x[i] =  xi_orig;

    return result;

} // approx_partial_derivative

void approx_g(double grad[], double x[])
{
    int i, j;

    for(i=0;i < vector_n; i++)
        grad[i] = approx_partial_derivative(f,i, x);

} // approx_g

void vector_mult_matrix(double res[], double x[], double Mat[3][3], int n) {
    int i, j;
    for(i=0; i<n; i++) {
        res[i] = 0;
        for(j=0; j<n; j++) {
            res[i] += x[j] * Mat[j][i];
        }
    }
}

double vector_mult_vector(double v1[], double v2[], int n) {
    double sum = 0;
    int i;
    for(i=0; i<n; i++) sum += v1[i] * v2[i];
    return sum;
}



void mutli_variable_optimization_schema(
    double xn[], int n, 
    FUN_PTR f, GRAD_FUN_PTR grad, double epsilon,
    VECTOR_CONVERGENCE_TEST v, 
    double lowb, double highb, int initial_m)
{
    int i, j, k;
    double h;
    double current_x[NMAX];
    double f_prev, f_curr, f_next;
    

    vector_n = n;
    objective_function = f;


    
    h = (highb - lowb) / (double)initial_m;
    double mid = (lowb + highb) / 2.0;

    for(i = 0; i < n; i++) {
        xn[i] = mid;
    }



    for(i = 0; i < 2; i++) {

        for(k = 0; k < n; k++) {
            

            copy_vector(current_x, xn, n);


            for(j = 0; j < initial_m - 1; j++) {
                

                current_x[k] = lowb + j * h;
                f_prev = f(current_x);

                current_x[k] = lowb + (j + 1) * h;
                f_curr = f(current_x);


                current_x[k] = lowb + (j + 2) * h;
                f_next = f(current_x);


                if (f_prev > f_curr && f_curr < f_next) {

                    xn[k] = lowb + (j + 1) * h;

                    break; 
                }
            }
        }
    }






    double x_prev_result[NMAX];
    double sum_diff, sum_prev;
    double ratio = 0.0;
    double current_epsilon = epsilon;
    

    double x0_for_steepest[NMAX];
    copy_vector(x0_for_steepest, xn, n);

    int iter_count = 0;

    do {
        iter_count++;

        if (iter_count > 1) {
            copy_vector(x_prev_result, xn, n);
        }



        steepest(xn, x0_for_steepest, n, f, grad, current_epsilon, v);


        copy_vector(x0_for_steepest, xn, n);


        if (iter_count > 1) {
            sum_diff = 0.0;
            sum_prev = 0.0;
            for(i = 0; i < n; i++) {
                sum_diff += fabs(xn[i] - x_prev_result[i]);
                sum_prev += fabs(x_prev_result[i]);
            }
            
            if (sum_prev == 0.0) sum_prev = 0.0000001;
            

            ratio = 1.0 - (sum_diff / sum_prev);

        }


        current_epsilon /= 2.0;


    } while (iter_count < 2 || ratio <= 0.995);

}


void find_initial_alphas(double (*falpha)(double), double *alpha_1, double *alpha_2) {
    int going_down_flag;
    double falpha1, falpha2, alpha1, alpha2, prev_alpha;

    falpha1 = (*falpha)(0.0);
    prev_alpha = alpha1 = 0.0;
    alpha2 = 0.0009765625; // 1/1024

    going_down_flag = 1;
    while(going_down_flag == 1) {
        prev_alpha = alpha1;
        falpha1 = (*falpha)(alpha1);
        falpha2 = (*falpha)(alpha2);
        
        if(falpha2 >= falpha1)
            going_down_flag = 0;
        else {
            alpha1 = alpha2; 
            alpha2 = 2.0 * alpha2; 
        }
    } 
    *alpha_1 = prev_alpha;
    *alpha_2 = alpha2;
}

double falpha(double alpha) {
    int i;
    for(i=0; i < vector_n; i++)
        gtemp[i] = xnm1[i] - alpha * grad_vector[i]; 
    return objective_function(gtemp);
}
double golden(double (*fp)(double), double x1, double x3, double eps)
{
    double x2, fx2,  fx3, x4, fx4;
    double phi = 1.618;
    double phi1 = 2.618;


    x2 = x1 + (x3-x1)/phi1;
    fx2 = (*fp)(x2);
    x4 = x1 + (x3-x1)/phi;
    fx4 = (*fp)(x4);


    do {
        if (fx2 > fx4)
        {

            x1 = x2;
            x2 = x4;
            fx2 = fx4;
            x4 = x1 + (x3-x1)/phi;
            fx4 = (*fp)(x4);
        }
        else
        {

            x3 = x4;
            x4 = x2;
            fx4 = fx2;
            x2 = x1 + (x3-x1)/phi1;
            fx2 = (*fp)(x2);
        }
    } while ( (x3 - x1) > eps);


    return ( (x1+x3)/2);
}


void steepest(double xn[], double x0[], int n,
FUN_PTR f, GRAD_FUN_PTR grad, double epsilon,
VECTOR_CONVERGENCE_TEST v)
{

    double temp, alpha_1, alpha_2, alpha_k;
    int i;

    vector_n = n;
    copy_vector(xn, x0, n);
    grad(grad_vector, x0);
    objective_function = f;
    copy_vector(xnm1, xn, n);
    grad(grad_vector, xnm1);
    copy_vector(diff, xn, n);

    while((v(diff, n, 0.001) == 0) || (v(grad_vector, n, 0.001) == 0) )
    {

        find_initial_alphas(falpha, &alpha_1, &alpha_2);

        alpha_k = golden(falpha, alpha_1, alpha_2, epsilon);

        //printf("alpha_k = %lf\n", alpha_k);

        for(i=0; i < n; i++)
            diff[i] =  alpha_k * grad_vector[i];
        for(i=0; i < n; i++)
            xn[i] = xnm1[i] - diff[i];

        //printf("xn:\n");
        for(i=0; i < n; i++)
        //    printf(" xn[%d] = %lf\n",  i, xn[i]);

        copy_vector(xnm1, xn, n);
        grad(grad_vector, xn);

    } // while


}  // steepest





int main()
{
    double xstar[10],  value;
    int i, j;


    mutli_variable_optimization_schema(xstar,  3,
       f, approx_g, 0.00001, vector_convergence_test,
       0.0, 1.0,  100);


    printf("\n\noptimal solution:\n xstar[0] = %lf\n,  xstar[1] = %lf,  xstar[2] = %lf\n",
                  xstar[0], xstar[1],  xstar[2] );

    printf("\n\nIn degrees: xstar[0] = %lf\n,  xstar[1] = %lf\n,  xstar[2] = %lf\n",
                  xstar[0]*180.0/M_PI, xstar[1]*180.0/M_PI, xstar[2]*180.0/M_PI);

    printf("\n\noptimal value  = %lf\n", f(xstar));


} // main
