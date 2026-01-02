#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NMAX 100
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Typedefs ---
typedef double (*FUN_PTR)(double[]); 
typedef void (*GRAD_FUN_PTR)(double grad[], double x[]); 
typedef int VECTOR_CONVERGENCE_TEST(double arr[], int n, double epsilon);

// --- Globals used by steepest/golden ---
FUN_PTR objective_function;
double grad_vector[NMAX];
double xnm1[NMAX];
double gtemp[NMAX];
int vector_n;

// --- Globals for the Example Problem (Quadratic) ---
double Q[3][3];
double b[3];
double ftemp_arr[3];

// --- Forward Declarations ---
double falpha(double alpha);
double golden(double (*fp)(double), double a, double b, double eps);
void steepest(double xn[], double x0[], int n, FUN_PTR f, GRAD_FUN_PTR grad, double epsilon, VECTOR_CONVERGENCE_TEST v);

// --- Helper Functions ---

void copy_vector(double dest[], double source[], int n) {
    int i;
    for(i=0; i < n; i++) dest[i] = source[i]; 
}

int vector_convergence_test(double arr[], int n, double epsilon) {
    int i;
    for(i=0; i < n; i++)
        if (fabs(arr[i]) > epsilon) return 0;
    return 1;
}

// Numerical Gradient Calculation (approx_g)
double approx_partial_derivative(double (*obj_f)(double x[]), int i, double x[]) {
    double temp1, temp2, xi_orig, result, h;
    double eps_const = 1048576.0; // 2^20

    xi_orig = x[i];
    h = (x[i] == 0) ? 0.000001 : x[i] / eps_const;

    x[i] = xi_orig + h;
    temp1 = (*obj_f)(x);

    x[i] = xi_orig - h;
    temp2 = (*obj_f)(x);

    result = (temp1 - temp2) / (2 * h);
    x[i] = xi_orig;

    return result;
}

void approx_g(double grad[], double x[]) {
    int i;
    // Assuming vector_n is set globally or passed. Using global vector_n safely if set, else assume from context (risky but standard for this course structure)
    // Better to use the 'n' known from context. In the schema, we set vector_n.
    for (i = 0; i < vector_n; i++) {
        grad[i] = approx_partial_derivative(objective_function, i, x);
    }
}

// Matrix/Vector helpers for the example function
void vector_matrix_mult(double res[], double x[], double Mat[3][3], int n) {
    int i, j;
    for(i=0; i<n; i++) {
        res[i] = 0;
        for(j=0; j<n; j++) {
            res[i] += x[j] * Mat[j][i]; // Note: usually x*Q means row vector x times matrix
        }
    }
}

double mult_vector_vector(double v1[], double v2[], int n) {
    double sum = 0;
    int i;
    for(i=0; i<n; i++) sum += v1[i] * v2[i];
    return sum;
}

// --- The Core Optimization Schema (The Assignment) ---

/*
 * This function implements the requested schema:
 * 1. Finds initial x0 using coordinate descent-like search.
 * 2. Refines the solution by repeatedly calling steepest descent with decreasing epsilon.
 */
void mutli_variable_optimization_schema(
    double xn[], int n, 
    FUN_PTR f, GRAD_FUN_PTR grad, double epsilon,
    VECTOR_CONVERGENCE_TEST v, 
    double lowb, double highb, int initial_m)
{
    int i, j, k, pass;
    double h;
    double current_x[NMAX]; // To verify valley condition
    double f_curr, f_next, f_prev;
    double val_curr;
    
    // Set globals needed for helpers
    vector_n = n;
    objective_function = f;

    // --- Part 1: Finding Initial x0 ---
    
    h = (highb - lowb) / (double)initial_m;
    double mid = (lowb + highb) / 2.0;

    // Initialize xm with midpoints
    for(i = 0; i < n; i++) {
        xn[i] = mid;
    }

    // Run 2 passes
    for(pass = 0; pass < 2; pass++) {
        // Iterate over each dimension
        for(k = 0; k < n; k++) {
            
            // Search along dimension k
            // We look for a triplet (j, j+1, j+2)
            // We need to construct vectors to evaluate f
            
            // Copy current best state to a temp vector for manipulation
            copy_vector(current_x, xn, n);

            for(j = 0; j < initial_m - 1; j++) {
                // Point 1: L + j*h
                current_x[k] = lowb + j * h;
                f_prev = f(current_x);

                // Point 2: L + (j+1)*h
                current_x[k] = lowb + (j + 1) * h;
                f_curr = f(current_x);

                // Point 3: L + (j+2)*h
                current_x[k] = lowb + (j + 2) * h;
                f_next = f(current_x);

                // Check Valley Condition: f(prev) > f(curr) < f(next)
                if (f_prev > f_curr && f_curr < f_next) {
                    // Found a valley! Update xn[k] to the middle point
                    xn[k] = lowb + (j + 1) * h;
                    // Move to next dimension
                    break; 
                }
            }
        }
    }

    printf("Initial x0 found: ");
    for(i=0; i<n; i++) printf("%lf ", xn[i]);
    printf("\n\n");

    // --- Part 2: Steepest Descent with Refinement ---

    double x_prev_result[NMAX];
    double sum_diff, sum_prev;
    double ratio = 0.0;
    double current_epsilon = epsilon;
    
    // Use the found xn as the starting point x0 for steepest
    double x0_for_steepest[NMAX];
    copy_vector(x0_for_steepest, xn, n);

    int iter_count = 0;

    do {
        iter_count++;
        // Keep result of previous run (if not first run)
        if (iter_count > 1) {
            copy_vector(x_prev_result, xn, n);
        }

        printf("--- Refinement Iteration %d (Epsilon: %lg) ---\n", iter_count, current_epsilon);
        
        // Run steepest descent
        // Note: xn is output, x0_for_steepest is input
        steepest(xn, x0_for_steepest, n, f, grad, current_epsilon, v);

        // Update x0 for next iteration to be the result of this one
        copy_vector(x0_for_steepest, xn, n);

        // Calculate ratio for stopping condition
        if (iter_count > 1) {
            sum_diff = 0.0;
            sum_prev = 0.0;
            for(i = 0; i < n; i++) {
                sum_diff += fabs(xn[i] - x_prev_result[i]);
                sum_prev += fabs(x_prev_result[i]);
            }
            
            // Avoid division by zero
            if (sum_prev == 0.0) sum_prev = 0.0000001; 
            
            // The logic: We want stability. 
            // If the change (sum_diff) is very small relative to the values (sum_prev), 
            // it means the ratio (1 - change_ratio) is close to 1.
            // Let's implement exactly what the prompt implies: similarity ratio.
            // Ratio = 1 - (diff / prev). If diff is small, Ratio -> 1.
            
            ratio = 1.0 - (sum_diff / (sum_prev + 0.000001)); // Add small epsilon to prevent overflow
            printf("Stability Ratio: %lf\n", ratio);
        }

        // Halve epsilon for next run
        current_epsilon /= 2.0;

    } while (iter_count < 2 || ratio <= 0.995); // Ensure at least 2 runs to calc ratio

}

// --- Steepest Descent Implementation (from previous lectures) ---

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