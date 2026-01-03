#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NMAX 100
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- הגדרות טיפוסים ---
typedef double (*FUN_PTR)(double[]); 
typedef void (*GRAD_FUN_PTR)(double grad[], double x[]); 
typedef int VECTOR_CONVERGENCE_TEST(double arr[], int n, double epsilon);

// --- משתנים גלובליים לשימוש ב-steepest/golden ---
FUN_PTR objective_function;
double grad_vector[NMAX];
double xnm1[NMAX];
double gtemp[NMAX];
int vector_n;

// --- משתנים גלובליים לדוגמה הספציפית (מטריצה Q ווקטור b) ---
double Q[3][3];
double b[3];
double ftemp_arr[3];

// --- הצהרות מוקדמות ---
double falpha(double alpha);
double golden(double (*fp)(double), double a, double b, double eps);
void steepest(double xn[], double x0[], int n, FUN_PTR f, GRAD_FUN_PTR grad, double epsilon, VECTOR_CONVERGENCE_TEST v);

// --- פונקציות עזר כלליות ---

void copy_vector(double dest[], double source[], int n) {
    int i;
    for(i=0; i < n; i++) dest[i] = source[i]; 
}

int vector_convergence_test(double arr[], int n, double epsilon) {
    int i;
    for(i=0; i < n; i++)
        if (fabs(arr[i]) > epsilon) return 0; // עדיין לא התכנסנו
    return 1; // התכנסנו
}

// --- חישוב נגזרת נומרית (Approx Gradient) ---

double approx_partial_derivative(double (*obj_f)(double x[]), int i, double x[]) {
    double temp1, temp2, xi_orig, result, h;
    double eps_const = 1048576.0; // 2^20

    xi_orig = x[i];
    // הגנה מפני חלוקה ב-0 וקביעת גודל צעד יחסי
    h = (fabs(x[i]) < 0.000001) ? 0.000001 : x[i] / eps_const;

    x[i] = xi_orig + h;
    temp1 = (*obj_f)(x);

    x[i] = xi_orig - h;
    temp2 = (*obj_f)(x);

    result = (temp1 - temp2) / (2 * h);
    x[i] = xi_orig; // החזרת הערך המקורי

    return result;
}

void approx_g(double grad[], double x[]) {
    int i;
    // משתמשים ב-vector_n הגלובלי שמוגדר בפונקציה הראשית
    for (i = 0; i < vector_n; i++) {
        grad[i] = approx_partial_derivative(objective_function, i, x);
    }
}

// --- פונקציות עזר למטריצות (עבור הדוגמה) ---

void vector_matrix_mult(double res[], double x[], double Mat[3][3], int n) {
    int i, j;
    for(i=0; i<n; i++) {
        res[i] = 0;
        for(j=0; j<n; j++) {
            res[i] += x[j] * Mat[j][i];
        }
    }
}

double mult_vector_vector(double v1[], double v2[], int n) {
    double sum = 0;
    int i;
    for(i=0; i<n; i++) sum += v1[i] * v2[i];
    return sum;
}

// --- הסכמה הראשית (הפרויקט המסכם) ---

/*
 * פונקציה: mutli_variable_optimization_schema
 * מבצעת שני שלבים:
 * 1. מציאת x0 התחלתי ע"י סריקה (Coordinate Descent).
 * 2. הרצת steepest descent בלולאה תוך הקטנת epsilon עד להתייצבות.
 */
void mutli_variable_optimization_schema(
    double xn[], int n, 
    FUN_PTR f, GRAD_FUN_PTR grad, double epsilon,
    VECTOR_CONVERGENCE_TEST v, 
    double lowb, double highb, int initial_m)
{
    int i, j, k, pass;
    double h;
    double current_x[NMAX]; // וקטור זמני לבדיקת הערכים
    double f_prev, f_curr, f_next;
    
    // הגדרת משתנים גלובליים נדרשים
    vector_n = n;
    objective_function = f;

    // === שלב 1: מציאת נקודת התחלה x0 ===
    
    h = (highb - lowb) / (double)initial_m;
    double mid = (lowb + highb) / 2.0;

    // אתחול וקטור התוצאה לנקודת האמצע
    for(i = 0; i < n; i++) {
        xn[i] = mid;
    }

    printf("Starting Initial Search (Range: %.2f - %.2f, Parts: %d)...\n", lowb, highb, initial_m);

    // ביצוע שני מעברים (Passes)
    for(pass = 0; pass < 2; pass++) {
        // מעבר על כל מימד (משתנה) בנפרד
        for(k = 0; k < n; k++) {
            
            // העתקת המצב הנוכחי לוקטור זמני
            copy_vector(current_x, xn, n);

            // סריקת הטווח במימד k
            // מחפשים שלישייה: f(j) > f(j+1) < f(j+2)
            for(j = 0; j < initial_m - 1; j++) {
                
                // נקודה 1
                current_x[k] = lowb + j * h;
                f_prev = f(current_x);

                // נקודה 2 (האמצעית)
                current_x[k] = lowb + (j + 1) * h;
                f_curr = f(current_x);

                // נקודה 3
                current_x[k] = lowb + (j + 2) * h;
                f_next = f(current_x);

                // בדיקת תנאי "עמק"
                if (f_prev > f_curr && f_curr < f_next) {
                    // מצאנו עמק! מעדכנים את הקואורדינטה k בוקטור הראשי
                    xn[k] = lowb + (j + 1) * h;
                    // עוברים למשתנה הבא (Break פנימי)
                    break; 
                }
            }
        }
    }

    printf("Initial x0 found: ");
    for(i=0; i<n; i++) printf("%lf ", xn[i]);
    printf("\n\n");

    // === שלב 2: Steepest Descent עם עידון (Refinement) ===

    double x_prev_result[NMAX]; // לשמירת התוצאה הקודמת להשוואה
    double sum_diff, sum_prev;
    double ratio = 0.0;
    double current_epsilon = epsilon;
    
    // xn כרגע מכיל את ההערכה הראשונית. נשתמש בו כ-x0 ל-steepest.
    double x0_for_steepest[NMAX];
    copy_vector(x0_for_steepest, xn, n);

    int iter_count = 0;

    do {
        iter_count++;
        // שמירת התוצאה הקודמת לפני ההרצה (החל מהאיטרציה השנייה)
        if (iter_count > 1) {
            copy_vector(x_prev_result, xn, n);
        }

        printf("--- Refinement Iteration %d (Epsilon: %lg) ---\n", iter_count, current_epsilon);
        
        // הרצת האלגוריתם
        // שים לב: x0_for_steepest הוא הקלט, xn הוא הפלט
        steepest(xn, x0_for_steepest, n, f, grad, current_epsilon, v);

        // עדכון x0 לאיטרציה הבאה להיות התוצאה של האיטרציה הנוכחית
        copy_vector(x0_for_steepest, xn, n);

        // חישוב יחס היציבות (Stability Ratio)
        if (iter_count > 1) {
            sum_diff = 0.0;
            sum_prev = 0.0;
            for(i = 0; i < n; i++) {
                sum_diff += fabs(xn[i] - x_prev_result[i]);
                sum_prev += fabs(x_prev_result[i]);
            }
            
            if (sum_prev == 0.0) sum_prev = 0.0000001; // מניעת חלוקה באפס
            
            // יחס הקרבה: 1 פחות (השינוי היחסי). ככל שהשינוי קטן, היחס מתקרב ל-1.
            ratio = 1.0 - (sum_diff / sum_prev);
            printf("Stability Ratio: %lf\n", ratio);
        }

        // הקטנת אפסילון פי 2
        current_epsilon /= 2.0;

    // תנאי העצירה: רצים לפחות פעמיים, ועוצרים כשהיחס גדול מ-0.995
    } while (iter_count < 2 || ratio <= 0.995);

}

// --- מימוש פונקציות steepest ו-golden (מהתרגילים הקודמים) ---

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

double golden(double (*fp)(double), double x1, double x3, double eps) {
    double x2, fx2, fx3, x4, fx4;
    double phi = 1.618033988;
    double phi1 = 2.618033988;

    x2 = x1 + (x3 - x1) / phi1;
    fx2 = (*fp)(x2);
    x4 = x1 + (x3 - x1) / phi;
    fx4 = (*fp)(x4);

    while ((x3 - x1) > eps) {
        if (fx2 > fx4) {
            x1 = x2;
            x2 = x4;
            fx2 = fx4;
            x4 = x1 + (x3 - x1) / phi;
            fx4 = (*fp)(x4);
        } else {
            x3 = x4;
            x4 = x2;
            fx4 = fx2;
            x2 = x1 + (x3 - x1) / phi1;
            fx2 = (*fp)(x2);
        }
    }
    return ((x1 + x3) / 2);
}

void steepest(double xn[], double x0[], int n, 
              FUN_PTR f, GRAD_FUN_PTR grad, double epsilon,
              VECTOR_CONVERGENCE_TEST v) 
{
    double alpha_1, alpha_2, alpha_k;
    int i;

    vector_n = n;
    copy_vector(xn, x0, n);
    grad(grad_vector, x0);
    objective_function = f;
    copy_vector(xnm1, xn, n);
    grad(grad_vector, xnm1);

    int max_iter = 1000; // הגנה מפני לולאה אינסופית
    int iter = 0;

    while(v(grad_vector, n, 0.001) == 0 && iter < max_iter) {
        find_initial_alphas(falpha, &alpha_1, &alpha_2);
        alpha_k = golden(falpha, alpha_1, alpha_2, epsilon);   
        
        for(i=0; i < n; i++)
            xn[i] = xnm1[i] - alpha_k * grad_vector[i];

        copy_vector(xnm1, xn, n);
        grad(grad_vector, xn);
        iter++;
    }
}

// --- Main: הרצת הדוגמה עם המטריצה Q ---

// פונקציית המטרה לדוגמה 1: 0.5 * x'Qx + b'x
double f_example1(double x[]) {
    double sum;
    vector_matrix_mult(ftemp_arr, x, Q, 3);
    sum = 0.5 * mult_vector_vector(ftemp_arr, x, 3);
    sum = sum + mult_vector_vector(b, x, 3);
    return sum;
}

int main() {
    double xstar[10]; // וקטור התוצאה
    
    // הגדרת המטריצה Q (סימטרית, מוגדרת חיובית)
    double qarr[9] = {2, -1, 0, -1, 2, -1, 0, -1, 2};
    int i, j;

    for(i=0; i < 3; i++)
        for(j=0; j < 3; j++)
            Q[i][j] = qarr[i*3 + j];

    printf("Q Matrix: \n");
    for(i=0; i < 3; i++) {
        for(j=0; j < 3; j++)  
            printf(" %6.2lf ", Q[i][j]);
        printf("\n");
    }

    // הגדרת הוקטור b
    b[0] = -1;
    b[1] = -2;
    b[2] = -3;
    
    // קריאה לפונקציה הראשית של הפרויקט
    // תחום חיפוש 0.0 עד 10.0, חלוקה ל-100 חלקים
    mutli_variable_optimization_schema(xstar, 3, 
       f_example1, approx_g,
       0.00001, vector_convergence_test,
       0.0, 10.0, 100);

    // הדפסת התוצאות הסופיות
    printf("\n\noptimal solution:\n xstar[0] = %lf\n,  xstar[1] = %lf,  xstar[2] = %lf\n",
           xstar[0], xstar[1], xstar[2] );

    printf("\n\nIn degrees (Approx): xstar[0] = %lf\n,  xstar[1] = %lf\n",
           xstar[0]*180.0/M_PI, xstar[1]*180.0/M_PI);

    printf("\n\noptimal value  = %lf\n", f_example1(xstar));

    return 0;
}