// steepest.c

#include <stdio.h>
#include <math.h>

#define NMAX 100


typedef double (*FUN_PTR)(double[]); // pointer to function 

typedef void (*GRAD_FUN_PTR)(double grad[], double x[] ); // pointer 
//to function 

typedef int VECTOR_CONVERGENCE_TEST(double arr[], int n, double epsilon);


FUN_PTR objective_function;

double grad_vector[NMAX];
double xnm1[NMAX];
double gtemp[NMAX];

int vector_n;

double falpha(double alpha)
{
  int i;

 printf("vector_n = %d, alpha = %lf\n", vector_n, alpha);

  for(i=0; i < vector_n; i++)
    gtemp[i] = xnm1[i] - alpha*grad_vector[i]; 

  for(i=0; i < vector_n; i++)
      printf(
    "gtemp[%d] = %lf,  xnm1[%d] = %lf, grad_vector[%d] = %lf\n",
      i,  gtemp[i], i, xnm1[i], i, grad_vector[i]);

  return  objective_function(gtemp);

} // falpha



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
  } /* else */  
} while ( (x3 - x1) > eps);

 return ( (x1+x3)/2);
} /* golden */




int vector_convergence_test(double arr[], int n, double epsilon)
{
  int i;

  for(i=0; i < n; i++)
    printf("arr[%d] = %lf\n", i, arr[i]);

  for(i=0; i < n; i++)
   if (fabs(arr[i]) > epsilon)
     return 0.0;
  return 1;

} // vector_convergence_test

void copy_vector(double dest[], double source[], int n)
{
  int i;
  
  for(i=0; i < n; i++)
   dest[i] = source [i];  

} // copy_vector


void find_initial_alphas(double (*falpha)(double), double *alpha_1, double *alpha_2)
{
  int going_down_flag;
  double falpha1, falpha2, alpha1, alpha2, prev_alpha;

  falpha1 = (*falpha)(0.0);
  prev_alpha = alpha1 = 0.0;

  alpha2 = 0.0009765625; // 1/1024

  going_down_flag =  1;

  while(going_down_flag == 1)
  {
    prev_alpha = alpha1;
    falpha1 = (*falpha)(alpha1);
    falpha2 = (*falpha)(alpha2);
    alpha1 = alpha2; 
    if(falpha2 >= falpha1)
       going_down_flag = 0;
    else
       alpha2 = 2.0*alpha2; 
   
  } // while

  *alpha_1 = prev_alpha;
  *alpha_2 = alpha2;
} // find_initial_alphas



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

  while(v(grad_vector, n, 0.001) == 0)
  {
      
    find_initial_alphas(falpha, &alpha_1, &alpha_2);
  
    alpha_k = golden(falpha, alpha_1, alpha_2, epsilon);   
    
    printf("alpha_k = %lf\n", alpha_k);
    
    for(i=0; i < n; i++)
      xn[i] = xnm1[i] - alpha_k * grad_vector[i];

    printf("xn:\n");
    for(i=0; i < n; i++)
      printf(" xn[%d] = %lf\n",  i, xn[i]);
    copy_vector(xnm1, xn, n);
    grad(grad_vector, xn);

  } // while 

 
}  // steepest


double f(double x[])
{
  return (-sin(x[0]+2*x[1]) - cos(3*x[0]+4*x[1])); 

} // f

void g(double grad[], double x[])
{
 
   grad[0] = -cos(x[0] + 2*x[1]) + 3*sin(3*x[0] + 4*x[1]); 
   grad[1] = -2*cos(x[0] + 2*x[1]) + 4*sin(3*x[0] + 4*x[1]); 

  printf(" grad[0] = %lf,  grad[1] = %lf\n", grad[0], grad[1]);

} // f


int main()
{
  double xstar[2], x0[2], value;

  x0[0] = 0.0;
  x0[1] = 0.0;

  steepest(xstar, x0, 2, 
   f, g, 0.0001, vector_convergence_test);

  printf("\noptimal solution:\n xstar[0] = %lf\n,  xstar[1] = %lf\n",
                xstar[0], xstar[1]);

  printf("\nIn Degrees:\n xstar[0] = %lf\n,  xstar[1] = %lf\n",
                xstar[0]*180.0/M_PI, xstar[1]*180.0/M_PI);


  printf("\noptimal value  = %lf\n", f(xstar));


} // main