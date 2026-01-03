// steepest_golden7.c - הגרסה המתקדמת עם נגזרת נומרית ותנאי עצירה משופר

#include <stdio.h>
#include <math.h>

#define NMAX 100

// --- הגדרות טיפוסים לפונקציות ---
typedef double (*FUN_PTR)(double[]); // מצביע לפונקציית המטרה
typedef void (*GRAD_FUN_PTR)(double grad[], double x[] ); // מצביע לפונקציית הגרדיאנט
typedef int VECTOR_CONVERGENCE_TEST(double arr[], int n, double epsilon); // מצביע לפונקציית בדיקת התכנסות

// --- משתנים גלובליים ---
FUN_PTR objective_function; // הפונקציה אותה אנו מנסים למזער
double grad_vector[NMAX];   // וקטור הגרדיאנט (השיפוע)
double xnm1[NMAX];          // המיקום הנוכחי (לפני הצעד)
double gtemp[NMAX];         // וקטור עזר לחישובים זמניים
double diff[NMAX];          // וקטור השינוי (ההפרש בין המיקום החדש לישן)
int vector_n;               // מספר המשתנים (ממד הבעיה)

// פונקציית עזר לחיסור וקטורים: dest = dest - source
void subtract_vector(double dest[], double source[], int n)
{
  int i;
  for(i=0; i < n; i++)
    dest[i] = dest[i] - source[i];
} 

/*
 * פונקציה: falpha
 * ----------------
 * מטרה: הופכת את הבעיה הרב-ממדית לחד-ממדית עבור "חיפוש קווי".
 * היא בודקת: מה יהיה ערך הפונקציה אם נזוז מרחק alpha בכיוון המינוס-גרדיאנט?
 */
double falpha(double alpha)
{
  int i;
  
  // הדפסת מידע לדיבוג
  printf("vector_n = %d, alpha = %lf\n", vector_n, alpha);

  // חישוב המיקום החדש הזמני: x_new = x_old - alpha * gradient
  for(i=0; i < vector_n; i++)
    gtemp[i] = xnm1[i] - alpha*grad_vector[i]; 

  // הדפסת המיקום החדש והגרדיאנט
  for(i=0; i < vector_n; i++)
      printf(
    "gtemp[%d] = %lf,  xnm1[%d] = %lf, grad_vector[%d] = %lf\n",
      i, i, i,
           gtemp[i],  xnm1[i], grad_vector[i]);

  // החזרת ערך הפונקציה בנקודה החדשה
  return  objective_function(gtemp);
} 


/*
 * פונקציה: golden
 * ----------------
 * מטרה: מציאת המינימום של פונקציה במשתנה אחד (alpha) בשיטת יחס הזהב.
 * קלט: פונקציה לבדיקה (fp), ותחום חיפוש [x1, x3].
 * פלט: ה-alpha האופטימלי (גודל הצעד הכי טוב).
 */
double golden(double (*fp)(double), double x1, double x3, double eps)
{
 double x2, fx2,  fx3, x4, fx4;
 double phi = 1.618;    // יחס הזהב
 double phi1 = 2.618;   // קבוע עזר

 // צעד ראשון: בחירת שתי נקודות פנימיות x2, x4 לפי יחס הזהב
 x2 = x1 + (x3-x1)/phi1;
 fx2 = (*fp)(x2);
 x4 = x1 + (x3-x1)/phi;
 fx4 = (*fp)(x4);

 // לולאה: הקטנת הטווח עד שההפרש קטן מ-eps
 do {
   if (fx2 > fx4)
   {
     // המינימום בחלק הימני -> זורקים את החלק השמאלי (x1 זז ל-x2)
     x1 = x2;
     x2 = x4;
     fx2 = fx4;
     x4 = x1 + (x3-x1)/phi; // חישוב נקודה חדשה
     fx4 = (*fp)(x4);
   }
   else 
   {
     // המינימום בחלק השמאלי -> זורקים את החלק הימני (x3 זז ל-x4)
     x3 = x4;
     x4 = x2;
     fx4 = fx2;
     x2 = x1 + (x3-x1)/phi1; // חישוב נקודה חדשה
     fx2 = (*fp)(x2);
   } 
 } while ( (x3 - x1) > eps);

 // החזרת אמצע הטווח שנשאר
 return ( (x1+x3)/2);
}


/*
 * פונקציה: vector_convergence_test
 * --------------------------------
 * מטרה: בדיקה האם וקטור מסוים "התאפס" (כל איבריו קרובים לאפס).
 * משמשת לבדיקת שני דברים:
 * 1. האם הגרדיאנט התאפס? (תנאי לקיצון)
 * 2. האם הצעד שעשינו (diff) התאפס? (תנאי להתכנסות)
 */
int vector_convergence_test(double arr[], int n, double epsilon)
{
  int i;
  for(i=0; i < n; i++)
    printf("arr[%d] = %lf\n", i, arr[i]); // הדפסה לדיבוג

  for(i=0; i < n; i++)
   if (fabs(arr[i]) > epsilon)
     return 0; // לא התכנס (נמצא ערך גדול מאפסילון)
  return 1;    // התכנס (כל הערכים קטנים מאפסילון)
}

// פונקציית עזר להעתקת וקטור
void copy_vector(double dest[], double source[], int n)
{
  int i;
  for(i=0; i < n; i++)
   dest[i] = source [i];  
}

/*
 * פונקציה: find_initial_alphas
 * -----------------------------
 * מטרה: מציאת טווח (Bracketing) לחיפוש ה-alpha.
 * מתחילה בצעד קטן ומכפילה אותו עד שהפונקציה מתחילה לעלות.
 * זה מבטיח שהמינימום נמצא איפשהו בין 0 לבין alpha2.
 */
void find_initial_alphas(double (*falpha)(double), double *alpha_1, double *alpha_2)
{
  int going_down_flag;
  double falpha1, falpha2, alpha1, alpha2, prev_alpha;

  falpha1 = (*falpha)(0.0);
  prev_alpha = alpha1 = 0.0;
  alpha2 = 0.0009765625; // התחלה מצעד קטן מאוד

  going_down_flag =  1;

  while(going_down_flag == 1)