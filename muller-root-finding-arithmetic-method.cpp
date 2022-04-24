/*
 *  Demonstration of the Müller's root-finding arithmetic method.
 *
 *  Copyright (C) 2009 Efstathios Chatzikyriakidis (stathis.chatzikyriakidis@gmail.com)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/*  We use the function: f(x) = x^6 - 2.
 *
 *  Assume we run the program.
 *
 *  Example 1.
 *
 *  Input.
 *
 *  x0 = 1, x1 = 2
 *  iterations = 20
 *  tolerance = 15
 * 
 *  Output.
 *
 *  Root x = +1.122462048309e+00
 *
 *  Example 2.
 *
 *  Input.
 *
 *  x0 = -1, x1 = -2
 *  iterations = 20
 *  tolerance = 15
 *
 *  Output.
 * 
 *  Root x = -1.122462048309e+00
 */

/* include standard c/c++ headers. */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
using namespace std;

/* max iterations and max tolerance digits. */
#define MAX_ITERATIONS 1000000
#define MAX_TOLERANCES 40

/* function that returns the y = f(x). */
double
function (double x)
{
  return (pow (x, 6) - 2);
}

/* function that returns the sign of a number. */
int
sign (double x)
{
  /* return the aproppriate sign value. */
  if (x > 0) return 1;
  if (x < 0) return -1;

  /* zero number does not has a sign. */
  return 0;
}

/* main function entry point. */
int
main (void)
{
  /* the initial edge points in x axes. */
  double x0, x1;

  /* tolerance digits, user
     iterations and current
     iteration number. */
  int n, usr_iters, iters;

  /* the tolerance in floating point value and
     the value of s formula in each iteration. */
  double tol, s;

  /* root find binary flag. */
  bool find_root = false;

  /* pointers to dynamic arrays (statistical tables). */
  double *x, *y, *d, *c;

  /* size of the statistical arrays. */
  int arr_size;

  /* init pointers to NULL. */
  x = y = d = c = NULL;

  /* print a welcome message.*/
  cout << "Muller's root-finding arithmetic method." << endl << endl;

  /* input the two edge points. */
  cout << "Input point x0: "; cin >> x0;
  cout << "Input point x1: "; cin >> x1;

  /* restriction: values x0, x1 ought to be different. */
  if (x0 == x1)
    {
      cout << endl << "Values x0, x1 should be different." << endl;
      return EXIT_FAILURE;
    }

  /* input number of iterations. */
  cout << "Iterations: "; cin >> usr_iters;

  /* restriction: iterations must be > 2 and <= MAX_ITERATIONS. */
  if (usr_iters <= 2 || usr_iters > MAX_ITERATIONS)
    {
      cout << endl << "Iterations value: 2<i<=" << MAX_ITERATIONS << endl;
      return EXIT_FAILURE;
    }

  /* input tolerance in digits. */
  cout << "Tolerance: "; cin >> n;

  /* restriction: tolerance must be > 0 and <= MAX_TOLERANCES. */
  if (n <= 0 || n > MAX_TOLERANCES)
    {
      cout << endl << "Tolerance value: 0<t<=" << MAX_TOLERANCES << endl;
      return EXIT_FAILURE;
    }

  /* the array size is the user input iterations plus one. */
  arr_size = usr_iters + 1;

  /* allocate enough memory for each array. we
     are doing this because we want to have the
     historic of all the statistical data. */

  /* x points on each iteration. */
  x = (double *) malloc (arr_size * sizeof (double));

  /* f(x) values on each iteration. */
  y = (double *) malloc (arr_size * sizeof (double));

  /* essential formulas on each iteration. */
  d = (double *) malloc (arr_size * sizeof (double));
  c = (double *) malloc (arr_size * sizeof (double));

  /* print a new line. */
  cout << endl;

  /* check for allocation problems. */
  if (x == NULL || y == NULL || d == NULL || c == NULL)
    {
      cout << "Dynamic memory allocation problem (";
      cout << arr_size * sizeof (double);
      cout << " bytes)." << endl;
      return EXIT_FAILURE;
    }

  /* init all arrays to zero. */
  for (int i = 0; i < arr_size; ++i)
    x[i] = y[i] = d[i] = c[i] = 0;

  /* solve the equation. */

  /* calculate tolerance in floating point. */
  tol = 0.5 * pow (10, -n);

  /* store initial edges points. */
  x[0] = x0;
  x[1] = x1;

  /* calculate and store the center point. */
  x[2] = (x[0] + x[1]) / 2;

  /* calculate and store initial f(x) values. */
  y[0] = function (x[0]);
  y[1] = function (x[1]);
  y[2] = function (x[2]);

  /* store the initial c formula. */
  c[0] = (y[1] - y[0]) / (x[1] - x[0]);

  /* loop iteration code - evaluate formulas. */
  for (int i = 2; i < usr_iters; ++i)
    {
      /* keep the current iteration. */
      iters = i;

      /* calculate the c and d formulas. */
      c[i-1] = (y[i-0] - y[i-1]) / (x[i] - x[i-1]);
      d[i-2] = (c[i-1] - c[i-2]) / (x[i] - x[i-2]);

      /* calculate the s formula. */
      s = c[i-1] + (x[i] - x[i-1]) * d[i-2];

      /* find next approching root value. */
      x[i+1] = x[i] - 2*y[i] / (s+sign (s) *
               sqrt (abs (pow (s, 2) - 4*y[i]*d[i])));

      /* find next f(x). */
      y[i+1] = function (x[i+1]);

      /* check the tolerance of the root. */
      if (abs (x[i+1] - x[i]) < tol)
        {
          /* we have found a root. */
          cout << "Method has converged to a root." << endl;
          find_root = true;
          break;
        }
    }

  /* check why the loop ends. */
  if (iters >= usr_iters-1)
    if (find_root == false)
      /* we haven't found a root. */
      cout << "Method it didn't reach the allowed tolerance." << endl;

  /* print a new line. */
  cout << endl;

  /* print header titles. */
  printf ("%-11s%-25s%-25s", "|step]", "|x]", "|f(x)]");
  printf ("%-25s%-s\n", "|d]", "|c]");

  /* print statistical data for all iterations. */
  for (int i = 0; i <= iters; ++i)
    {
      printf ("| %08d | %+022.12e | %+022.12e ", i+1, x[i], y[i]);
      printf ("| %+022.12e | %+022.12e\n", d[i], c[i]);
    }

  /* free the dynamic allocated memory. */
  free (x); free (y);
  free (d); free (c);

  /* successfully return. */
  return EXIT_SUCCESS;
}
