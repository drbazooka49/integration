#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

using namespace std;
 
//TODO: relative error and absolute error

//Midpoint rule

//integral function
double function(double x)
{
    return (1 / ( 1 + ( x * x)));
}

void MidpointRule(double a, double b, int n)
{
    double xi;
    double INTGRL = 0.0; //running sum
    

    double width = double((b-a)/n); 
 
    for (double i = a; i <= n ; i++)
    {
        xi = a + (i - 1) * width;
        INTGRL += (width * function(xi)); // width * height of ith rectangle 
    }
    
    cout << "Midpoint rule result: " << INTGRL << "\n";
}

//Trapezoidal rule
//TODO if n == 0
void TrapezoidalRule(double a, double b, int n)
{
  double sum = 0.0;
  double step;
  int i;
 
  step = (b - a) / (1.0 * n);
  for ( i = 1 ; i < n ; ++i ) {
    sum += function(a + i * step);
  }
  sum += (function(a) + function(b)) / 2;
  sum *= step; 
  cout << "Trapezoidal rule result: " << sum << "\n";
}

//Simpson's Rule
void SimpsonsRule(double a, double b, int n)
{
    double sum = 0.0;
    double k;
    int i;
    double step = (b - a) / n;

    sum = function(a) + function(b);
    for(i = 1; i < n; i++)
    {
        k = a + i * step;
        if(i % 2 == 0)
        {
            sum += (2 * function(k));
        }
        else
        {
            sum += (4 * function(k));
        }
    }

    sum *= step / 3;

    cout << "Simpsons rule result: " << sum << "\n";
}

//Bode formula
void BodeFormula(double a, double b, double eps)
{   
    double sum = 0.0;
    double temp_sum;
    double step;
    double temp_step;
    int n = 1;
    int i;
 
    n = 1;
    do
    {
        temp_sum = sum; 
        sum = 0.0; 
        n *= 2;
        step = (b - a) / n;
        temp_step = step / 4;
 
        for (i = 1; i <= n; ++i)
            sum += 7 * function(a + i * step - 4 * temp_step) + 32 * function(a + i * step - 3 * temp_step) + 12 * function(a + i * step - 2 * temp_step) +
                32 * function(a + i * step - temp_step) + 7 * function(a + i * step);
 
        sum *= 2 * temp_step / 45;
 
    } while (fabs(temp_sum - sum) < eps);
 
    cout << "Bode formula result: " << sum << "\n";
}

//Gauss-Legendre
void GaussMethod(double a ,double b,double n)
{
    double step = (b - a) / n;
    double sum = 0.0;
    int i;
    double temp_step, x1, x2;
    
    temp_step = a + step / 2;
    //x12 = 
    x1 = temp_step - (step / 2) * 0.5773502692;
    x2 = temp_step + (step / 2) * 0.5773502692;
    for(i = 1; i <= n; i++)
    {
        sum += function(x1) + function(x2);
        temp_step += step;
        x1 = temp_step - (step / 2) * 0.5773502692;
        x2 = temp_step + (step / 2) * 0.5773502692;
    }
    sum *= (step / 2);
    cout << "Gauss result: " << sum << "\n";
}

int main()
{
    double beginInterval;
    double endInterval;
    double epsilon = 0.001;
    int n;


    cout << "Beginning of interval: ";
    cin >> beginInterval;

    cout << "End of interval: ";
    cin >> endInterval;
    
    cout << "Number of iterations: ";
    cin >> n;

    if(n == 0)
    {
        cout << "0 intervals\n ";
        return 1;
    }
    else
    {
        MidpointRule(beginInterval, endInterval, n);
        TrapezoidalRule(beginInterval, endInterval, n);
        SimpsonsRule(beginInterval, endInterval, n);
        BodeFormula(beginInterval, endInterval, epsilon);
        GaussMethod(beginInterval, endInterval, n);
        return 0;
    }
}
    