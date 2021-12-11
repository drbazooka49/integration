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
  //if (0 == n)
  //{
    //cout << "Trapezoidal rule result: " << sum << "\n";
  //}
 
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

int main()
{
    double beginInterval;
    double endInterval; 
    int n; 

    cout << "Beginning of interval: ";
    cin >> beginInterval;

    cout << "End of interval: ";
    cin >> endInterval;
    
    cout << "Number of iterations: ";
    cin >> n;

    MidpointRule(beginInterval, endInterval, n);
    TrapezoidalRule(beginInterval, endInterval, n);
    SimpsonsRule(beginInterval, endInterval, n);
    return 0;
}

/*
int main()
{
    int a, b, n, k, m = 1;
 
    cout << "a: "; //Начальная точка интеграла
    cin >> a;
    cout << "b: "; //Конечная точка интеграла
    cin >> b;
    cout << "n: "; //Точка деления интервала
    cin >> n;
    cout << "k: "; //Кол-во раз, которое следует вычислять интеграл
    cin >> k;
    cout << endl;
 
    float arg_x, sol; //arg_x это аргумент x, sol это шаг
 
    // внешний цикл, чтобы вычислить n, 2n, 4n, ... , 512n
    for (int s = 1; s <= k; s++)
    {
        // вычисление шага
        sol = (float)(b - a) / (float)(2 * n * m);
 
        float rez = 0;
        float pag1 = 0; //pag == временная переменная (temp)
        float pag2 = 0;
        float pag3 = 0;
        float j = (b - a) / sol; // используется, чтобы контролировать границу цикла и указатель на индекс последнего члена массива
        float* y = new float[j];
 
        for (int i = 0; i <= j; i++)
        {
            arg_x = (a + i * sol); // вычисление аргумента x
 
            if (arg_x < 1)
                y[i] = 1 / (arg_x * arg_x + 4 * arg_x + 5);
            else
                y[i] = (arg_x * arg_x) / 10;
 
            if (i == 0 || i == j) // pag1 хранит y[0] и сумму [j]
                pag1 += y[i];
            else if (i % 2 != 0)
                pag2 += 4 * y[i]; // pag2 вычисляет непарные y[i]
            else
                pag3 += 2 * y[i]; // pag3 вычисляет парные y[i]
        }
        // rez хранит значения интеграла, которые вычислены при помощи метода Симпсона
        rez = (pag1 + pag2 + pag3) * (sol / 3);
        cout << "Кол-во делений интервала = " << j << "; шаг = " << sol << "; интеграл = " << setprecision(10) << rez << endl;
        cout << endl;
        m *= 2;
    }
    _getch();
}

//метод Гаусса

double X[100],Y[100],dY[100],DY[100],d2Y[100],D2Y[100],pog1[100],pog2[100],h,I;
int i,a,b,k,r,v,m;
double f(double x)
{
double t;
t=sqrt(x)-pow(cos(x),2);
return(t);
}
double Gauss2(double A,double B,double M)
{
double h,x12,x1,x2,s=0,integral;
int i;
h=(B-A)/M;
x12=A+h/2;
x1=x12-(h/2)*0.5773502692;
x2=x12+(h/2)* 0.5773502692;
for(i=1;i<=m;i++)
{
s+=f(x1)+  f(x2);
x12=x12+h;
x1=x12-(h/2)* 0.5773502692;
x2=x12+(h/2)* 0.5773502692;
}
integral=(h/2)*s;
return (integral);
}
void main(void)
{
cout<<"a=";
cin>>a;
cout<<"b=";
cin>>b;
cout<<"Vibirete shag:\n"<<"1 h=0.2\n"<<"2 h=0.1\n"<<"3 h=0.05\n";
cin>>r;
switch(r)
{
case 1:
h=0.2;
break;
case 2:
h=0.1;
break;
case 3:
h=0.05;
break;
default:
cout<<"Vibrano ne dopustimoe zna4enie";
}
cout<<"Viberete m:\n1 m=10\n2 m=20\n3 m=40\n";
cin>>v;
switch(v)
{
case 1:
m=10;
break;
case 2:
m=20;
break;
case 3:
m=40;
break;
default:
cout<<"Vibrano ne dopustimoe zna4enie";
}
k=(b-a)/h;
for(i=1;i<=k+1;i++)
{
X[i]=a+(i-1)*h;
Y[i]=f(X[i]);
}
//первая производная
dY[1]=-(3*Y[1]-4*Y[2]+Y[3])/(2*h);
dY[k+1]=(Y[k-1]-4*Y[k]+3*Y[k+1])/(2*h);
for(i=2;i<=k;i++) dY[i]=(Y[i+1]-Y[i-1])/(2*h);
for(i=1;i<=k+1;i++) DY[i]=2*X[i]-5*sin(X[i]);
//вторая производная
for(i=1;i<=k+1;i++) D2Y[i]=2-5*cos(X[i]);
for(i=2;i<=k;i++) d2Y[i]=(Y[i-1]-2*Y[i]+Y[i+1])/(h*h);
for(i=1;i<=k+1;i++)
{
pog1[i]=fabs(fabs(dY[i])-fabs(DY[i]));
pog2[i]=fabs(fabs(d2Y[i])-fabs(D2Y[i]));
}
printf("X\tY\tdY\tDY\tpogr\td2Y\tD2Y\tpogr\n");
for(i=1;i<=k+1;i++) printf("%2.1f\t%4.3f\t%4.3f\t%4.3f\t%5.4f\t%4.3f\t%4.3f\t%5.4f\n",X[i],Y[i],dY[i],DY[i],pog1[i],d2Y[i],D2Y[i],pog2[i]);
I=Gauss2(a,b,m);
cout<<endl<<"Integral="<<I<<endl;
}

//формула Боде
using namespace std;
 
double Bode(double a, double b, double eps)
{
    double h, s, s1, z;
    int i, n;
 
    s = 0; n = 1;
    do
    {
        s1 = s; s = 0; n *= 2; h = (b - a) / n; z = h / 4;
 
        for (i = 1; i <= n; ++n)
            s += 7 * pow(a + i * h - 4 * z, 2) + 32 * pow(a + i * h - 3 * z, 2) + 12 * pow(A + i * h - 2 * z, 2) +
                32 * pow(a + i * h - z, 2) + 7 * pow(a + i * h, 2);
 
        s *= 2 * z / 45;
 
    } while (fabs(s1 - s) >= eps);
 
    return s;
}
 
int main()
{
    cout << Bode(1, 2, 0.001) << endl;
}
*/
