#include <iostream>
#include <cstdio>
#include <cstdlib>
#include<fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

#define n 16
typedef struct ee{
     double the_error = 1000000;
     int index = 0;
} theE;
double p[n+1][n+1]={0},x[n+1],y[n+1],sy[n+1]={0},a[n+1]={0};
theE E[n+1];
double sigma_x(int m)   //m is degree n is number of points
{
    double sum = 0.0;
    int k;
    for(k = 1 ; k<=n ;k++)
    {
        sum+=pow(x[k],m);
    }
    return sum;
}

void makematrix_A()
{
    int i,j;
    for(i=1;i<=n;i++)
        for(j=1;j<=n;j++)
            p[i][j]=sigma_x(i+j-2);
}


double sigma_y(int i)
{
    double sum =0.0;
    int k;
    for(k=1;k<=n;k++)
    {
        sum+=(y[k]*pow(x[k],i-1));
    }
    return sum;
}

void makematrix_B()
{
    int i;
    for(i=1;i<=n;i++)
    {
        sy[i] = sigma_y(i);
    }

}

void calculate_a(int m)
{
    makematrix_A();
    makematrix_B();
    int i,j,k;
    double d;
    for(k = 1;k<m+1;k++){
        for(i = k+1; i<=m+1;i++)
        {
            d = p[i][k]/p[k][k];
            for(j = k;j <=m+1;j++)
            {
                p[i][j] = p[i][j] - p[k][j] * d;
            }
          sy[i] = sy[i] - sy[k]*d;
        }
    }

    for(i = m+1; i>=1;i--)
    {
        double sum = 0.0;
        for(j = i+1 ; j <= m+1;j ++)
        {
            sum += p[i][j]*a[j];
        }
         a[i] = (sy[i] - sum)/p[i][i];
    }

}
double calculate_fx(int m,int i)
{
    int k;
    double sum =0.0;

    for(k=1;k<=m+1;k++)
    {

        sum= sum+ a[k]*pow(x[i],k-1);
    }

    return sum;
}

double calculate_E(int m)
{
    int i;
    double sum=0.0;
    for(i=1;i<=n;i++)
    {
        sum += pow((y[i] - calculate_fx(m,i)),2);

    }
    return sum;
}
bool compare(const theE &x, const theE &y)
{
    return x.the_error < y.the_error;
}
int main()
{
    int i ,j ,m=1;
    freopen("HW3.txt","r",stdin);
    for(i = 1; i <= n; i ++)
    {
        cin>>x[i]>>y[i];
    }
    while(m<n)
    {
        calculate_a(m);
        E[m].the_error = calculate_E(m);
        E[m].index = m;
        m++;
    }
    sort(E+1,E+n,compare);
    cout<<" m  | "<<"error"<<endl;
    cout<<"-------------"<<endl;
    for(i=1;i<n;i++)
    {
        cout<<setw(2)<<E[i].index<<"  | "<<E[i].the_error<<endl;

    }
    printf("The best choice is %d. \nThe error of %d is %f",E[1].index,E[1].index,E[1].the_error);
    return 0;
}
