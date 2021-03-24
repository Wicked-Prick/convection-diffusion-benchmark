#include "var.h"
#include "dif.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <time.h>
#include <omp.h>
#include <fstream>

using namespace std;
using namespace Diff;
using namespace Var;

ostream & operator << (ostream & os, const vector<double> & v)
{
    for(auto i = 0; i < v.size(); ++i)
    {
        os << "vector["<< i <<"]" << "=" << v[i] << endl;
    }
    return os;
}

vector<double> operator *(const vector< double> &x,const  vector<double> &y)
{
    vector<double> res = x;
    res.resize(x.size());
    for (int i = 0; i < res.size() && i < y.size(); i++) {
        res[i] *= y[i];
    }
    
    return res;
}

vector<double> operator -(const vector< double> &x,const  vector<double> &y)
{
    vector<double> res = x;
    res.resize(x.size());
    for (int i = 0; i < res.size() && i < y.size(); i++) {
        res[i] -= y[i];
    }
    
    return res;
}

vector<double> operator -(const vector< double> &x, double y)
{
    vector<double> res = x;
    res.resize(x.size());
    for (int i = 0; i < res.size(); i++) {
        res[i] -= y;
    }
    
    return res;
}

double dmax1(double x, double y)
{
    return (x > y) ? x : y;
}

vector <double> fan(const vector <double>  &x, const  vector <double>  &y,  double z)
{
    vector<double> tmp(L1);
    tmp = x*x - y*y - z ;
    return tmp;
}

void start()
{
    XL = 0., XLR = 1.;
    YL = 0., YLR = 1.;
    TIME = 0., DT = 0.001;
    
    L2 = L1 - 1;
    M2 = M1 - 1;    
    
    grid1(L1, XL, XLR, XU, X);
    grid1(M1, YL, YLR, YV, Y);

    for(int i = 1; i < L1; ++i)
    {
        XDIF[i] = X[i] - X[i - 1];       
    }   
    
    for(int i = 1; i < L2; ++i)
    {
        XCV[i] = XU[i + 1] - XU[i];
    }
    
   
    for (int j = 1; j < M1; ++j)
    {
        YDIF[j] = Y[j] - Y[j - 1];
    }
  
    
    for (int j = 1; j < M2; ++j)
    {
        YCV[j] = YV[j + 1] - YV[j];
    }

    for(int j = 0; j < M1; ++j)
    {
        for(int i = 0; i < L1; ++i)
        {
            T0[i][j] = 0.;
            T[i][j] = T0[i][j];
        }
    }
    for(int j = 0; j < M1; ++j)
    {
        T[1][j] = 100.;
    }
    for(int i = 0; i < L1; ++i)
    {
        T[i][M1] = 100.;
    }

    for(int j = 0; j < M1; ++j)
    {
        for(int i = 0; i < L1; ++i)
        {
            Ta[i][j] = 100.;
        }
    }
    
    for(int j = 0; j < M1; ++j)
    {
        for(int i = 0; i < L1; ++i)
        {
            if(i > j) Ta[i][j] = 0.;
        }
    }

    for(int j = 0; j < M1; ++j)
    {
        for(int i = 0; i < L1; ++i)
        {
            V[i][j] = 2.;
            U[i][j] = 2.;
            
        }
    }

    /*for(int j = 0; j < M1; ++j)
    {
        for(int i = 0; i < L1; ++i)
        {
            T0[i][j] = X[i] * X[i] - Y[j] * Y[j];
            T[i][j] = T0[i][j];
        }
    }
    
    for(int j = 0; j < M1; ++j)
    {
        for(int i = 0; i < L1; ++i)
        {
            V[i][j] = X[i];
            U[i][j] = Y[j];
            
        }
    }
    
    for(int j = 0; j < M1; ++j)
    {
        T[L1 - 1][j] = T[L2][j] + 2. * (X[L1] - X[L2]);
    }
      
    for(int i = 0; i < L1; ++i)
    {
        T[i][M1 - 1] = T[i][M2] - 2. * (Y[M1] - Y[M2]);
    } */ 
}

void grid1(int n,  double a, double b, vector <double> &cv, vector <double> &node)
{
    DX = b / double(n - 2);
    cv[1] = a;
    
    for(int i = 2; i < n; ++i)
    {
        cv[i] = cv[i - 1] + DX;
    }
    
    node[0] = cv[1];

    for(int i = 1;i < n - 1; ++i)
    {
        node[i] = 0.5*(cv[i + 1] + cv[i]);
    }
    
    node[n] = cv[n];  
}



void output()
{
    delT = 0.;
    delT1 = 0.;
    int k = 0;
    for(int j = 0; j < M1; ++j)
    {
        for(int i = 0; i < L1; ++i)
        {
            //RET = fan(X,Y,TIME);   
            delT = dmax1(delT, fabs(T[i][j] -  Ta[i][j]));
            delT1 = dmax1(delT1, 100. * fabs((T[i][j] -  Ta[i][j]) /  Ta[i][j]));
        }
    }
    
    cout.setf(ios::scientific);
    cout <<"TIME \t\t " " T(1,1) \t\t " " Tan(1,1) \t\t " " delt \t\t " " delt1" << "\n" 
    << TIME << "\t"<< T[1][1] << " \t\t " <<  Ta[1][1] << "\t\t" << delT << "\t" << delT1 << endl;
}

void rt()
{
    RMAX = 0.;
    
    for(int j = 1; j < M2; ++j)
    {
        for(int i = 1; i < L2; ++i)
        {
            RMAX = dmax1(RMAX, fabs(1. - T[i][j] / T1[i][j]));
        }
    }
    std::cout << RMAX << endl;
}

void allfile()
{
    ofstream ost;
    ost.setf(ios::scientific);
    ost.open("all.dat");
    ost << "VARIABLES = X, Y, T, Ta, delT" << endl;
    ost << "zone i = " << L1 << "\t" << ", j = " << M1 << "\t" << ", f=point " << endl;
    
    /*for(int j = 1; j < M2; ++j)
    {
        for(int i = 1; i < L2; ++i)
        {  
            RET[i] = X[i]*X[i] - Y[j]*Y[j] - TIME;
            ost << X[i] << "\t" << Y[j] << "\t" << T[i][j] << "\t\t" << 
            RET[i]  << "\t\t" << fabs(T[i][j] - RET[i] ) << endl; 
        }
    }*/
    int j = M2;
    for(int i = 1; i < L2; ++i)
    { 
            ost << X[i] << "\t" << Y[j] << "\t" << T[i][j] << "\t\t" << 
            Ta[i][j]  << "\t\t" << fabs(T[i][j] - Ta[i][j] ) << endl; 
            --j;
    }
    ost.close();

}

int main()
{
    double ftime, etime, exec_time;
    power_law pw;
    QUICK qk;
    TVD tvd;
    start();
    ftime = omp_get_wtime(); 

    while(TIME <= (alltime - 0.5*DT))
    {
        TIME += DT;
        RMAX = 1.;
        
        while (RMAX > eps)
        {    
           T1 = T; 
           pw.coeff_of_conduct();
           pw.neighbor_coeff();
           pw.gamsor();
           pw.solve();
           rt();
        }
        output();

        T0 = T;       
    }
   
    etime = omp_get_wtime();
    exec_time = etime - ftime; 
    cout << "\n" << "All taken time " << exec_time << endl;
    allfile();
    
    return 0;
}