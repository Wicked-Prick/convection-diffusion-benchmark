#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <time.h>
#include <omp.h>
#include <fstream>

using namespace std;

const int L1 = 24 , M1 = 24, Lmax = 24;
const double eps = 1.e-5, alltime = 10.;
int L2, M2;
double TIME, DT, DX, XL, XLR, delT, delT1, RMAX, YL, YLR, DENOM, AP0;
vector<double>  X(L1), XU(L1), Y(M1), YV(M1), XDIF(L1), YDIF(M1), 
                XCV(L1), YCV(M1), PT(Lmax), QT(Lmax), RET(L1);

vector< vector <double> > U(L1, vector<double> (M1)), V(L1, vector<double> (M1)), T(L1, vector<double> (M1)),
        Ta(L1, vector<double> (M1)), T1(L1, vector<double> (M1)), T0(L1, vector<double> (M1)), 
        GAMI(L1, vector<double> (M1)), GAMJ(L1, vector<double> (M1)), CON(L1, vector<double> (M1)) , 
        APS(L1, vector<double> (M1)), AIP(L1, vector<double> (M1)), AIM(L1, vector<double> (M1)), 
        AJM(L1, vector<double> (M1)), AJP(L1, vector<double> (M1)), AP(L1, vector<double> (M1)),
        TEMPS(L1, vector<double> (M1)), RHO(L1, vector<double> (M1)), 
        B(L1, vector<double> (M1)), AIPP(L1, vector<double> (M1)) , AJPP(L1, vector<double> (M1)); 

void output();
void rt();
vector <double> fan( vector<double> const &, vector<double> const &, double); 
double dmax1(double, double);
void allfile();
void start();
void grid1(int,  double, double, vector <double> &, vector <double> &); 

class dif
{         
    protected:
        double Fe, Fw, Fs, Fn, De, Dw, Dn, Ds;
        double DENSe, DENSw, DENSn, DENSs;
    public:       
        virtual void coeff_of_conduct()
        {
            int i, j;
            //#pragma omp parallel for private(i, Dw, De, Ds, Dn, Fw, Fe, Fs, Fn) default(shared)
            for( j = 1; j < M2; ++j)
            {
                for( i = 1; i < L2; ++i)
                {
                    Dw = GAMI[i][j] * YCV[j] / XDIF[i];
                    De = GAMI[i + 1][j] * YCV[j] / XDIF[i + 1];
                    Ds = GAMJ[i][j] * XCV[i] / YDIF[j];
                    Dn = GAMJ[i][j + 1] * XCV[i] / YDIF[j + 1];

                    DENSw = RHO[i - 1][j] + (RHO[i][j] - RHO[i - 1][j]) * (XU[i] - X[i - 1]) / XDIF[i];
                    DENSe = RHO[i][j] + (RHO[i + 1][j] - RHO[i][j]) * (XU[i + 1] - X[i]) / XDIF[i + 1];
                    DENSs = RHO[i][j - 1] + (RHO[i][j] - RHO[i][j - 1]) * (YV[j] - Y[j - 1]) / YDIF[j];
                    DENSn = RHO[i][j] + (RHO[i][j + 1] - RHO[i][j]) * (YV[j + 1] - Y[j]) / YDIF[j + 1];

                    Fw = DENSw * U[i][j] * YCV[j];
                    Fe = DENSe * U[i][j] * YCV[j];
                    Fs = DENSs * V[i][j] * XCV[i];
                    Fn = DENSn * V[i][j] * XCV[i];         
                }
            } 
        }
        virtual void neighbor_coeff() = 0;
        virtual void solve()
        {
            int i, j;
            //---------------- X_dir ------------->
            //#pragma omp parallel for private(i)
            for(j = 1; j < M2; ++j)
            {
                for(i = 1;i < L2; ++i)
                {
                    TEMPS[i][j] = B[i][j] + AJM[i][j] * T[i][j - 1] + AJP[i][j] * T[i][j + 1];
                }
            }
            //#pragma omp parallel for private(i, DENOM, PT, QT) 
            for(j = 1; j < M2; ++j)
            {
                PT[0] = 0.;
                QT[0] = T[0][j];
                for(i = 1; i < L2; ++i)
                {
                    DENOM = AP[i][j] - PT[i - 1] * AIM[i][j];
                    PT[i] = AIP[i][j] / (DENOM + 1.e-30);
                    QT[i] = (TEMPS[i][j] + AIM[i][j] * QT[i - 1]) / (DENOM + 1.e-30);
                }
                for(i = L2 - 1; i > 1; --i)
                {
                    T[i][j] = T[i + 1][j] * PT[i] + QT[i];
                }
            }

            //----------------Y_dir ------------->
            //#pragma omp parallel for private(i)
            for(j = 1; j < M2; ++j)
            {
                for( i = 1;i < L2; ++i)
                {
                    TEMPS[i][j] = B[i][j] + AIM[i][j] * T[i - 1][j] + AIP[i][j] * T[i + 1][j];
                }
            }
            //#pragma omp parallel for private(j, DENOM, PT, QT) 
            for(i = 1; i < L2; ++i)
            {
                PT[0] = 0.;
                QT[0] = T[i][0];
                for( j = 1; j < M2; ++j)
                {
                    DENOM = AP[i][j] - PT[j - 1] * AJM[i][j];
                    PT[j] = AIP[i][j] / (DENOM + 1.e-30);
                    QT[j] = (TEMPS[i][j] + AJM[i][j] * QT[j - 1]) / (DENOM + 1.e-30);
                }
                for( j = M2 - 1; j > 1; --j)
                {
                    T[i][j] = T[i][j + 1] * PT[j] + QT[j];
                }
            }
        }   
        virtual double alpha_param(double f)
        {
            return f > 0. ? 1. : 0.;
        }  
        void gamsor();
      
       
};

class power_law : virtual public dif
{
    private:
        double max_Fw, max_Fe, max_Fs, max_Fn;
        double acof_w, acof_e, acof_s, acof_n;
    public:
        void neighbor_coeff();
        double acof(double x, double y);
};

class QUICK : public virtual dif
{
    private: 
        double alpha_w, alpha_e, alpha_s, alpha_n; 
    public:
        void neighbor_coeff();  
};

class TVD : public virtual dif
{
    private:
        double max_Fw, max_Fe, max_Fs, max_Fn;
        double grad_wp, grad_wm, grad_ep, grad_em,
               grad_sm, grad_sp, grad_nm, grad_np;
        double alpha_w, alpha_e, alpha_s, alpha_n;
        double psy_wp, psy_wm, psy_ep, psy_em,
               psy_sm, psy_sp, psy_nm, psy_np;         
    public:
        void neighbor_coeff(); 
        double limiter_param(double grad);
};

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
    for(int i = 2;i < n; ++i)
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

void dif::gamsor()
{
    for(int j = 1; j < M2; ++j)
    {
        for(int i = 1; i < L2; ++i)
        {
            CON[i][j] = 0.;
            APS[i][j] = 0.;
            RHO[i][j] = 1.;           
        }
    }
    for(int j = 1; j < M1; ++j)
    {
        for(int i = 1; i < L1; ++i)
        {
            GAMI[i][j] = 0.;
            GAMJ[i][j] = 0.; 
        }
    }
    /*for(int j = 1; j < M2; ++j)
    {
        for(int i = 1; i < L2; ++i)
        {
            CON[i][j] = - XCV[i] * YCV[j];
            APS[i][j] = 0.;
            RHO[i][j] = 1.;           
        }
    }

    for(int j = 0; j < M1; ++j)
    {
        T[1][j] = -(TIME + Y[j] * Y[j]);
    }
    
    for(int i = 0; i < L1; ++i)
    {
        T[i][1] = X[i] * X[i] - TIME;
    }

    for(int j = 1; j < M1; ++j)
    {
        for(int i = 1; i < L1; ++i)
        {
            GAMI[i][j] = 1.;
            GAMJ[i][j] = 1.;
            
        }
    }
   
    for(int i = 1; i < L1; ++i)
    {
        GAMJ[i][M1-1] = 0.;
    }

    for(int j = 1; j < M1; ++j)
    {
        GAMI[L1-1][j] = 0.;
    }
     
    for(int j = 1; j < M2; ++j)
    {
        CON[L2][j] = CON[L2][j] + 2. * YCV[j];
    } 

    for(int i = 1; i < L2; ++i)
    {
        CON[i][M2] = CON[i][M2] - 2. * XCV[i];
    } */
}

double power_law::acof(double x, double y)
{
    return (x == 0) ? 0 : dmax1(0., pow(1. - 0.1 * fabs(y / x), 5)); 
}

void power_law::neighbor_coeff()
{
    #pragma omp parallel for private(AP0, DT) default(shared)
    for(int j = 1; j < M2 ; ++j)
    {
        for(int i = 1; i < L2 ; ++i)
        {
            max_Fw = dmax1(0., Fw);
            max_Fe = dmax1(0., -Fe);
            max_Fs = dmax1(0., Fs);
            max_Fn = dmax1(0., -Fn);

            acof_w = acof(Dw, Fw);
            acof_e = acof(De, Fe);
            acof_s = acof(Ds, Fs);
            acof_n = acof(Dn, Fn);

            AIM[i][j] = Dw * acof_w + max_Fw;
            AIP[i][j] = De * acof_e + max_Fe;
            AJM[i][j] = Ds * acof_s + max_Fs;
            AJP[i][j] = Dn * acof_n + max_Fn;          
        }
    }
    #pragma omp parallel for  
    for(int j = 1; j < M2; ++j)
    {
        for(int i = 1; i < L2; ++i)
        {
            AP0 = RHO[i][j] * XCV[i] * YCV[j] / DT;
            AP[i][j] = - APS[i][j] + AIM[i][j] + AIP[i][j]
                        + AJM[i][j] + AJP[i][j] + AP0; 
            B[i][j] = CON[i][j] + AP0 * T0[i][j];             
        }
    }  

}

void QUICK::neighbor_coeff()
{
    for(int j = 1; j < M2 ; ++j)
    {
        for(int i = 1; i < L2 ; ++i)
        {
            alpha_w = alpha_param(Fw);
            alpha_e = alpha_param(Fe);
            alpha_s = alpha_param(Fs);
            alpha_n = alpha_param(Fn);

            AIM[i][j] = Dw + 6./8. * alpha_w * Fw + 1./8. * alpha_e*Fe +
                        3./8. * (1. - alpha_w) * Fw; 
            AIP[i][j] = De - 3./8. * alpha_e * Fe - 6./8. *(1. - alpha_e)*Fe -
                        1./8. * (1. - alpha_w)*Fw;
            AJM[i][j] = Ds + 6./8. * alpha_s * Fs + 1./8. * alpha_n*Fn +
                        3./8. * (1. - alpha_s) * Fs;
            AJP[i][j] = Dn - 3./8. * alpha_n * Fn - 6./8. *(1. - alpha_n)*Fn -
                        1./8. * (1. - alpha_s)*Fs;
            
            AIPP[i][j] = 1./8. * (1. - alpha_e)*Fe;
            AJPP[i][j] = 1./8. * (1. - alpha_n)*Fn;                                                  
        }
    }
    for(int j = 1; j < M2 ; ++j)
    {
        for(int i = 1; i < L2 ; ++i)
        {
            AP0 = RHO[i][j] * XCV[i] * YCV[j] / DT;
            AP[i][j] = - APS[i][j] + AIM[i][j] + AIP[i][j]
                        + AJM[i][j] + AJP[i][j] + AIPP[i][j] + AJPP[i][j] + AP0; 
            B[i][j] = CON[i][j] + AP0  * T0[i][j];             
        }
    }   
}

double TVD::limiter_param(double grad)
{
    double psy;
    
    if(grad <= 0.) psy = 0.;
    else if(grad > 0. && grad <= (1. / (2.+1.)) ) psy = (2. * grad) / 2.;
    else if(grad > (1./(2.+1.)) && grad <= 0.5) psy = 0.5 - grad / 2.;
    else if(grad > 0.5 && grad <= (2. / (2. + 1.)) ) psy = grad / 2.;
    else if(grad > (2. / (2. + 1.)) && grad <= 1.) psy = ((1.- grad) * 2.) / 2.;
    else if(grad >= 1.) psy = 0; 

    return psy;
}

void TVD::neighbor_coeff()
{

    for( int j = 1; j < M2 ; ++j)
    {
        for( int i = 1; i < L2 ; ++i)
        {
            max_Fw = dmax1(0., Fw);
            max_Fe = dmax1(0., -Fe);
            max_Fs = dmax1(0., Fs);
            max_Fn = dmax1(0., -Fn);
            
            alpha_w = alpha_param(Fw);
            alpha_e = alpha_param(Fe);
            alpha_s = alpha_param(Fs);
            alpha_n = alpha_param(Fn);
            
            AIM[i][j] = Dw + max_Fw;
            AIP[i][j] = De + max_Fe;
            AJM[i][j] = Ds + max_Fs;
            AJP[i][j] = Dn + max_Fn;

            /*if((i - 2) < 0) T[i][j] = 0.;
            else grad_wp = (T.at(i - 1).at(j) - T.at(i - 2).at(j)) / (T.at(i).at(j) - T.at(i - 1).at(j));
            if((i + 2) > L2) T[i][j] = 0.;
            else grad_em = (T.at(i + 2).at(j) - T.at(i + 1).at(j)) / (T.at(i + 1).at(j) - T.at(i).at(j));
            if((j - 2) < 0) T[i][j] = 0.;
            else grad_sp = (T.at(i).at(j - 1) - T.at(i).at(j - 2)) / (T.at(i).at(j) - T.at(i).at(j - 1));
            if((j + 2) > M2) T[i][j] = 0.;
            else grad_nm = (T[i][j + 2] - T[i][j + 1]) / (T[i][j + 1] - T[i][j]); */
            
            grad_nm = (T[i][j + 1] - T[i][j + 1]) / (T[i][j + 1] - T[i][j]); 
            grad_sp = (T.at(i).at(j - 1) - T.at(i).at(j - 1)) / (T.at(i).at(j) - T.at(i).at(j - 1));
            grad_em = (T.at(i + 1).at(j) - T.at(i + 1).at(j)) / (T.at(i + 1).at(j) - T.at(i).at(j));
            grad_wp = (T.at(i - 1).at(j) - T.at(i - 1).at(j)) / (T.at(i).at(j) - T.at(i - 1).at(j));
            grad_wm = (T.at(i + 1).at(j) - T.at(i).at(j)) / (T.at(i).at(j) - T.at(i - 1).at(j));
            grad_ep = (T.at(i).at(j) - T.at(i - 1).at(j)) / (T.at(i + 1).at(j) - T.at(i).at(j));
            grad_sm = (T[i][j - 1] - T[i][j]) / (T[i][j] - T[i][j - 1]);
            grad_np = (T[i][j] - T[i][j - 1]) / (T[i][j + 1] - T[i][j]);
           

            psy_wp = limiter_param(grad_wp);
            psy_wm = limiter_param(grad_wm);
            psy_ep = limiter_param(grad_ep);
            psy_em = limiter_param(grad_em);
            psy_sm = limiter_param(grad_sm);
            psy_sp = limiter_param(grad_sp);
            psy_nm = limiter_param(grad_nm);
            psy_np = limiter_param(grad_np);

            CON[i][j] = 0.5 * Fe * ((1. - alpha_e) * psy_em - alpha_e * psy_ep) * 
                        (T[i + 1][j] - T[i][j]) + 
                        0.5 * Fw * (alpha_w * psy_wp - (1. - alpha_w) * psy_wm)  *
                        (T[i][j] - T[i - 1][j]) +
                        0.5 * Fn * ((1. - alpha_n) * psy_nm - alpha_n * psy_np) * 
                        (T[i][j + 1] - T[i][j]) + 
                        0.5 * Fs * (alpha_s * psy_sp - (1. - alpha_s) * psy_sm)  *
                        (T[i][j] - T[i][j - 1]);

        }
    }
    
    for(int j = 1; j < M2; ++j)
    {
        for(int i = 1; i < L2; ++i)
        {
            AP0 = RHO[i][j] * XCV[i] * YCV[j] / DT;
            AP[i][j] = - APS[i][j] + AIM[i][j] + AIP[i][j]
                        + AJM[i][j] + AJP[i][j] + AP0; 
            B[i][j] = CON[i][j] + AP0 * T0[i][j];             
        }
    }  
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