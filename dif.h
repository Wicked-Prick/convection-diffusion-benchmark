#ifndef DIF_H
#define DIF_H
#include "var.h"

using namespace Var;

namespace Diff
{
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
}

#endif