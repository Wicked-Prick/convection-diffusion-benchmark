#include "dif.h"
#include "var.h"
#include <cmath>

using namespace Diff;
using namespace Var;

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
