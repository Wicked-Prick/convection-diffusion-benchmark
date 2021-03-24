#ifndef VAR_H
#define VAR_H
#include <vector>

namespace Var
{
        const int L1 = 24 , M1 = 24, Lmax = 24;
        const double eps = 1.e-5, alltime = 10.;
        static int L2, M2;
        static double TIME, DT, DX, XL, XLR, delT, delT1, RMAX, YL, YLR, DENOM, AP0;
        static std::vector<double>  X(L1), XU(L1), Y(M1), YV(M1), XDIF(L1), YDIF(M1), 
                        XCV(L1), YCV(M1), PT(Lmax), QT(Lmax), RET(L1);
        static std::vector< std::vector <double> > U(L1, std::vector<double> (M1)), V(L1, std::vector<double> (M1)), T(L1, std::vector<double> (M1)),
                Ta(L1, std::vector<double> (M1)), T1(L1, std::vector<double> (M1)), T0(L1, std::vector<double> (M1)), 
                GAMI(L1, std::vector<double> (M1)), GAMJ(L1, std::vector<double> (M1)), CON(L1, std::vector<double> (M1)) , 
                APS(L1, std::vector<double> (M1)), AIP(L1, std::vector<double> (M1)), AIM(L1, std::vector<double> (M1)), 
                AJM(L1, std::vector<double> (M1)), AJP(L1, std::vector<double> (M1)), AP(L1, std::vector<double> (M1)),
                TEMPS(L1, std::vector<double> (M1)), RHO(L1, std::vector<double> (M1)), 
                B(L1, std::vector<double> (M1)), AIPP(L1, std::vector<double> (M1)) , AJPP(L1, std::vector<double> (M1)); 

}

void output();
void rt();
std::vector <double> fan( std::vector<double> const &, std::vector<double> const &, double); 
double dmax1(double, double);
void allfile();
void start();
void grid1(int,  double, double, std::vector <double> &, std::vector <double> &); 

#endif