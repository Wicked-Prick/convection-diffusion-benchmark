#ifndef _DIF_H_
#define _DIF_H_
#include "var.h"

using namespace var;

namespace dif
{
    class Dif
    {         
        protected:
            double Fe, Fw, Fs, Fn, De, Dw, Dn, Ds;
            double DENSe, DENSw, DENSn, DENSs;
        public:       
            virtual void coeff_of_conduct();
            virtual void neighbor_coeff() = 0;
            virtual void solve();
            virtual double alpha_param(double );
            void gamsor(); 
            Dif(){};  
    };

    class Power_law : virtual public Dif
    {
        private:
            double max_Fw, max_Fe, max_Fs, max_Fn;
            double acof_w, acof_e, acof_s, acof_n;
        public:
            void neighbor_coeff();
            double acof(double x, double y);
            Power_law(): Dif(){};
    };

    class QUICK : public virtual Dif
    {
        private: 
            double alpha_w, alpha_e, alpha_s, alpha_n; 
        public:
            void neighbor_coeff(); 
            QUICK() : Dif(){};
    };

    class TVD : public virtual Dif
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
            TVD() : Dif(){};
    };    
}

#endif