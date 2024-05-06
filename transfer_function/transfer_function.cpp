#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "iostream"
#include "stdlib.h"
#include "time.h"
#include <complex>
#include "cmath"
#include "cstdlib"

using namespace std;
// Subroutin : written in 1 August 2015 by Farbod Hassani; All of the transfer functions as the refrences suggest is added to this version.
// Refrences: 1. BBKS transfer function in appendix G : http://adsabs.harvard.edu/abs/1986ApJ...304...15B
//2. Eisenstein and Hu transfer function (Baryonic feature) : http://lanl.arxiv.org/abs/astro-ph/9709112
//3. Eisenstein and Hu transfer function( Neutrino+ other cosmological effect added) :http://lanl.arxiv.org/abs/astro-ph/9710252
//4. Smith et al. halo fit transfer function (low resolution fitting formula): http://arxiv.org/abs/astro-ph/0207664v2
// 5: Takahashi et al. revising the halo fit model (high resolution fitting formula)  http://arxiv.org/abs/1208.2701v2


//------------------------------ Constant Definition ------------------------------------------------

#define PI 3.141592

//------------------------------------Sub-Routines/ Functions-------------------------------------------------

//********************************************************************************************************


// ***************MATHEMATHICAL FUNCTIONS *******************

//--------Sub-Routine for squared ----------------

double square (double x )
{
    double y;
    y=x*x;
    
    return (y);
}

double Inverse (double x )
{
    double y;
    y=1.0/x;
    
    return (y);
}

//---------------------------------------------------

//--------Sub-Routine for cube ---------

double cube(double x )
{
    double y;
    y=x*x*x;
    
    return (y);
}

//********************************************************************************************************

// ************* Cosmological Parameters *******

double E(double om_m,double om_L,double om_r, double x)   // The hubble parameter defined for LCDM - model
{
    double EE;   //the Hubble parameter
    double om_k;
    om_k=1-om_m-om_L-om_r;
    EE=sqrt(om_m*cube(1.0+x)+om_L+om_r*pow(1.0+x,4.0)+om_k*pow(1.0+x,2.0));
    return(EE);
    
    
}


//********************************************************************************************************

//*****************Comoving Distance *****************

double  Chi(double om_m, double om_L,double om_r, double h, double x)   // The comoving distance 
// The output is Mpc/h (Note the h)

{
    double om_k;
    om_k=1-om_m-om_L-om_r;
    
    int    x_top;     // The # loops over redshift
    double dx;        // the differential in redshift
    double Chii;      // The comoving distance normalized to cH_0^-1
    
    
    x_top=5000;        // Set the # of loops over redshift
    dx=(x-0.0)/x_top; //Set the differential magnitude
    Chii=0.0;         //Reset the value
    x=0.0;
    
    
    for(int i=1; i<x_top; ++i)    //integration loop over redshift
        
    {
        
        Chii=Chii+0.5*(3000)*dx*((1.0/E(om_m,om_L,om_r,x))+(1.0/E(om_m,om_L,om_r, x+dx)));  //comoving integral, main term
        
        
        x=x+dx;
    }
    
    return(Chii);
}





//********************************************************************************************************
//********************************************************************************************************
//********************************************************************************************************
//********************************************************************************************************

double Growth(double z, double om_m, double  om_L, double om_r) //A4 equation in Eisenstein and Hu article  http://lanl.arxiv.org/abs/astro-ph/9709112

{
    
    double om_k;
    om_k=1-om_m-om_L-om_r;
    double Dz;
    double om_La,om_ma;
    om_La=om_L/square(E(om_m, om_L, om_r,   z));
    om_ma=(om_m*cube(1+z))/square(E(om_m, om_L, om_r, z));
    
    Dz=2.5*(1.0/(1+z))*om_ma*Inverse(pow(om_ma,4.0/7.0)-om_La+(1.0+om_ma/2.0)*(1.0+om_La/70.0));
    
    
    
    return(Dz);
}

// BBKS Transfer function according to  function in appendix G : http://adsabs.harvard.edu/abs/1986ApJ...304...15B

// Fitting for just CDM universe (G3 case 1 in the article)
double Tk_BBKS_CDM(double k,double om_r, double om_CDM, double  om_b , double  om_L, double h)

{
    
    
    // --------parameter definition ----
    
    double Tk;     // Transfer function (BBKS);    /
    double  q;   //q parameter in the article (G3)
    double om_k,om_m;
    om_m=om_b+om_CDM;
    om_k=1-om_m-om_L-om_r;
    //q=k/(om_CDM*square(h)); // In the case of three  flavors of relativistic neutrino
    q=k/(om_CDM*square(h)); // [Mpc]^(-1)
    Tk=log(1.0+2.34*q)*pow(1+3.89*q+pow(16.1*q,2.0)+pow(5.46*q,3.0)+pow(6.71*q,4.0),-0.25)/(2.34*q);
    
    return(Tk);
}

/* double Tk_BBKS_CDM_h(double k,double om_r, double om_m, double om_CDM, double  om_b , double  om_L, double h) // K is calculated in Mpc^(-1)h

{
    double Tk;
    Tk=Tk_BBKS_CDM( k*h, om_r,  om_m,  om_CDM,   om_b ,   om_L,  h);
    return(Tk);
}  */

// BBKS Power spectrum as is defined in http://lanl.arxiv.org/abs/astro-ph/9709112according in A1
// Fitting for just CDM universe (G3 case 1 in the article)
double Pk_BBKS_CDM(double k,double om_r, double om_CDM, double  om_b , double  om_L ,double A, double h, double ns ,double z)

{
    
    
    // --------parameter definition ----
    
    double Pk;     // Transfer function (BBKS);    /
    double om_k;
    double om_m;
    om_m=om_b+om_CDM;
    om_k=1-om_m-om_L-om_r;
   // A=1.94*pow(10,-5)*pow(om_m,-0.785-0.05*log(om_m) )*exp(-0.95*(ns-1)-0.169*square(ns-1));
   

//** Pk=(2*PI*PI)*square(A)*pow(3000*Inverse(h),ns+3)*pow(k,ns)*square(Growth( z,  om_m,   om_L,  om_r))*square(Tk_BBKS_CDM( k, om_r,   om_CDM,   om_b ,   om_L,   h)); //*pow(k/0.05, 0.5 * alpha * log(k/0.05));     //LINEAR Power-Spectrum defined with BBKS transfer function
  

  //Pk=A*pow(k,ns)*square(Growth( z,  om_m,   om_L,  om_r))*square(Tk_BBKS_CDM( k, om_r,  om_m,  om_CDM,   om_b ,   om_L,   h));

 Pk=A*pow(k,ns)*square(Growth( z,  om_m,   om_L,  om_r))*square(Tk_BBKS_CDM( k, om_r,   om_CDM,   om_b ,   om_L,   h));



    return(Pk);


}




//Fitting in the presence of Baryon
double Tk_BBKS_Baryon(double k,double om_r, double om_CDM, double  om_b , double  om_L , double h)

{
    
    // --------parameter definition ----
    
    double Tk;     // Transfer function (BBKS);    /
    double  q;   //q parameter in the article (G3)
    double Rj; // Dark matter density
    double om_k;
    double om_m;
    om_m=om_b+om_CDM;
    om_k=1-om_m-om_L-om_r;
    Rj=1.6*pow(om_CDM*square(h),-0.5)*pow(10,-3);// Because Rj is the scale which baryon supress the CDM power and it is in Kpc not Mpc, so I have multiplied by 10^(-3)
    
    
    Tk=Tk_BBKS_CDM( k, om_r,   om_CDM,  om_b ,   om_L  , h )*Inverse(1+square(k*Rj)/2);
    
    return(Tk);
}


// BBKS Power spectrum as is defined in http://lanl.arxiv.org/abs/astro-ph/9709112according in A1
// Fitting for just CDM universe (G3 case 1 in the article)
double Pk_BBKS_Baryon(double k,double om_r, double om_CDM, double  om_b , double  om_L ,double A, double h, double ns ,double z)

{
    
    
    // --------parameter definition ----
    
    double Pk;     // Transfer function (BBKS);    /
    double om_k;
    double om_m;
    om_m=om_b+om_CDM;
    om_k=1-om_m-om_L-om_r;
    //A=1.94*pow(10,-5)*pow(om_m,-0.785-0.05*log(om_m) )*exp(-0.95*(ns-1)-0.169*square(ns-1));
    Pk=(2*PI*PI)*square(A)*pow(3000*Inverse(h),ns+3)*pow(k,ns)*square(Growth( z,  om_m,   om_L,  om_r))*square(Tk_BBKS_Baryon( k, om_r,    om_CDM,   om_b ,   om_L ,  h)); //*pow(k/0.05, 0.5 * alpha * log(k/0.05));     //LINEAR Power-Spectrum defined with BBKS transfer function
    return(Pk);
}






               //***************************************** EH transfer function*****

               
               double Tk_EH_full(double k,double om_r , double om_CDM, double  om_b , double  om_L,  double h,  double CMB_temp)  // CMB_temp is the temprature of CMB
               
    {
        
        
        //*****************************************
        // --------parameter definition ----
        double om_m;
        om_m=om_b+om_CDM;
        double z_drag, z_drag_b1, z_drag_b2, R, theta_cmb,z_eq,k_eq, om_bh,om_mh,om_CDMh, sound_horizon_fit,T_EH, T_b, T_c ,T;
        theta_cmb=(CMB_temp)/(2.7);
        om_mh=om_m*square(h);
        om_bh=om_b*square(h);
        om_CDMh=om_CDM*square(h);
        z_eq=2.5*pow(10,4)*om_mh*pow(theta_cmb,-4);// ~
        k_eq= 0.0746*om_mh/square(theta_cmb);  // ~
        z_drag_b1=0.313*pow(om_mh,-0.419)*(1+0.607*pow(om_mh,0.674));       //fitting parameter / ~
        
        z_drag_b2= 0.238*pow(om_mh,0.223);       //fitting parameter / ~
        
        z_drag = 1291*pow( om_mh,0.251)*(1.0+z_drag_b1*pow( om_bh,z_drag_b2))/(1.0+0.659*pow( om_mh,0.828));         // the redshift for drag epoch ~
        
        sound_horizon_fit = 44.5*log(9.83/om_mh)/(pow(1.0+10.0*pow(om_bh,0.75),0.5)); //  ~fitting formula for sound horizon in the Mpc unit, equation 26 in the ref.
        double om_k;
        om_k=1-om_m-om_L-om_r;
        
        //*****************************************

        // CDM part T_c in Eq.16
        double qq, alpha_c, a_1, a_2, beta_c, b_1, b_2, ff, CC_1,CC_2, T0_sup_alpha,T0_sup;

        qq=k/(13.41*k_eq);  //  ~
        a_1=pow(46.9*om_mh,0.670)*(1+pow(32.1*om_mh,-0.532));
        
        a_2=pow(12.0*om_mh,0.424)*(1+pow(45.0*om_mh,-0.582));
        
        alpha_c=pow(a_1,-(om_b/om_m))*pow(a_2,-pow((om_b/om_m),3));
        
        b_1=0.944*Inverse(1+pow(458*om_mh,-0.708));
        
        b_2=pow(0.395*om_mh,-0.0266);
        
        beta_c=Inverse(1+b_1*(pow(om_CDM/om_m,b_2)-1)); //
        
        
        ff=Inverse(1+pow(k*sound_horizon_fit/5.4,4)); //
        CC_1 = (14.2/alpha_c)+386/(1+69.9*pow(qq,1.08)); //For alpha_c
        CC_2 = 14.2+386/(1+69.9*pow(qq,1.08)); // For when we want alpha_c=1
        
        T0_sup_alpha = log(2.71828+1.8*beta_c*qq)/(log(2.71828+1.8*beta_c*qq)+CC_1*square(qq));  //in eq. 17 in the referance T_0 for alpha_c
        
        T0_sup = log(2.71828+1.8*beta_c*qq)/(log(2.71828+1.8*beta_c*qq)+CC_2*square(qq));  //in eq. 17 in the referance T_0 for alpha_c=1
        T_c=ff*T0_sup+(1-ff)*T0_sup_alpha;  //

        //*****************************************
        
        //Baryon part T_b in Eq. 16 of the reference
        double T0_sup_baryon, k_silk, alpha_b,beta_node,s_tilde, beta_b, R_drag, GG, bessel_0,yy ;
        k_silk=1.6*pow(om_bh,0.52)*pow(om_mh,0.73)*(1+pow(10.4*om_mh,-0.95));
        
        R_drag=31.5*om_bh*pow(theta_cmb,-4)*1000/z_drag;
        yy=(1+z_eq)/(1+z_drag);
        
        GG=yy*(-6*sqrt(1+yy)+(2+3*yy)*log((sqrt(1+yy)+1)/(sqrt(1+yy)-1))); //
        
        T0_sup_baryon= log(2.71828+1.8*qq)/(log(2.71828+1.8*qq)+CC_2*square(qq));
        
        alpha_b=2.07*k_eq*sound_horizon_fit*pow(1+R_drag,-0.75)*GG;
        
        beta_b= 0.5+(om_b/om_m)+(3-2*(om_b/om_m))*(sqrt(1+square(17.2*om_mh)));
        
        beta_node=8.41*pow(om_mh,0.435);
        
        s_tilde=sound_horizon_fit/(pow(1+pow(beta_node/(k*sound_horizon_fit),3),1/3)); //
        
        bessel_0=sin(k*s_tilde)/(k*s_tilde);
        //T_b=(T0_sup_baryon/(1+square(k*sound_horizon_fit))+alpha_b*exp(-pow(k/k_silk,1.4))/(1+pow(beta_b/(k*sound_horizon_fit),3)))*bessel_0;
        T_b=bessel_0*((T0_sup_baryon/(1+square(k*sound_horizon_fit/5.2)))+alpha_b*exp(-pow(k/k_silk,1.4))/(1+pow(beta_b/(k*sound_horizon_fit),3)));
        // The Baryon + CMD effects. Complete fitting formula.
        T=(om_b/om_m)*T_b+(om_CDM/om_m)*T_c;
        
        
        
        return(T);
    }

//************ EH Power spectrum ********************
double Pk_EH(double k,double om_r,  double om_CDM, double  om_b , double  om_L,double A,double CMB_temp, double h, double ns ,double z)

{
    
    
    // --------parameter definition ----
    
    double Pk;     // Transfer function (BBKS);    /
    double om_k;
    double om_m;
    om_m=om_b+om_CDM;
    om_k=1-om_m-om_L-om_r;
    //A=1.94*pow(10,-5)*pow(om_m,-0.785-0.05*log(om_m) )*exp(-0.95*(ns-1)-0.169*square(ns-1));
    //Pk=(2*PI*PI)*square(A)*pow(3000*Inverse(h),ns+3)*pow(k,ns)*square(Growth( z,  om_m,   om_L,  om_r))*square(Tk_EH_full( k, om_r,  om_m,  om_CDM,   om_b ,   om_L ,  h, CMB_temp)); //*pow(k/0.05, 0.5 * alpha * log(k/0.05));     //LINEAR Power-Spectrum defined with BBKS transfer function
 
 // **  Pk=(2*PI*PI)*square(A)*pow(3000*Inverse(h),ns+3)*pow(k,ns)*square(Growth( z,  om_m,   om_L,  om_r))*square(Tk_EH_full( k, om_r,    om_CDM,   om_b ,   om_L ,  h, CMB_temp)); //*pow(k/0.05, 0.5 * alpha * log(k/0.05));     //LINEAR Power-Spectrum defined with BBKS transfer function

  Pk=A*pow(k,ns)*square(Growth( z,  om_m,   om_L,  om_r))*square(Tk_EH_full( k, om_r,    om_CDM,   om_b ,   om_L ,  h, CMB_temp));   
    
    return(Pk);
}



// Halo model transfer function according to Smith et al. halo fit transfer function (low resolution fitting formula): http://arxiv.org/abs/astro-ph/0207664v2
// ****** NON-LINEAR Power Spectrum(HALO MODEL) ****************************
// Based on halo-model fit by Smith, R.E., Peacock, J.A., Jerkins, A. et al. Mon. Not. R. Astron. Soc. 341, 1311 (2003)

 
double P_NL_halomodel(double k, double Pk_linear, double Omega_m, double ksigma,double CC,  double neff )  // ksigma is the wavenumber that matter variance is equal to unity, there is a function to calculate this, neff=-1.550 Effective spectral index in LCDM  // CC=0.384; Effective second derivative of power-spectrum
// n is defined as a power of matter power spectrum as Eq 1 in reference.
 {
 
 // PARAMETERS FROM Smith et al. paper MNRAS, 2002 for halo model
 
 double P_NL_Hprime;  // Nonlinear Halo Model -parameters
 double P_NL_Q;    // Nonlinear matter-spectrum in quasi linear regime Q-term
 double P_NL_H;     // Nonlinear matter-power spectrum in halo-scale
 double P_NL_halomodel; //the final nonlinear halo model power spectrum,
 double Delta_L;   //the dimensionless power spectrum of matter
double Del_NL_Hprime;  // Nonlinear Halo Model -parameters
double Del_NL_Q;    // Nonlinear matter-spectrum in quasi linear regime Q-term
double Del_NL_H;     // Nonlinear matter-power spectrum in halo-scale
double Del_NL_halomodel; //the final nonlinear halo model power spectrum,
     
     
     
 double y;      //  parameter defined as y=k/ksigma
 
 
 double alphan;  //free parameter defined in halo-model
 double betan;  //free parameter defined in halo-model
 double gamman;  //free parameter defined in halo-model
 double mun, nun;  //free parameters in Halo-fit model
 double an,bn,cn;   //free parameters in Halo-fit model
 double f1,f2,f3;   //free parameters in Halo-fit model
 
 
 //free parameters based on NBODY Simulation
 
     alphan=1.3884+0.3700*neff-0.1452*square(neff);
     betan =0.8291+0.9854*neff+0.3401*square(neff);
     gamman=0.8649+0.2989*neff+0.1631*CC;
     
     an=pow(10,+1.4861+1.8369*neff+1.6762*square(neff)+0.7940*cube(neff)+0.1670*pow(neff,4.0)-0.6206*CC);
     bn=pow(10,+0.9463+0.9466*neff+0.3084*square(neff)-0.9400*CC);
     cn=pow(10,-0.2807+0.6669*neff+0.3214*square(neff)-0.0793*CC);
     
     mun=pow(10,-3.5442+0.1908*neff);
     nun=pow(10,+0.9589+1.2857*neff);
 

 //parametrs for LCDM
 
 f1=pow(Omega_m,-0.0307);
 f2=pow(Omega_m,-0.0585);
 f3=pow(Omega_m,+0.0743);
 
 
 y=k/ksigma;   //normalized wavenumber to the wavenumber of nonlinearity
 
 Delta_L=k*k*k*Pk_linear/(2.0*3.1415*3.1415);  //Dimensionless Power-Spectrum of matter with EH transfer function
 
 Del_NL_Hprime=an*pow(y,3*f1)/(1+bn*pow(y,f2)+pow(cn*f3*y,3-gamman));  //Auxilary parameter for Halo-fit model
 
 Del_NL_Q=Delta_L*(pow(1+Delta_L,betan)/(1+alphan*Delta_L))*exp(-(y/4.0+y*y/8.0));   // Quasi-Linear Power-Spectrum
 // P_NL_H=(2*3.1415*3.1415/(k*k*k))*P_NL_Hprime/(1.0+mun/y+nun/(y*y));            //Halo Model Power spectrum -SB
 Del_NL_H=Del_NL_Hprime/(1.0+mun/y+nun/(y*y));    //Halo Model power spectrum-FH
 Del_NL_halomodel=Del_NL_H+Del_NL_Q;   //The final halo model
     P_NL_halomodel=(2*3.1415*3.1415/(k*k*k))*Del_NL_halomodel;
 
 return(P_NL_halomodel);
 
 }

//************************************* Halo fit-high resolution- Takahashi et al.// // Halo fit Power spectrum (high resolution) according to Takahashi et al. revising the halo fit model (high resolution fitting formula)  http://arxiv.org/abs/1208.2701v2

double P_NL_halomodel_Takahashi(double k, double Pk_linear, double Omega_m,double Omega_w, double ksigma,double CC,double neff, double w_DE )  // ksigma is the wavenumber that matter variance is equal to unity, there is a function to calculate this, neff=-1.550 Effective spectral index in LCDM  // CC=0.384; Effective second derivative of power-spectrum
// n is defined as a power of matter power spectrum as Eq 1 in reference. w_DE is w for Dark energy which is well fitted by -1. Omega_w, Dark energy density
{
    
    // PARAMETERS FROM Smith et al. paper MNRAS, 2002 for halo model
    
    double P_NL_Hprime;  // Nonlinear Halo Model -parameters
    double P_NL_Q;    // Nonlinear matter-spectrum in quasi linear regime Q-term
    double P_NL_H;     // Nonlinear matter-power spectrum in halo-scale
    double P_NL_halomodel; //the final nonlinear halo model power spectrum,
    double Delta_L;   //the dimensionless power spectrum of matter
    double Del_NL_Hprime;  // Nonlinear Halo Model -parameters
    double Del_NL_Q;    // Nonlinear matter-spectrum in quasi linear regime Q-term
    double Del_NL_H;     // Nonlinear matter-power spectrum in halo-scale
    double Del_NL_halomodel; //the final nonlinear halo model power spectrum,
    
    
    
    double y;      //  parameter defined as y=k/ksigma
    
    
    double alphan;  //free parameter defined in halo-model
    double betan;  //free parameter defined in halo-model
    double gamman;  //free parameter defined in halo-model
    double mun, nun;  //free parameters in Halo-fit model
    double an,bn,cn;   //free parameters in Halo-fit model
    double f1,f2,f3;   //free parameters in Halo-fit model
    
    
    //free parameters based on NBODY Simulation
    
    alphan=abs(6.0835+1.3373*neff-0.1959*square(neff)-5.5274*CC);
    betan =2.0379-0.7354*neff+0.3157*square(neff)+1.2490*pow(neff,3)+0.3980*pow(neff,4)-0.1682*CC;
    gamman=0.1971-0.0843*neff+0.8460*CC;
    
    an=pow(10,1.5222+2.8553*neff+2.3706*square(neff)+0.9903*cube(neff)+0.2250*pow(neff,4.0)-0.6038*CC)+0.1749*Omega_w*(1+w_DE); // Omega_w can be as a function of redshift
    bn=pow(10,-0.5642+0.5864*neff+0.5716*square(neff)-1.5474*CC+0.2279*Omega_w*(1+w_DE));
    cn=pow(10,0.3698+2.0404*neff+0.8161*square(neff)+0.5869*CC);
    
    mun=0;
    nun=pow(10,5.2105+3.6902*neff);
    
    
    //parametrs for LCDM
    
    f1=pow(Omega_m,-0.0307);
    f2=pow(Omega_m,-0.0585);
    f3=pow(Omega_m,+0.0743);  //where Omega_m( Ωm ) is the matter density parameter at redshift z.
    
    
  /*  y=k/ksigma;   //normalized wavenumber to the wavenumber of nonlinearity
    
    Delta_L=k*k*k*Pk_linear/(2.0*3.1415*3.1415);  //Dimensionless Power-Spectrum of matter with EH transfer function
    
    Del_NL_Hprime=an*pow(y,3*f1)/(1+bn*pow(y,f2)+pow(cn*f3*y,3-gamman));  //Auxilary parameter for Halo-fit model
    
    Del_NL_Q=Delta_L*(pow(1+Delta_L,betan)/(1+alphan*Delta_L))*exp(-(y/4.0+y*y/8.0));   // Quasi-Linear Power-Spectrum
    // P_NL_H=(2*3.1415*3.1415/(k*k*k))*P_NL_Hprime/(1.0+mun/y+nun/(y*y));            //Halo Model Power spectrum -SB
    Del_NL_H=Del_NL_Hprime/(1.0+mun/y+nun/(y*y));    //Halo Model power spectrum-FH
    Del_NL_halomodel=Del_NL_H+Del_NL_Q;   //The final halo model
    P_NL_halomodel=(2*3.1415*3.1415/(k*k*k))*Del_NL_halomodel; */
    
    y=k/ksigma;   //normalized wavenumber to the wavenumber of nonlinearity
    
    Delta_L=k*k*k*Pk_linear/(2.0*3.1415*3.1415);  //Dimensionless Power-Spectrum of matter with EH transfer function
    
    P_NL_Hprime=an*pow(y,3*f1)/(1+bn*pow(y,f2)+pow(cn*f3*y,3-gamman));  //Auxilary parameter for Halo-fit model
    
    P_NL_Q=Pk_linear*(pow(1+Delta_L,betan)/(1+alphan*Delta_L))*exp(-(y/4.0+y*y/8.0));   // Quasi-Linear Power-Spectrum
    P_NL_H=(2*3.1415*3.1415/(k*k*k))*P_NL_Hprime/(1.0+mun/y+nun/(y*y));            //Halo Model Power spectrum
    P_NL_halomodel=P_NL_H+P_NL_Q;   //The final halo modela
    
    return(P_NL_halomodel);
    
}





double Sigmaa(double R1,double Omega_r, double Omega_cdm, double  Omega_b , double  Omega_L ,double A, double h, double ns, double z)
{
    double kk; //wavenumbber
    double k_ini;  //initial and final value of k
    double k_fin; //initial and final value of k
    double dk;
    double k_top;
    
    double Omega_m;   // density parameter of matter(cold dark matter + baryonic matter)
    Omega_m=Omega_cdm+Omega_b;
    
    k_ini=1e-10;
    k_fin=1e10;
    k_top=8e5;
    dk=(log10(k_fin)-log10(k_ini))/k_top;
    kk=k_ini;
    double sigmam;
    double Win;
    double x;
    sigmam=0.0;
    
    
    for(int i=1; i<k_top+1; ++i)
    {
        
        x=kk*R1;
        
        Win = 3.0*(sin(x) -x*cos(x))/square(x);
        
        
        
        // sigmam=sigmam+1*(kk*log(10)*dk)*Pk_BBKS_CDM(kk, Omega_r,Omega_cdm,  Omega_b , Omega_L , Apresent, h, ns ,z)*square(Win)*square(kk)/(2*square(PI));
      //  sigmam=sigmam+1*(kk*(pow(10,dk)-1))*Pk_BBKS_CDM(kk, Omega_r,Omega_cdm,  Omega_b , Omega_L , A, h, ns ,z)*square(Win)*square(kk)/(2*square(PI));

// Variance for BBKS
   //  sigmam=sigmam+1*(kk*log(10)*dk)*Pk_BBKS_CDM(kk, Omega_r,Omega_cdm,  Omega_b , Omega_L , A, h, ns ,z)*square(Win)*square(kk)/(2*square(PI));
   //     kk=kk*pow(10.0,dk);

// Variance for EH

sigmam=sigmam+1*(kk*(pow(10,dk)-1))*Pk_EH(kk, Omega_r,  Omega_cdm,  Omega_b ,  Omega_L,  A,  2.7, h,  ns , z)*square(Win)*square(kk)/(2*square(PI));
        kk=kk*pow(10.0,dk);

// Variance for TAKAHASHI



    }
    
    return(sigmam);
    
}





// *******************Linear-EH Power Spectrum ***************
// Effective power spectrum from Eisenstein and Hu article

double Pk_linEH(double k,double om_r ,double om_m, double  om_b , double  om_L , double h, double A, double ns,double z)

{
    
    
    // --------parameter definition ----
    
    double Tk;     // Transfer function (BBKS)
    double Pk;    // matter power spectrum, with BBKS transfer function
    double TEH;     // Eisnstein-Hu transfer function
    double PkEH;     // matter power spectrum, with Eisnstein-Hu transfer function
    double Delta_L; //linear dimensionless power-spectrum
    double D;  // growth fiunction
    
    
    double    ag,sg;             // Eisnstein-Hu transfer function parameters
    double   gamma,gamma1,gammaeff;   // Eisnstein-Hu transfer function parameters
    double   L0EH,C0EH;         // Eisnstein-Hu transfer function parameters
    
    
    D=1.0;   //set the growth function in present time (z=0.0)
    
    // ******TRANSFER FUNCTION
    
    
    //A) BBKS TRANSFER FUNCTION
    gamma=om_m*h*h;
    Tk=((log(1.0+2.34*k/gamma))/(2.34*k/gamma))*pow(1+3.89*k/gamma+pow(16.2*k/gamma,2.0)+pow(5.47*k/gamma,3.0)+pow(16.71*k/gamma,4.0),-0.25); //BBKS
    
    //B)EISNSTEIN-HU TRANSFER FUNCTION
    
    ag=1.0-0.328*log(431.0*om_m*h*h)*(om_b/om_m)+0.38*log(22.3*gamma)*(om_b/om_m)*(om_b/om_m);
    sg=(44.5*log(9.83/gamma))/(sqrt(1+10*pow(om_b*h*h,3.0/4.0)));
    gammaeff=gamma*(ag+(1.0-ag)/(1+pow(0.43*sg*k,4.0)));    // Introducing the effect of baryons
    gamma1=gammaeff;
    
    L0EH=log(2*exp(1.0)+1.8*k/gamma1);
    C0EH=14.2+731/(1.0+62.5*k/gamma1);
    
    
    TEH=L0EH/(L0EH+C0EH*(k/gamma)*(k/gamma));   //EH- Transfer function
    
   // Pk=  A*pow(k,ns)*square(Growth(z,  om_m,  om_L, om_r))*square(Tk);     //LINEAR Power-Spectrum defined with BBKS transfer function
    PkEH=(2*PI*PI)*square(A)*pow(3000*Inverse(h),ns+3)*pow(k,ns)*square(Growth(z,  om_m,  om_L, om_r))*square(TEH); //LINEAR Power-Spectrum defined with EH transfer function
    
    
    return(PkEH);
}

// Transfer function that is defined in Cosmology weinberg
double Pk_linW(double k, double om_r, double om_m, double  om_b , double  om_L , double h, double A, double ns, double z)

{
    // --------parameter definition ----
    
    double Tk;     // Transfer function (BBKS)
    double PkW;    // matter power spectrum, with BBKS transfer function
    double TkW;     // Eisnstein-Hu transfer function
    double Delta_L; //linear dimensionless power-spectrum
    double D;  // growth fiunction
    
    // ******TRANSFER FUNCTION
    
    double   kappa,alpha,s;   // Weinberg 2008 transfer function parameters
    
    s=44.5*log(9.83/(om_m*square(h)))/sqrt(1.0+10.0*pow(om_b*square(h),3.0/4.0));
    alpha=1-0.328*log(431*om_m*square(h))*(om_b/om_m)+0.38*log(22.3*om_m*square(h))*square(om_b/om_m);
    kappa=(k*sqrt(om_r))/(om_m*h/3000.0)*pow(alpha+(1-alpha)/(1.0+pow(0.43*k*s,4.0)),-1.0);
    
    TkW=((log(1.0+square(0.124*kappa)))/(square(0.124*kappa)))*(sqrt((1+square(1.257*kappa)+pow(0.4452*kappa,4.0)+pow(0.2197*kappa,6.0))/(1+square(1.606*kappa)+pow(0.8568*kappa,4.0)+pow(0.3927*kappa,6.0))));
    D=1.0;   //set the growth function in present time (z=0.0)
    
    
    
    
    
    
    PkW=(2*PI*PI)*square(A)*pow(3000*Inverse(h),ns+3)*pow(k,ns)*square(Growth(z, om_m,  om_L,om_r ))*square(TkW);     //LINEAR Power-Spectrum defined with BBKS transfer function
    
    return(PkW);
}




