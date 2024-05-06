module TransferFunction

#cpp equivalent functions for ease of transfering to julia
#---------------------------------------------------------#
function pow(a,b)
    return a^b
end

function square(a)
    return a^2
end

function Inverse(a)
    return 1/a
end
#---------------------------------------------------------#
function Tk_EH_full(k, om_r , om_CDM, om_b, om_L, h, CMB_temp)  # CMB_temp is the temprature of CMB           
    
    # --------parameter definition ----
    om_m=om_b+om_CDM;
    theta_cmb=(CMB_temp)/(2.7);
    om_mh=om_m*square(h);
    om_bh=om_b*square(h);
    om_CDMh=om_CDM*square(h);
    z_eq=2.5*pow(10,4)*om_mh*pow(theta_cmb,-4);# ~
    k_eq= 0.0746*om_mh/square(theta_cmb);  # ~
    z_drag_b1=0.313*pow(om_mh,-0.419)*(1+0.607*pow(om_mh,0.674));       #fitting parameter / ~
    
    z_drag_b2= 0.238*pow(om_mh,0.223);       #fitting parameter / ~
    
    z_drag = 1291*pow( om_mh,0.251)*(1.0+z_drag_b1*pow( om_bh,z_drag_b2))/(1.0+0.659*pow( om_mh,0.828));         # the redshift for drag epoch ~
    
    sound_horizon_fit = 44.5*log(9.83/om_mh)/(pow(1.0+10.0*pow(om_bh,0.75),0.5)); #  ~fitting formula for sound horizon in the Mpc unit, equation 26 in the ref.
    om_k=1-om_m-om_L-om_r;
    
    #*****************************************
    # CDM part T_c in Eq.16
    qq=k/(13.41*k_eq);  
    a_1=pow(46.9*om_mh,0.670)*(1+pow(32.1*om_mh,-0.532));
    
    a_2=pow(12.0*om_mh,0.424)*(1+pow(45.0*om_mh,-0.582));
    
    alpha_c=pow(a_1,-(om_b/om_m))*pow(a_2,-pow((om_b/om_m),3));
    
    b_1=0.944*Inverse(1+pow(458*om_mh,-0.708));
    
    b_2=pow(0.395*om_mh,-0.0266);
    
    beta_c=Inverse(1+b_1*(pow(om_CDM/om_m,b_2)-1)); 
    
    
    ff=Inverse(1+pow(k*sound_horizon_fit/5.4,4)); 
    CC_1 = (14.2/alpha_c)+386/(1+69.9*pow(qq,1.08)); #For alpha_c
    CC_2 = 14.2+386/(1+69.9*pow(qq,1.08)); # For when we want alpha_c=1
    
    T0_sup_alpha = log(2.71828+1.8*beta_c*qq)/(log(2.71828+1.8*beta_c*qq)+CC_1*square(qq));  #in eq. 17 in the referance T_0 for alpha_c
    
    T0_sup = log(2.71828+1.8*beta_c*qq)/(log(2.71828+1.8*beta_c*qq)+CC_2*square(qq));  #in eq. 17 in the referance T_0 for alpha_c=1
    T_c=ff*T0_sup+(1-ff)*T0_sup_alpha;  #
    #*****************************************
    
    #Baryon part T_b in Eq. 16 of the reference
    k_silk=1.6*pow(om_bh,0.52)*pow(om_mh,0.73)*(1+pow(10.4*om_mh,-0.95));
    
    R_drag=31.5*om_bh*pow(theta_cmb,-4)*1000/z_drag;
    yy=(1+z_eq)/(1+z_drag);
    
    GG=yy*(-6*sqrt(1+yy)+(2+3*yy)*log((sqrt(1+yy)+1)/(sqrt(1+yy)-1))); 
    
    T0_sup_baryon= log(2.71828+1.8*qq)/(log(2.71828+1.8*qq)+CC_2*square(qq));
    
    alpha_b=2.07*k_eq*sound_horizon_fit*pow(1+R_drag,-0.75)*GG;
    
    beta_b= 0.5+(om_b/om_m)+(3-2*(om_b/om_m))*(sqrt(1+square(17.2*om_mh)));
    
    beta_node=8.41*pow(om_mh,0.435);
    
    s_tilde=sound_horizon_fit/(pow(1+pow(beta_node/(k*sound_horizon_fit),3),1/3)); 
    
    bessel_0=sin(k*s_tilde)/(k*s_tilde);
    #T_b=(T0_sup_baryon/(1+square(k*sound_horizon_fit))+alpha_b*exp(-pow(k/k_silk,1.4))/(1+pow(beta_b/(k*sound_horizon_fit),3)))*bessel_0;
    T_b=bessel_0*((T0_sup_baryon/(1+square(k*sound_horizon_fit/5.2)))+alpha_b*exp(-pow(k/k_silk,1.4))/(1+pow(beta_b/(k*sound_horizon_fit),3)));
    # The Baryon + CMD effects. Complete fitting formula.
    T=(om_b/om_m)*T_b+(om_CDM/om_m)*T_c;
    
    return(T);
end
end