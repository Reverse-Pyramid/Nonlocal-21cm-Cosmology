module ParamsFunctions

#---------------------parameters---------------------

#General Params
G = 6.67 * 10^(-11) 
Omega_m = 0.317
Omega_l = 1 - Omega_m

#ΛCDM_Params
H0_CDM = 67.4
h_CDM = 0.674
deltaH0_CDM = 0.5
OmegaBh2_CDM = 0.0224
deltaOmegeBh2_CDM = 0.0001
OmegaCh2_CDM = 0.12
deltaOmegaCh2_CDM = 0.001
wB_CDM = 0
wC_CDM = 0
wLambda_CDM = -1.03
deltawLambda_CDM = 0.03
A_COBE = 2.089 * 10^(-9)
n_s = 0.9667

#Nonlocal Params
alpha_s = 0.0034
beta_s = -1.0
H0 = 68.42
h = 0.6842
deltaH0 = 0.5
OmegaBh2 = 0.02217
deltaOmegeBh2 = 0.0001
OmegaCh2 = 0.1176
deltaOmegaCh2 = 0.001
wB = 0
w = 0 #CDM
wC = 0
wlambda = -1.0
dP_dRho = 0
deltawLambda = 0.03

#Other Params
h = 0.6842
OmegaB = 0.02217 / h^2
OmegaC = 0.1176 / h^2
Omega0 = 1.0
thetaCMB = 2.7
a1 = (46.9 * Omega0 * h^2)^0.670 * (1 + (32.1 * Omega0 * h^2)^(-0.532))
a2 = (12.0 * Omega0 * h^2)^0.424 * (1 + (45.0 * Omega0 * h^2)^(-0.582))
alpha_c = a1^(-OmegaB/Omega0) * a2^(-(OmegaB/Omega0)^3)
b1 = 0.944 / (1 + (458 * Omega0 * h^2)^(-0.708))
b2 = (0.395 * Omega0 * h^2)^(-0.0266)
beta_c = 1 / (1 + (b1 * (OmegaC / Omega0)^b2 - 1))
keq = 7.46 * 10^(-2) * Omega0 * h^2 * thetaCMB^(-2)

#---------------------functions---------------------

#Nonlocal theory functions
Gamma(z) = (1 + S(0)) / (1 + S(z))

Gamma_prime_z(z) = - Sprime_z(z) * (1 + S(0)) / (1 + S(z)) ^ 2

H_nonlocal(z) = H0 * Gamma(z) * sqrt(OmegaBh2 / h^2 / (1 + S(0)) * (1 + z) ^ (3 * (1 + wB)) + OmegaCh2 / h^2 / (1 + S(0)) * (1 + z) ^(3 * (1 + wC)) + (1 - (OmegaBh2 + OmegaCh2) / h^2 / (1 + S(0))) * (1 + z) ^ (3 * (1 + wlambda)))

H_nl(z) = 1 / (1 + z) * H_nonlocal(z)

H_nonlocal_prime_z(z) = H_nonlocal(z) * ((OmegaBh2 / h^2 / (1 + S(0)) * (3 * (1 + wB)) * (1 + z) ^ (3 * (1 + wB) - 1) + OmegaCh2 / h^2 / (1 + S(0)) * (3 * (1 + wC)) * (1 + z) ^(3 * (1 + wC) - 1) + (1 - (OmegaBh2 + OmegaCh2) / h^2 / (1 + S(0))) * (3 * (1 + wlambda)) * (1 + z) ^ (3 * (1 + wlambda) - 1)) / 2 / (OmegaBh2 / h^2 / (1 + S(0)) * (1 + z) ^ (3 * (1 + wB)) + OmegaCh2 / h^2 / (1 + S(0)) * (1 + z) ^(3 * (1 + wC)) + (1 - (OmegaBh2 + OmegaCh2) / h^2 / (1 + S(0))) * (1 + z) ^ (3 * (1 + wlambda))) + Gamma_prime_z(z) / Gamma(z))

S(z) = alpha_s * (1 + z) ^ beta_s

Sprime_z(z) = alpha_s * beta_s * (1 + z) ^ (beta_s - 1)

Ssecond_z(z) = alpha_s * beta_s * (beta_s - 1) * (1 + z) ^ (beta_s - 2)

beta(z) = - H_nonlocal(z) * Sprime_z(z) / (1 + S(z))

betaprime(z) = H_nonlocal_prime_z(z) * H_nonlocal(z) * Sprime_z(z) / (1 + S(z)) + H_nonlocal(z) ^ 2 * Ssecond_z(z) / (1 + S(z)) - beta(z) ^ 2

gothic_R(z, k) = 2 * k^2 - 3 * beta(z) * H_nl(z)

A1(z, k) = - (18 * (1 + w) * H_nl(z)^3 + 4 * k^2 * beta(z) + 3 * H_nl(z)^2 * (7 - 3 * w - 6 * dP_dRho) * beta(z) + 12 * k^2 * H_nl(z) * (w - dP_dRho) - 6 * H_nl(z) * betaprime(z)) / (2 * (1 + S(z)) * gothic_R(z, k))

A2(z, k) = -1 / (1 + z) * (-4 * k^4 + 54 * (1 + w) * H_nl(z)^4 + 12 * H_nl(z) * beta(z) * (2 * k^2 + 3 * H_nl(z)^2) + 9 * H_nl(z)^2 * (2 * k^2 * (1 + w) - beta(z)^2 - 2 * betaprime(z))) / (2 * k^2 * (1 + S(z)) * gothic_R(z, k) * (1 + 2)^(-1))

B1(z, k) = - (2 + 3 * w) * H_nl(z) - beta(z) + (3 * H_nl(z) * (2 * k^2 + 3 * (1 + w) * H_nl(z)^2)) / gothic_R(z, k)

B2(z, k) = k^2 * (-2 * k^2 * dP_dRho + 3 * (1 + w) * H_nl(z)^2 + 3 * (1 + dP_dRho) * beta(z) * H_nl(z)) / ((1 + w) * gothic_R(z, k))



end