module ParamsFunctions

using DifferentialEquations
using Plots
using DataFrames
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

H_nl_prime(z) = H_nl(z)^2 - H_nonlocal(z) * H_nonlocal_prime_z(z) / (1 + z)

S(z) = alpha_s * (1 + z) ^ beta_s

Sprime_z(z) = alpha_s * beta_s * (1 + z) ^ (beta_s - 1)

Ssecond_z(z) = alpha_s * beta_s * (beta_s - 1) * (1 + z) ^ (beta_s - 2)

Sprime(z) = -H_nonlocal(z) * Sprime_z(z)

Ssecond(z) = H_nonlocal(z) * H_nonlocal_prime_z(z) * Sprime_z(z) + H_nonlocal(z)^2 * Ssecond_z(z)

beta(z) = - H_nonlocal(z) * Sprime_z(z) / (1 + S(z))

betaprime(z) = H_nonlocal_prime_z(z) * H_nonlocal(z) * Sprime_z(z) / (1 + S(z)) + H_nonlocal(z) ^ 2 * Ssecond_z(z) / (1 + S(z)) - beta(z) ^ 2

gothic_R(z, k) = 2 * k^2 - 3 * beta(z) * H_nl(z)

Omega_tilde_nl(z) = (OmegaBh2 / h^2 / (1 + S(0)) * (1 + z) ^ (3 * (1 + wB)) + OmegaCh2 / h^2 / (1 + S(0)) * (1 + z) ^(3 * (1 + wC))) / (OmegaBh2 / h^2 / (1 + S(0)) * (1 + z) ^ (3 * (1 + wB)) + OmegaCh2 / h^2 / (1 + S(0)) * (1 + z) ^(3 * (1 + wC)) + (1 - (OmegaBh2 + OmegaCh2) / h^2 / (1 + S(0))) * (1 + z) ^ (3 * (1 + wlambda)))

Omega_tilde_nl_prime(z) = H_nonlocal(z) * (Omega_tilde_nl(z)^2 * (1 - (OmegaBh2 + OmegaCh2) / h^2 / (1 + S(0))) / ((OmegaBh2 + OmegaCh2) / h^2 / (1 + S(0)))) * (-3) / (1 + z)^4

function initial_conditions_phi(D_100, D_prime_100, k)
    A1(k) = 2 * H_nl(100) + 0.5 * beta(100) + 2 * Sprime(100) / (1 + S(100)) - k^2 / 3 / H_nl(100)
    A2(k) = -2 * beta(100)^2 + 2 * Ssecond(100) / (1 + S(100)) - 2 * Sprime(100) / (1 + S(100)) * H_nl_prime(100) / H_nl(100) + k^2 * H_nl_prime(100) / 3 / H_nl(100)^2 - H_nl_prime(100) - 0.5 * betaprime(100)
    B1() = -0.5 * (H_nl_prime(100) * Omega_tilde_nl(100) + H_nl(100) * Omega_tilde_nl_prime(100))
    B2() = -0.5 * H_nl(100) * Omega_tilde_nl(100)
    C2(k) = k^2 / 3 / H_nl(100) + H_nl(100) + 0.5 * beta(100)
    D1() = 0.5 * H_nl(100) * Omega_tilde_nl(100)
    phi_100 = 1 / (A2(k) / A1(k) - C2(k)) * ((B1() / A1(k) - D1()) * D_100 + B2() / A1(k) * D_prime_100)
    phi_prime_100 = D1() * D_100 - C2(k) * phi_100
    return phi_100, phi_prime_100
end

function phi_solve(k)
    phi_100, phi_prime_100 = initial_conditions_phi(1 / 101, -1 / (101)^2 * (-H_nonlocal(100)), k)
    M1(z) = 3 * H_nl(z) + beta(z) + 2 * Sprime(z) / (1 + S(z))
    M2(z) = - 2 * beta(z)^2 + 2 * Ssecond(z) / (1 + S(z)) - 2 * Sprime(z) / (1 + S(z)) * H_nl_prime(z) / H_nl(z)
    u0 = [phi_100, phi_prime_100]
    zspan = (100.0, 0.0)
    function phi_ODE!(du, u, p, t)
        phi, phi_dot = u
        du[2] = dphi_doubledot = -(H_nonlocal_prime_z(t) - M1(t)) / H_nonlocal(t) * phi_dot - M2(t) / H_nonlocal(t)^2 * phi
        du[1] = dphi = phi_dot
    end
    prob = ODEProblem(phi_ODE!, u0, zspan)
    sol = solve(prob)
    return sol
end

function plot_phi(k)
    sol = phi_solve(k)
    plot(sol, layout = (2,1), xaxis = "z", label = ["ϕ, k = $k" "dϕ/dz, k = $k"])
end

function plot_phi!(k)
    sol = phi_solve(k)
    plot!(sol, layout = (2,1), xaxis = "z", label = ["ϕ, k = $k" "dϕ/dz, k = $k"])
end

function d_solve(k)
    sol = phi_solve(k)
    df = DataFrame(sol)
    z = df[!, 1]
    phi = df[!, 2]
    phi_dot = df[!, 3]
    d_sol = Float64[]
    num = length(z)
    for i in 1:num
        d = 2 / (H_nl(z[i]) * Omega_tilde_nl(z[i])) * (phi_dot[i] + (k^2 / 3 / H_nl(z[i]) + H_nl(z[i]) + 0.5 * beta(z[i])) * phi[i])
        push!(d_sol, d)
    end
    return z, d_sol
end

function plot_d(k)
    z_arr, d_arr = d_solve(k)
    plot(z_arr, d_arr, xaxis = "z", yaxis = "D", label = "k = $k")
end

function plot_d!(k)
    z_arr, d_arr = d_solve(k)
    plot!(z_arr, d_arr, xaxis = "z", yaxis = "D", label = "k = $k")
end

end