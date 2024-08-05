module ParamsFunctions
include("./TransferFunction.jl")
import .TransferFunction as TF
using DifferentialEquations
using Integrals
using Plots
using DataFrames
using CSV
#---------------------Config-------------------------

save_step_length = 0.01

#---------------------parameters---------------------

#General Params
G = 6.67 * 10^(-11)
c_kms = 3 * 10 ^ 5
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
#n_s = 0.9667
n_s = 0.96
CMB_temp = 2.726
Omega_r = 0

#Nonlocal Params
alpha_s = 0.0034
#alpha_s = 0
beta_s = -1.0
#H0 = 68.42
#h = 0.6842
H0 = 70
h = 0.7
deltaH0 = 0.5
#OmegaBh2 = 0.02217
OmegaBh2 = 0.0226
deltaOmegeBh2 = 0.0001
#OmegaCh2 = 0.1176
OmegaCh2 = 0.112
deltaOmegaCh2 = 0.001
wB = 0
w = 0 #CDM
wC = 0
wlambda = -1.0
dP_dRho = 0
deltawLambda = 0.03

#Other Params
OmegaB = OmegaBh2 / h^2
OmegaC = OmegaCh2 / h^2
Omega0 = 1.0
thetaCMB = 2.7
a1 = (46.9 * Omega0 * h^2)^0.670 * (1 + (32.1 * Omega0 * h^2)^(-0.532))
a2 = (12.0 * Omega0 * h^2)^0.424 * (1 + (45.0 * Omega0 * h^2)^(-0.582))
alpha_c = a1^(-OmegaB / Omega0) * a2^(-(OmegaB / Omega0)^3)
b1 = 0.944 / (1 + (458 * Omega0 * h^2)^(-0.708))
b2 = (0.395 * Omega0 * h^2)^(-0.0266)
beta_c = 1 / (1 + (b1 * (OmegaC / Omega0)^b2 - 1))
keq = 7.46 * 10^(-2) * Omega0 * h^2 * thetaCMB^(-2)

#---------------------functions---------------------

#Nonlocal theory functions
Gamma(z) = (1 + S(0)) / (1 + S(z))

Gamma_prime_z(z) = -Sprime_z(z) * (1 + S(0)) / (1 + S(z))^2

H_nonlocal(z) = H0 * Gamma(z) * sqrt(OmegaBh2 / h^2 / (1 + S(0)) * (1 + z)^(3 * (1 + wB)) + OmegaCh2 / h^2 / (1 + S(0)) * (1 + z)^(3 * (1 + wC)) + (1 - (OmegaBh2 + OmegaCh2) / h^2 / (1 + S(0))) * (1 + z)^(3 * (1 + wlambda)))

H_nl(z) = 1 / (1 + z) * H_nonlocal(z)

H_nonlocal_prime_z(z) = H_nonlocal(z) * ((OmegaBh2 / h^2 / (1 + S(0)) * (3 * (1 + wB)) * (1 + z)^(3 * (1 + wB) - 1) + OmegaCh2 / h^2 / (1 + S(0)) * (3 * (1 + wC)) * (1 + z)^(3 * (1 + wC) - 1) + (1 - (OmegaBh2 + OmegaCh2) / h^2 / (1 + S(0))) * (3 * (1 + wlambda)) * (1 + z)^(3 * (1 + wlambda) - 1)) / 2 / (OmegaBh2 / h^2 / (1 + S(0)) * (1 + z)^(3 * (1 + wB)) + OmegaCh2 / h^2 / (1 + S(0)) * (1 + z)^(3 * (1 + wC)) + (1 - (OmegaBh2 + OmegaCh2) / h^2 / (1 + S(0))) * (1 + z)^(3 * (1 + wlambda))) + Gamma_prime_z(z) / Gamma(z))

H_nl_prime(z) = H_nl(z)^2 - H_nonlocal(z) * H_nonlocal_prime_z(z) / (1 + z)

S(z) = alpha_s * (1 + z)^beta_s

Sprime_z(z) = alpha_s * beta_s * (1 + z)^(beta_s - 1)

Ssecond_z(z) = alpha_s * beta_s * (beta_s - 1) * (1 + z)^(beta_s - 2)

Sprime(z) = -H_nonlocal(z) * Sprime_z(z)

Ssecond(z) = H_nonlocal(z) * H_nonlocal_prime_z(z) * Sprime_z(z) + H_nonlocal(z)^2 * Ssecond_z(z)

beta(z) = -H_nonlocal(z) * Sprime_z(z) / (1 + S(z))

betaprime(z) = H_nonlocal_prime_z(z) * H_nonlocal(z) * Sprime_z(z) / (1 + S(z)) + H_nonlocal(z)^2 * Ssecond_z(z) / (1 + S(z)) - beta(z)^2

gothic_R(z, k) = 2 * k^2 - 3 * beta(z) * H_nl(z)

Omega_tilde_nl(z) = (OmegaBh2 / h^2 / (1 + S(0)) * (1 + z)^(3 * (1 + wB)) + OmegaCh2 / h^2 / (1 + S(0)) * (1 + z)^(3 * (1 + wC))) / (OmegaBh2 / h^2 / (1 + S(0)) * (1 + z)^(3 * (1 + wB)) + OmegaCh2 / h^2 / (1 + S(0)) * (1 + z)^(3 * (1 + wC)) + (1 - (OmegaBh2 + OmegaCh2) / h^2 / (1 + S(0))) * (1 + z)^(3 * (1 + wlambda)))

Omega_tilde_nl_prime(z) = H_nonlocal(z) * (Omega_tilde_nl(z)^2 * (1 - (OmegaBh2 + OmegaCh2) / h^2 / (1 + S(0))) / ((OmegaBh2 + OmegaCh2) / h^2 / (1 + S(0)))) * (-3) / (1 + z)^4

q(k) = k / 13.41 / keq
C(k) = 14.2 / alpha_c + 386.0 / (1 + 69.9 * q(k)^1.08)
T0(k) = log(exp(1.0) + 1.8 * beta_c * q(k)) / (log(exp(1.0) + 1.8 * beta_c * q(k)) + C(k) * q(k)^2.0) #Basic Transfer function
W2(x) = (3 * (sin(x) - x * cos(x)) / x^3)^2

function initial_conditions_phi(D_100, D_prime_100, k)
    A1(k) = 2 * H_nl(100) + 0.5 * beta(100) + 2 * Sprime(100) / (1 + S(100)) - k^2 / 3 / H_nl(100)
    A2(k) = -2 * beta(100)^2 + 2 * Ssecond(100) / (1 + S(100)) - 2 * Sprime(100) / (1 + S(100)) * H_nl_prime(100) / H_nl(100) + k^2 * H_nl_prime(100) / 3 / H_nl(100)^2 - H_nl_prime(100) - 0.5 * betaprime(100) + 3 * (1 - Omega_tilde_nl(100)) * H_nl(100) * (H_nl(100) + beta(100))
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
    M2(z) = 3 * (1 - Omega_tilde_nl(z)) * H_nl(z) * (H_nl(z) + beta(z)) - 2 * beta(z)^2 + 2 * Ssecond(z) / (1 + S(z)) - 2 * Sprime(z) / (1 + S(z)) * H_nl_prime(z) / H_nl(z)
    u0 = [phi_100, -phi_prime_100 / H_nonlocal(100)]
    zspan = (100.0, 0.0)
    function phi_ODE!(du, u, p, t)
        phi, phi_dot = u
        du[2] = dphi_doubledot = -(H_nonlocal_prime_z(t) - M1(t)) / H_nonlocal(t) * phi_dot - M2(t) / H_nonlocal(t)^2 * phi
        du[1] = dphi = phi_dot
    end
    prob = ODEProblem(phi_ODE!, u0, zspan)
    sol = solve(prob, saveat = save_step_length)
    return sol
end

function phi_solve(k, z_target)
    phi_100, phi_prime_100 = initial_conditions_phi(1 / 101, -1 / (101)^2 * (-H_nonlocal(100)), k)
    M1(z) = 3 * H_nl(z) + beta(z) + 2 * Sprime(z) / (1 + S(z))
    M2(z) = 3 * (1 - Omega_tilde_nl(z)) * H_nl(z) * (H_nl(z) + beta(z)) - 2 * beta(z)^2 + 2 * Ssecond(z) / (1 + S(z)) - 2 * Sprime(z) / (1 + S(z)) * H_nl_prime(z) / H_nl(z)
    u0 = [phi_100, -phi_prime_100 / H_nonlocal(100)]
    zspan = (100.0, z_target)
    function phi_ODE!(du, u, p, t)
        phi, phi_dot = u
        du[2] = dphi_doubledot = -(H_nonlocal_prime_z(t) - M1(t)) / H_nonlocal(t) * phi_dot - M2(t) / H_nonlocal(t)^2 * phi
        du[1] = dphi = phi_dot
    end
    prob = ODEProblem(phi_ODE!, u0, zspan)
    sol = solve(prob, saveat = save_step_length)
    return sol
end

function plot_phi(k)
    sol = phi_solve(k)
    plot(sol, layout=(2, 1),title = ["ϕ(z)" "dϕ/dz(z)"], xlabel="z", label=["ϕ, k = $k" "dϕ/dz, k = $k"])
end

function plot_phi!(k)
    sol = phi_solve(k)
    plot!(sol, layout=(2, 1), xlabel="z", label=["ϕ, k = $k" "dϕ/dz, k = $k"])
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
        d = 2 / (H_nl(z[i]) * Omega_tilde_nl(z[i])) * (-H_nonlocal(z[i]) * phi_dot[i] + (k^2 / 3 / H_nl(z[i]) + H_nl(z[i]) + 0.5 * beta(z[i])) * phi[i])
        push!(d_sol, d)
    end
    return z, d_sol
end

function d_solve(k, z_target)
    sol = phi_solve(k, z_target)
    df = DataFrame(sol)
    z = df[!, 1]
    phi = df[!, 2]
    phi_dot = df[!, 3]
    d_sol = Float64[]
    num = length(z)
    for i in 1:num
        d = 2 / (H_nl(z[i]) * Omega_tilde_nl(z[i])) * (-H_nonlocal(z[i]) * phi_dot[i] + (k^2 / 3 / H_nl(z[i]) + H_nl(z[i]) + 0.5 * beta(z[i])) * phi[i])
        push!(d_sol, d)
    end
    return z, d_sol
end

function plot_d(k)
    z_arr, d_arr = d_solve(k)
    plot(z_arr, d_arr, title = "D(z)",xlabel="z", ylabel="D", label="k = $k - Nonlocal")
end

function plot_d!(k)
    z_arr, d_arr = d_solve(k)
    plot!(z_arr, d_arr, xlabel="z", ylabel="D", label="k = $k - Nonlocal")
end

function power_spectrum_solve(k_order_min, k_order_max, z)
    data_count = 1000
    k_range = 10 .^ range(k_order_min, stop = k_order_max, length=data_count)
    k_range = k_range .* H_nonlocal(0) ./ (3 * 10^5)
    ps_arr = Float64[]
    for k in k_range
        z_arr, d_arr = d_solve(k, z)
        ps = TF.Tk_EH_full(k, Omega_r , OmegaC, OmegaB, Omega_l, h, CMB_temp)^2 * d_arr[end] ^ 2 * k ^ n_s
        push!(ps_arr, ps)
    end
    ps_arr = ps_arr .* (4 * pi * (3 * 10^5 / H_nonlocal(0))^4 / Omega_m^2 * A_COBE)
    return k_range ./ h, ps_arr
end

function power_spectrum_plot(k_order_min, k_order_max, z)
    k_arr, ps_arr = power_spectrum_solve(k_order_min, k_order_max, z)
    plot(k_arr, ps_arr, title = "Nonlocal Power Spectrum", xlabel="k (1/Mpc)", ylabel="P(k)", label="z = $z", xaxis = :log, yaxis = :log)
end

function power_spectrum_plot!(k_order_min, k_order_max, z)
    k_arr, ps_arr = power_spectrum_solve(k_order_min, k_order_max, z)
    plot!(k_arr, ps_arr, xlabel="k (1/Mpc)", ylabel="P(k)", label="z = $z", xaxis = :log, yaxis = :log)
end

function save_power_spectrum(title, k_order_min, k_order_max, z)
    K, Ps = power_spectrum_solve(k_order_min, k_order_max, z)
    d =DataFrame(k_over_h = K, P = Ps)
    CSV.write(title, d)
end

function CAMB_PS()
    df = DataFrame(CSV.File("../../Data/test_matterpower.csv", delim="   "))
    k_h_CAMB = df[:,:"k/h"]
    P_CAMB = df[:,:"P"]
    return k_h_CAMB, P_CAMB
end

function PS_ratio_plot()
    df = DataFrame(CSV.File("../../Data/LCDM_matterpower.csv", delim="   "))
    k_h_CAMB = df[:,:"k/h"]
    P_CAMB = df[:,:"P"]
    k_range = k_h_CAMB .* h
    ps_arr = Float64[]
    for k in k_range
        z_arr, d_arr = d_solve(k, 0)
        ps = TF.Tk_EH_full(k, Omega_r , OmegaC, OmegaB, Omega_l, h, CMB_temp)^2 * k ^ n_s * d_arr[end] ^ 2 
        push!(ps_arr, ps)
    end
    ps_arr = ps_arr .* (8 * pi^2 / 25 * (3 * 10^5 / H_nonlocal(0))^4 / Omega_m^2 * A_COBE)
    p_ratio = ps_arr ./ P_CAMB .- 1
    normalize = ps_arr[250] / P_CAMB[250]
    ps_arr = ps_arr ./ normalize
    #plot(k_h_CAMB, p_ratio, title = "Power spectrum ratio with CAMB", xlabel = "k/h", ylabel = "P_nl / P_CAMB - 1", xaxis = :log)
    plot(k_h_CAMB, P_CAMB, label = "CAMB", xlabel = "k/h", ylabel = "P", xaxis = :log, yaxis = :log)
    plot!(k_h_CAMB, ps_arr, label = "Nonlocal", xaxis = :log, yaxis = :log )
end

function PS_ratio_ratio_plot()
    df = DataFrame(CSV.File("../../Data/LCDM_matterpower.csv", delim="   "))
    k_h_CAMB = df[:,:"k/h"]
    P_CAMB = df[:,:"P"]
    k_range = k_h_CAMB .* h
    ps_arr = Float64[]
    for k in k_range
        z_arr, d_arr = d_solve(k, 0)
        ps = TF.Tk_EH_full(k, Omega_r , OmegaC, OmegaB, Omega_l, h, CMB_temp)^2 * k ^ n_s * d_arr[end] ^ 2 
        push!(ps_arr, ps)
    end
    ps_arr = ps_arr .* (8 * pi^2 / 25 * (3 * 10^5 / H_nonlocal(0))^4 / Omega_m^2 * A_COBE)
    normalize = ps_arr[250] / P_CAMB[250]
    ps_arr = ps_arr ./ normalize
    p_ratio = ps_arr ./ P_CAMB .- 1
    plot(k_h_CAMB, p_ratio, title = "Power spectrum ratio with CAMB", xlabel = "k/h", ylabel = "P_nl / P_CAMB - 1", xaxis = :log)
end

function PS_integrate(title, R)
    df = DataFrame(CSV.File(title, delim=","))
    k_h_0 = df[:,:"k_over_h"]
    k_h_0 = k_h_0 .* h
    P_0 = df[:,:"P"]
    #index = findfirst(x -> x>=1, k_h_0)
    index = 1
    k_h = k_h_0[index:end]
    P = P_0[index:end]
    integrable = (1/(2 * pi^2)) .* (k_h .^ 2) .* (P) .* W2.(k_h .* (R * h))
    prob = SampledIntegralProblem(integrable, k_h)
    method = TrapezoidalRule()
    sigma = solve(prob, method)
    return sigma.u
end

function PS_integrate_CDM(title, R)
    df = DataFrame(CSV.File(title, delim="   "))
    k_h_0 = df[:,:"k/h"]
    k_h_0 = k_h_0 .* h
    P_0 = df[:,:"P"]
    #index = findfirst(x -> x>=0.01, k_h_0)
    index = 1
    k_h = k_h_0[index:end]
    P = P_0[index:end]
    integrable = (1/(2 * pi^2)) .* (k_h .^ 2) .* (P) .* W2.(k_h .* (R * h))
    prob = SampledIntegralProblem(integrable, k_h)
    method = TrapezoidalRule()
    sigma = solve(prob, method)
    return sigma.u
end

function PS_integrate_plot(title)
    data_count = 100
    R_range = 10 .^ range(-1, stop = 3, length=data_count)
    sigmaR = sqrt.(PS_integrate.(title, R_range)) #./ PS_integrate(title, 8)) .* 0.8
    plot(R_range, sigmaR, xaxis = :log, yaxis = :log, title = "σ8 vs R", xlabel = "R(Mpc/h)", ylabel = "σ(R)", legend = false)
end

function PS_integrate_plot_CDM(title)
    data_count = 100
    R_range = 10 .^ range(-1, stop = 3, length=data_count)
    sigmaR = sqrt.(PS_integrate_CDM.(title, R_range)) #./ PS_integrate(title, 8)) .* 0.8
    plot(R_range, sigmaR, xaxis = :log, yaxis = :log, title = "σ8 vs R", xlabel = "R(Mpc/h)", ylabel = "σ(R)", legend = false)
end

function PS_fix_sigma8(title)
    df = DataFrame(CSV.File(title, delim=","))
    k_h = df[:,:"k_over_h"]
    P_0 = df[:,:"P"]
    sigma_8_ratio = (0.8^2 / PS_integrate("../../Data/Nonlocal_Linear_Powerspectrum.csv", 8))
    Ps = P_0 .* sigma_8_ratio
    d =DataFrame(k_over_h = k_h, P = Ps)
    CSV.write("../../Data/Nonlocal_powerspectrum_sigma8.csv", d)
end

function PS_fixed_sigma8_plot(title)
    df = DataFrame(CSV.File(title, delim=","))
    k_h = df[:,:"k_over_h"]
    P = df[:,:"P"]
    df = DataFrame(CSV.File("../../Data/LCDM_matterpower.csv", delim="   "))
    k_h_CDM = df[:,:"k/h"]
    P_CDM = df[:,:"P"]
    plot(k_h, P, title = "Nonlocal Power Spectrum Corrected \n with σ8", xlabel = "k/h", ylabel = "P", label = "Nonlocal fixed", xaxis = :log, yaxis = :log)
    plot!(k_h_CDM, P_CDM, label = "LCDM")
end

function power_spectrum_solve_linear(k_order_min, k_order_max, z)
    data_count = 1000
    k_range = range(10.0^k_order_min, stop = 10.0^k_order_max, length=data_count)
    k_range = k_range .* H_nonlocal(0) ./ (3 * 10^5) ./ h
    ps_arr = Float64[]
    for k in k_range
        k *= h
        z_arr, d_arr = d_solve(k, z)
        ps = TF.Tk_EH_full(k, Omega_r , OmegaC, OmegaB, Omega_l, h, CMB_temp)^2 * d_arr[end] ^ 2 * k ^ n_s
        push!(ps_arr, ps)
    end
    ps_arr = ps_arr .* (4 * pi * (3 * 10^5 / H_nonlocal(0))^4 / Omega_m^2 * A_COBE)
    return k_range, ps_arr
end

function save_power_spectrum_linear(title, k_order_min, k_order_max, z)
    K, Ps = power_spectrum_solve_linear(k_order_min, k_order_max, z)
    d =DataFrame(k_over_h = K, P = Ps)
    CSV.write(title, d)
end
end