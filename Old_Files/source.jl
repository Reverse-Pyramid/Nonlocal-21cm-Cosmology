using Plots, Colors
#Parameters
z_eq = 3400
#Functions
y(z) = (1 + z_eq) / (1 + z)
D1(z) = 1 + 3 / 2 * y(z)
D2(z) = (1 + 3 / 2 * y(z)) * log1p((sqrt(1 + y(z)) + 1) / (sqrt(1 + y(z)) - 1) - 1) - 3 * sqrt(1 + y(z))

#plots
plot(D1, title = "D(z) Normal Plot", xlabel = "z", ylabel = "D(z)", xlims = (1, 20), label = "D1", dpi = 500)
plot!(D2, label = "D2")

plot(D1, title = "D(z) Log-Log Plot", xlabel = "z", ylabel = "D(z)", xlims = (1, 20), label = "D1", xaxis = :log, yaxis = :log, dpi = 500)
plot!(D2, label = "D2")

#Parameters
Omega_m = 0.317
Omega_l = 1 - Omega_m
#a(z) = 1 / (1 + z)
H(a) = sqrt(Omega_m * a^(-3) + Omega_l)
steps = 100000
#Functions
f(a) = 1 / (a * H(a)) ^ 3
function integral(f, c1, c2, steps)
    delta = (c2 - c1) / steps
    global results = 0
    for i in 1:steps
        global results += f(c1 + i * delta) * delta
    end
    return results
end

Dplus(a) = 5 * Omega_m / 2 * H(a) * integral(f, 0, a, steps) 

#plots
plot(Dplus, title = "D+(a)/a Log Normal Plot", xlabel = "Log(a)", ylabel = "D+(a)/a", xlims = (0.01, 1), label = "Omega_l = 0.69" , dpi = 500)
Omega_m = 0.28
plot!(Dplus, label = "Omega_l = 0.72")

#more parameters
Omega_m = 0.317
alpha = 0.01
beta_s = -1.0
H0 = 67.4
h = 0.674
deltaH0 = 0.5
OmegaBh2 = 0.0224
deltaOmegeBh2 = 0.0001
OmegaCh2 = 0.12
deltaOmegaCh2 = 0.001
wB = 0
wC = 0
wLambda = -1.03
deltawLambda = 0.03
G = 6.67 * 10^(-11) 

#nonlocal functions
S(a) = alpha / a ^ beta_s
Gamma(a) = (1 + S(1)) / (1 + S(a))
H_NonLocal(a) = Gamma(a) * sqrt(OmegaBh2 / h^2 / (1 + S(1)) / a ^ (3 * (1 + wB)) + OmegaCh2 / h^2 / (1 + S(1)) / a ^(3 * (1 + wC)) + (1 - (OmegaBh2 + OmegaCh2) / h^2 / (1 + S(1))) / a ^ (3 * (1 + wLambda)))

f_NonLocal(a) = 1 / (a * H_NonLocal(a)) ^ 3
Dplus_NonLocal_Scaled(a) = 5 * Omega_m / 2 * H_NonLocal(a) * integral(f_NonLocal, 0, a, steps) / a * 1.1305729371069084 #constant to normalize the nonlocal D

#plots
plot(Dplus, title = "D+(a)/a for normal and nonlocal - Log Normal Plot", xlabel = "Log(a)", ylabel = "D+(a)/a", xlims = (0.01, 1), label = "Normal", xaxis = :log , dpi = 500)
Omega_m = 0.28
plot!(Dplus_NonLocal_Scaled, label = "nonlocal scaled")


G = 6.67 * 10^(-11) 
H(z) = sqrt(Omega_m * (1 + z)^(3) + Omega_l)
function VelocityVerlet(j, D00, D10)
    f1(z, k) = 2 * H(z)
    f2(z, k) = -3 / 2 * H(z)^2 * Omega_m * (1 + z)^(3)
    D2(z, k, D1, D0) = f1(z, k) * D1 + f2(z, k) * D0
    steps = 20 * 10^j
    h = 1 / 10^j
    D0s = Float64[D00]
    D1s = Float64[D10]
    z = Float64[0]
    for i in 1:Int(steps)
        push!(z, i * h)
        push!(D0s, D0s[i] + D1s[i] * h + 1 / 2 * D2(z[i], 0, D1s[i], D0s[i]) * h^2)
        push!(D1s, D1s[i] + 1 / 2 * (D2(z[i], 0, D1s[i], D0s[i]) + D2(z[i + 1], 0, D1s[i], D0s[i + 1])) * h)
    end
plot(z, D0s, background = :black, dpi = 500, xlabel = "z", ylabel = "D(z)", title = "", lw = 2, xlims = (1, 15), ylims = (-1, 1))
end


G = 6.67 * 10^(-11) 
H(a) = sqrt(Omega_m * (a)^(-3) + Omega_l)
tilde_Omega_m(a) = Omega_m / (Omega_m + Omega_l * a^3)
function Euler(j, D00, D10)
    f1(a, k) = -3 / a * (1 - 0.5 * tilde_Omega_m(a))
    f2(a, k) =  3 / 2 * tilde_Omega_m(a) * a^(-2)
    D2(a, k, D1, D0) = f1(a, k) * D1 + f2(a, k) * D0
    steps = 0.99 * 10^j
    h =  1 / 10^j
    D0s = Float64[D00]
    D1s = Float64[D10]
    a = Float64[0.01]
    for i in 1:Int(steps)
        push!(a, 0.01 + i * h)
        push!(D1s, D1s[i] + D2(a[i], 0, D1s[i], D0s[i]) * h)
        push!(D0s, D0s[i] + D1s[i] * h)
    end
    #a = 1 ./ a .- 1
    D0a = D0s ./ a 
plot(a, D0a, background = :black, dpi = 500, xlabel = "a", ylabel = "D(a)", title = "", lw = 2 , xaxis = :log)
end

Euler(5, 0.01, 1)


G = 6.67 * 10^(-11) 
Hz(z) = sqrt(Omega_m * (1 + z)^(3) + Omega_l)
tilde_Omega_mz(z) = Omega_m / (Omega_m + Omega_l * (1 + z) ^ (-3))
function Eulerz(j, D00, D10)
    f1(z, k) = (1 - 1.5 * tilde_Omega_mz(z)) / (1 + z)
    f2(z, k) =  3 / 2 * tilde_Omega_mz(z) * (1 + z)^(-2)
    D2(z, k, D1, D0) = f1(z, k) * D1 + f2(z, k) * D0
    steps = 100 * 10^j
    h =  -1 / 10^j
    D0s = Float64[D00]
    D1s = Float64[D10]
    z = Float64[100]
    for i in 1:Int(steps)
        push!(z,  100 + i * h)
        push!(D1s, D1s[i] + D2(z[i], 0, D1s[i], D0s[i]) * h)
        push!(D0s, D0s[i] + D1s[i] * h)
    end
plot(z, D0s, background = :black, dpi = 500, xlabel = "z", ylabel = "D(z)", title = "", lw = 2)
end

Eulerz(4, 1 / 101, -1 / (101)^2)

plot(tilde_Omega_mz, xlims = (0, 100))

alpha_s = 0.0034
beta_s = -4.0
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
G = 6.67 * 10^(-11) 

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

function nonlocal_delta_plot(k_dimless)
    j = 3
    k = k_dimless * H_nonlocal(0)
    steps = 100 * 10^j
    steplen =  -1 / 10^j
    delta_init = 1 / (101)
    deltaprime_init = -1 / (101)^2
    theta_init = (deltaprime_init / H_nonlocal(100) - A1(100, k) * delta_init) / A2(100, k)
    #thetaprime_init = - B1(100, k) * theta_init - B2(100, k) * delta_init
    thetaz = Float64[theta_init]
    #thetaprimez = Float64[thetaprime_init]
    deltaz = Float64[delta_init]
    deltaprimez = Float64[deltaprime_init]
    Z = Float64[100]
    thetaprime(z, k, theta, delta) =  (B1(z, k) * theta + B2(z, k) * delta) / H_nonlocal(z)
    deltaprime(z, k, theta, delta) =  (A1(z, k) * delta + A2(z, k) * theta) / H_nonlocal(z)
    for i in 1:Int(steps)
        push!(Z,  100 + i * steplen)
        push!(deltaz, deltaz[i] + deltaprime(Z[i], k, thetaz[i], deltaz[i]) * steplen)
        push!(thetaz, thetaz[i] + thetaprime(Z[i], k, thetaz[i], deltaz[i]) * steplen)
    end
    plot!(Z, deltaz, background = :black, dpi = 500, xlabel = "z", ylabel = "δ(z)", lw = 2, label = "k = $k_dimless")
end


Eulerz(4, 1 / 101, -1 / (101)^2)
nonlocal_delta_plot(10.0)
nonlocal_delta_plot(7.0)
nonlocal_delta_plot(6.0)

Eulerz(4, 1 / 101, -1 / (101)^2)
nonlocal_delta_plot(20.0)

k_chosen = 5.0
Z = []
Rs = []
for i in 1:10000
    push!(Z, 0.01 * i)
    push!(Rs, gothic_R(0.01 * i, k_chosen))
end
plot(Z, Rs)

