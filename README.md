# Nonlocal-21cm-Cosmology
Our work regarding numerical simulation of nonlocal gravity and it's effect on the 21cm power spectrum.

---
jupyter:
  kernelspec:
    display_name: Julia 1.10.0
    language: julia
    name: julia-1.10
  language_info:
    file_extension: .jl
    mimetype: application/julia
    name: julia
    version: 1.10.0
  nbformat: 4
  nbformat_minor: 2
---

::: {.cell .code}
``` julia
include("../src/ParamsFunctions.jl")
include("../src/ParamsFunctionsLCDM.jl")
include("../src/TransferFunction.jl")
import .ParamsFunctions as PF
import .ParamsFunctionsLCDM as PFL
import .TransferFunction as TF
using Plots
```
:::

::: {.cell .markdown slideshow="{\"slide_type\":\"slide\"}"}
## \# Non-local 21cm Julia package {#-non-local-21cm-julia-package}

In this presentation we will showcase our code related to our work on
21cm power spectrum of LLNL cosmology.`<br>`{=html} Our package includes
numerical solvers for the field equation, power spectrum, transfer
function, sigma 8, halo and 21cm bias, and density
contrast.`<br>`{=html} Parameters can be easily modified to obtain
different models and adjust the strength of non-local effects.
:::

::: {.cell .markdown slideshow="{\"slide_type\":\"slide\"}"}
## \### Parameters {#-parameters}

We chose our parameters from \"Planck 2018\" and \"Tabatabaei 2023 -
LLNL: cosmological perturbations\". `<br>`{=html} Nonlocal Parameters:
`<br>`{=html} $\alpha = 0.0034$, $\beta = -1$ and
$S(z) = \alpha (1 + z)^\beta$
:::

::: {.cell .markdown slideshow="{\"slide_type\":\"slide\"}"}
## \### Field Equations {#-field-equations}

$$
 \phi' + (\frac{k^2}{3\mathcal{H}}+\mathcal{H}+\frac{1}{2}\mathcal{B})\phi = \frac{\kappa a^2}{6(1+\bar{S})\mathcal{H}}\delta\rho = \frac{\mathcal{H} \rho_m}{2 \rho} \delta 
$$ $$
 \phi'' + (3\mathcal{H}+\mathcal{B})\phi' - (3w\mathcal{H}^2 + 3w\mathcal{H}\mathcal{B} +2\mathcal{B}^2)\phi + 2\mathcal{H}\frac{\delta S'}{1+\bar{S}} =0
$$ The second equation is solved by DifferentialEquations.jl ODE solver
in reverse from $z=100$ to $z=0$ (or any other z). `<br>`{=html} Since
non-local effects are miniscule at $z=100$ the initial conditions of
$\phi$ are calculated using $D$ and $D'$ (density contrast for
$\Lambda \text{CDM}$) at $z=100$.`<br>`{=html} Different values of $k$
can also be given to the function.
:::

::: {.cell .code execution_count="4" slideshow="{\"slide_type\":\"slide\"}"}
``` julia
PF.plot_phi(100)
PF.plot_phi!(1000)
```

::: {.output .display_data}
![](presentation/b22e132cba5ca679e11779f3263996580e45f276.png)
:::
:::

::: {.cell .markdown slideshow="{\"slide_type\":\"slide\"}"}
## \### Density Contrast {#-density-contrast}

Using $\phi$ and $\phi'$ we can calculate the density contrast of both
non-local and $\Lambda \text{CDM}$ models.
:::

::: {.cell .code execution_count="6" slideshow="{\"slide_type\":\"slide\"}"}
``` julia
PF.plot_d(1)
PFL.plot_d_LCDM!(1)
PF.plot_d!(100)
PFL.plot_d_LCDM!(100)
```

::: {.output .display_data}
![](presentation/7817637da8f4cba19a118aeeba6a29b4dc331040.png)
:::
:::

::: {.cell .markdown slideshow="{\"slide_type\":\"slide\"}"}
The ratio of density contrast of the two models can also be plotted
:::

::: {.cell .code execution_count="14" slideshow="{\"slide_type\":\"slide\"}"}
``` julia
k = 1
z_nl, d_nl = PF.d_solve(k)
z_LCDM, d_LCDM = PFL.d_solve_LCDM(k)
d_ratio = (d_nl ./ d_LCDM) .^ 2 .- 1
plot(z_nl, d_ratio, title = "(D_nl / D_LCDM)^2 - 1", label = "k = $k", xlabel = "z", ylabel = "(D_nl / D_LCDM)^2 - 1")
k = 100
z_nl, d_nl = PF.d_solve(k)
z_LCDM, d_LCDM = PFL.d_solve_LCDM(k)
d_ratio = (d_nl ./ d_LCDM) .^ 2 .- 1
plot!(z_nl, d_ratio, title = "(D_nl / D_LCDM)^2 - 1", label = "k = $k", xlabel = "z", ylabel = "(D_nl / D_LCDM)^2 - 1")
```

::: {.output .display_data}
![](presentation/db324009753c987d4f6315a273c89a240d7fb4fb.png)
:::
:::

::: {.cell .markdown}
## \### Transfer Function {#-transfer-function}

We used the Eisenstein-Hu transfer function throughout our code.
`<br>`{=html} This module is a direct port of a c++ module by Farbod
Hassani.
:::

::: {.cell .code execution_count="15"}
``` julia
Omega_r = 0
h = 0.7
Omega_c = 0.112 / h ^ 2
Omega_b = 0.0226 / h ^ 2
Omega_l = 1 - 0.317
CMB_temp = 2.726

k_order_min = -1.0
k_order_max = 4.0
data_count = 1000
k_range = 10 .^ range(k_order_min, stop = k_order_max, length=data_count)
k_range = k_range .* PF.H_nonlocal(0) ./ (3 * 10^5)

TF_array = TF.Tk_EH_full.(k_range, Omega_r , Omega_c, Omega_b, Omega_l, h, CMB_temp)
plot(k_range, TF_array, xaxis = :log, yaxis = :log, xlabel = "k (1/Mpc)", ylabel = "TF", title = "EH Transfer Function", legend = false)
```

::: {.output .display_data}
![](presentation/853cc77b4b694c2a0837f25c3e849968f9e4896d.png)
:::
:::

::: {.cell .markdown}
## \## Power Spectrum {#-power-spectrum}

Using the growth function and EH transfer function we calculated, we can
obtain the non-local power spectrum `<br>`{=html} using the equation
below: $$
LLNL \quad PS = A \; {T_{EH}}(k) ^ 2 \; D(k,z) ^ 2 
$$
:::

::: {.cell .code execution_count="25"}
``` julia
PF.power_spectrum_plot(-1.0, 4.0, 0)
PF.power_spectrum_plot!(-1.0, 4.0, 0.5)
PF.power_spectrum_plot!(-1.0, 4.0, 1)
PF.power_spectrum_plot!(-1.0, 4.0, 2)
PF.power_spectrum_plot!(-1.0, 4.0, 10)
```

::: {.output .display_data}
![](presentation/64d8e75d75c70d6001ca5469fc9f0c882c3f83bb.png)
:::
:::

::: {.cell .markdown}
The resulting power spectrum can be multiplied by a constant value so we
fix this value `<br>`{=html} by calculating $\sigma$ with regards to
$R$. `<br>`{=html} But first we can fix the max value of the non-local
power spectrum, with the one obtained from CAMB `<br>`{=html} to
illustrate the non-local effects.
:::

::: {.cell .code execution_count="26"}
``` julia
PF.PS_ratio_plot()
```

::: {.output .display_data}
![](presentation/79d0c8bf0fc2a72bdbfd8923606fdc7d39d3be70.png)
:::
:::

::: {.cell .code execution_count="27"}
``` julia
PF.PS_ratio_ratio_plot()
```

::: {.output .display_data}
![](presentation/b4d17dfb4ba5dddcfe65c61aed4c1ec221138e6a.png)
:::
:::

::: {.cell .markdown}
## \## $\sigma_8$ Fixing {#-sigma_8-fixing}

First $\sigma_8$ is calculated with regards to $R$ using the formula
below: (at z=0) $$
\sigma^2(R) = \int \frac{1}{2 \pi ^2} (\frac{k}{h})^2 P(k) W^2(k \cdot R) dk
$$ where $W(x)$ is a window function defined as: $$
W(x) = 3 \frac{\sin(x) - x \cos(x)}{x^3}
$$ The integration is done numerically using the Integrals.jl package
and the trapezoidal rule method.`<br>`{=html} An important consideration
is that until now our power spectrum was calculated with equal
logarithmic distribution `<br>`{=html} which causes issues with
numerical integration methods. Thus, another power spectrum should be
calculated with linear `<br>`{=html} distribution of scale which will
integrated to plot $\sigma(R)$. `<br>`{=html} `<br>`{=html} After
calculating $\sigma$ for the non-local power spectrum, we fix $\sigma_8$
with the $\Lambda \text{CDM}$ value calculated in the `<br>`{=html} same
manner. The scaled power spectrum will have the correct order of
magnitude.
:::

::: {.cell .code execution_count="30"}
``` julia
PF.PS_integrate_plot("../../Data/Nonlocal Power Spectrum.csv")
```

::: {.output .display_data}
![](presentation/4857cec8a8218b572ba232a972ff1bb215c49f3d.png)
:::
:::

::: {.cell .code execution_count="31"}
``` julia
PF.PS_fix_sigma8("../../Data/Nonlocal Power Spectrum.csv")
PF.PS_fixed_sigma8_plot("../../Data/Nonlocal_powerspectrum_sigma8.csv")
```

::: {.output .display_data}
![](presentation/89331a0f4707d9b561183a210c9a11d0179620e8.png)
:::
:::

::: {.cell .markdown}
## \## Bias {#-bias}

As was discussed in a previous research group meeting, the
Press-Schechter bias is defined as: $$
\begin{align*}
b &= 1 + \frac{\nu^2 - 1}{\delta_c} \\
\nu &= \frac{\delta_c}{\sigma(R)}
\end{align*}
$$ $R$ can be calculated by knowing tracer mass and density. $$
R = \left(  \frac{3 M}{4 \pi \rho}  \right) ^ {(\frac{1}{3})}
$$ The biased power spectrum is calculated as below: $$
P_{Halo} = b^2 \; P_M
$$
:::

::: {.cell .code execution_count="46"}
``` julia
PF.PS_bias_applied_compare("../../Data/Nonlocal_powerspectrum_sigma8.csv", 0, 10.0^13)
```

::: {.output .display_data}
![](presentation/37c1cc55f333e07f1dd0957fe0214cef7a0965df.png)
:::
:::

::: {.cell .markdown slideshow="{\"slide_type\":\"slide\"}"}
## \## Summary {#-summary}

Through out this project we developed a piece of code that provides
useful checks and numerical solutions for the non-local gravity model.
`<br>`{=html} In this presentation we displayed our work for LLNL. The
parameters and functions can be easily modified to generate results for
new models.
:::
