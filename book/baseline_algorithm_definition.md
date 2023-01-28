# Baseline Algorithm Definition
## Introduction

In the following the main algorithm for the multi parameter retrieval is
described in detail. The algorithm is based on the works of
{cite}`Pedersen1991,Scarlat2017,Scarlat2018,Scarlat2020`. The algorithm is
divided into several steps, which are described in the following.
```{mermaid}
graph TD
	subgraph Input data
		subgraph CIMR L1b
			TBs[TBs]
			TBe[TB error]
		end
		subgraph External
			TEC[TEC]
			ECMWF[ECMWF analysis]
			ERA5["Historical ERA5"]
		end
	end
	subgraph Inversion
		F[Forward model]
		cost[Cost function]
		apriori[A priori] 
		OS[Optimal state]

	end

	subgraph Output data L1R
		geo[Geophysical variables]
		unc[Geophsysical uncertainties]
		outtb[Grightness tempreature residuals]
	end

	resampling[Resampling processor L1R]
	icorr[Ionospheric correction]

	TBs --> icorr
	TBe --> icorr
	TEC --> icorr


	icorr --> resampling
	ECMWF --> resampling
	resampling -- Measurement vector --> cost
	resampling -- Measurement uncertainty --> cost
	resampling -- State vector --> apriori 
	resampling -- State uncertainty --> apriori

	ERA5 --Covariance estimate--> apriori
	apriori --> cost
	F-->cost
	cost-->F

	cost --> OS
	OS --> geo
	OS --> unc 
	OS --> outtb
```

# Retrieval Method

## retrieval definition
The retrieval is following a typical scheme with the objective to minimize

$$
\begin{align}
χ^2(\mathbf y,\mathbf x,\mathbf S_e,\mathbf S_a,\mathbf x_a) = \left(\mathbf y-F(\mathbf x)\right)^\mathbf{T}\mathbf S_e^{-1}(\mathbf
y-F(\mathbf x))+(\mathbf x_a-\mathbf x)^{\mathbf T}\mathbf S_a^{-1}(\mathbf x_a
- \mathbf x) 
\end{align}
$$ (eq:chi2)


where 

```{math}
:label: eqxy
\begin{align}
\mathbf y= \begin{bmatrix}
T_{b,h,1.4}\\
T_{b,v,1.4}\\
T_{b,h,6.9}\\
T_{b,v,6.9}\\
T_{b,h,10.7}\\
T_{b,v,10.7}\\
T_{b,h,18.7}\\
T_{b,v,18.7}\\
T_{b,h,36.5}\\
T_{b,v,36.5}\\
\end{bmatrix}, \quad
\mathbf x= \begin{bmatrix}
\text{WSP}\\
\text{TWV}\\
\text{CLW}\\
\text{SST}\\
\text{IST}\\
\text{SIC}\\
\text{MYIF}\\
\text{SIT}\\
\text{SSS}\\
\end{bmatrix}
\end{align}
```

i.e., $\mathbf{y}$ is the vector of input brightness temperatures (measurement vector) and
$\mathbf x$ is the vector with the corresponding physical quantities (state
vector), $\mathbf S_e$ is the error covariance matrix of the input brightness
temperatures, $\mathbf S_a$ is the covariance matrix of the a priori values,
$\textbf x_a$ is the a priori value, and $F$ is the forward operator, which is
the compositional forward model described in {ref}`fw-model`.  The optimal
solution $\hat {\mathbf y}$ which minimizes equation {eq}`eq:chi2` is found by
iteration. The uncertainty of $\mathbf y$ can then be obtained as 

```{math}
:label: eq:error
\mathbf {\hat S}_a= (\mathbf S_a^{-1}+\mathbf{ \hat M}^\mathbf T \mathbf S_e^{-1}\mathbf {\hat M})^{-1}, 
```

with $\mathbf {\hat M}$ is the Jacobian of the optimal solution of
{eq}`eq:chi2`. The uncertainty of the individual parameters is then given by the
diagonal elements of $\mathbf {\hat S}_a$. 



(fw-model)= 
## Forward Model
The Forward Model is the compositional forward model. It consists of individual
components, namely the ocean surface, the sea ice, and the atmosphere. At the
low frequency channels of the {term}`CIMR` satellite, the sensitivity to
atmospheric parameters is relatively small, but nevertheless the atmosphere
needs to be considered.
The forward model for Ocean and Atmosphere for frequencies 6.9, 10.7, 18.7, and 36.5 GHz is used from the AMSR2 ATBD from {cite}`Wentz2000`. 
The surface contribution to the brightness temperature is given by
```{math}
:label: eq:surface
T_{b,s} = C_{\text{ow}}ε_{\text{ow}}T_{\text{ow}}+C_{\text{fyi}}ε_{\text{fyi}}T_{\text{fyi}}+C_{\text{myi}}ε_{\text{myi}}T_{\text{myi}},
```

where $C_{\text{ow}}$, $C_{\text{fyi}}$, and $C_{\text{myi}}$ are the area
fraction of ocean water, first year ice, and multi year ice, respectively.
$ε_{\text{ow}}$, $ε_{\text{fyi}}$, and $ε_{\text{myi}}$ are the emissivity of
ocean water, first year ice, and multi year ice, respectively. $T_{\text{ow}}$,
$T_{\text{fyi}}$, and $T_{\text{myi}}$ are the brightness temperature of ocean
water, first year ice, and multi year ice, respectively. $C_{\text{ow}}$,
$C_{\text{fyi}}$, and $C_{\text{myi}}$ are adding up to one. With the fruequency dependent emissivities
for first year ice and multi year ice derived by {cite}`Mathew2009`. 
```{math}
:label: eq:Nizy
U_{T, t, p}=(a_{t}T_{C}+b_{t}+273.15)ε_{t, p}
```
with $U_{t,h}$ being the temperature corrected upwelling brightness temperature for the polarization $p$
at the air temperature at the surface $T_{C}$ (in °C), $a_{t}$ and $b_{t}$ are the frequency
dependent coefficients from {cite}`Mathew2009` (see {numref}`tab:emtemp`), and
$ε_{t,p}$ is the frequency dependent emissivity for the polarization $p$ from {numref}`tab:c_ice`

The ice thickness dependence at all frequencies is introduced via modification of the emission from the ice surface.
The ice emissivity is modified by the ice thickness $\text{SIT}$ according to
```{math}
:label: eq:ice_thickness
T_{b,p} = a_{p}-(a_{p}-b_{p})\exp\left(-\frac{\text{SIT}}{c_{p}}\right)
```
with the index $p$ indicate polarization ($h$ or $v$) the coefficients $a_{p}$,
$b_{p}$, and $c_{p}$ from {cite}`Scarlat2020` (see {numref}`tab:fy_thick`). To
combine ice temperatuer and ice thickness dependence, the ice emissivity is
modified by the ice thickness  by substituting $a_{p}$ with the ice temperature
dependent $U_{T, t, p}$ from {eq}`eq:Nizy` for $t=\text{FYI}$. This was not performed in
{cite}`Scarlat2020` but is essential for the minimization of the cost function
to not introduce discontinuities in the forward model. The MYI emissivity is not
affected by the ice thickness in this forward model, as it was not part of the
thickness sensitivity study in {cite}`Scarlat2020`. 

To get back to emissivity in order to account for the atmospheric contribution
to the brightness temperatures at surface level, the brightness temperature is just devided by the ice surface temperature
```{math}
:label: eq:ist
ε_p = \frac{T_{b,p}}{\text{IST}}
```

With
$C_{\text{myi}}+C_{\text{fyi}}=\text{SIC}$ and $\text{SIC}*C_\text{myi} =
\text{MYIF}$ being part of the state vector {eq}`eqxy`, the state defines
the surface area fraction of all three considered surface types.

The reflectivity of
the surface, $R_{\text{surf}}$, is given by 
```{math}
:label: eq:surfref
R_{\text{surf}} = 1 -
ε_{\text{ow}}C_{\text{ow}} - ε_{\text{fyi}}C_{\text{fyi}} -
ε_{\text{myi}}C_{\text{myi}}.
```

The atmospheric contribution to the brightness temperature is calculated from
the parametrization from {cite}`Wentz2000`. They fitted the donwelling and
upwelling effective temperature using a least squares fit to the {term}`TWV`.
The fit is given by 
```{math}
:label: eq:atm
\begin{align}
T_D=b_0+b_1V+b_2V^2+b_3V^3+b_4V^4+b_5\zeta(T_s-T_v)\\
D_U=T_D+b_6+b7V\\
\end{align}
```
where $T_v = 273.16+0.8337 V - 3.029\cdot 10^{-5}V^{3.33}$ for V<48 and
$T_v=301.16$ for V>48, $\zeta(x)=1.05x(1-x^2)/1200$ for $|x|<20\ \text{K}$ and
$\zeta(x)=\text{sign}(x)*14\ \text{K}$ for $|x|>20\ \text{K}$. 

The absorption by oxygen is given by 
```{math}
:label: eq:abs_oxy
A_o = a_{O1}+a_{O2}(T_{D}-270)
```

The vapor absorption is given by
```{math}
:label: eq:abs_vap
A_V = a_{V1}V+a_{V2}V^2
```

The liquid water absorption is given by
```{math}
:label: eq:abs_liq
A_L = a_{L1}(x-a_{L2}(T_L - 283))L
```

with $L$ being the liquid water path.

The total atmospheric attenuation is given by a combination of the individual terms as from {eq}`eq:abs_oxy`, {eq}`eq:abs_vap`, and {eq}`eq:abs_liq`:

```{math}
:label: eq:abs_tot
\tau = \exp\left(-\frac{A_o+A_V+A_L}{\cos\theta}\right)
```
with $\theta$ being the incidence angle.

## Emission from water
While for the sea ice an empirical model is used for the emissivity, for the
emission from the ocean the model uses the Fresnell reflection coefficient as a
bases, which relies on the dielectric constant of the sea water. The emission from calm sea water after {cite}`Wentz2012` is given by
```{math}
:label: eq:em_ocean
\begin{align}
E_{0p} &= 1-|r_p|^2\\ 
r_v &= \frac{ε\cos(θ_i)-\sqrt{ε-\sin^2(θ_i)}}{ε\cos(θ_i)+\sqrt{ε-\sin^2(θ_i)}}\\
r_h &= \frac{\cos(θ_i)-\sqrt{ε-\sin^2(θ_i)}}{\cos(θ_i)+\sqrt{ε-\sin^2(θ_i)}}
\end{align}
```

with $ε$ being the dielectric constant of the sea water.
To account for the roughness and other disturbances on the ocean surface, the power eflectivity at each polariztaion $R_{0p}=|r_p|^2$. The contribution of different ocean surface types are modeled by different parameters again following {cite}`Wentz2000`. The reflectivity is a composition of foam covered ocean and clear ocean. With a reduction of the ocean surface emission through the foam by a factor κ, the composite reflectivity is given by
```{math} 
R = (1-f_{\text{foam}})\cdot R_{\text{clear}}+f_{\text{foam}}\cdot κ\cdot R_{\text{clear}}
``` 
with $f_{\text{foam}}$ being the fraction of the ocean surface covered by foam and $R_{\text{clear}}$ being the reflectivity of the clear ocean. With a small loss from diffraction a term $β$ we can express the reflectivity as 
```{math}
R_{\text{clear}} = (1-\beta)R_{\text{geo}}
```
with $R_{\text{geo}}$ being the reflectivity from a standart geometric optics model{cite}`Wentz75`. The combination can then be combining the foam and the diffraction term into one quantity $F=f_{\text{foam}} + \beta -f_{foam}\cdot β - f_{foam} κβ$, which is a monotonic function of wind speed and is addressed by {cite}`Wentz2000` as *catch-all* term. They determined $F$ empirically from experiments with various radiometeres. A fit for F is given by 
```{math}
:label: eq:catchall
\begin{aligned}
F & = m_1W  &(W<W_1)\\
F & = m_1W + \frac{(m_2-m_1)(W-W_1)^2}{2(W_2-W_1)}  &(W_1 \leq W \leq W_2)\\
F & = m_2W - (m_2-m_1)(W_2+W_1)  &(W>W_2),
\end{aligned}
```
a quadratic spline with knots at $W_1=3~\text{m/s}$ and $W_2=12~\text{m/s}$ for v- and $W_1=7~\text{m/s}$ and $W_2=12~\text{m/s}$ for h-polarization. The coefficients $m_1$ and $m_2$ are given in {numref}`tab:mc_m` in the appendix.





### ocean surface roughness

The ocean surface roughness is used as modeled by {cite}`Wentz2000` as
```{math}
R_{geo} = R_0 - (r_0 + r1(θ_i-53) + r_2(T_S-288) + r_3(θ_i-53)(T_S-288))W
```
Where $R_0$ is the specular reflection, θ$_i$ is the incidence angle, T$_S$ is
the sea surface temperature, and W is the wind speed. The coefficients for each
polarization and frequency are in the appendix in {numref}`tab:mc_geo`.



### Dielectric constant of sea water
The dielectric constant of sea water depends on salinity and temperature {cite}`Wentz2004` as 
```{math}
:label: eq:eps_seawater
\begin{aligned}
\varepsilon(T, S)=\frac{\varepsilon_{\mathrm{S}}(T, S)-\varepsilon_1(T, S)}{1+i \nu / \nu_1(T, S)} & +\frac{\varepsilon_1(T, S)-\varepsilon_{\infty}(T, S)}{1+i \nu / \nu_2(T, S)} \\
& +\varepsilon_{\infty}(T, S)-i \frac{\sigma(T, S)}{\left(2 \pi \varepsilon_0\right) \nu}
\end{aligned}
```
with ε$_1$ as the intermediate frequency dielectric constant, ε$_S$ as the static dielectric constant, ε$_\infty$ as the high frequency dielectric constant, $\nu_1$ and $\nu_2$ as first and second Debye relaxition frequencies, respectively and σ as the conductivity. Fits for S=0 are given by {cite}`Wentz2004` as

```{math}
:label: eq:eps_fresh
\begin{aligned}
\varepsilon_{\mathrm{S}}(T, S=0)&=\frac{3.70886 \cdot 10^4-8.2168 \cdot 10^1
T}{4.21854 \cdot 10^2+T}\\
 \varepsilon_1(T, S=0)&=a_0+a_1 T+a_2 T^2 \\
 \nu_1(T, S=0)&=\frac{45+T}{a_3+a_4 T+a_5 T^2} \\
 \varepsilon_{\infty}(T, S=0)&=a_6+a_7 T \\
 \nu_2(T, S=0)&=\frac{45+T}{a_8+a_9 T+a_{10} T^2} \\
\end{aligned}
```
where the salinity dependence is modeled {cite}`Wentz2004` as
```{math}
:label: eq:eps_sal
\begin{aligned}
\varepsilon_{\mathrm{S}}(T, S) & =\varepsilon_{\mathrm{S}}(T, S=0) \cdot \exp \left[b_0 S+b_1 S^2+b_2 T S\right] \\
\nu_1(T, S) & =\nu_1(T, S=0) \cdot\left[1+S \cdot\left(b_3+b_4 T+b_5 T^2\right)\right] \\
\varepsilon_1(T, S) & =\varepsilon_1(T, S=0) \cdot \exp \left[b_6 S+b_7 S^2+b_8 T S\right] \\
\nu_2(T, S) & =\nu_2(T, S=0) \cdot\left[1+S \cdot\left(b_9+b_{10} T\right)\right] \\
\varepsilon_{\infty}(T, S) & =\varepsilon_{\infty}(T, S=0) \cdot\left[1+S \cdot\left(b_{11}+b_{12} T\right)\right]
\end{aligned}
```

However, Wentz et al. {cite}`Wentz2012` published an update to the salinity dependence of ε$_S$ and ν$_1$ as
```{math}
:label: eq:eps_water_sal_update
\begin{aligned}
ε_S(T_S,S) & = ε_S(T_S,S=0) \cdot \exp [b_0S +b_1S^2+b_2T_SS]\\
\nu_1\left(T_S, S\right) & =\nu_1\left(T_S, S=0\right) \cdot  {\left[1+S \cdot\left(d_0+d_1 T_S+d_2 T_S^2+d_3 T_S^3+d_4 t_s^4\right)\right] }
\end{aligned}
```

For conductivity we follow {cite}`Wentz2004` and {cite}`Wentz2012` who used the earlier work by {cite}`Stogryn1995` as
```{math}
:label: eq:sigma
\sigma(T,S) = \sigma(T,S=35) \cdot R_{15}(S)\frac{R_T(S)}{R_{15}(S)},
```

with

```{math}
:label: eq:sigma_support
\begin{aligned}
 \sigma(T, S=35) & = 2.903602+8.607 \cdot 10^{-2} \cdot T+4.738817 \cdot 10^{-4} \cdot T^2\\
& -2.991 \cdot 10^{-6} \cdot T^3+4.3047 \cdot 10^{-9} \cdot T^4 \\
 R_{15}(S) & = S \cdot \frac{37.5109+5.45216 \cdot S+1.4409 \cdot 10^{-2} \cdot S^2}{1004.75+182.283 \cdot S+S^2} \\
 \frac{R_T(S)}{R_{15}(S)} & =1+\frac{\alpha_0(T-15)}{\alpha_1+T} \\
 \alpha_0 & = \frac{6.9431+3.2841 \cdot S-9.9486 \cdot 10^{-2} \cdot S^2}{84.850+69.024 \cdot S+S^2} \\
 \alpha_1 & = 49.843-0.2276 \cdot S+0.198 \cdot 10^{-2} \cdot S^2
\end{aligned}
```





## L-band forward model
The addition for 1.4GHz was done by {cite}`Scarlat2020` and is based on {cite}`Ruf2003`. 

The atmospheric attennuation for L-band is
```{math}
:label: eq:tau_l
\tau = \exp\left(-\frac{0.009364 + 0.0000024127V }{\cos(\theta)}\right),
```
with $V$ being the {term}`TWV` in mm and $\theta$ being the incidence angle.

The up- and downwelling brightness temperature for L-band is given by
```{math}
:label: eq:tbud_l
\begin{align}
T_{b,\text{up}} = (1-\tau)(\text{ST}+258.15)\\
T_{b,\text{down}} = (1-\tau)(\text{ST}+263.15)
\end{align}
``` 
with ST being the surface temperature in K. Over open ocean, this is
{term}`SST` and over ice it is {term}`IST`. The composite upwelling and
downwelling brightness temperature is thenn
calculated with the weight of the ice concentration.


The effect of wind roughenning for L-band over open ocean is given by
```{math}
:label: eq:rough_L
\begin{align}
\epsilon_h = \epsilon_{h,0} + u(0.0007 + 0.000015\theta) \\
\epsilon_v = \epsilon_{v,0} + 0.0007u, 
\end{align}
```
with $u$ being the wind speed in m/s.

The brightness temperature which would be observed by an instrument in orbit is then given by
```{math}
T_{b,p} = T_{b,\text{up}} + \left((T_{c} 𝜏 + T_{b,\text{down}})(1 − \epsilon_p) + T_{b,s}\right) \tau,
```
with $p$ being the polarization and $T_c$ being the cosmic background temperature which, for L-band, is given as 6&nbsp;K  by {cite}`Ruf2003`. $T_{b,s}$ is the surface emitted brightness temperature (from {eq}`eq:surface`).

## from template
### CIMR Level-1b re-sampling approach

Subsection Text


### Algorithm Assumptions and Simplifications

Subsection Text

### Level-2 end to end algorithm functional flow diagram

Subsection Text

### Functional description of each Algorithm step

Subsection Text

### Mathematical description

SubSubsection Text
### Input data

SubSubsection Text

### Output data

SubSubsection Text

### Auxiliary data

SubSubsection Text

### Ancillary data

SubSubsection Text

### Validation process

SubSubsection Text


