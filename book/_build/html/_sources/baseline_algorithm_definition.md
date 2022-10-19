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
		subgraph external
			TEC[TEC]
			ECMWF[ECMWF analysis]
			ERA5["historical ERA5 (historical)"]
		end
	end
	subgraph inversion
		F[Forward model]
		cost[cost function]
		apriori[A priori] 
		OS[optimal state]

	end

	subgraph Output data L1R
		geo[geophysical variables]
		unc[geophsysical uncertainties]
		outtb[brightness tempreature residuals]
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
%\text{SSS}\\
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
T_{b,s} = C_{\text{ow}}*ε_{\text{ow}}*T_{\text{ow}}+C_{\text{fyi}}*ε_{\text{fyi}}*T_{\text{fyi}}+C_{\text{myi}}*ε_{\text{myi}}*T_{\text{myi}},
```

where $C_{\text{ow}}$, $C_{\text{fyi}}$, and $C_{\text{myi}}$ are the area
fraction of ocean water, first year ice, and multi year ice, respectively.
$ε_{\text{ow}}$, $ε_{\text{fyi}}$, and $ε_{\text{myi}}$ are the emissivity of
ocean water, first year ice, and multi year ice, respectively. $T_{\text{ow}}$,
$T_{\text{fyi}}$, and $T_{\text{myi}}$ are the brightness temperature of ocean
water, first year ice, and multi year ice, respectively. $C_{\text{ow}}$,
$C_{\text{fyi}}$, and $C_{\text{myi}}$ are adding up to one. With the fruequency dependent emissivities
for first year ice and multi year ice derived by {cite}`Mathew2009`.  With
$C_{\text{myi}}+C_{\text{fyi}}=\text{SIC}$ and $\text{SIC}*C_\text{myi} =
\text{MYIF}$ being part of the state vector {eq}`eqxy`, the state defines
the surface area fraction of all three considered surface types.

The reflectivity of
the surface, $R_{\text{surf}}$, is given by 
```{math}
:label: eq:surfref
R_{\text{surf}} = 1 -
ε_{\text{ow}}*C_{\text{ow}} - ε_{\text{fyi}}*C_{\text{fyi}} -
ε_{\text{myi}}*C_{\text{myi}}.
```

## L-band forward model
The addition for 1.4GHz was done by {cite}`Scarlat2020` and is based on {cite}`Ruf2003`. 






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


