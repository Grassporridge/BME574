# Current antibiotic consumer-resource model

This file describes the current `model.py` with **balanced feeding categories** and an explicit **fitness cost of resistance**.

For species `i = 1,...,25` and resources `j in {R0, M1, M2}`:

## Balanced resource-use categories

There are seven non-empty feeding categories:

1. `R0`
2. `M1`
3. `R0+M1`
4. `M2`
5. `R0+M2`
6. `M1+M2`
7. `R0+M1+M2`

The default generator repeats categories `1..7` and shuffles their species assignments. Therefore category counts differ by at most one. With 25 species the category counts are four species in four categories and three species in the remaining three categories.

## Continuous resistance

Each species has a resistance score

\[
\rho_i \in [0,1],
\]

mapped to IC50 by

\[
IC50_i=IC50_{min}\left(\frac{IC50_{max}}{IC50_{min}}\right)^{\rho_i}.
\]

Higher `rho_i` therefore gives higher `IC50_i` and weaker antibiotic killing at the same antibiotic concentration.

## Fitness cost of resistance

Resistance also reduces maximum resource-supported growth. Define the maximum resistance cost parameter `c_max = resistance_cost_max`. The species-specific growth multiplier is

\[
f_i^{growth}=1-c_{max}\rho_i.
\]

The effective maximum growth rate is

\[
\mu_{i,max}^{eff}=\mu_{i,max}^{intrinsic}(1-c_{max}\rho_i).
\]

The default is

\[
c_{max}=0.20,
\]

so a species with `rho_i=1` has a 20% reduction in maximum growth capacity, a species with `rho_i=0.5` has a 10% reduction, and a species with `rho_i=0` has no resistance cost.

The uptake capacities `V_ij` are normalized so that growth at saturating permitted resources equals this cost-adjusted `mu_max`.

## Resource uptake

\[
u_{ij}(R_j)=V_{ij}\frac{R_j}{K_{ij}+R_j}.
\]

`V_ij` is the cost-adjusted maximum uptake capacity and `K_ij` is the half-saturation constant. If a species cannot consume a resource, its corresponding `V_ij` is zero.

## Gross growth

\[
g_i=(1-\ell_i)u_{i0}+u_{i1}+u_{i2}.
\]

The leakage fraction `l_i` applies to R0 uptake. Leaked R0 energy is secreted as M1/M2 rather than retained for growth.

## Antibiotic killing

\[
d_i(A)=k_{max}\frac{A^h}{IC50_i^h+A^h}.
\]

There is no basal or starvation mortality in the current ODE.

## Species biomass equation

\[
\frac{dN_i}{dt}=N_i\left[g_i-d_i(A)\right].
\]

Expanded:

\[
\frac{dN_i}{dt}=N_i\left[(1-\ell_i)V_{i0}\frac{R_0}{K_{i0}+R_0}
+V_{i1}\frac{M_1}{K_{i1}+M_1}
+V_{i2}\frac{M_2}{K_{i2}+M_2}
-k_{max}\frac{A^h}{IC50_i^h+A^h}\right].
\]

The fitness cost is already contained in the cost-adjusted `V_ij` values through `mu_max`.

## Resource equations

\[
\frac{dR_0}{dt}=-\sum_i N_i u_{i0}.
\]

If `q_i1 + q_i2 = 1` are secretion fractions:

\[
\frac{dM_1}{dt}=-\sum_i N_i u_{i1}+\sum_i \ell_i q_{i1}N_i u_{i0},
\]

\[
\frac{dM_2}{dt}=-\sum_i N_i u_{i2}+\sum_i \ell_i q_{i2}N_i u_{i0}.
\]

## Antibiotic decay

\[
\frac{dA}{dt}=-k_{deg}A,
\]

so

\[
A(t)=A_0 e^{-k_{deg}t}.
\]

## Default initial conditions

- Total initial bacterial biomass: `0.01`
- 25 species with equal initial biomass in the default non-structured case: `N_i(0)=0.01/25=0.0004`
- `R0(0)=1`
- `M1(0)=0`
- `M2(0)=0`
- `A(0)=initial_antibiotic` (default 1 for treatment; 0 for control analyses)
- `resistance_cost_max=0.20`

## Interpretation

Without antibiotic, resistance can still be disadvantageous because high-resistance species have reduced maximum growth capacity. Under antibiotic exposure, the same resistance trait reduces antibiotic killing. This creates an explicit tradeoff between ecological growth performance and antibiotic protection.
