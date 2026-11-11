# Equations for the aligned no-antibiotic control

For species i and foods j = 0, 1, 2 (R0, M1, M2):

\[
\mu_i^{intrinsic}=s_i\mu_{ref},\qquad
\mu_i^{effective}=\mu_i^{intrinsic}(1-c\rho_i).
\]

Resistance rho is sampled uniformly on [0,1]. The default cost c is 0.20.
IC50 is retained for matching the antibiotic model, but does not affect these equations.

Food preferences w are nonnegative, zero for forbidden foods, and sum to one.
Let f_i0 = 1-leak_i and f_i1 = f_i2 = 1. Uptake capacity is:

\[
V_{ij}=\frac{\mu_i^{effective}w_{ij}}{\sum_k w_{ik}f_{ik}},
\qquad u_{ij}=V_{ij}\frac{R_j}{K_{ij}+R_j}.
\]

Gross growth and living biomass obey:

\[
g_i=\sum_j f_{ij}u_{ij},\qquad \frac{dN_i}{dt}=N_i g_i.
\]

There is no baseline, starvation, or antibiotic death. Thus N cannot decrease
under these ODEs. Growth approaches zero as usable foods vanish.

\[
\frac{dR_0}{dt}=-\sum_iN_i u_{i0},
\]

\[
\frac{dM_j}{dt}=-\sum_iN_i u_{ij}
+\sum_i leak_i q_{ij}N_i u_{i0},\quad j=1,2.
\]

Secretion fractions q_i1 + q_i2 = 1. Only R0 consumption produces metabolites.
Resources are substitutable; there is no replenishment, dilution, or necromass recycling.
Food categories are balanced across the seven nonempty food subsets.

Defaults: 25 species, total initial biomass 0.01, R0=1, M1=M2=0,
equal initial shares unless structured_community is enabled, time 0–300.
Time and concentrations use normalized units. The accounting budget is initial
resource plus initial biomass; loss counters remain zero in normal control runs.
