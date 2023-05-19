# Monte Carlo Simulation of 2D Random Bond Ising Model 
Originally Created 3/12/2022 Zhaoyi Li
   
## Technical Details
   
Hamiltonian: The Hamiltonian is given by $H=J_{ij}\sum_{<i,j>}\sigma_i\sigma_j$, where each spin $\sigma_i$ takes on the value $\{-1,1\}.$
### Assuming constant coupling 

$$J_{ij}=J$$

In this case, the thermodynamic variablesare given by:
$$\chi\sim\Delta M, C\sim\Delta E$$
Want to compute $\langle E\rangle$, $\langle E^2\rangle$, $\langle M\rangle$, $\langle M^2\rangle$ with thermodynamic average: $$\langle X\rangle = \frac{1}{Z}\sum_{\{\sigma\}}X(\{\sigma\})e^{-\beta H(\{\sigma\})}$$


### Adding randomness in $J_{\langle ij\rangle}$

$$P[J_{ij}] = (1-p)\delta(J_{ij}-J)+p\delta(J_{ij}+J)$$

Here $\langle X\rangle = \mathbb{E}[\frac{1}{Z}\sum_{\{\sigma\}}X(\{\sigma\})e^{-\beta H(\{\sigma\})}]$
![Phase.png](readme_pictures/Phase.png)

## Background

$\cdot$ spin glass <br>
![1920px-Silica.svg.png](attachment:1920px-Silica.svg.png)
$\cdot$ quantum error correction <br>

![The-vertex-Av-and-plaquette-Bp-operators-of-the-toric-code-as-defined-in-2-Edges.png](attachment:The-vertex-Av-and-plaquette-Bp-operators-of-the-toric-code-as-defined-in-2-Edges.png)
![Data-obtained-from-numerical-simulations-of-the-toric-code-failure-rate-close-to.png](attachment:Data-obtained-from-numerical-simulations-of-the-toric-code-failure-rate-close-to.png)

## Methodology
Metropolis algorithm:

$$w(a\leftarrow b) = \min{(1,e^{\beta(E_a-E_b)})}$$
therefore satisfying the condition:
$$\frac{w(a\leftarrow b)}{w(b\leftarrow a)}=\frac{e^{-\beta E_a}}{e^{-\beta E_b}}$$