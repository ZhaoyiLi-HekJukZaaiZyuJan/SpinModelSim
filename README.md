# Monte Carlo Simulation of 2D Random Bond Ising Model
   Originally Created 3/12/2022 Zhaoyi Li
   
   
   
   
   Hamiltonian: $H=J_{ij}\sum_{<i,j>}\sigma_i\sigma_j$ <br>

### Assuming constant coupling$J_{ij}=J$

Thermodynamic variables:
$$\chi\sim\Delta M, C\sim\Delta E$$
Want to compute $\langle E\rangle$, $\langle E^2\rangle$, $\langle M\rangle$, $\langle M^2\rangle$ with thermodynamic average: $$\langle X\rangle = \frac{1}{Z}\sum_{\{\sigma\}}X(\{\sigma\})e^{-\beta H(\{\sigma\})}$$

### Adding randomness in $J_{\langle ij\rangle}$

$P[J_{ij}] = (1-p)\delta(J_{ij}-J)+p\delta(J_{ij}+J)$ <br>

Here $\langle X\rangle = \mathbb{E}[\frac{1}{Z}\sum_{\{\sigma\}}X(\{\sigma\})e^{-\beta H(\{\sigma\})}]$
![Phase.png](readme_pictures/Phase.png)

## Background

$\cdot$ spin glass <br>
![1920px-Silica.svg.png](attachment:1920px-Silica.svg.png)
$\cdot$ quantum error correction <br>

![The-vertex-Av-and-plaquette-Bp-operators-of-the-toric-code-as-defined-in-2-Edges.png](attachment:The-vertex-Av-and-plaquette-Bp-operators-of-the-toric-code-as-defined-in-2-Edges.png)
![Data-obtained-from-numerical-simulations-of-the-toric-code-failure-rate-close-to.png](attachment:Data-obtained-from-numerical-simulations-of-the-toric-code-failure-rate-close-to.png)