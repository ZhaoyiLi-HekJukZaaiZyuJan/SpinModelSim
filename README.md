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

<section>
  <h2>Variables</h2>
<p>This program accepts command-line arguments, which can be used to customize its behavior. The available options for each subprogram are listed in the following table, along with any options that are not available. Hover your mouse over to see more details. </p>
<table>
  <tr>
    <th>Option</th>
    <th>Description</th>
    <th>Type</th>
    <th>Default Value</th>

  </tr>
  <tr>
    <td>-f, --fname</td>
    <td title="">Output filename, the default file name is in the <code>../data/</code> folder</td>
    <td>string</td>
    <td>N/A</td>
  </tr>
  <tr>
    <td>--out</td>
    <td>output in this directory as an .out file</td>
    <td>bool</td>
    <td>0</td>
  </tr>
  <tr>
    <td>--waitSweep</td>
    <td title="">number of sweeps to wit before sampling</td>
    <td>int</td>
    <td>10</td>
  </tr>
    <tr>
    <td>--rptSweep</td>
    <td title="">number of sweeps sampled</td>
    <td>int</td>
    <td>10000</td>
  </tr>
  <tr>
    <td>-n</td>
    <td>Size of Lattice</td>
    <td>int</td>
    <td>3</td>
  </tr>
  <tr>
    <td>-J</td>
    <td>interaction strength</td>
    <td>float between 0-1</td>
    <td>0</td>
  </tr>
  <tr>
    <td>-p</td>
    <td title="randomness for the RBIM">randomness</td>
    <td>float between 0-1</td>
    <td>0</td>
  </tr>
  <tr>
    <td>--tmin</td>
    <td>minimal temperature</td>
    <td>float</td>
    <td>0.1</td>
  </tr>
  <tr>
    <td>--tmax</td>
    <td>maximal temperature</td>
    <td>float</td>
	<td>2.1</td>
  </tr>
  <tr>
    <td>--nt</td>
    <td>number of bins (data points -1)</td>
    <td>nt</td>
	<td>10</td>
  </tr>
</table>