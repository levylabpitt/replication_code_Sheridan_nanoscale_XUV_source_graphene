---
title: "Sternheimer "
series: "Tutorials"
authors: "David Strubbe"
---

### Preparing the ground state

In this tutorial, we are going to calculate the dielectric constant of diamond-structure silicon, using the linear-response Sternheimer equation [^footnote-4] (aka density-functional perturbation theory [^footnote-2]) and the quantum theory of polarization [^footnote-3]. Like in the tutorial on dynamic polarizabilities of a molecule, we will use an electric field as a perturbation, but unlike in the case of a finite system, we will need to use a different form for the electric field due to the peculiar properties of electric fields in periodic systems. An electric field $\vec{\mathcal{E}}$ gives rise to a term in the Hamiltonian $\vec{\mathcal{E}} \cdot \vec{r}$, but such a perturbation is not periodic. If we use such a perturbation, we lose the ability to use Bloch's Theorem
( $\psi_{\vec{k}} \left( \vec{r} \right) = e^{i \vec{k} \cdot \vec{r}} u_{\vec{k}} \left( \vec{r} \right)$ ) 
and all the other convenient related formalism for periodic systems. Instead, we can express the Hamiltonian term as 
$\vec{\mathcal{E}} \cdot i \vec{\nabla}_{\vec{k}}$, which is now no longer a local potential but instead an operator that applies to the periodic part of the wavefunction $u_{\vec{k}} \left( \vec{r} \right)$. The general theory is laid out in Refs. [^footnote-4] and [^footnote-5].

While some codes such as Quantum ESPRESSO use this formalism to calculate dielectric constants, or Raman intensities (related to $\frac{\partial \epsilon}{\partial R}$ for an atomic position $R$), they evaluate the derivative with respect to $k$ by finite differences, i.e. 
$\frac{\partial}{\partial k_i} u_{\vec{k}} \approx \left( u_{\vec{k} + \Delta \vec{k}_i} - u_{\vec{k}} \right) / \Delta k_i $. 
In Octopus, instead we use linear response to evaluate that derivative, via the same Sternheimer equation we would use for the electric field perturbation. The perturbation in this case, to the effective Bloch Hamiltonian $H_{\vec{k}} = \frac{\hbar^2}{2m} \left( i \nabla + \vec{k} \right) + V \left( \vec{r}, \vec{r}' \right)$ is $\frac{\partial H_{\vec k}}{\partial k_i} = -i \frac{\partial}{\partial k_i} + k_i + \left[ V_{\rm nl}, r_i \right]$. The main contribution is from the kinetic energy term, but the potential can contribute as well insofar as it is nonlocal, in which case it is effectively a function of $k$ too. The Hartree potential is local. The exchange-correlation potential is local in the Kohn-Sham framework, i.e. LDA, GGA, and meta-GGA functionals, but it is nonlocal when hybrid functionals are used. Moreover, when pseudopotentials are used, as we do in Octopus, these are almost always non-local (angular-momentum-dependent) and therefore there is a non-zero commutator to include. This perturbation is referred to as $k \cdot p$ in Octopus, as it is the same one used in the classic $k \cdot p$ perturbation theory for semiconductor effective masses. Indeed, from our $k \cdot p$ run, we will be able to obtain group velocities and effective masses -- albeit with some complication in both cases whenever there is a degeneracy, in which case a degenerate perturbation theory needs to be used to consider properly the way that bands split when the $k$-point changes from the symmetric point. One simplification of using this perturbation compared to other ones is that since this is a non-physical perturbation, i.e. it does not constitute a change of the physical system, but only of an internal parameter we use to describe states, it should not be solved self-consistently, and so the calculation is faster.

Once we have obtained $\frac{u_{n \vec{k}}}{\partial k_i}$ for each band $n$, k-point $\vec{k}$, and Cartesian direction $i$, we can now calculate the electric field response as we would for a finite system, and obtain a polarizability $\alpha_{ij}$. We generally prefer to think in periodic systems of the intensive quantity $\epsilon_{ij} = 1 + \frac{4 \pi}{V} \alpha_{ij}$. In the simplest case, we can find the static dielectric constant, i.e. the response at frequency $\omega = 0$. Note that this includes only the electronic contribution, not any ionic contributions, and thus is referred to as $\epsilon_\infty$ since the frequency is infinite with respect to phonons (but is very small as far as the electrons are concerned). We can also extend this approach straightforwardly to the dielectric function at nonzero frequencies $\omega$, in which case we are indeed using TDDFT. To prevent divergences at resonant frequencies, and allow us to obtain the imaginary part $\epsilon_2$ (related to optical absorption) as well as the real part $\epsilon_1$, we introduce a small imaginary part of the frequency $i \eta$, which will broaden out these resonances into a Lorentzian lineshape.

Now that we have laid out the theory, let's do some calculations. We begin with a ground-state SCF run as usual, with the input file below. We ask for somewhat higher precision than usual, which helps with Sternheimer calculations' convergence later. We use {{<variable "FilterPotentials">}}{{<code "= filter_none">}} as this is sometimes helpful for convergence in Sternheimer.

```text
#include_input doc/tutorials/periodic_systems/sternheimer/1.gs/inp
```

### Performing a k dot p calculation

Next, we perform the kdotp run, based on the ground-state results. We need the {{<variable "CalculationMode">}}{{<code "= kdotp">}} . The {{<variable ExtraStates>}} is not needed but we will be able to see the velocity and effective mass of the lowest conduction band in this way.

```text
#include_input doc/tutorials/periodic_systems/sternheimer/2.kdotp/inp
```

Look at the file {{<file "kdotp/velocity">}} to find the band velocities for each band and k-point. This is not the typical approach to finding velocities. How else might you obtain them? They play a key role in Boltzmann transport calculations of electric conductivity for example, which you can do in the {{<file "BoltzTraP">}} code.

Then, look at the effective masses in the kdotp/kpoint_1_1 etc. This is a $3 \times 3$ tensor defined from the band energies $\epsilon_{n \vec{k}}$. The inverse effective mass is $m^{-1}_{ij} \left( n, \vec{k} \right) = \frac{1}{\hbar^2} \frac{\partial^2 \epsilon_{n \vec{k}}}{\partial k_i \partial k_j}$. The effective mass tensor itself is obtained by matrix inversion of $m^{-1}$. Note that the complication of degenerate perturbation theory is not dealt with currently in the code (as pointed out in a warning), so focus on either the non-degenerate states, or look at invariants of a degenerate subspace, in which you trace over the degenerate states. The end of the main output file classifies the degenerate subspaces for your convenience.

### Calculating the linear reponse

Finally, we do a em_resp run.

```text
#include_input doc/tutorials/periodic_systems/sternheimer/3.em_resp/inp
```

Look in {{<file "em_resp/freq_0.000/epsilon">}} to find the static dielectric constant. You should obtain a diagonal isotropic tensor ($\epsilon_{xx} = \epsilon_{yy} = \epsilon_{zz}$) as required by the $T_d$ point group of silicon. The value obtained is somewhat overestimated compared to experiment. Sometimes you will hear people attribute this discrepancy to the "band-gap" problem of Kohn-Sham DFT, thinking of a simple perturbation theory expression for the dielectric constant in which a too-small energy denominator will give a too-large result. This reasoning is a misunderstanding however, since the eigenvalues were not used in the calculation of the dielectric constant, and moreover the dielectric constant is a ground-state property which should be calculable according to the Hohenberg-Kohn Theorems regardless of whether electronic excitation energies are computed correctly. We should look instead for the source of any discrepancies in other deficiencies of LDA. We can also look at the variation of $\epsilon$ with frequency by checking {{<file "em_resp/freq_0.1000/epsilon">}}, which is the result for frequency $\hbar \omega$ = 0.1 Ha = 2.72 eV. Is symmetry still preserved? How big are the real and imaginary parts? You will see in the iterations that some lines list a negative number for the state index, i.e. "-2": this means the negative frequency for state 2, as we need both $\omega$ and $-\omega$ for the dielectric constant calculations. As an exercise, you can calculate values from 0 to 0.2 Ha, with 10 points in between (using the {{<variable "EMFreqs">}} block), and make plots of the real and imaginary parts. Compare to results from time-propagation if you have done that tutorial -- the results are identical in principle, but can differ due to the use of $i \eta$ here, and the fact that that time-propagation is not linear response and depends somewhat on the strength of the perturbation used.

Two final notes: We can study partially periodic systems (e.g. monolayer hexagonal boron nitride) through this same approach, where we use the finite scheme for finite directions, and this $k \cdot p$ scheme for periodic directions. We can also extend this scheme to the calculation of magnetic susceptibilities, and even to nonlinear properties like magneto-optical susceptibilities [^footnote-6].

### References

[^footnote-1]: X. Andrade, S. Botti, M. A. L. Marques, and A. Rubio, "Time-dependent density functional theory scheme for efficient calculations of dynamic (hyper)polarizabilities," J. Chem. Phys. 126, 184106 (2007)
[^footnote-2]: S. Baroni, S. de Gironcoli, A. Dal Corso, and P. Gianozzi, "Phonons and related crystal properties from density-functional perturbation theory," Rev. Mod. Phys. 73, 515 (2001).
[^footnote-3]: R. Resta, "Macroscopic polarization in crystalline dielectrics: the geometric phase approach," Rev. Mod. Phys. 66, 899 (1994)
[^footnote-4]: D. A. Strubbe, L. Lehtovaara, A. Rubio, M. A. L Marques, and S. G. Louie, "Response Functions in TDDFT: Concepts and Implementation," in Fundamentals of Time-Dependent Density Functional Theory, Lecture Notes in Physics, Springer (2012).
[^footnote-5]: X. Andrade, D. A. Strubbe, U. De Giovannini, A. H. Larsen, M. J. T. Oliveira, J. Alberdi-Rodríguez, A. Varas, I. Theophilou, N. Helbig, M. Verstraete, L. Stella, F. Nogueira, A. Aspuru-Guzik, A. Castro, M. A. L. Marques, and Á. Rubio, "Real-space grids and the Octopus code as tools for the development of new simulation approaches for electronic systems," Phys. Chem. Chem. Phys. 17, 31371-31396 (2015).
[^footnote-6]: I. V. Lebedeva, D. A. Strubbe, I. V. Tokatly, and A. Rubio, “Orbital magneto-optical response of periodic insulators from first principles,” npj Comput. Mater. 5, 32 (2019).


{{<tutorial-footer>}}