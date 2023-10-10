---
title: "Use of symmetries in optical spectra from time-propagation"
#tags: ["Advanced", "Time-dependent", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "oct-propagation_spectrum"]
series: "Tutorials"
tutorials: "Optical Response"
difficulties: "advanced"
theories: "DFT"
calculation_modes: "Time-dependent"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Optical Absorption"
utilities: "oct-propagation_spectrum"
weight: 6
description: "Reduce the number of time-propagations required to compute the absorption cross-section."
---


In this tutorial we will see how the spatial symmetries of a molecule can be used to reduce the number of time-propagations required to compute the absorption cross-section.

## Introduction

The dynamic polarizability (related to optical absorption cross-section via $\sigma = \frac{4 \pi \omega}{c} \mathrm{Im} \alpha $) is, in its most general form, a 3x3 tensor. The reason is that we can shine light on the system polarized in any of the three Cartesian axes, and for each of these three cases measure how the dipole of the molecule oscillates along the three Cartesian axes. This usually means that to obtain the full dynamic polarizability of the molecule we usually need to apply 3 different perturbations along $x, y, z\,$, by setting {{< variable "TDPolarizationDirection" >}} to 1, 2, or 3.

However, if the molecule has some symmetries, it is in general possible to reduce the total number of calculations required to obtain the tensor from 3 to 2, or even 1.[^footnote-1]

To use this formalism in {{< octopus >}}, you need to supply some extra information. The most important thing that the code requires is the information about equivalent axes, that is, directions that are related through some symmetry transformation. Using these axes, we construct a reference frame and specify it with the {{< variable "TDPolarization" >}} block. Note that these axes need not be orthogonal, but they must be linearly-independent. The {{< variable "TDPolarizationEquivAxes" >}} tells the code how many equivalent axes there are. Ideally, the reference frame should be chosen to maximize the number of equivalent axes. When using three equivalent axes, an extra input variable, {{< variable "TDPolarizationWprime" >}}, is also required.

Let us give a couple of examples, which should make all these concepts easier to understand.

## Methane

As seen in previous tutorials, the methane molecule has $T_d$ symmetry. This means that it is trivial to find three linearly-independent equivalent axes such that only one time-propagation is needed to obtain the whole tensor. As it happens, we can use the usual $x$, $y$, and $z$ directions, with all of them being equivalent (check the [symmetry operations](https://en.wikipedia.org/wiki/Symmetry_operation) of the $T_d$ [point group](https://en.wikipedia.org/wiki/Molecular_symmetry) if you are not convinced). Therefore we can perform just one propagation with the perturbation along the $x$ direction adding the following to the input file used previously: 

{{< code-block >}}
#include_input_snippet doc/tutorials/optical_response/use_of_symmetries_in_optical_spectra_from_time-propagation/2.td/inp polarization
{{< /code-block >}}

{{% expand "complete input file" %}}
{{< code-block >}}
#include_input doc/tutorials/optical_response/use_of_symmetries_in_optical_spectra_from_time-propagation/2.td/inp
{{< /code-block >}}
{{% /expand %}}


Note that we had omitted the blocks {{< variable "TDPolarization" >}} and {{< variable "TDPolarizationWprime" >}} in the previous tutorials, as these are their default 

Once the time-propagation is finished, you will find, as usual, a {{< file "td.general/multipoles" >}} file. This time the file contains all the necessary information for {{< octopus >}} to compute the full tensor, so running the {{< file "oct-propagation_spectrum" >}} utility will produce two files: {{< file "cross_section_vector.1" >}} and {{< file "cross_section_tensor" >}}. The former should look like this (assuming you used the converged grid parameters from the
{{< tutorial "Convergence of the optical spectra" >}} tutorial):

{{< code-block >}}
#include_file doc/tutorials/optical_response/use_of_symmetries_in_optical_spectra_from_time-propagation/3.spectrum/cross_section_tensor-head.txt
...
{{< /code-block >}}

Try comparing the spectrum for each component of the $\sigma$ tensor.

## Linear molecule

Now let us look at a linear molecule. In this case, you might think that we need two calculations to obtain the whole tensor, one for the direction along the axis of the molecule, and another for the axis perpendicular to the molecule. The fact is that we need only one, in a specially chosen direction, so that our field has components both along the axis of the molecule and perpendicular to it. Let us assume that the axis of the molecule is oriented along the $x\,$ axis. Then we can use

{{< code-block >}}
 %{{< variable "TDPolarization" >}}
  1/sqrt(2) | -1/sqrt(2) | 0
  1/sqrt(2) |  1/sqrt(2) | 0
  1/sqrt(2) |  0         | 1/sqrt(2)
 %
 {{< variable "TDPolarizationDirection" >}} = 1
 {{< variable "TDPolarizationEquivAxes" >}} = 3
 %{{< variable "TDPolarizationWprime" >}}
  1/sqrt(2) |  0         | 1/sqrt(2)
 %
{{< /code-block >}}

You should try to convince yourself that the three axes are indeed equivalent and linearly independent. The first and second axes are connected through a simple reflection in the $xz$ plane, transforming the $y$ coordinate from $-1/\sqrt{2}$ into $1/\sqrt{2}$. {{< variable "TDPolarizationWprime" >}} should be set to the result obtained by applying the inverse symmetry operation to the third axis. This actually leaves the third axis unchanged.

## Planar molecule

Finally, let us look at a general planar molecule (in the $xy$ plane). In principle we need only two calculations (that is reduced to one if more symmetries are present like, ''e.g.'', in benzene). In this case we chose one of the polarization axes on the plane and the other two rotated 45 degrees:

{{< code-block >}}
 %{{< variable "TDPolarization" >}}
  1/sqrt(2) | 0 | 1/sqrt(2)
  1/sqrt(2) | 0 |-1/sqrt(2)
  0         | 1 | 0
 %
 
 {{< variable "TDPolarizationEquivAxes" >}} = 2
{{< /code-block >}}

In this case, we need two runs, one for {{< variable "TDPolarizationDirection" >}} equal to 1, and another for it equal to 3. Note that if there are less than 3 equivalent axes, {{< variable "TDPolarizationWprime" >}} is irrelevant.

[^footnote-1]: {{< article title="On the use of Neumann's principle for the calculation of the polarizability tensor of nanostructures" authors="M.J.T. Oliveira, A. Castro, M.A.L. Marques, and A. Rubio" journal="J. Nanoscience and Nanotechnology" volume="8" pages="1-7" year="2008" doi="10.1166/jnn.2008.142" arxiv="0710.2624v1" >}}


{{< tutorial-footer >}}
