---
title: "Visualization"
#series: "Manual"
Weight: 102
---


Every given number of time iterations, or after ground-state calculations, some of the functions that characterise the system may be written to disk so that they may be analized. Files are written within {{< file "static/" >}} output directory after the self-consistent field, or within {{< file "td.x/" >}} directories, during evolution, where “x” stands for the iteration number at which each write is done. 

The function that you want to plot is selected by the {{< variable "Output" >}} variable and the output format is chosen by the 
{{< variable "OutputFormat" >}}.

### dx

This is an OpenDX network, aimed at the visualization of wavefunctions. To be able to use it, you need to have properly installed the [OpenDX](https://opendx.org) program, as well as the Chemistry extensions developed at the Cornell Theory Center. Unfortunately, since this software is old and unmaintained, you are likely to have trouble finding and installing it.

Once these are working, you may follow the {{< versioned-link "tutorial/basics/visualization/#benzene-molecule" "tutorial for the benzene molecule" >}}
### XCrySDen

Atomic coordinates (finite or periodic), forces, and functions on a grid can be plotted with the free program [https://www.xcrysden.org/ XCrySDen]. Its XSF format also can be read by [V_sim](https://inac.cea.fr/sp2m/L_Sim/V_Sim/index.en.html) and [Vesta](https://www.geocities.jp/kmo_mma/crystal/en/vesta.html). Beware, these all probably assume that your output is in Angstrom units (according to the [specification](https://www.xcrysden.org/doc/XSF.html)), so use UnitsOutput = eV_Angstrom, or your data will be misinterpreted by the visualization software.

### CUBE

The Gaussian cube format (see https://gaussian.com/cubegen/, https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html and https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html for a more detailed description of the format) can be output, and can be read by VMD, XCrysDen, Avogadro, and other software. Note that CUBE files are always in atomic units, so the UnitsOutput input option will be ignored.

### PDB

Everything is supposed to be in Angstroms: https://deposit.rcsb.org/adit/docs/pdb_atom_format.html

### XYZ

Generally considered to be in Angstroms: https://openbabel.org/wiki/XYZ_%28format%29, https://en.wikipedia.org/wiki/XYZ_file_format, https://www.molpro.net/info/2012.1/doc/manual/node100.html, https://departments.icmab.es/leem/siesta/Documentation/Manuals/siesta-3.1-manual/node32.html

{{< manual-foot prev="Manual:Geometry Optimization" next="Manual:Advanced ways of running Octopus" >}}
---------------------------------------------
