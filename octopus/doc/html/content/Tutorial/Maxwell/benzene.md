---
Title: " Coupled Maxwell-TDDFT propagation"
series: "Tutorials"
tutorials: "Maxwell"
Weight: 20
draft: yes
---

### Coupled Maxwell-TDDFT propagation 

{{< notice warning >}}
This needs to be removed as Maxwell-matter coupling is not yet working.{{< /notice >}}

#### Ground state of benzene


To simulate a coupled Maxwell-Kohn-Sham system, we have to obtain first an initial Kohn-Sham system. Here in our case a groundstate of benzene. 

```
# ----- Calculation node and parallelization ----------------------------------------------------

CalculationMode = gs
ParDomains = auto
ParStates = no

# ----- Matter box variables --------------------------------------------------------------------

lsize_ma = 10.0
dx_ma = 0.5

Dimensions = 3
BoxShape = parallelepiped

%Lsize
lsize_ma | lsize_ma | lsize_ma
%

%Spacing
dx_ma | dx_ma | dx_ma
%

# ----- Species variables -----------------------------------------------------------------------

XYZCoordinates = "geometry_z_dir.bohr.xyz"

# ----- Hatter calculation variables

Extrastates = 2
XCFunctional = lda_x + lda_c_gl
EigenSolver = cg
EigenSolverTolerance = 1.00e-12
EigensolverMaxIter = 50
CoanelDens = 1.00e-11
MaximumIter = 1000
ConvForce = 0.0
SmearingFunction = fermi_dirac
Smearing = 0.001
MixingScheme = broyden
Mixing = 0.02
LCAOStart = lcao_states
```
#### Coupled dynamics of benzene and a cosinoidal EM pulse 

** Folder [/02_maxwell_ks_propagation/01_ed_coupling/01_no_backreaction/](./02_maxwell_ks_propagation/01_ed_coupling/01_no_backreaction/) **

A cosinoidal laser pulse hits a molecule in groundstate usind electric dipole approximation without any matter to Maxwell back-reaction. 

```
# ----- Calculation mode and parallelization ----------------------------------------------------

CalculationMode = maxwell_ks
ParDomains = auto
ParStates = no
MaxwellParDomains = auto
MaxwellParStates = no

# ----- Matter box variables --------------------------------------------------------------------

lsize_ma = 10.0
dx_ma = 0.5
Dimensions = 3
BoxShape = parallelepiped

%Lsize
lsize_ma | lsize_ma | lsize_ma
%

%Spacing
dx_ma | dx_ma | dx_ma
%

# ----- Maxwell box variables -------------------------------------------------------------------

# free maxwell box limit of 10.0 plus for the 5.0 absorbing pml boundaries 
# plus 2.0 for the incident wave boundaries with der_order = 4 times dx_mx
lsize_mx = 17.0 
dx_mx = 0.5
MaxwellDimensions = 3
MaxwellBoxShape = parallelepiped

%MaxwellLsize
lsize_mx | lsize_mx | lsize_mx
%

%MaxwellSpacing
dx_mx | dx_mx | dx_mx
%

# ----- Species variables -----------------------------------------------------------------------

XYZCoordinates = "geometry_z_dir.bohr.xyz"

# ----- Matter calculation variables ------------------------------------------------------------

Extrastates = 2
XCFunctional = lda_x + lda_c_gl
EigenSolver = cg
EigenSolverTolerance = 1.00e-12
EigensolverMaxIter = 50
CoanelDens = 1.00e-12
MaximumIter = 1000
DerivativesStencil = stencil_starplus
DerivativesOrder = 4
TDExpOrder = 4
AbsorbingBoundaries = not_absorbing

# ----- Maxwell calculation variables -----------------------------------------------------------

MaxwellHamiltonianOperator = faraday_ampere
MaxwellTDOperatorMethod = maxwell_op_fd

MaxwellTDPropagator = maxwell_etrs
MatterToMaxwellCoupling = no
MaxwellToMatterCoupling = yes

MaxwellCouplingOrder = electric_dipole_coupling
MaxwellTransFieldCalculationMethod = trans_field_poisson
MaxwellPoissonSolver = isf
MaxwellPoissonSolverBoundaries = multipole

MaxwellDerivativesStencil = stencil_starplus
MaxwellDerivativesOrder = 4
MaxwellTDExpOrder = 4

%MaxwellBoundaryConditions
maxwell_plane_waves | maxwell_plane_waves | maxwell_plane_waves
%

%MaxwellAbsorbingBoundaries
cpml | cpml | cpml
%

MaxwellABPMLWidth = 5.0
MaxwellABPMLKappaMax = 1.0
MaxwellABPMLAlphaMax = 1.0 
MaxwellABPMLPower = 2.0
MaxwellABPMLReflectionError = 1.0e-16

# ----- Output variables ------------------------------------------------------------------------

OutputFormat = plane_x + plane_y + plane_z + vtk + xyz + axis_x

# ----- Matter output variables -----------------------------------------------------------------

Output = potential + density + current + geometry + forces + elf

OutputInterval = 50
TDOutput = energy + multipoles + laser + geometry

# ----- Maxwell output variables ----------------------------------------------------------------

MaxwellOutput = maxwell electric field + maxwell magnetic field + maxwell_energy_density + maxwell_trans_§lectric_field #_(has to be written in one line)

MaxwellOutputInterval = 1

MaxwellTDOutput = maxwell_energy + maxwell_fields

%MaxwellFieldsCoordinate
0.00 | 0.00 | 0.00
%

# ----- Time step variables ---------------------------------------------------------------------
TDTimeStep = 0.002
TDMaxSteps = 200

TDEnergyUpdateIter = 1
MaxwellTDIntervalSteps = 1
TDMaxwellTDRelaxationSteps = 0
TDMaxwellKSRelaxationSteps = 0

MaxwellTDETRSApprox = no
CurrentPropagationTest = no

# ----- Maxwell field variables -----------------------------------------------------------------

lambda = 10.0
omega = 2 * pi * c / lambda
kx = omega / c
E2 = 0.05
pw = 10.0
ps = - 5 * 5.0

%UserDefinedMaxwellIncidentWaves
plane_waves_mx_function | 0 | 0 | E2 | "plane_waves_function" | plane_wave
%

%MaxwellFunctions
"plane_waves_function" | mxf_cosinoidal_wave | kx | 0 | 0 | ps | 0 | 0 | pw
%

# cosinoidal pulse
%UserDefinedInitialMaxwellStates
3 | formula | electric_field | “ Ez*cos(kx*(x-ps))*cos(pi/2*(x-ps-2*pw)/pw) * step(pw-abs((kx*x-kx*ps)/kx^2)) “ # (in one line)
2 | formula | magnetic_field | “ -1/c*Ez*cos(kx*(x-ps))*cos(pi/2*(x-ps-2*pw)/pw) * step(pw-abs((kx*x-kx*ps)/kx^2)) " # (in one line)
%
```

To explore the effect of the back reaction of the field (considering the fully coupled Maxwell-matter interaction) and the effect of higher order coupling terms (magnetic dipole and electric quadrupole) run the input files in folders {{< file "02_maxwell_ks_propagation/01_ed_coupling/02_fully_coupled/" >}} , {{< file "2_maxwell_ks_propagation/02_ed_md_eq_coupling/01_no_backreaction/" >}} and {{< file "02_maxwell_ks_propagation/02_ed_md_eq_coupling/01_no_backreaction/" >}} . 

{{< tutorial-footer >}}

