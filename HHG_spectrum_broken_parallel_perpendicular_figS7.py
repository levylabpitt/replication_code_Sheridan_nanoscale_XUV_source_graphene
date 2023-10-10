#!/usr/bin/env python3
# -*- coding: utf-8 mode: python -*-
"""
HHG_spectrum_broken_parallel_perpendicular.py - 
"""
# Copyright (c) 2023, Davis Welakuh
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
# =============================================================================


import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# -----------------------------------------------------------------------------
# INPUT DATA
# -----------------------------------------------------------------------------

au_to_eV = 27.21138386
au_to_fs = 0.024189

max_Energy = 100
dd_Energy = 0.01

# -----------------------------------------------------------------------------
# LOAD THE REQUIRED DATA
# -----------------------------------------------------------------------------

# 001 --> 0.1 TW/cm^2 ; 002 --> 3.0 TW/cm^2


data_broken_ramp = np.loadtxt("data/01.parallel/02.broken-symmetry/"
                              "002.ramp_Ab_0.50/td.general/total_current")
data_broken_par_001 = np.loadtxt("data/01.parallel/02.broken-symmetry/"
                                 "402.IL_0.10_TWcm2_ramp_Ab_0.50/"
                                 "td.general/total_current")
data_broken_par_002 = np.loadtxt("data/01.parallel/02.broken-symmetry/"
                                 "602.IL_3.00_TWcm2_ramp_Ab_0.50/"
                                 "td.general/total_current")

time_t_broken = data_broken_par_001[:, 1]
time_t_fs = time_t_broken*au_to_fs

current_xy_broken_par_001 = \
    (data_broken_par_001[:, 2] - data_broken_ramp[:, 2]) + \
    (data_broken_par_001[:, 3] - data_broken_ramp[:, 3])
current_xy_broken_par_002 = \
    (data_broken_par_002[:, 2] - data_broken_ramp[:, 2]) + \
    (data_broken_par_002[:, 3] - data_broken_ramp[:, 3])


data_broken_per_001 = np.loadtxt("data/02.perperdicular/02.broken-symmetry/"
                                 "402.IL_0.10_TWcm2_ramp_Ab_0.50/"
                                 "td.general/total_current")
data_broken_per_002 = np.loadtxt("data/02.perperdicular/02.broken-symmetry/"
                                 "602.IL_3.00_TWcm2_ramp_Ab_0.50/"
                                 "td.general/total_current")

current_xy_broken_per_001 = \
    (data_broken_per_001[:, 2] - data_broken_ramp[:, 2]) + \
    (data_broken_per_001[:, 3] - data_broken_ramp[:, 3])
current_xy_broken_per_002 = \
    (data_broken_per_002[:, 2] - data_broken_ramp[:, 2]) + \
    (data_broken_per_002[:, 3] - data_broken_ramp[:, 3])

# -----------------------------------------------------------------------------
# DEFINE FUNCTIONS FOR POTENTIAL AND ABSORPTION SPECTRUM
# -----------------------------------------------------------------------------


def integration_trap(time, current, max_Energy, d_energy):
    
    au_to_fs = 0.024189; au_to_eV = 27.2113845
    time_dt = time[1] - time[0]
    e_step = min(4.13566733/((max(time))*au_to_fs)/au_to_eV, d_energy/au_to_eV)
    
    Energy_max = 4.13566733/(time_dt*au_to_fs)/au_to_eV/2
    Energy_max = min(Energy_max, max_Energy/au_to_eV)
    N_energy = int(np.ceil(Energy_max/e_step))
    num_iter = len(current)

    zz = np.zeros(N_energy)
    En = np.zeros(N_energy)
    # Re = np.zeros(N_energy)
    # Im = np.zeros(N_energy)

    tmp = np.arange(0, num_iter, 1)*time_dt
    
    for ie in range(N_energy-1):
        ee = ie*e_step
        
        etmp_sin = time_dt*sum(np.sin(tmp*ee)*current)
        etmp_sin = etmp_sin - time_dt*np.sin((num_iter-1)*ee*time_dt)*current[num_iter-1]/2
        
        etmp_cos = time_dt*sum(np.cos(tmp*ee)*current)
        etmp_cos = etmp_cos - time_dt*current[0]/2
        etmp_cos = etmp_cos - time_dt*np.cos((num_iter-1)*ee*time_dt)*current[num_iter-1]/2
         
        zz[ie] = etmp_sin**2 + etmp_cos**2
        En[ie] = au_to_eV*ee
        # Re[ie+1] = etmp_cos
        # Im[ie+1] = -etmp_sin

    return En, zz #, Re, Im


def integrate(energy, function):
    
    d_energy = energy[1] - energy[0]
    integral_result = function.sum()*d_energy
    
    return integral_result
    

#------------------------------------------------------------------------------
# COMPUTE HARMONIC SPECTRA
#------------------------------------------------------------------------------

mask_broken = np.hanning(len(current_xy_broken_par_002))

current_xy_broken_par_001 = current_xy_broken_par_001*mask_broken
current_xy_broken_par_002 = current_xy_broken_par_002*mask_broken

current_xy_broken_per_001 = current_xy_broken_per_001*mask_broken
current_xy_broken_per_002 = current_xy_broken_per_002*mask_broken

print()
print("... start computing HHG spectra")

energy, HHG_xy_broken_par_001 = \
    integration_trap(time_t_broken, current_xy_broken_par_001,
                     max_Energy, dd_Energy)
energy, HHG_xy_broken_par_002 = \
    integration_trap(time_t_broken, current_xy_broken_par_002,
                     max_Energy, dd_Energy)
    
energy, HHG_xy_broken_per_001 = \
    integration_trap(time_t_broken, current_xy_broken_per_001,
                     max_Energy, dd_Energy)
energy, HHG_xy_broken_per_002 = \
    integration_trap(time_t_broken, current_xy_broken_per_002,
                     max_Energy, dd_Energy)

# Time derivative    
HHG2_xy_broken_par_001 = HHG_xy_broken_par_001*energy**2
HHG2_xy_broken_par_002 = HHG_xy_broken_par_002*energy**2

HHG2_xy_broken_per_001 = HHG_xy_broken_per_001*energy**2
HHG2_xy_broken_per_002 = HHG_xy_broken_per_002*energy**2

print("...finished computing HHG spectra")

# equivalent to semilogy
log10_HHG2_xy_broken_par_001 = np.log10(HHG2_xy_broken_par_001)
log10_HHG2_xy_broken_par_002 = np.log10(HHG2_xy_broken_par_002)

log10_HHG2_xy_broken_per_001 = np.log10(HHG2_xy_broken_per_001)
log10_HHG2_xy_broken_per_002 = np.log10(HHG2_xy_broken_per_002)

# -----------------------------------------------------------------------------
# Data for plots
# -----------------------------------------------------------------------------

n_H = 11
omega = 1.54980
xposition = np.fromiter(((2.0*nn+1)*omega for nn in range(0, n_H)),
                        dtype=float)

top_ticks = xposition
top_label = np.fromiter((xposition[nn]/omega for nn in range(0, n_H)),
                        dtype=int)

bottom_label = np.arange(1, 16, 1)

# extract the maximum HHG peak intensities
min_int = 100
max_int = 235
max_shift = 150 # max_int
max_vals_HHG_par_001 = np.zeros(len(bottom_label))
max_vals_HHG_per_001 = np.zeros(len(bottom_label))
max_vals_HHG_par_002 = np.zeros(len(bottom_label))
max_vals_HHG_per_002 = np.zeros(len(bottom_label))


for ii in range(len(bottom_label)):
    max_vals_HHG_par_001[ii] = \
        integrate(energy[min_int:max_int],
                  log10_HHG2_xy_broken_par_001[min_int:max_int])
    max_vals_HHG_per_001[ii] = \
        integrate(energy[min_int:max_int],
                  log10_HHG2_xy_broken_per_001[min_int:max_int])
    max_vals_HHG_par_002[ii] = \
        integrate(energy[min_int:max_int],
                  log10_HHG2_xy_broken_par_002[min_int:max_int])
    max_vals_HHG_per_002[ii] = \
        integrate(energy[min_int:max_int],
                  log10_HHG2_xy_broken_per_002[min_int:max_int])
    min_int = max_int
    max_int = max_int + max_shift


# -----------------------------------------------------------------------------
# PLOTS START HERE
# -----------------------------------------------------------------------------


fig = plt.figure(figsize=(13, 10))
gs0 = gridspec.GridSpec(3, 2)
gs0.update(hspace = 0.40, wspace = 0.15, left = 0.14,
            right = 0.84, bottom = 0.10, top = 0.99)

gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[0:2, :],
                                        hspace=0.05)
# FIRST SUPLOT:
ax_00 = fig.add_subplot(gs00[0, :])
ax_00.plot(time_t_fs, current_xy_broken_par_001, "b", linestyle='-')
ax_00.plot(time_t_fs, current_xy_broken_per_001, "r", linestyle='--')
ax_00.legend(["Parallel", "Perpendicular"], loc='upper center', ncol=2, 
           bbox_to_anchor=(0.5, 1.20), prop={'size': 13})
ax_00.margins(x=0)
ax_00.set_ylabel(r'current, (a.u.)', fontsize=15)
ax_00.text(1.0,-0.004, r"$I_{L}=100 $GW/cm$^2$", fontsize=15)
ax_00.text(1.0, 0.004, "(a)", fontsize=14)
ax_00.set_ylim(-0.006, 0.006)
ax_00.xaxis.set_visible(False)
ax_00.get_yaxis().set_label_coords(-0.08, 0.5)
ax_00.tick_params(axis='both', which='major', labelsize=13)
# SECOND SUPLOT:
ax_10 = fig.add_subplot(gs00[1, :])
ax_10.plot(time_t_fs, current_xy_broken_par_002, "b", linestyle='-')
ax_10.plot(time_t_fs, current_xy_broken_per_002, "r", linestyle='--')
ax_10.margins(x=0)
ax_10.set_xlabel('Time (fs)', fontsize=15)
ax_10.set_ylabel(r'current, (a.u.)', fontsize=15)
ax_10.text(1.0,-0.08, r"$I_{L}=3 $TW/cm$^2$", fontsize=15)
ax_10.text(1.0, 0.08, "(b)", fontsize=14)
ax_10.set_ylim(-0.1, 0.1)
ax_10.get_yaxis().set_label_coords(-0.08, 0.5)
ax_10.tick_params(axis='both', which='major', labelsize=13)

gs01 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[2, :],
                                        hspace=0.5, wspace=0.15)
# THIRD SUPLOT:
ax_01 = fig.add_subplot(gs01[0, 0])
ax_01.semilogy(energy, HHG2_xy_broken_par_001, color='b')
ax_01.semilogy(energy, HHG2_xy_broken_per_001, color='r')
for xc in xposition:
    ax_01.axvline(x=xc, color='gray', linestyle='--', linewidth=0.9)
ax_01.set_xlabel(r'Photon energy (eV)', fontsize=15)
ax_01.set_ylabel(r"log(Intensity) (a.u.)", fontsize=15)
ax_01.set_xlim(-0.5, 35.0)
ax_01.text(30.0, 1e0, "(c)", fontsize=14)
ax_01.text(12.0, 1e0, r"$I_{L}=100 $GW/cm$^2$", fontsize=12)
ax_01.get_yaxis().set_label_coords(-0.17, 0.5)
ax_01.tick_params(axis='both', which='major', labelsize=13)
# top axis:
ax_01_1 = ax_01.twiny()
ax_01_1.set_xlim(ax_01.get_xlim())
ax_01_1.set_xticks(top_ticks)
ax_01_1.set_xticklabels(top_label)
ax_01_1.set_xlabel("Harmonic order", fontsize=13)
ax_01_1.tick_params(axis='both', which='major', labelsize=13)
# FOURTH SUPLOT:
ax_11 = fig.add_subplot(gs01[0, 1])
ax_11.semilogy(energy, HHG2_xy_broken_par_002, color='b')
ax_11.semilogy(energy, HHG2_xy_broken_per_002, color='r')
for xc in xposition:
    ax_11.axvline(x=xc, color='gray', linestyle='--', linewidth=0.9)
ax_11.set_xlabel(r'Photon energy (eV)', fontsize=15)
ax_11.set_xlim(-0.5, 35.0)
ax_11.text(30.0, 1e2, "(d)", fontsize=15)
ax_11.text(15.0, 1e2, r"$I_{L}=3 $TW/cm$^2$", fontsize=12)
ax_11.get_yaxis().set_label_coords(-0.17, 0.5)
ax_11.tick_params(axis='both', which='major', labelsize=13)
# top axis:
ax_11_1 = ax_11.twiny()
ax_11_1.set_xlim(ax_01.get_xlim())
ax_11_1.set_xticks(top_ticks)
ax_11_1.set_xticklabels(top_label)
ax_11_1.set_xlabel("Harmonic order", fontsize=13)
ax_11_1.tick_params(axis='both', which='major', labelsize=13)
plt.savefig('Figures/current_harm_spectrum_par_per.png', bbox_inches="tight", dpi=200)
plt.show()




fig = plt.figure(figsize=(10, 14))
gs0 = gridspec.GridSpec(4, 1)
gs0.update(hspace = 0.6, wspace = 0.5, left = 0.14,
           right = 0.84, bottom = 0.10, top = 0.93)

gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[0:2, :],
                                        hspace=0.05)
# FIRST SUPLOT:
ax_00 = fig.add_subplot(gs00[0, :])
ax_00.plot(time_t_fs, current_xy_broken_par_001, "b", linestyle='-')
ax_00.plot(time_t_fs, current_xy_broken_per_001, "r", linestyle='--')
ax_00.legend(["Parallel", "Perpendicular"], loc='upper center', ncol=2, 
           bbox_to_anchor=(0.5, 1.20), prop={'size': 13})
ax_00.margins(x=0)
ax_00.set_xlabel('Time (fs)', fontsize=15)
ax_00.set_ylabel(r'current, (a.u.)', fontsize=15)
ax_00.text(1.0,-0.004, r"$I_{L}=100 $GW/cm$^2$", fontsize=15)
ax_00.text(1.0, 0.004, "(a)", fontsize=14)
# ax_00.set_xlim(-0.1, 40.0)
ax_00.set_ylim(-0.006, 0.006)
ax_00.xaxis.set_visible(False)
ax_00.get_yaxis().set_label_coords(-0.12, 0.5)
ax_00.tick_params(axis='both', which='major', labelsize=13)
# SECOND SUPLOT:
ax_10 = fig.add_subplot(gs00[1, :])
ax_10.plot(time_t_fs, current_xy_broken_par_002, "b", linestyle='-')
ax_10.plot(time_t_fs, current_xy_broken_per_002, "r", linestyle='--')
ax_10.margins(x=0)
ax_10.set_xlabel('Time (fs)', fontsize=15)
ax_10.set_ylabel(r'current, (a.u.)', fontsize=15)
ax_10.text(1.0,-0.08, r"$I_{L}=3 $TW/cm$^2$", fontsize=15)
ax_10.text(1.0, 0.08, "(b)", fontsize=14)
# ax_10.set_xlim(-0.1, 40.0)
ax_10.set_ylim(-0.1, 0.1)
ax_10.get_yaxis().set_label_coords(-0.12, 0.5)
ax_10.tick_params(axis='both', which='major', labelsize=13)

gs01 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[2:, :],
                                        hspace=0.05, wspace=0.15)
# THIRD SUPLOT:
ax_01 = fig.add_subplot(gs01[0, :])
ax_01.semilogy(energy, HHG2_xy_broken_par_001, color='b')
ax_01.semilogy(energy, HHG2_xy_broken_per_001, color='r')
for xc in xposition:
    ax_01.axvline(x=xc, color='gray', linestyle='--', linewidth=0.9)
ax_01.set_xlabel(r'Photon energy (eV)', fontsize=15)
ax_01.set_ylabel(r"log(Intensity) (a.u.)", fontsize=15)
ax_01.set_xlim(-0.5, 35.0)
ax_01.text(30.0, 1e0, "(c)", fontsize=14)
ax_01.text(12.0, 1e0, r"$I_{L}=100 $GW/cm$^2$", fontsize=12)
ax_01.xaxis.set_visible(False)
ax_01.get_yaxis().set_label_coords(-0.1, 0.5)
ax_01.tick_params(axis='both', which='major', labelsize=13)
# top axis:
ax_01_1 = ax_01.twiny()
ax_01_1.set_xlim(ax_01.get_xlim())
ax_01_1.set_xticks(top_ticks)
ax_01_1.set_xticklabels(top_label)
ax_01_1.set_xlabel("Harmonic order", fontsize=13)
ax_01_1.tick_params(axis='both', which='major', labelsize=13)
# FOURTH SUPLOT:
ax_11 = fig.add_subplot(gs01[1, :])
ax_11.semilogy(energy, HHG2_xy_broken_par_002, color='b')
ax_11.semilogy(energy, HHG2_xy_broken_per_002, color='r')
for xc in xposition:
    ax_11.axvline(x=xc, color='gray', linestyle='--', linewidth=0.9)
ax_11.set_xlabel(r'Photon energy (eV)', fontsize=15)
ax_11.set_xlim(-0.5, 35.0)
ax_11.set_ylabel(r"log(Intensity) (a.u.)", fontsize=15)
ax_11.text(30.0, 1e2, "(d)", fontsize=15)
ax_11.text(15.0, 1e2, r"$I_{L}=3 $TW/cm$^2$", fontsize=12)
ax_11.get_yaxis().set_label_coords(-0.1, 0.5)
ax_11.tick_params(axis='both', which='major', labelsize=13)
plt.savefig('Figures/current_harm_spectrum_par_per_1.png', bbox_inches="tight", dpi=200)
plt.show()




# defining the subplots
gs = gridspec.GridSpec(2, 1)
gs.update(hspace=0.05)

fig = plt.figure(figsize=(7, 6))
# First plot: 
ax = fig.add_subplot(gs[0])
ax.plot(bottom_label, max_vals_HHG_par_002, color='blue', marker='s')
ax.plot(bottom_label, max_vals_HHG_per_002, color='red', marker='o')
ax.legend(["Parallel", "Perpendicular"], loc='upper center', ncol=2, 
           bbox_to_anchor=(0.5, 1.25), prop={'size': 13})
ax.set_xticks(bottom_label)
ax.text(10.0, 0.0, r"$I_{L}=3$ TW/cm$^2$", fontsize=15)
ax.xaxis.set_visible(False)
ax.get_yaxis().set_label_coords(-0.1, 0.5)
ax.tick_params(axis='both', which='major', labelsize=14)
# Second plot: 
ax = fig.add_subplot(gs[1])
ax.plot(bottom_label, max_vals_HHG_par_001, color='blue', marker='s')
ax.plot(bottom_label, max_vals_HHG_per_001, color='red', marker='o')
ax.set_xticks(bottom_label)
ax.set_xlabel("Harmonic order", fontsize=13)
ax.text(10.0, -3.5, r"$I_{L}=100$ GW/cm$^2$", fontsize=15)
ax.get_yaxis().set_label_coords(-0.1, 0.5)
ax.tick_params(axis='both', which='major', labelsize=14)
fig.text(0.01, 0.5, r"Integrated Intensity (a.u.)",
         va="center", rotation="vertical", fontsize=18)
plt.savefig('Figures/integrated_HHG_intensity_par_per.png', bbox_inches="tight", dpi=200)
plt.show()


