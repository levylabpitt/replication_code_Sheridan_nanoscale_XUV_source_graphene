#!/usr/bin/env python3
# -*- coding: utf-8 mode: python -*-
"""
HHG_spectrum_intact_parallel_figS5.py - 
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

data_intact_001 = np.loadtxt("data/01.parallel/01.intact-symmetry/"
                             "004.IL_0.10_TWcm2/td.general/total_current")
data_intact_002 = np.loadtxt("data/01.parallel/01.intact-symmetry/"
                             "009.IL_3.00_TWcm2/td.general/total_current")

data_001_laser = np.loadtxt("data/01.parallel/01.intact-symmetry/"
                            "004.IL_0.10_TWcm2/td.general/laser")
data_002_laser = np.loadtxt("data/01.parallel/01.intact-symmetry/"
                            "009.IL_3.00_TWcm2/td.general/laser")

time_t_intact = data_intact_001[:, 1]
time_t_fs = time_t_intact*au_to_fs

current_xy_intact_001 = data_intact_001[:, 2] + data_intact_001[:, 3]
current_xy_intact_002 = data_intact_002[:, 2] + data_intact_002[:, 3]


A_field_001 = data_001_laser[:, 2]
A_field_002 = data_002_laser[:, 2] 


# -----------------------------------------------------------------------------
# DEFINE FUNCTIONS FOR POTENTIAL AND ABSORPTION SPECTRUM
# -----------------------------------------------------------------------------

def my_mask_fucntion(signal, last_vals):
    mask_func = np.ones(len(signal))
    aa = np.arange(0, last_vals+1, 1)
    bb = np.flipud((np.cos(aa/(last_vals)*np.pi/2))**2)
    mask_func[-len(bb):] = np.flipud(bb)
    return mask_func


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

#------------------------------------------------------------------------------
# COMPUTE HARMONIC SPECTRa
#------------------------------------------------------------------------------

mask = np.hanning(len(current_xy_intact_001))

current_xy_intact_001 = current_xy_intact_001*mask
current_xy_intact_002 = current_xy_intact_002*mask

print()
print("... start computing HHG spectra")

energy, HHG_xy_intact_001 = \
    integration_trap(time_t_intact, current_xy_intact_001,
                     max_Energy, dd_Energy)
energy, HHG_xy_intact_002 = \
    integration_trap(time_t_intact, current_xy_intact_002,
                     max_Energy, dd_Energy)

# Time derivative
HHG2_xy_intact_001 = HHG_xy_intact_001*energy**2
HHG2_xy_intact_002 = HHG_xy_intact_002*energy**2

print("...finished computing HHG spectra")

# -----------------------------------------------------------------------------
# Data for plots
# -----------------------------------------------------------------------------

n_H = 11
omega = 1.54980
xposition = np.fromiter(((2.0*nn+1)*omega for nn in range(0, n_H)),  dtype=float)

top_ticks = xposition
top_label = np.fromiter((xposition[nn]/omega for nn in range(0, n_H)),  dtype=int)

# -----------------------------------------------------------------------------
# PLOTS START HERE
# -----------------------------------------------------------------------------


fig = plt.figure(figsize=(13, 8))
gs1 = gridspec.GridSpec(2, 2)
gs1.update(hspace = 0.05, wspace = 0.5, left = 0.14,
            right = 0.84, bottom = 0.10, top = 0.93)

# FIRST SUPLOT:
ax_00 = plt.subplot(gs1[0, 0])
ax_00.plot(time_t_fs, A_field_001, 'r', linestyle='-')
ax_00.margins(x=0)
ax_00.set_xlabel('Time (fs)', fontsize=15)
ax_00.set_ylabel(r'$A_{L}(t)$ (a.u.)', fontsize=15, color="r")
ax_00.text(22.0, 4.1, r"$I_{L}=100 $GW/cm$^2$", color='b', fontsize=12)
ax_00.text(1.0, 4.1, "(a)", fontsize=14)
ax_00.set_xlim(-0.1, 40.0)
ax_00.set_ylim(-5.0, 5.0)
ax_00.xaxis.set_visible(False)
ax_00.get_yaxis().set_label_coords(-0.13, 0.5)
ax_00.tick_params(axis='both', which='major', labelsize=13)
# First figure: right-hand side
ax_00_1 = ax_00.twinx()
ax_00_1.plot(time_t_fs, current_xy_intact_001, "b", linestyle='-')
ax_00_1.margins(x=0)
ax_00_1.set_ylabel('current (a.u.)', fontsize=15, color="b")
# SECOND SUPLOT:
ax_01 = plt.subplot(gs1[0, 1])
ax_01.semilogy(energy, HHG2_xy_intact_001, color='blue')
for xc in xposition:
    ax_01.axvline(x=xc, color='gray', linestyle='--', linewidth=0.9)
ax_01.set_xlabel(r'Photon energy (eV)', fontsize=15)
ax_01.set_ylabel(r"log(Intensity) (a.u.)", fontsize=15)
ax_01.set_xlim(-0.5, 35.0)
ax_01.text(30.0, 1e0, "(b)", fontsize=14)
ax_01.text(12.0, 1e0, r"$I_{L}=100 $GW/cm$^2$", color='b', fontsize=12)
ax_01.xaxis.set_visible(False)
ax_01.get_yaxis().set_label_coords(-0.17, 0.5)
ax_01.tick_params(axis='both', which='major', labelsize=13)
# top axis:
ax_01_1 = ax_01.twiny()
ax_01_1.set_xlim(ax_01.get_xlim())
ax_01_1.set_xticks(top_ticks)
ax_01_1.set_xticklabels(top_label)
ax_01_1.set_xlabel("Harmonic order", fontsize=13)
ax_01_1.tick_params(axis='both', which='major', labelsize=13)
# THIRD SUPLOT:
ax_10 = plt.subplot(gs1[1, 0])
ax_10.plot(time_t_fs, A_field_002, 'r', linestyle='-')
ax_10.margins(x=0)
ax_10.set_xlabel('Time (fs)', fontsize=15)
ax_10.set_ylabel(r'$A_{L}(t)$ (a.u.)', fontsize=15, color="r")
ax_10.text(25.0, 21.0, r"$I_{L}=3 $TW/cm$^2$", color='b', fontsize=12)
ax_10.text(1.0, 21.0, "(c)", fontsize=14)
ax_10.set_xlim(-0.1, 40.0)
ax_10.set_ylim(-25.0, 25.0)
ax_10.get_yaxis().set_label_coords(-0.13, 0.5)
ax_10.tick_params(axis='both', which='major', labelsize=13)
# First figure: right-hand side
ax_10_1 = ax_10.twinx()
ax_10_1.plot(time_t_fs, current_xy_intact_002, "b", linestyle='-')
ax_10_1.margins(x=0)
ax_10_1.set_ylabel('current (a.u.)', fontsize=15, color="b")
# SECOND SUPLOT:
ax_11 = plt.subplot(gs1[1, 1])
ax_11.semilogy(energy, HHG2_xy_intact_002, color='blue')
for xc in xposition:
    ax_11.axvline(x=xc, color='gray', linestyle='--', linewidth=0.9)
ax_11.set_xlabel(r'Photon energy (eV)', fontsize=15)
ax_11.set_ylabel(r"log(Intensity) (a.u.)", fontsize=15)
ax_11.set_xlim(-0.5, 35.0)
ax_11.text(30.0, 1e2, "(d)", fontsize=15)
ax_11.text(15.0, 1e2, r"$I_{L}=3 $TW/cm$^2$", color='b', fontsize=12)
ax_11.get_yaxis().set_label_coords(-0.17, 0.5)
ax_11.tick_params(axis='both', which='major', labelsize=13)
plt.savefig('Figures/current_harm_spectrum_odd.png', bbox_inches="tight", dpi=200)
plt.show()



fig = plt.figure(figsize=(10, 16))
gs0 = gridspec.GridSpec(4, 1)
gs0.update(hspace = 0.6, wspace = 0.5, left = 0.14,
           right = 0.84, bottom = 0.10, top = 0.93)

gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[0:2, :],
                                        hspace=0.05)
# FIRST SUPLOT:
ax_00 = plt.subplot(gs00[0, :])
ax_00.plot(time_t_fs, A_field_001, 'r', linestyle='-')
ax_00.margins(x=0)
ax_00.set_xlabel('Time (fs)', fontsize=15)
ax_00.set_ylabel(r'$A_{L}(t)$ (a.u.)', fontsize=15, color="r")
ax_00.text(30.0, 3.9, r"$I_{L}=100 $GW/cm$^2$", color='b', fontsize=14)
ax_00.text(1.0, 3.9, "(a)", fontsize=14)
ax_00.set_xlim(-0.1, 40.0)
ax_00.set_ylim(-5.0, 5.0)
ax_00.xaxis.set_visible(False)
ax_00.get_yaxis().set_label_coords(-0.1, 0.5)
ax_00.tick_params(axis='both', which='major', labelsize=13)
# First figure: right-hand side
ax_00_1 = ax_00.twinx()
ax_00_1.plot(time_t_fs, current_xy_intact_001, "b", linestyle='-')
ax_00_1.margins(x=0)
ax_00_1.set_ylabel('current (a.u.)', fontsize=15, color="b")
# SECOND SUPLOT:
ax_10 = plt.subplot(gs00[1, :])
ax_10.plot(time_t_fs, A_field_002, 'r', linestyle='-')
ax_10.margins(x=0)
ax_10.set_xlabel('Time (fs)', fontsize=15)
ax_10.set_ylabel(r'$A_{L}(t)$ (a.u.)', fontsize=15, color="r")
ax_10.text(30.0, 19.0, r"$I_{L}=3 $TW/cm$^2$", color='b', fontsize=14)
ax_10.text(1.0, 19.0, "(b)", fontsize=14)
ax_10.set_xlim(-0.1, 40.0)
ax_10.set_ylim(-25.0, 25.0)
ax_10.get_yaxis().set_label_coords(-0.1, 0.5)
ax_10.tick_params(axis='both', which='major', labelsize=13)
# Second figure: right-hand side
ax_10_1 = ax_10.twinx()
ax_10_1.plot(time_t_fs, current_xy_intact_002, "b", linestyle='-')
ax_10_1.margins(x=0)
ax_10_1.set_ylabel('current (a.u.)', fontsize=15, color="b")

gs01 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[2:, :],
                                        hspace=0.05, wspace=0.15)
# THIRD SUPLOT:
ax_20 = plt.subplot(gs01[0, :])
ax_20.semilogy(energy, HHG2_xy_intact_001, color='blue')
for xc in xposition:
    ax_20.axvline(x=xc, color='gray', linestyle='--', linewidth=0.9)
ax_20.set_xlabel(r'Photon energy (eV)', fontsize=15)
ax_20.set_ylabel(r"log(Intensity) (a.u.)", fontsize=15)
ax_20.set_xlim(-0.5, 35.0)
ax_20.text(30.0, 1e0, "(c)", fontsize=14)
ax_20.text(15.0, 1e0, r"$I_{L}=100 $GW/cm$^2$", color='b', fontsize=14)
ax_20.xaxis.set_visible(False)
ax_20.get_yaxis().set_label_coords(-0.13, 0.5)
ax_20.tick_params(axis='both', which='major', labelsize=13)
# top axis:
ax_20_1 = ax_20.twiny()
ax_20_1.set_xlim(ax_20.get_xlim())
ax_20_1.set_xticks(top_ticks)
ax_20_1.set_xticklabels(top_label)
ax_20_1.set_xlabel("Harmonic order", fontsize=13)
ax_20_1.tick_params(axis='both', which='major', labelsize=13)
# FOURTH SUPLOT:
ax_30 = plt.subplot(gs01[1, :])
ax_30.semilogy(energy, HHG2_xy_intact_002, color='blue')
for xc in xposition:
    ax_30.axvline(x=xc, color='gray', linestyle='--', linewidth=0.9)
ax_30.set_xlabel(r'Photon energy (eV)', fontsize=15)
ax_30.set_ylabel(r"log(Intensity) (a.u.)", fontsize=15)
ax_30.set_xlim(-0.5, 35.0)
ax_30.text(30.0, 1e2, "(d)", fontsize=15)
ax_30.text(15.0, 1e2, r"$I_{L}=3 $TW/cm$^2$", color='b', fontsize=14)
ax_30.get_yaxis().set_label_coords(-0.13, 0.5)
ax_30.tick_params(axis='both', which='major', labelsize=13)
plt.savefig('Figures/current_harm_spectrum_odd_1.png', bbox_inches="tight", dpi=200)
plt.show()


