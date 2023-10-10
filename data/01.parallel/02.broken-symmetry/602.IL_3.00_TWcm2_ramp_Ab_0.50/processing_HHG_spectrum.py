#!/usr/bin/env python3
# -*- coding: utf-8 mode: python -*-
"""
processing_HHG_spectrum.py - 
"""
# Copyright (c) 2022, Davis Welakuh
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


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


#------------------------------------------------------------------------------
# INPUT DATA
#------------------------------------------------------------------------------


last_vals = 150
au_to_fs = 0.024189
au_to_eV = 27.2113845

max_Energy = 100
dd_Energy = 0.01

mn_val = 2000
mx_val = 11000

#mn_val = 0
#mx_val = 13238

#------------------------------------------------------------------------------
# LOAD DATA
#------------------------------------------------------------------------------

data_current = np.loadtxt("td.general/total_current")
data_current_1 = np.loadtxt("../402.I0-3TWcm2-ramp-A0-00.50/td.general/total_current")

time_t = data_current[:, 1]
current_x = data_current[:, 2] # - data_current_1[:, 2]
current_y = data_current[:, 3] # - data_current_1[:, 3]
current_z = data_current[:, 4] # - data_current_1[:, 4]

current = current_x + current_y + current_z

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

mask = 1.0
#mask = my_mask_fucntion(time_t, last_vals)
mask = np.hanning(len(current_x[mn_val:mx_val]))

current_x = current_x[mn_val:mx_val]*mask
current_y = current_y[mn_val:mx_val]*mask
current_z = current_z[mn_val:mx_val]*mask

energy, HHG_x = integration_trap(time_t, current_x, max_Energy, dd_Energy)
energy, HHG_y = integration_trap(time_t, current_y, max_Energy, dd_Energy)
energy, HHG_z = integration_trap(time_t, current_z, max_Energy, dd_Energy)


# Time derivative
HHG2_x = HHG_x*energy**2
HHG2_y = HHG_y*energy**2
HHG2_z = HHG_z*energy**2

HHG2 = HHG2_x + HHG2_y + HHG2_z

#np.savetxt("data_HHG/HHG_spectrum.txt", np.c_[energy, HHG2_x])


# -----------------------------------------------------------------------------
# PLOTS START HERE
# -----------------------------------------------------------------------------


n_H = 11
omega = 1.54980
xposition = np.fromiter(((2.0*nn+1)*omega for nn in range(0, n_H)),  dtype=float)

top_ticks = xposition
top_label = np.fromiter((xposition[nn]/omega for nn in range(0, n_H)),  dtype=int)


# defining the subplots
gs = gridspec.GridSpec(2, 1)
gs.update(hspace=0.10)

fig = plt.figure(figsize=(6,7))
# bottom plot:
ax = fig.add_subplot(gs[0])
ax.semilogy(energy, HHG2, color='dodgerblue')
for xc in xposition:
    ax.axvline(x=xc, color='gray', linestyle='--', linewidth=0.9)
ax.margins(x=0)
ax.set_xlabel(r'Photon energy (eV)', fontsize=15)
ax.set_ylabel(r"log(Intensity) (a.u.)", fontsize=15)
ax.set_xlim(-0.5, 35.0)
#ax.set_ylim(-0.5, 5.0)
#ax.xaxis.set_visible(False)
ax.get_yaxis().set_label_coords(-0.12, 0.5)
ax.tick_params(axis='both', which='major', labelsize=13)
# top plot:
ax1 = ax.twiny()
ax1.set_xlim(ax.get_xlim())
ax1.set_xticks(top_ticks)
ax1.set_xticklabels(top_label)
ax1.set_xlabel("Harmonic order", fontsize=13)
ax1.tick_params(axis='both', which='major', labelsize=13)
plt.savefig('graphene_HHG_spectrum.png', bbox_inches="tight", dpi=200)
plt.show()



