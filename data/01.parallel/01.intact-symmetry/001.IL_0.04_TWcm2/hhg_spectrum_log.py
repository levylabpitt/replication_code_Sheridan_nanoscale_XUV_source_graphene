#!/usr/bin/env python
# -*- coding: utf-8 mode: python -*-
"""
hhg_spectrum.py - 
"""
# Copyright (c) 2021, Davis Welakuh
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
from scipy import sparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# -----------------------------------------------------------------------------
# INPUT DATA
# -----------------------------------------------------------------------------

au_to_eV = 27.21138386
au_to_fs = 0.024189

ww_min = 0.0
ww_max = 1.5
dd_ww = 0.000367493
ww_omega = np.arange(ww_min, ww_max, dd_ww)

# -----------------------------------------------------------------------------
# LOAD THE REQUIRED DATA
# -----------------------------------------------------------------------------

data_laser = np.loadtxt("td.general/laser")
data_current = np.loadtxt("td.general/total_current")

time = data_laser[:, 1]
time_fs = time*au_to_fs
time_dt = time[1] - time[0]
A_field = data_laser[:, 2]

time_t = data_current[:, 1]
current_x = data_current[:, 2]
current_y = data_current[:, 3]
current_z = data_current[:, 4]

current = current_x + current_y + current_z

# -----------------------------------------------------------------------------
# DEFINE REQUIRED FUNCTIONS
# -----------------------------------------------------------------------------


def mask_function(time):
    Tcut = (1.0)*time[-1] #Mess with this
    tind = np.argmin(np.abs(time - Tcut))
    tcut = time[tind]
    mask = 1.0 - (time[:tind]/tcut)**2 + 2*(time[:tind]/tcut)**3 #window function
    return tind, mask

def my_mask_function(signal, last_vals):
    mask_func = np.ones(len(signal))
    aa = np.arange(0, last_vals+1, 1)
    bb = np.flipud((np.cos(aa/(last_vals)*np.pi/2))**2)
    mask_func[-len(bb):] = np.flipud(bb)
    return mask_func


def first_derivative_current(time_dt, current):
    d_dt_current = np.zeros(len(current))
    for ii in range(1, len(current)-1):
        d_dt_current[ii] = (current[ii+1] - current[ii-1])/(2.0*time_dt)
    return d_dt_current


def fourier_transform(signal, time, energy):
    ft_signal = np.zeros(len(energy)).astype(np.complex128)
    de = energy[1]-energy[0]
    for ei, e in enumerate(energy):
        ft_signal[ei] = np.sum(np.exp(1j*e*time[:tind])*signal)*de
    return ft_signal


def my_fourier_transform(signal, time_dt, time, ww_omega):
    ft_signal = np.zeros(len(ww_omega))
    for idx in range(len(ww_omega)):
        ft_signal[idx] = \
            abs((signal*np.exp(1j*ww_omega[idx]*time)).sum()*time_dt)**2
    return ft_signal


def get_HHG_energy(signal, time_step):
    sample_rate_fs = 1.0/time_step
    length_signal = len(signal) 
    interval = np.arange(length_signal)
    time = length_signal/sample_rate_fs
    # two sides frequency range
    frequency = interval/time
    # one side frequency range
    frequency = frequency[range(int(length_signal/2))]
    omega = 2.0*np.pi*frequency    
    return omega


def fourier_transform_FT(signal, time_step):
    length_signal = len(signal)
    normalized_fft =  np.fft.fft(signal)/length_signal
    normalized_fft = normalized_fft[range(int(length_signal/2))]
    ft_signal = abs(normalized_fft)
    return ft_signal


# -----------------------------------------------------------------------------
# COMPUTE HARMONIC SPECTRUM
# -----------------------------------------------------------------------------

last_vals = 50

mask_func = 1.0
#mask_func = np.hanning(len(current))
#mask_func = my_mask_function(current, last_vals)

current_x_mask = current_x*mask_func
current_y_mask = current_y*mask_func
current_z_mask = current_z*mask_func

ww_omega = get_HHG_energy(time_t, time_dt)

ft_current_x = my_fourier_transform(current_x_mask, time_dt, time_t, ww_omega)
ft_current_y = my_fourier_transform(current_y_mask, time_dt, time_t, ww_omega)
ft_current_z = my_fourier_transform(current_z_mask, time_dt, time_t, ww_omega)

hs_current_x = abs(1j*ft_current_x*ww_omega)**2
hs_current_y = abs(1j*ft_current_y*ww_omega)**2
hs_current_z = abs(1j*ft_current_z*ww_omega)**2

hs_current = hs_current_x + hs_current_y + hs_current_z

ww_omega = ww_omega*au_to_eV

# -----------------------------------------------------------------------------
# PLOTS START HERE
# -----------------------------------------------------------------------------


# defining the subplots
gs = gridspec.GridSpec(2, 1)
gs.update(hspace=0.55)


fig = plt.figure(figsize=(6,7))
# left-hand side
ax = fig.add_subplot(gs[0])
ax.plot(time_fs, A_field, 'r', linestyle='-')
ax.margins(x=0)
ax.set_xlabel('Time (fs)', fontsize=15)
ax.set_ylabel('Vector potential (a.u.)', fontsize=14, color="r")
# ax.set_xlim(23.5, 26.0)
# ax.set_ylim(-0.2, 52.0)
ax.get_yaxis().set_label_coords(-0.1, 0.5)
ax.tick_params(axis='both', which='major', labelsize=13)
# right-hand side
ax1 = ax.twinx()
ax1.plot(time_fs, current, "b", linestyle='-')
ax1.margins(x=0)
ax1.set_ylabel('current (a.u.)', fontsize=14, color="b")
#ax1.tick_params(axis ='y', labelcolor = "b")
plt.savefig('graphene_current_A_field.png', bbox_inches="tight", dpi=200)
plt.show()



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
ax.plot(ww_omega, hs_current, color='dodgerblue')
for xc in xposition:
    ax.axvline(x=xc, color='gray', linestyle='--', linewidth=0.9)
ax.margins(x=0)
ax.set_xlabel(r'Photon energy (eV)', fontsize=15)
ax.set_ylabel(r"log(Intensity) (a.u.)", fontsize=15)
ax.set_yscale('log') 
ax.set_xlim(0.0, 35.0)
#ax.set_ylim(1e-16, 1e+2)
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
plt.savefig('graphene_spectrum.png', bbox_inches="tight", dpi=200)
plt.show()
