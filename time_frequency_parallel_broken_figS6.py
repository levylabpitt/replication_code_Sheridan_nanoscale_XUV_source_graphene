#!/usr/bin/env python3
# -*- coding: utf-8 mode: python -*-
"""
time_frequency_parallel_broken.py - generates figure S6
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


import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec


# -----------------------------------------------------------------------------
# INPUT DATA
# -----------------------------------------------------------------------------

lmbda = 800 #in nm
omega = 4.13566733e-15*299792458/(lmbda*1.0e-9)*0.0367493
period = 2.0*np.pi/omega
Tpulse = 12*period

au_to_fs = 0.024189
au_to_eV = 27.21138386

NN = 300

ww_min = 0.0
ww_max = 1.105
dd_ww = 0.00367493
#ww_omega = np.arange(ww_min, ww_max, dd_ww)
ww_omega = np.linspace(ww_min, ww_max, NN)

tau_min = 0.0
tau_max = 2647.6 # 1434.2
tau_step = 0.5
#time_tau = np.arange(tau_min, tau_max, tau_step)
time_tau = np.linspace(tau_min, tau_max, NN)


sigma_001 = 10.335
sigma_002 = 10.335 # (0.25 fs)
#amplitude = 1.0/(np.sqrt(2.0*np.pi)*sigma)
amplitude_001 = 50.5
amplitude_002 = 20.5

Intensity_001 = np.zeros((len(time_tau), len(ww_omega)), dtype=complex)
Intensity_002 = np.zeros((len(time_tau), len(ww_omega)), dtype=complex)

print()
print("... the period is τ =", period/8)
print("... FWHM pulse duration τ =", Tpulse*au_to_fs)
print("... time step d tau =", time_tau[1] - time_tau[0])
print("... omega step d omega =", ww_omega[1] - ww_omega[0])
print("... the shape of Intensity is", np.shape(Intensity_001))

#sys.exit()

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------

# 001 --> 0.1 TW/cm^2 ; 002 --> 3.0 TW/cm^2

data_broken_ramp = np.loadtxt("data/01.parallel/02.broken-symmetry/"
                              "002.ramp_Ab_0.50/td.general/total_current")
data_broken_001 = np.loadtxt("data/01.parallel/02.broken-symmetry/"
                             "402.IL_0.10_TWcm2_ramp_Ab_0.50/"
                             "td.general/total_current")
data_broken_001_laser = np.loadtxt("data/01.parallel/02.broken-symmetry/"
                                   "402.IL_0.10_TWcm2_ramp_Ab_0.50/"
                                   "td.general/laser")
data_broken_002 = np.loadtxt("data/01.parallel/02.broken-symmetry/"
                             "602.IL_3.00_TWcm2_ramp_Ab_0.50/"
                             "td.general/total_current")
data_broken_002_laser = np.loadtxt("data/01.parallel/02.broken-symmetry/"
                                   "602.IL_3.00_TWcm2_ramp_Ab_0.50/"
                                   "td.general/laser")

# -----------------------------------------------------------------------------
# EXTRACT DATA
# -----------------------------------------------------------------------------

time_t = data_broken_001[:, 1]
time_fs = time_t*au_to_fs
time_dt = time_t[1] - time_t[0]

A_field_001 = data_broken_001_laser[:, 2] + data_broken_001_laser[:, 5]
A_field_002 = data_broken_002_laser[:, 2] + data_broken_002_laser[:, 5]

current_xy_broken_001 = (data_broken_001[:, 2] - data_broken_ramp[:, 2]) + \
    (data_broken_001[:, 3] - data_broken_ramp[:, 3])
current_xy_broken_002 = (data_broken_002[:, 2] - data_broken_ramp[:, 2]) + \
    (data_broken_002[:, 3] - data_broken_ramp[:, 3])

# print(time_t[-1])
# sys.exit(1)

# -----------------------------------------------------------------------------
# DEFINE FUNCTIONS
# -----------------------------------------------------------------------------


def power_spectrum(time_dt, amplitude, sigma, time_t, time_tau,
                   current_t, ww_omega, Intensity):

    for ii in range(len(time_tau)):
        for jj in range(len(ww_omega)):
            Intensity[ii, jj] = \
                (amplitude*np.exp(-0.5*(time_t - time_tau[jj])**2/sigma**2) * \
                 np.hanning(len(current_t))*current_t * \
                 np.exp(1j*ww_omega[ii]*time_t)).sum()*time_dt
            
    return Intensity


# -----------------------------------------------------------------------------
# COMPUTE HARMONIC SPECTRUM
# -----------------------------------------------------------------------------


filename_001 = "data/01.parallel/02.broken-symmetry/" \
    "graphene_spectrogram_001.npy"
filename_002 = "data/01.parallel/02.broken-symmetry/" \
    "graphene_spectrogram_002.npy"

if os.path.isfile(filename_001) and os.path.isfile(filename_002):
    Intensity_001 = np.load(filename_001)
    Intensity_002 = np.load(filename_002)
else:
    Intensity_001 = power_spectrum(time_dt, amplitude_001, sigma_001, time_t,
                                   time_tau, current_xy_broken_001, ww_omega,
                                   Intensity_001)
    Intensity_002 = power_spectrum(time_dt, amplitude_002, sigma_002, time_t,
                                   time_tau, current_xy_broken_002, ww_omega,
                                   Intensity_002)

    Intensity_001 = abs(ww_omega*Intensity_001)**2
    Intensity_002 = abs(ww_omega*Intensity_002)**2
    np.save(filename_001, Intensity_001)
    np.save(filename_002, Intensity_002)


ww_omega = ww_omega*au_to_eV
time_tau = time_tau*au_to_fs
time_t = time_t*au_to_fs

X_axis, Y_axis = np.meshgrid(time_tau, ww_omega)
Z_axis_001 = np.clip(np.log10(Intensity_001), -6, 2)
Z_axis_002 = np.clip(np.log10(Intensity_002), -6, 2)


# -----------------------------------------------------------------------------
# PLOT DATA
# -----------------------------------------------------------------------------

aa = 1.0
bb = 0.4

fig = plt.figure(figsize=(6,4))
gs1 = gridspec.GridSpec(2, 2)
gs1.update(hspace = 0.1, wspace = 0.3, left = 0.14,
            right = 0.84, bottom = 0.10, top = 0.93)

# FIRST SUPLOT:
ax1 = plt.subplot(gs1[0, :])
CP1 = ax1.contourf(X_axis, Y_axis, Z_axis_001, 100, cmap=plt.cm.jet)
cbar1 = fig.colorbar(CP1, ax=ax1, ticks=[0, -2.5, -5.0])
cbar1.set_label('Log(Intensity) (Arb. u.)', fontsize=10)
ax1.plot(time_t, aa*A_field_001+10, 'w', linestyle='-')
ax1.xaxis.set_visible(False)
ax1.set_xlabel(r'Time (fs)', fontsize=12)
# ax1.set_ylabel('Photon energy (eV)', fontsize=12)
ax1.text(2.0, 27.0, r"$I_{L}=100$GW/cm$^2$", color='w', fontsize=12)
ax1.text(45.0, 27.0, r"$A_{b}(t) + A_{L}(t)$", color='w', fontsize=12)
# SECOND SUPLOT:
ax2 = plt.subplot(gs1[1, :])
CP2 = ax2.contourf(X_axis, Y_axis, Z_axis_002, 100, cmap=plt.cm.jet)
cbar2 = fig.colorbar(CP2, ax=ax2, ticks=[1.75, 0, -1.75, -3.5, -5.25])
cbar2.set_label('Log(Intensity) (Arb. u.)', fontsize=10)
ax2.plot(time_t, bb*A_field_002+10, 'w', linestyle='-')
ax2.set_xlabel(r'Time (fs)', fontsize=12)
# ax2.set_ylabel('Photon energy (eV)', fontsize=12)
ax2.text(2.0, 27.0, r"$I_{L}=3$TW/cm$^2$", color='w', fontsize=12)
ax2.text(45.0, 27.0, r"$A_{b}(t) + A_{L}(t)$", color='w', fontsize=12)
fig.text(0.04, 0.5, "Photon energy (eV)", va="center", rotation="vertical", fontsize=15)
plt.savefig('Figures/graphene_spectrogram_odd_even.png', bbox_inches="tight", dpi=200)
plt.show()
