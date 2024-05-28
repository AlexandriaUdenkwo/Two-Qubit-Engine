#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 12:04:37 2024

@author: katermurch
"""



import numpy as np
import matplotlib.pyplot as plt

#savedir = r"/Users/katermurch/Documents/Python/entanglement_engine_data"

# Constants
g = 10
delta = 13
theta = np.arctan(g / delta)
omega = np.sqrt(delta**2 + g**2)

# Time array
t = np.linspace(0, 2 * np.pi / omega, 1000)

# Define the functions
a = np.cos(theta / 2) ** 2 * np.exp(1j * omega * t / 2) + np.sin(theta / 2) ** 2 * np.exp(-1j * omega * t / 2)
astar = np.cos(theta / 2) ** 2 * np.exp(-1j * omega * t / 2) + np.sin(theta / 2) ** 2 * np.exp(1j * omega * t / 2)
b = -np.cos(theta / 2) * np.sin(theta / 2) * (np.exp(1j * omega * t / 2) - np.exp(-1j * omega * t / 2))
bstar = -np.cos(theta / 2) * np.sin(theta / 2) * (np.exp(-1j * omega * t / 2) - np.exp(1j * omega * t / 2))

# First plot: (10 astar a + (10 + delta) bstar b)
y1 = (10 * astar * a + (10 + delta) * bstar * b).real
# Second plot: (g / 2 * (astar b + bstar a))
y2 = (g / 2 * (astar * b + bstar * a)).real

y0 = y1[0]

fig, ax = plt.subplots(figsize=(3, 2))

ax.plot(t, y1, label=r'$\langle H_\mathrm{loc}\rangle$', color='red')
ax.plot(t, y2, label=r'$\langle V\rangle$', color='blue')
ax.plot(y0, label=r'$\langle H\rangle$', color='black', linestyle='--')
ax.set_ylim(-5.5, 17)

# Set custom x-axis ticks and labels
ticks = [0, np.pi / omega, 2 * np.pi / omega]
dic = {np.pi / omega: r"$\pi/\Omega$", 2 * np.pi / omega: r"$2\pi/\Omega$"}
labels = [dic.get(t, t) for t in ticks]
ax.set_xticks(ticks)
ax.set_xticklabels(labels)

# Set custom y-axis ticks and labels
ticks2 = [0, 10]
dic2 = {10: r"$\omega_\mathrm{A}$"}
labels2 = [dic2.get(t, t) for t in ticks2]
ax.set_yticks(ticks2)
ax.set_yticklabels(labels2)

# Set axis labels
ax.set_xlabel(r'$t_0$')
ax.set_ylabel('Energy')

# Add legend
ax.legend()

# Save and show the plot
plt.savefig('first_plot.pdf')
plt.show()



# First plot: (10 astar a + (10 + delta) bstar b)
y3 = (10 * astar * a + (10 + delta) * bstar * b).real * np.heaviside(-t+np.pi/omega,1) + np.heaviside(t-np.pi/omega,1)*y1[500]
# Second plot: (g / 2 * (astar b + bstar a))
y4 = (g / 2 * (astar * b + bstar * a)).real * np.heaviside(-t+np.pi/omega,1)
y5 = y1[0] * np.heaviside(-t+np.pi/omega,1) + np.heaviside(t-np.pi/omega,1)*y1[500]


fig, ax = plt.subplots(figsize=(3, 2))

ax.plot(t, y3, label=r'$\langle H_\mathrm{loc}\rangle$', color='red')
ax.plot(t, y4, label=r'$\langle V\rangle$', color='blue')
ax.plot(t, y5, label=r'$\langle H\rangle$', color='black', linestyle='--')
ax.set_ylim(-5.5, 17)

# Set custom x-axis ticks and labels
ticks = [0, np.pi / omega, 2 * np.pi / omega]
dic = {np.pi / omega: r"$\pi/\Omega$", 2 * np.pi / omega: r"$2\pi/\Omega$"}
labels = [dic.get(t, t) for t in ticks]
ax.set_xticks(ticks)
ax.set_xticklabels(labels)

# Set custom y-axis ticks and labels
ticks2 = [0, 10]
dic2 = {10: r"$\omega_\mathrm{A}$"}
labels2 = [dic2.get(t, t) for t in ticks2]
ax.set_yticks(ticks2)
ax.set_yticklabels(labels2)

# Set axis labels
ax.set_xlabel(r'$t_0$')
ax.set_ylabel('Energy')

# Add legend
ax.legend(loc='center right')

# Save and show the plot
plt.savefig('second_plot.pdf')
plt.show()
