# This file is part of Poroelasticity1d.
# Poroelasticity1d is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Poroelasticity1d is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Poroelasticity1d. If not, see <https://www.gnu.org/licenses/>.
# The full text of the license can be found in the file LICENSE.md.

# 1D poroelasticity system using Biot's equations.
# Mixed boundary conditions implementation for mechanics [M] and hydrological flow [H] implemneted as [MMHH]
# [M, H]: 0 for Dirichlet and 1 for Neumann
# units [m hr MPa]

# import libraries/ modules
import sys
import os
import numpy as np
from scipy.sparse import csc_matrix

# for plotting/ printing
import matplotlib.pyplot as plt
np.set_printoptions(precision = 15)

# get current directory
cwd = os.getcwd()
# define external source folder path
sources_path = cwd + "/ext_sources"
sys.path.insert(0, sources_path)
# define analytical solution folder path
analytical_sol_path = cwd + "/analytical_solutions"
sys.path.insert(0, analytical_sol_path)
# check if there is an output folder in current directory
# output folder path
out_path = cwd + "/Output"
isExist = os.path.exists(out_path)
# if no output folder, create one
if not isExist:
  print("Creating Output folder: ", out_path)
  os.makedirs(out_path)

# import necessary files
from parameters import  *
from grid import *
from source_f import *
from source_g import *
from exact_u import *
from exact_du import *
from exact_p import *
from exact_dp import *

def driver(Nxcells, Ntcells, example_number):

  Nxcells = int(Nxcells)
  Ntcells = int(Ntcells)
  example_number = int(example_number)
  # spatial/temporal grid parameters
  #  construct/ extract grid
  # a, b, h, M, xcc, xn, tau, t = grid(np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]), 3)
  a, b, h, M, xcc, xn, tau, t, Tend = grid(Nxcells, Ntcells, example_number)
  
  # set boundary flags [M, M, H, H]
  BC = np.array([1, 0, 0, 1])
  # boundary flags
  Neum_ua = BC[0]
  Neum_ub = BC[1]
  Neum_pa = BC[2]
  Neum_pb = BC[3]
  # get physical parameters
  lam_const, mu_const, betaf, phi, kappa, viscosity, alpha, rhof, rhos, G = parameters(xcc, t[0], example_number)
  # calculate average density
  rho_avg = np.zeros((M,1))
  for j in range(0,M-1):
    rho_avg[j] = phi[j] * rhof + (1 - phi[j]) * rhos[j]
  # compute FE matrices
  # displacement stiffness matrix
  Auu = np.zeros((M+1, M+1)) # change to sparse (and all other matrices too)
  for j in range(1,M):
    Auu[j,j-1] = (-1/h[j-1]) * (lam_const[j-1] + 2*mu_const[j-1])
    Auu[j,j] = (1/h[j-1]) * (lam_const[j-1] + 2*mu_const[j-1]) + (1/h[j]) * (lam_const[j] + 2*mu_const[j])
    Auu[j,j+1] = -(1/h[j]) * (lam_const[j] + 2*mu_const[j])

  Auu[0,0] = (1/h[0]) * (lam_const[0] + 2*mu_const[0])
  Auu[0,1] = (-1/h[0]) * (lam_const[0] + 2*mu_const[0])
  Auu[M,M-1] = (-1/h[M-1]) * (lam_const[M-1] + 2*mu_const[M-1])
  Auu[M,M] = (1/h[M-1]) * (lam_const[M-1] + 2*mu_const[M-1])
  # pressure-displacement stiffness matrix
  Apu = np.zeros((M+1,M))
  for j in range(1,M):
    Apu[j,j-1] = 1
    Apu[j,j] = -1
  
  Apu[0,0] = -1
  Apu[M,M-1] = 1
  # displacement-pressure stiffness matrix
  Aup = Apu.transpose()
  # pressure mass matrix
  Mpp = np.zeros((M,M))
  for j in range(0,M):
    Mpp[j,j] = betaf * phi[j] * h[j]
  # transmissibilities
  T = np.zeros((M+1,1))
  for j in range(1,M):
    T[j] = ( (h[j-1]/2)*(kappa[j-1]/viscosity)**(-1) + (h[j]/2)*(kappa[j]/viscosity)**(-1) )**(-1)

  T[0] = ( (h[0]/2)*(kappa[0]/viscosity)**(-1) )**(-1)
  T[M] = ( (h[M-1]/2)*(kappa[M-1]/viscosity)**(-1) )**(-1)
  # flux mass matrix
  Mqfqf = np.zeros((M+1,M+1))
  for j in range(0,M+1):
    Mqfqf[j,j] = T[j]**(-1)
  # pressure-flux stiffness matrix
  Apqf = np.zeros((M+1,M))
  for j in range(1,M):
    Apqf[j,j-1] = 1
    Apqf[j,j] = -1
  
  Apqf[0,0] = -1
  Apqf[M,M-1] = 1
  # flux-pressure stiffness matrix
  Aqfp = Apqf.transpose()
  # gravity contribution
  # Mechanics
  G_vector_M = np.zeros((M+1,1))
  for j in range(1,M):
    G_vector_M[j] = rho_avg[j] * (1/2) * (h[j-1] + h[j]) * G

  G_vector_M[0] = rho_avg[0] * (1/2) * h[0] * G
  G_vector_M[M] = rho_avg[M-1] * (1/2) * h[M-1] * G
  # Flow
  G_vector_H = np.zeros((M+1,1))
  for j in range(1,M):
    G_vector_H[j] = rhof * (1/2) * (h[j-1] + h[j]) * G

  G_vector_H[0] = rhof * (1/2) * h[0] * G
  G_vector_H[M] = rhof * (1/2) * h[M-1] * G
  # boundary conditions, eliminate nodes from matrices
  # displacement boundary conditions
  if example_number <=3:
    Ua = exact_u(a,t,example_number)
    Ub = exact_u(b,t,example_number)
    tNa = ( (lam_const[0] + 2*mu_const[0])*exact_du(a,t,example_number) - alpha*exact_p(a,t,example_number) )* (-1)
    tNb = ( (lam_const[M-1] + 2*mu_const[M-1])*exact_du(b,t,example_number) - alpha*exact_p(b,t,example_number) ) * (1)

  else:
    Ua = 0.0 * np.ones((t.size,1))
    Ub = 0.0 * np.ones((t.size,1))
    tNa = 1e-1 * np.ones((t.size,1))
    tNb = 0.0 * np.ones((t.size,1))
    
  if Neum_ua == 0 and Neum_ub == 0:
    free_nodes_u = np.arange(1,M,1)
  elif Neum_ua == 0 and Neum_ub == 1:
    free_nodes_u = np.arange(1,M+1,1)
  elif Neum_ua == 1 and Neum_ub == 0:
    free_nodes_u = np.arange(0,M,1)
  else:
    print("Cannot have only Neumann BC for displacement... stopping")
    sys.exit()
  # eliminate nodes from displacement matrices
  Auu = Auu[np.ix_(free_nodes_u,free_nodes_u)]
  Apu = Apu[np.ix_(free_nodes_u)]
  Aup = Apu.transpose()
  G_vector_M = G_vector_M[np.ix_(free_nodes_u)]
  # pressure, flux boundary conditions
  Pa = exact_p(a,t,example_number)
  Pb = exact_p(b,t,example_number)
  qfa = -(kappa[0]/viscosity)*exact_dp(a,t,example_number) * (-1)
  qfb = -(kappa[M-1]/viscosity)*exact_dp(b,t,example_number) * (1)

  if Neum_pa == 0 and Neum_pb == 0:
    free_nodes_qf = np.arange(0,M+1,1)
  elif Neum_pa == 0 and Neum_pb == 1:
    free_nodes_qf = np.arange(0,M,1)
  elif Neum_pa == 1 and Neum_pb == 0:
    free_nodes_qf = np.arange(1,M+1,1)
  else:
    free_nodes_qf = np.arange(1,M,1)
  # eliminate nodes from flux matrix
  Mqfqf = Mqfqf[np.ix_(free_nodes_qf,free_nodes_qf)]
  Apqf = Apqf[np.ix_(free_nodes_qf)]
  Aqfp = Apqf.transpose()
  G_vector_H = G_vector_H[np.ix_(free_nodes_qf)]
  # block matrix
  Mqfqfinv = np.linalg.inv(Mqfqf)
  A = np.block([[Auu, -alpha*Apu], [-alpha*Aup, -Mpp - tau*Aqfp.dot(Mqfqfinv).dot(Apqf)]])
  # initialize values/ vectors
  # initial fluid content
  etaf = np.multiply(betaf*phi, exact_p(xcc,t[0], example_number)) + alpha*exact_du(xcc, t[0], example_number)
  # boundary/ source terms
  pressure_boundary = np.zeros((free_nodes_qf.size,1))
  # settlement
  settlement = np.zeros((t.size,1))
  # time loop
  for n in range(1,t.size):
    # mechanical sources and boundary conditions
    # source contribution
    rhs_f_vector = np.zeros((xn.size,1))
    rhs_f_val = source_f(xn,t[n],example_number)
    for j in range (1,M):
      rhs_f_vector[j] = (h[j-1] + h[j])/2 * rhs_f_val[j]
    
    rhs_f_vector[0] = (h[0]/2) * rhs_f_val[0]
    rhs_f_vector[M] = (h[M-1]/2) * rhs_f_val[M]
    # eliminate boundary nodes
    rhs_f_vector = rhs_f_vector[np.ix_(free_nodes_u)]
    # displacement boundary conditions
    rhs_f_vector[0] += (1 - Neum_ua) * (1/h[0]) * (lam_const[0] + 2*mu_const[0]) * Ua[n] + Neum_ua * tNa[n]
    rhs_f_vector[-1] += (1 - Neum_ub) * (1/h[M-1]) * (lam_const[M-1] + 2*mu_const[M-1]) * Ub[n] + Neum_ub * tNb[n]
    # gravity contribution
    rhs_f_vector += G_vector_M
    # flow sources and boundary conditions
    # source contribution
    rhs_g_vector = np.zeros((xcc.size,1))
    rhs_g_vector = tau * np.multiply(h, source_g(xcc,t[n],example_number))
    # previous time step contribution
    rhs_g_vector += np.multiply(h, etaf)
    # displacement boundary contribution
    rhs_g_vector[0] += alpha * (1 - Neum_ua) * Ua[n]
    rhs_g_vector[M-1] += alpha * (1 - Neum_ub) * (-Ub[n])
    # pressure boundary contribution
    pressure_boundary[0] = (1 - Neum_pa) * Pa[n]
    pressure_boundary[-1] = -(1 - Neum_pb) * Pb[n]
  
    rhs_g_vector += -tau * Aqfp.dot(Mqfqfinv).dot(pressure_boundary + G_vector_H)
    # flux boundary contribution
    rhs_g_vector[0] += -Neum_pa * tau * qfa[n]
    rhs_g_vector[-1] += -Neum_pb * tau * qfb[n]
  
    # solve system
    UP = np.linalg.solve(A, np.concatenate([rhs_f_vector, -rhs_g_vector]))
    # extract displacement and pressure
    U = UP[0:free_nodes_u.size]
    P = UP[free_nodes_u.size:free_nodes_u.size + M + 1]
    # post process
    # add Dirichlet boundary conditions to displacement vector
    U_BC = np.concatenate([Ua[n]*np.ones((1,1)), np.zeros((M-1,1)), Ub[n]*np.ones((1,1))])
    U_BC[np.ix_(free_nodes_u)] = U
    U = U_BC
    # compute flux
    Q = np.linalg.solve(Mqfqf, Apqf.dot(P) + pressure_boundary + G_vector_H)
    # add flux boundary conditions
    Q_BC = np.concatenate([-qfa[n]*np.ones((1,1)), np.zeros((M-1,1)), qfb[n]*np.ones((1,1))])
    Q_BC[np.ix_(free_nodes_qf)] = Q
    Q = Q_BC
    # update fluid content
    etaf = betaf * np.multiply(phi,P);
    for  j in range(0,M):
      etaf[j] += alpha * (U[j+1] - U[j])/h[j]
    # settlement values
    settlement[n] = U[0]
  # time loop ends

  # plot solution
  # plot displacement
  fig, ax = plt.subplots()

  if example_number <=3:
    ax.plot(xn, exact_u(xn, Tend, example_number), ls = '-', color = 'black', linewidth = 2, label = 'Exact')
    ax.plot(xn, U, ls = '--', marker = '.', markerfacecolor = 'red', markersize = 10, color = 'red', linewidth = 3, label = 'Numerical')
  else:
    plt.plot(xn, U, ls = '--', marker = '.', markerfacecolor = 'red', markersize = 10, color = 'red', linewidth = 3, label = 'Numerical')
  
  plt.xticks([a, (a+b)/2, b])
  plt.xticks(fontsize = 16)
  plt.yticks(fontsize = 16)
  ax = plt.gca
  plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
  plt.ylabel('Displacement $u$ [m]', fontsize = 18)
  plt.xlabel('$x$ [m]', fontsize = 18)
  plt.title('t = '+str(Tend)+', M = '+str(M)+ r', $\tau$ = '+str(tau), fontsize = 18)
  plt.legend(loc = 'upper left', fontsize = 14)
  plt.legend(frameon = True)
  plt.tight_layout()
  plt.savefig('Output/displacement.png', dpi = 200)
  print("displacement outputted")
  # plot pressure
  fig, ax = plt.subplots()

  if example_number <=3:
    ax.plot(xcc, exact_p(xcc, Tend, example_number), ls = '-', color = 'black', linewidth = 2, label = 'Exact')
    ax.plot(xcc, P, ls = '', marker = 'o', markerfacecolor = 'none', markeredgecolor = 'blue', markersize = 10, color = 'blue', linewidth = 3, label = 'Numerical')
  else:
    plt.plot(xcc, P, ls = '', marker = 'o', markerfacecolor = 'none', markeredgecolor = 'blue', markersize = 10, color = 'skyblue', linewidth = 3, label = 'Numerical')
  
  plt.xticks([a, (a+b)/2, b])
  plt.xticks(fontsize = 16)
  plt.yticks(fontsize = 16)
  ax = plt.gca
  plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
  plt.ylabel('Pressure $p$ [MPa]', fontsize = 18)
  plt.xlabel('$x$ [m]', fontsize = 18)
  plt.title('t = '+str(Tend)+', M = '+str(M)+ r', $\tau$ = '+str(tau), fontsize = 18)
  plt.legend(loc = 'upper left', fontsize = 14)
  plt.legend(frameon = True)
  plt.tight_layout()
  plt.savefig('Output/pressure.png', dpi = 200)
  print("pressure outputted")
  # plot flux
  fig, ax = plt.subplots()

  plt.plot(xn, Q, ls = '--', marker = '.', markerfacecolor = 'red', markeredgecolor = 'red', markersize = 10, color = 'red', linewidth = 3, label = 'Numerical')
  
  plt.xticks([a, (a+b)/2, b])
  plt.xticks(fontsize = 16)
  plt.yticks(fontsize = 16)
  ax = plt.gca
  plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
  plt.ylabel('Flux $q_f$ [m/hr]', fontsize = 18)
  plt.xlabel('$x$ [m]', fontsize = 18)
  plt.title('t = '+str(Tend)+', M = '+str(M)+ r', $\tau$ = '+str(tau), fontsize = 18)
  plt.legend(loc = 'upper left', fontsize = 14)
  plt.legend(frameon = True)
  plt.tight_layout()
  plt.savefig('Output/flux.png', dpi = 200)
  print("flux outputted")
  # plot settlement
  if example_number >= 5:
    print('Maximum settlement: '+str(settlement[-1]))
    fig, ax = plt.subplots()

    plt.plot(t, settlement, ls = '-', marker = '.', markerfacecolor = 'blue', markeredgecolor = 'blue', markersize = 10, color = 'blue', linewidth = 2, label = 'Numerical')
  
    plt.xticks([0, Tend])
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    ax = plt.gca
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.ylabel('Settlement $s$ [m]', fontsize = 18)
    plt.xlabel('$t$ [hr]', fontsize = 18)
    plt.title(' M = '+str(M)+ r', $\tau$ = '+str(tau), fontsize = 18)
    #plt.box(False)
    plt.tight_layout()
    plt.savefig('Output/settlement.png', dpi = 200)
    print("settlement outputted")

if __name__ == '__main__':
  args = sys.argv[1:]
  globals()[args[0]](args[1], args[2], args[3])
