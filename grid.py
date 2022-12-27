# This file is part of Poroelasticity1d.
# Poroelasticity1d is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Poroelasticity1d is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Poroelasticity1d. If not, see <https://www.gnu.org/licenses/>.
# The full text of the license can be found in the file LICENSE.md.

# Spatial (can be non-uniform) and temporal grid
# if M >= 1, then construct uniform spatial grid with number of cells M otherwise return manually constructed non-uniform spatial grid.
# only uniform temporal grid for now using Mt cells

import numpy as np
from scipy.sparse import csc_matrix
from math import *

# function definiton
def grid(Nxcells, Ntcells, example_number):
  # domain specification (end points) based on example_number
  if example_number <= 3:
    a = 0
    b = 1.0
    Tend = 1.0
  elif example_number == 4:
    a = 0
    b = 0.1
    Tend = 24
  elif example_number == 5:
    a = 0
    b = 1
    Tend = 24
  else:
    # custom values
    a = 0
    b = 1
    Tend = 1
  # define grid using Nxcells
  if type(Nxcells) is int:
    M  = Nxcells
    h = (b-a)/M
    xn = np.arange(a,b+h,h)[:, np.newaxis]        # nodal grid
    xcc = np.arange(a + h/2,b,h)[:,np.newaxis]    # cell centered grid
    h = (b-a)/M * np.ones((M,1))
  else:
    # define nodal grid including endpoints
    el_nodes = Nxcells
    M = el_nodes.size - 1;
    h = np.zeros((M,1))
    for j in range(0,el_nodes.size - 1):
      h[j] = el_nodes[j+1] - el_nodes[j]
     
    xn = el_nodes
    xcc = np.zeros((M,1))
    for j in range(0,el_nodes.size - 1):
      xcc[j] = (el_nodes[j] + el_nodes[j+1])/2
    #
    
  # define time step tau using Ntcells
  tau = (Tend - 0)/Ntcells
  t = np.arange(0,Tend+tau,tau)[:,np.newaxis]
  
  return a, b, h, M, xcc, xn, tau, t, Tend
  
  


























