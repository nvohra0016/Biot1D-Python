# This file is part of Poroelasticity1d.
# Poroelasticity1d is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Poroelasticity1d is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Poroelasticity1d. If not, see <https://www.gnu.org/licenses/>.
# The full text of the license can be found in the file LICENSE.md.

# List of physical parameters for different examples.
# Units are [meter, hour, megapascal] ([m, hr, MPa])
import numpy as np
from math import *

# function definiton
def parameters(x, t, example_number): # x is cell centered grid
  if example_number <= 3:
    lam_const = np.ones((x.size,1))
    mu_const = np.ones((x.size,1))
    betaf = 1.0
    phi = np.ones((x.size,1))
    kappa = np.ones((x.size,1))
    viscosity = 1.0
    alpha = 1.0
    rhof = 1.0
    rhos = 1.0 * np.ones((x.size,1))
    G = 0.0
  elif example_number == 4:
    # clay parameters
    E = 20
    nu = 0.30
    lam_const = (E*nu)/((1 + nu)*(1-2*nu)) * np.ones((x.size,1))
    mu_const = E/(2*(1+nu)) * np.ones((x.size,1))
    kappa = 1e-17 * np.ones((x.size,1))
    phi = 0.50 * np.ones((x.size,1))
    betaf = 4.16e-4
    viscosity = 2.7822e-13
    alpha = 1.0
    rhof = 998.21 * (1/3600) * (1/3600) * 1e-6
    rhos = 2700 * (1/3600) * (1/3600) * 1e-6 * np.ones((x.size,1))
    G = 1.27290528 * 1e8
  elif example_number == 5:
    # heterogeneous case sand/clay
    E_sand = 15
    nu_sand = 0.25
    E_clay = 20
    nu_clay = 0.30
    kappa_sand = 1e-12
    kappa_clay = 1e-17
    phi_sand = 0.30
    phi_clay = 0.50
    
    lam_const = np.ones((x.size,1))
    mu_const = np.ones((x.size,1))
    kappa = np.ones((x.size,1))
    phi = np.ones((x.size,1))
    rhos = 1.0 * np.ones((x.size,1))
    for j in range(0,x.size):
      if x[j] >= 0.5:
        lam_const[j] = (E_sand*nu_sand)/((1 + nu_sand)*(1-2*nu_sand))
        mu_const[j] = E_sand/(2*(1+nu_sand))
        kappa[j] = kappa_sand
        phi[j] = phi_sand
        rhos[j] = 2650 * (1/3600) * (1/3600) * 1e-6
      else:
        lam_const[j] = (E_clay*nu_clay)/((1 + nu_clay)*(1-2*nu_clay))
        mu_const[j] = E_clay/(2*(1+nu_clay))
        kappa[j] = kappa_clay
        phi[j] = phi_clay
        rhos[j] = 2700 * (1/3600) * (1/3600) * 1e-6
        
    betaf = 4.16e-4
    viscosity = 2.7822e-13
    alpha = 1.0
    rhof = 998.21 * (1e-6 * (1/3600)*(1/3600))
    G = 1.27290528e8 * 0
    
  else:
    print("Example not implemented. Input parameter/ grid values")
    
  return lam_const, mu_const, betaf, phi, kappa, viscosity, alpha, rhof, rhos, G


























