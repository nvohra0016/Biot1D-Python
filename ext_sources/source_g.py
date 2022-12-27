# This file is part of Biot1D-Python.
# Biot1D-Python is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Biot1D-Python is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Biot1D-Python. If not, see <https://www.gnu.org/licenses/>.
# The full text of the license can be found in the file LICENSE.md.

# Mass source rate "g".
# Units to use are [MPa hr / m^2] (SI units [kg / m^3 s])

import numpy as np
from scipy.sparse import csc_matrix
from parameters import *
from math import *

# function definiton
def source_g(x,t,example_number):
  lam_const, mu_const, betaf, phi, kappa, viscosity, alpha, rhof, rhos, G = parameters(x, t, example_number)
  
  if example_number == 1:
    return np.multiply(( (pi/2)*np.cos(pi*t/2)*(betaf * phi[0] + alpha) + (1/viscosity)*np.sin(pi*t/2)*pi*pi*kappa[0] ), np.sin(pi*x))
  elif example_number == 2:
    return np.multiply( (-betaf*phi[0] - alpha*(pi/2) + kappa[0]*(pi*pi/(viscosity*4)) )*np.cos(pi*x/2), np.exp(-t))
  elif example_number == 3:
    return 0 + 0*x + 0*t
  else:
    return 0*x + 0*t


























