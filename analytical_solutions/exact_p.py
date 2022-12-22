# This file is part of Poroelasticity1d.
# Poroelasticity1d is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Poroelasticity1d is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Poroelasticity1d. If not, see <https://www.gnu.org/licenses/>.
# The full text of the license can be found in the file LICENSE.md.

# Analytical pressure expression.

import numpy as np
from scipy.sparse import csc_matrix
from math import *

# function definiton
def exact_p(x,t,example_number):

  if example_number == 1:
    return np.multiply(np.sin(pi*t/2),np.sin(pi*x))
  elif example_number == 2:
    return np.multiply(np.cos(pi*x/2),np.exp(-t))
  elif example_number == 3:
    return 1 + x + 0*t
  else:
    return 0*x + 0*t

























