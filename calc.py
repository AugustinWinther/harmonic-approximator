"""Functions calculating energies of electron in H.O.

"""

# Third party
import numpy as np
from numpy.typing import NDArray

# Physical constants in atomic units
PI = np.pi
HBAR = 1.0

# Physical properties in atomic units
M = 1.0
W = 1.0

# Numerical constants
PREC = 19.0  # How many x-values between 0 and L to use for numerical integration.
             # Larger L need larger PREC value. 19 is good value for L=11

def T(n_row: NDArray, L: int) -> NDArray:
    """Vector containing kinetic energies corresponding to each n value in n_row
    
    """
    # Using index 0 on n_row to make output 1D array instead of 2D row vector
    return (0.5/M)*(n_row[0]*PI*HBAR/L)**2

def V(x_col: NDArray, n_row: NDArray, dx: float, L: int) -> NDArray:
    """Matrix containing potential energies V_mn
    
    """
    S = (x_col - 0.5*L) * np.sin((PI/L)*x_col*n_row)
    return (W**2)/(L*M) * (S.T @ S) * dx

def numerical_energies(N: int, L: int) -> NDArray:
    """Returns N first numerically approximated energies of H.O. 
    
    """
    # Create column vector with all x values used in integration
    # We use midpoint Riemann summation, hence the + 0.5
    dx = L/PREC
    x_col = ((np.arange(PREC, dtype=np.float64) + 0.5)*dx)[:, np.newaxis]

    # Create row-vector with all n values
    n_row = np.arange(1, N + 1, dtype=np.float64)[np.newaxis, :]

    # First add potential to total energy matrix
    H = V(x_col, n_row, dx, L)

    # Then add kinetic along the diagonal of the total energy matrix
    H[np.diag_indices(N)] += T(n_row, L)

    # Find eigenvalues of total energy matrix (eigenenergies)
    approx_energies = np.linalg.eigvalsh(H)

    return approx_energies

def analytical_energies(N: int) -> NDArray:
    """Returns N first analytically calculated energies of H.O. 
    
    """
    return W*(np.arange(0, N, dtype=np.float64) + 0.5)
