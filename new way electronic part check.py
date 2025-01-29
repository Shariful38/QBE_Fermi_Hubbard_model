import numpy as np
from scipy.linalg import eigh

# Model parameters
L = 6  # Number of lattice sites
t = 1.0  # Hopping parameter
U = 4.0  # On-site interaction
tol = 1e-6  # Convergence tolerance
max_iter = 100  # Max iterations
N_elec = L  # Total number of electrons (half-filling)

# Initialize occupation numbers randomly (between 0 and 1)
n_up = np.random.rand(L)
n_down = np.random.rand(L)
mu = 0.0  # Initial guess for chemical potential

def construct_mean_field_hamiltonian(n_up, n_down, mu):
    """Constructs the mean-field Hamiltonian for the Fermi-Hubbard model with chemical potential."""
    H = np.zeros((2*L, 2*L), dtype=complex)  # 2L because we consider both spin-up and spin-down
    
    # Hopping term (-t sum c†_i c_j)
    for i in range(L - 1):
        H[i, i + 1] = H[i + 1, i] = -t  # Spin-up hopping
        H[L + i, L + i + 1] = H[L + i + 1, L + i] = -t  # Spin-down hopping

    # Mean-field interaction term (U * n_opp) and chemical potential term (-mu * n)
    for i in range(L):
        H[i, i] = U * n_down[i] - mu  # Spin-up sees spin-down density
        H[L + i, L + i] = U * n_up[i] - mu  # Spin-down sees spin-up density

    return H

def fermi_dirac_occupations(eigvals, mu, beta=1e6):
    """Compute occupation numbers using the Fermi-Dirac distribution at T=0 (beta -> infinity)."""
    return (eigvals - mu < 0).astype(float)  # Step function at T=0

# Self-consistent loop
for iteration in range(max_iter):
    # Construct mean-field Hamiltonian
    H_mf = construct_mean_field_hamiltonian(n_up, n_down, mu)
    
    # Solve for eigenvalues and eigenvectors
    eigvals, eigvecs = eigh(H_mf)

    # Compute new occupation numbers using Fermi-Dirac filling
    occupations = fermi_dirac_occupations(eigvals, mu)
    
    # Compute new density matrices
    n_up_new = np.sum(np.abs(eigvecs[:L, :])**2 * occupations, axis=1)
    n_down_new = np.sum(np.abs(eigvecs[L:, :])**2 * occupations, axis=1)

    # Adjust chemical potential to enforce electron number conservation
    N_current = np.sum(n_up_new) + np.sum(n_down_new)
    mu += (N_current - N_elec) * 0.1  # Newton-Raphson step to correct total electron number

    # Check convergence
    if np.linalg.norm(n_up_new - n_up) < tol and np.linalg.norm(n_down_new - n_down) < tol:
        print(f"Converged in {iteration+1} iterations.")
        break
    
    # Update densities
    n_up, n_down = n_up_new, n_down_new

# Print final results
print("Final spin-up occupation numbers:", n_up)
print("Final spin-down occupation numbers:", n_down)
print("Final chemical potential:", mu)
print("Final eigenvalues of H_mf:", eigvals)
print("Final eigenvectors of H_mf:", eigvecs)
