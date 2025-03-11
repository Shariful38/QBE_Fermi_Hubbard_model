#comparison between pyscf and our code
from pyscf import scf, gto
import numpy as np

n = 4  # Number of sites
mol_hub = gto.M()
mol_hub.nelectron = n // 2
mol_hub.incore_anyway = True

# Define the hopping term (kinetic energy)
h1 = np.zeros([n, n], dtype=np.float64)
for i in range(n - 1):
    h1[i, i + 1] = h1[i + 1, i] = -1.0
h1[n - 1, 0] = h1[0, n - 1] = -1.0  # Periodic boundary condition

# Define the two-electron interaction tensor (set to zero to exclude U)
eri = np.zeros([n] * 4, dtype=np.float64)

# Set up the Hartree-Fock calculation
rhf_hub = scf.RHF(mol_hub)
rhf_hub.get_hcore = lambda *args: h1
rhf_hub.get_ovlp = lambda *args: np.eye(n)
rhf_hub._eri = eri  # No interaction term
rhf_hub.init_guess = '1e'
rhf_hub.kernel()


## our code

import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt

########## Electron Hamiltonian Formation
# Define convergence parameters
#max_iterations = 100  # Maximum number of SCF iterations
tolerance = 1e-6  # Convergence tolerance for expectation values

# Define system parameters
num_sites = 4  # Number of lattice sites

num_spin = 2  # Number of spin states (spin-up and spin-down)
num_electrons = M = num_sites  # Number of electrons (assuming half-filling)

t = 1.0  # Hopping parameter
#U = 2.0  # On-site interaction strength
#mu = 1.0  # Chemical potential

# Random initial wavefunctions for spin-up and spin-down states (normalized)
psi_up = np.random.rand(num_sites, num_electrons // 2)
psi_down = np.random.rand(num_sites, num_electrons // 2)

# Normalize the wavefunctions (columns of Slater determinant matrix)
def normalize(psi):
    for i in range(psi.shape[1]):
        psi[:, i] /= np.linalg.norm(psi[:, i])
    return psi

psi_up = normalize(psi_up)
psi_down = normalize(psi_down)

# Slater determinant matrices
Phi_up = psi_up.copy()
Phi_down = psi_down.copy()

delta_threshold = 1e-6  # Convergence threshold
iteration=0
E_prev = 0
energy_values = []  # To store energy values for plotting

# SCF loop
while True:

    # Initialize the Hamiltonian matrices for spin-up and spin-down
    H_up = np.zeros((num_sites, num_sites), dtype=complex)
    H_down = np.zeros((num_sites, num_sites), dtype=complex)

    # Add hopping terms and interactions for both spin-up and spin-down
    for i in range(num_sites-1):  # Loop over lattice sites

        # For spin-up
         H_up[i, i+1] = H_up[i+1, i] = -t

        # For spin-down
         H_down[i, i+1] = H_down[i+1, i] = -1
    H_up[num_sites-1, 0] = H_up[0, i-1] = -1
    H_down[num_sites-1, 0] = H_down[0, i-1] = -1


     # Diagonalize the Hamiltonian separately for spin-up and spin-down
    eigenvalues_up, eigenvectors_up = np.linalg.eigh(H_up)
    eigenvalues_down, eigenvectors_down = np.linalg.eigh(H_down)

    # Update the wavefunctions for spin-up and spin-down
    # For spin-up: select the eigenvectors corresponding to the lowest eigenvalues (occupied states)
    num_occupied_up = num_electrons // 2
    num_occupied_down = num_electrons // 2

    new_psi_up = normalize(eigenvectors_up[:, :num_occupied_up])
    new_psi_down = normalize(eigenvectors_down[:, :num_occupied_down])

    # Construct the Slater determinant using the tensor product
    slater_matrix = np.kron(new_psi_up, new_psi_down)  # ⊗ Tensor product
    H_total = np.kron(H_up, np.eye(M)) + np.kron(np.eye(M), H_down)

   # Compute energy using the Hamiltonian
    E_curr = np.trace(slater_matrix.T.conj() @ H_total @ slater_matrix)
      # Store energy for plotting
    energy_values.append(E_curr)
    # Print the energy at each iteration
    print(f"our code Iteration {iteration + 1}: Energy = {E_curr}")

    # Convergence check
    if abs(E_curr - E_prev) < delta_threshold:
        print(f"SCF Converged after {iteration + 1} iterations.")
        break



    # Update wavefunctions for next iteration
    psi_up = new_psi_up
    psi_down = new_psi_down
    E_prev = E_curr
    iteration += 1



# Final results
#print("Final Slater Determinant (Slater Matrix):\n", slater_matrix)
#print("Final expectation values (n_up):", np.real(np.diag(new_psi_up @ new_psi_up.T.conj())))
#print("Final expectation values (n_down):", np.real(np.diag(new_psi_down @ new_psi_down.T.conj())))
#print("Eigenvalues (spin-up):", eigenvalues_up)
#print("Eigenvalues (spin-down):", eigenvalues_down)