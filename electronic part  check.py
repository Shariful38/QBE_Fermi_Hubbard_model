import numpy as np
########## Electron hamiltonian formation
# Define convergence parameters
max_iterations = 100  # Maximum number of SCF iterations
tolerance = 1e-6  # Convergence tolerance for expectation values

# Define parameters
num_sites = 4  # Number of lattice sites
num_k = num_sites  # Number of electron momentum states (equal to num_sites for periodic boundary)
num_spin = 2  # Number of spin states (spin-up and spin-down)
num_branches = 3  # Number of phonon branches (lambda_1, lambda_2, lambda_3)
num_momenta = num_sites  # Number of momentum points for phonons (q-values)

t = 1.0  # Hopping parameter
U = 2.0  # On-site interaction strength
mu = 1.0  # Chemical potential
V = 1.0  # Normalization constant for electron-phonon interaction

# Random initial wavefunctions for spin-up and spin-down states (normalized)
psi_up = np.random.rand(num_k)
psi_down = np.random.rand(num_k)
psi_up /= np.linalg.norm(psi_up)
psi_down /= np.linalg.norm(psi_down)

# Create number operators for spin-up and spin-down in momentum space
n_up = np.zeros((num_k, num_k))
n_down = np.zeros((num_k, num_k))

for i in range(num_k):
    n_up[i, i] = 1  # Diagonal number operator for spin-up
    n_down[i, i] = 1  # Diagonal number operator for spin-down
# SCF loop
for iteration in range(max_iterations):

# Initialize the Hamiltonian in momentum basis
    H_k = np.zeros((num_k * num_spin, num_k * num_spin), dtype=complex)

# Function to compute flattened index for spin and momentum
    def index(k, spin):
        return k * num_spin + spin

# Term 1: Hopping term (-t * sum_k,σ c_k^† c_k)
# Define the dispersion relation for the lattice (e.g., 1D chain, nearest neighbor)
# Here, we assume a 1D lattice with spacing d=1 for simplicity
    epsilon_k = -2 * t * np.cos(2 * np.pi * np.arange(num_k) / num_k)

# Fill the momentum-space Hamiltonian
    for k in range(num_k):
        for sigma in range(num_spin):
            idx = k * num_spin + sigma  # Index for spin and momentum
            H_k[idx, idx] += epsilon_k[k]  # Add diagonal term for momentum k and spin sigma

# Term 2: On-site interaction (U/2 * sum_k (n_up_k * n_down_k))





# Calculate the expectation values <n_up> and <n_down> for each momentum state
    expectation_n_up = np.zeros(num_k)
    expectation_n_down = np.zeros(num_k)

    for k in range(num_k):
    # Expectation value of the number operator for spin-up (n_up)
        expectation_n_up[k] = np.real(np.dot(psi_up.conj()[k], np.dot(n_up[k, k], psi_up[k])))

    # Expectation value of the number operator for spin-down (n_down)
        expectation_n_down[k] = np.real(np.dot(psi_down.conj()[k], np.dot(n_down[k, k], psi_down[k])))

    for k in range(num_k):
        idx_up = index(k, 0)  # Spin-up index
        idx_down = index(k, 1)  # Spin-down index
        H_k[idx_up, idx_up] += U / 2 * expectation_n_down[k]  # Interaction term for spin-up
        H_k[idx_down, idx_down] += U / 2 * expectation_n_up[k]  # Interaction term for spin-down

# Term 3: Chemical potential (-mu * sum_k,σ c_k^† c_k)
    for k in range(num_k):
        for sigma in range(num_spin):
            idx_k = index(k, sigma)
            H_k[idx_k, idx_k] += -mu

     # Diagonalize the Hamiltonian to obtain eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(H_k)

    # Extract the new wavefunctions for spin-up and spin-down
    new_psi_up = eigenvectors[:, 0][:num_k]  # First eigenvector for spin-up
    new_psi_down = eigenvectors[:, 1][:num_k]  # Second eigenvector for spin-down

    # Normalize the new wavefunctions
    new_psi_up /= np.linalg.norm(new_psi_up)
    new_psi_down /= np.linalg.norm(new_psi_down)

    # Check for convergence
    delta_up = np.linalg.norm(new_psi_up - psi_up)
    delta_down = np.linalg.norm(new_psi_down - psi_down)

    # Print the current iteration number and differences (optional)
    print(f"Iteration {iteration + 1}: Δψ_up = {delta_up}, Δψ_down = {delta_down}")

    if delta_up < tolerance and delta_down < tolerance:
        print(f"SCF converged after {iteration + 1} iterations.")
        break

    # Update wavefunctions for the next iteration
    psi_up = new_psi_up
    psi_down = new_psi_down
else:
    print("SCF did not converge within the maximum number of iterations.")

# Final results
print("Final expectation values for spin-up:", expectation_n_up)
print("Final expectation values for spin-down:", expectation_n_down)
print("Eigenvalues of the Hamiltonian:", eigenvalues)
print("Eigenvectors", eigenvectors)