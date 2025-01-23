
# ########### In position basis
# import numpy as np

# # Parameters
# t = 1.0  # Hopping parameter
# L = 4    # Number of lattice sites
# n_spin = 2  # Two spin states: 0 (up) and 1 (down)
# U = 2.0  # On-site interaction strength
# mu = 0.5  # Example chemical potential


# # Total number of states
# n_states = L * n_spin

# # Initialize the Hamiltonian matrix
# H_el = np.zeros((n_states, n_states), dtype=complex)

# # Function to get the index for a given site and spin
# def index(site, spin):
#     return site * n_spin + spin

# # Add hopping terms to the Hamiltonian
# for i in range(L):  # Loop over lattice sites
#     for j in range(L):  # Loop over neighboring sites
#         if i != j:  # Avoid diagonal elements (on-site terms)
#             for sigma in range(n_spin):  # Loop over spin states (0: spin-up, 1: spin-down)
#                 # Get indices for the current site and spin
#                 idx_i = index(i, sigma)  # State index for site i and spin sigma
#                 idx_j = index(j, sigma)  # State index for site j and spin sigma
                
#                 # Add hopping term and its Hermitian conjugate
#                 H_el[idx_i, idx_j] += -t  # c_i^\dagger c_j (hopping from j to i)
#                 H_el[idx_j, idx_i] += -t  # Hermitian conjugate: c_j^\dagger c_i (hopping from i to j)



# #### Onsite interaction term calculation
# # Define the number of sites and initial state (this could be a more complex quantum state)
# #num_sites = 4  # Example for 10 sites
# psi_up = np.random.rand(n_states)  # Random initial state for spin-up (normalized)
# psi_down = np.random.rand(n_states)  # Random initial state for spin-down (normalized)

# # Normalize the states (important in quantum mechanics)
# psi_up /= np.linalg.norm(psi_up)
# psi_down /= np.linalg.norm(psi_down)

# # Create the number operator for spin-up and spin-down (example: diagonal matrices)
# n_up = np.zeros((n_states, n_states))
# n_down = np.zeros((n_states, n_states))

# for i in range(L):
#     n_up[i, i] = 1  # Spin-up number operator (1 at each site for spin-up)
#     n_down[i, i] = 1  # Spin-down number operator (1 at each site for spin-down)

# # Calculate the expectation values <n_up> and <n_down>
# expectation_n_up = np.real(np.dot(psi_up.conj(), np.dot(n_up, psi_up)))
# expectation_n_down = np.real(np.dot(psi_down.conj(), np.dot(n_down, psi_down)))

# # On-site interaction term
# for i in range(L):  # Loop over all sites
#     idx_up = index(i, 0)      # Spin-up index at site i
#     idx_down = index(i, 1)    # Spin-down index at site i

#     # Interaction term: U * n_up * n_down
#     H_el[idx_up, idx_up] += U / 2 * expectation_n_up  # Contribution from spin-down
#     H_el[idx_down, idx_down] += U / 2 * expectation_n_down  # Contribution from spin-up

# # Chemical potential term calculation
# chemical_potential_term = -mu * (n_up + n_down)

# #### electron-phonon term 

# # Define parameters
# num_sites=4
# num_branches = 3  # Number of phonon branches (lambda_1, lambda_2, lambda_3)
# num_momenta = 5  # Number of momentum points for each phonon branch (q-values)
# num_k = 5  # Number of electron momentum states
# num_spin = 2  # Number of spin states (spin-up and spin-down)

# # Define the phonon state (could be thermal, coherent, etc.)
# # Here, we assume a random initial state and normalize it.
# psi_phonon = np.random.rand(num_sites, num_branches, num_momenta)
# psi_phonon /= np.linalg.norm(psi_phonon)

# # Define the electron state (random state for simplicity, normalized)
# psi_electron = np.random.rand(2*num_sites)  # Random electron state
# psi_electron /= np.linalg.norm(psi_electron)

# # Define the coupling constant g_lambda(k, q)
# g_lambda = np.random.rand(num_branches, num_momenta, num_k)  # Coupling constants for each branch, q, k

# # Define the electron creation and annihilation operators
# c = np.zeros((num_k, num_k, num_spin))  # Electron annihilation operator matrix for each spin
# c_dag = np.zeros((num_k, num_k, num_spin))  # Electron creation operator matrix for each spin

# # Fill the creation and annihilation matrices for each spin
# for k in range(num_k):
#     for sigma in range(num_spin):
#         c[k, k, sigma] = 1  # Annihilation operator at momentum k and spin sigma
#         c_dag[k, k, sigma] = 1  # Creation operator at momentum k and spin sigma

# # Define the phonon annihilation and creation operators
# b = np.zeros((num_sites, num_branches, num_momenta))  # Boson annihilation operators
# b_dag = np.zeros((num_sites, num_branches, num_momenta))  # Boson creation operators

# # For simplicity, assume diagonal representations
# for alpha in range(num_branches):
#     for q in range(num_momenta):
#         for i in range(num_sites):
#             b[i, alpha, q] = 1  # Boson annihilation operator at site i, branch alpha, q
#             b_dag[i, alpha, q] = 1  # Boson creation operator at site i, branch alpha, q

# # Calculate the expectation value of b_q,lambda + b_q,lambda^dagger
# # Compute the expectation value for each branch and momentum
# expectation_B = np.zeros((num_branches, num_momenta))

# for alpha in range(num_branches):
#     for q in range(num_momenta):
#         # Compute expectation value <psi_phonon | B | psi_phonon>
#         expectation_B[alpha, q] = np.real(
#             np.dot(psi_phonon[:, alpha, q].conj(), np.dot(B[:, alpha, q], psi_phonon[:, alpha, q]))
#         )

# # Now calculate the full interaction Hamiltonian with sum over spin sigma
# V = 1  # Normalization constant or volume factor (example)
# interaction_H = 0

# # Sum over all branches, q, k, and spin states
# for alpha in range(num_branches):
#     for q in range(num_momenta):
#         for k in range(num_k):
#             for sigma in range(num_spin):
#                 # Calculate the electron-phonon interaction term, summing over spin sigma
#                 interaction_H += V * g_lambda[alpha, q, k] * np.real(np.dot(psi_electron.conj(), np.dot(c_dag[k + q, k, sigma], psi_electron))) * expectation_b_bdag[alpha, q]





############ In momentum basis

import numpy as np
######## Electronic Hamiltonian formation
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

# # Random initial values for the expectation values and states
# n_up = np.random.rand(num_k)  # Expectation value of spin-up density (momentum basis)
# n_down = np.random.rand(num_k)  # Expectation value of spin-down density (momentum basis)
# g_lambda = np.random.rand(num_branches, num_momenta, num_k)  # Coupling constants
# expectation_B = np.random.rand(num_branches, num_momenta)  # Expectation value of b_q + b_q^dagger



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

# Term 4: Electron-phonon interaction

# Random initial phonon wavefunction (normalized)
psi_phonon = np.random.rand(num_sites, num_branches, num_momenta)
psi_phonon /= np.linalg.norm(psi_phonon)

# Create the phonon annihilation (b) and creation (b_dag) operators in momentum space
b = np.zeros((num_branches, num_momenta))
b_dag = np.zeros((num_branches, num_momenta))

for p in range(num_branches):
    for q in range(num_momenta):
        b[p, q] = 1  # Diagonal bosonic annihilation operator
        b_dag[p, q] = 1  # Diagonal bosonic creation operator
#### create another loop for number of branches
# Compute the expectation values <b_q + b_q^dagger>
expectation_B = np.zeros((num_branches, num_momenta))

for alpha in range(num_branches):
    for q in range(num_momenta):
        B_operator = b + b_dag  # Sum of annihilation and creation operators
        expectation_B[alpha, q] = np.real(
            np.dot(psi_phonon[:, alpha, q].conj(), np.dot(B_operator[alpha,q], psi_phonon[:, alpha, q]))
        )

g_lambda = np.random.rand(num_branches, num_momenta, num_k)  # Random coupling constants

for alpha in range(num_branches):  # Loop over phonon branches
    for q in range(num_momenta):  # Loop over phonon momentum transfer
        for k in range(num_k):  # Loop over electron momentum states
            for sigma in range(num_spin):  # Loop over spin states
                idx_k = k * num_spin + sigma  # Index for momentum state k and spin σ
                idx_k_plus_q = ((k + q) % num_k) * num_spin + sigma  # Index for k+q, spin σ

                # Add electron-phonon interaction term
                H_k[idx_k_plus_q, idx_k] += (
                    (1 / V)
                    * g_lambda[alpha, q, k]
                    * expectation_B[alpha, q]
                )

                # Add Hermitian conjugate
                H_k[idx_k, idx_k_plus_q] += (
                    (1 / V)
                    * g_lambda[alpha, q, k]
                    * expectation_B[alpha, q]
                )

# Print the resulting Hamiltonian
print("Full electronic Hamiltonian in Momentum Basis:\n", H_k)

