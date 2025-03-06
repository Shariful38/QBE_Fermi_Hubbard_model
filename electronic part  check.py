##### Total electronic hamiltonian

import numpy as np
import scipy.linalg
import qutip as qt

# Define system parameters
num_sites = 4  # Number of lattice sites
num_spin = 2  # Number of spin states (spin-up and spin-down)
num_electrons = M = num_sites  # Half-filling
N = 5  # Truncated phonon Hilbert space per site
V = num_sites  # Volume factor for normalization
t = 1.0  # Hopping parameter
U = 2.0  # On-site interaction strength
mu = 1.0  # Chemical potential
alpha = 2.0  # Coherent state displacement parameter
delta_threshold = 1e-6  # Convergence criterion

# Initialize electron wavefunctions (normalized Slater determinants)
def normalize(psi):
    for i in range(psi.shape[1]):
        psi[:, i] /= np.linalg.norm(psi[:, i])
    return psi

psi_up = normalize(np.random.rand(num_sites, num_electrons // 2))
psi_down = normalize(np.random.rand(num_sites, num_electrons // 2))

# Initialize bosonic state in a coherent state
coherent_state = qt.coherent(N, alpha)
bosonic_state_real = qt.tensor([coherent_state for _ in range(num_sites)])

# Define electron and phonon operators in real space
c_list = [qt.tensor([qt.destroy(2) if i == j else qt.qeye(2) for i in range(num_sites)]) for j in range(num_sites)]
c_dag_list = [c.dag() for c in c_list]
b_list = [qt.tensor([qt.destroy(N) if i == j else qt.qeye(N) for i in range(num_sites)]) for j in range(num_sites)]
b_dag_list = [b.dag() for b in b_list]

# Define discrete momentum indices
k_indices = range(L)
q_indices = range(L)

# Fourier transform to momentum space for electron operators
c_k = [(1/np.sqrt(L)) * sum(np.exp(-1j * (2 * np.pi * k_idx / L) * j) * c_list[j] for j in range(L)) for k_idx in k_indices]
c_k_dag = [c.dag() for c in c_k]

# Fourier transform to momentum space for phonon operators
b_q = [(1/np.sqrt(L)) * sum(np.exp(-1j * (2 * np.pi * q_idx / L) * j) * b_list[j] for j in range(L)) for q_idx in q_indices]
b_q_dag = [b.dag() for b in b_q]

# Define electron-phonon coupling function g_lambda_kq
g_lambda_kq = lambda k_idx, q_idx: np.exp(-((2 * np.pi * k_idx / L) - (2 * np.pi * q_idx / L))**2)  # Example function

# Self-consistent field (SCF) loop
iteration = 0
E_prev = 0
while True:
    # Initialize Hamiltonian matrices
    H_up = np.zeros((num_sites, num_sites), dtype=complex)
    H_down = np.zeros((num_sites, num_sites), dtype=complex)



    # Add hopping terms
    for i in range(num_sites-1):
        H_up[i, i+1] = H_up[i+1, i] = -t
        H_down[i, i+1] = H_down[i+1, i] = -t
    H_up[-1, 0] = H_up[0, -1] = -t  # Periodic boundary
    H_down[-1, 0] = H_down[0, -1] = -t

    # Compute Green's function for interaction terms
    def compute_green_function(Phi):
        Phi_dagger = Phi.conj().T
        U_matrix = np.random.rand(num_sites, num_sites)
        U_matrix = (U_matrix + U_matrix.T) / 2  # Make U_matrix Hermitian
        B_matrix = scipy.linalg.expm(U_matrix)
        Phi_prime = B_matrix @ Phi
        intermediate = np.linalg.inv(Phi_dagger @ Phi_prime)
        return Phi_prime @ intermediate @ Phi_dagger

    G_up = compute_green_function(psi_up)
    G_down = compute_green_function(psi_down)

    # Compute expectation values of number operators
    expectation_n_up = np.real(np.diag(G_up))
    expectation_n_down = np.real(np.diag(G_down))

    # Add interaction terms
    for i in range(num_sites):
        H_up[i, i] += U / 2 * expectation_n_down[i]
        H_down[i, i] += U / 2 * expectation_n_up[i]
    for i in range(num_sites):
      H_up[i, i] += -mu  # Add chemical potential for spin-up states
      H_down[i, i] += -mu  # Add chemical potential for spin-down states


    # Compute expectation value of bosonic field
    boson_sum_real = sum(b_list[j] + b_dag_list[-j] for j in range(num_sites))
    expectation_value_real = qt.expect(boson_sum_real, bosonic_state_real)
    expectation_value_momentum = (1 / np.sqrt(L)) * sum(
        np.exp(-1j * (2 * np.pi * q_idx / L) * j) * expectation_value_real
        for j in range(L) for q_idx in q_indices)
    H_eph = sum(
          (1/V) * g_lambda_kq(k_idx, q_idx) / 2 * (
          c_k_dag[(k_idx + q_idx) % L] * c_k[k_idx] * expectation_value_momentum
            )
    for k_idx in k_indices for q_idx in q_indices
      )
    H_eph = (H_eph + H_eph.dag()) / 2  # Ensure Hermiticity


    #print("H_up shape:", H_up.shape)
    #print("H_eph shape:", H_eph.shape)


    # Total Hamiltonian


    identity_electron = np.eye(4)  # Identity matrix for the electron part
    H_up_extended = np.kron(H_up, identity_electron)  # Expand H_up to the full Hilbert space (16, 16)
    H_down_extended = np.kron(H_down, identity_electron)  # Expand H_down to the full Hilbert space (16, 16)

    # Now you can add H_eph, as their dimensions match
    H_total_up = H_up_extended + H_eph.full()
    H_total_down = H_down_extended + H_eph.full()
    H_total_up = (H_total_up + H_total_up.conj().T) / 2  # Ensure Hermiticity
    H_total_down = (H_total_down + H_total_down.conj().T) / 2 # Ensure Hermiticity




    #print("H_up", H_total_up)
    #print("H_down", H_total_down)




    # Diagonalize the Hamiltonian
    eigenvalues_up, eigenvectors_up = np.linalg.eigh(H_total_up)
    eigenvalues_down, eigenvectors_down = np.linalg.eigh(H_total_down)

    # Update wavefunctions
    num_occupied = num_electrons // 2
    # Extract only the electronic part of the eigenvectors (first num_sites = 4 components)
    eigenvectors_up_projected = eigenvectors_up[:num_sites, :num_occupied]

    # Normalize and update psi_up
    psi_up = normalize(eigenvectors_up_projected)

    # similiarly for spin_down
    eigenvectors_down_projected = eigenvectors_down[:num_sites, :num_occupied]

    # Normalize and update psi_down
    psi_down = normalize(eigenvectors_down_projected)

    # Compute total energy
    E_curr = np.real(np.trace(H_up @ G_up) + np.trace(H_down @ G_down))
    print(f"Iteration {iteration + 1}: Energy = {E_curr}")

    # Check convergence
    if abs(E_curr - E_prev) < delta_threshold:
        print(f"SCF Converged after {iteration + 1} iterations.")
        break

    E_prev = np.real(E_curr)
    iteration += 1


