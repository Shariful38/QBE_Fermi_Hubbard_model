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
U = 2.0  # On-site interaction strength
mu = 1.0  # Chemical potential

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
    H_up[i-1, 0] = H_up[0, i-1] = -1
    H_down[i-1, 0] = H_down[0, i-1] = -1
    

    # Term 2: On-site interaction using Green's function
    # Compute Green's function G = I - Phi' (Phi† Phi')^-1 Phi†
    def compute_green_function(Phi):
        Phi_dagger = Phi.conj().T
        U_matrix=np.random.rand(num_sites, num_sites)
        B_matrix = scipy.linalg.expm(U_matrix)
        Phi_prime = B_matrix @ Phi
        intermediate = np.linalg.inv(Phi_dagger @ Phi_prime)
        #return np.eye(num_sites) - Phi_prime @ intermediate @ Phi_dagger
        return Phi_prime @ intermediate @ Phi_dagger
    G_up = compute_green_function(Phi_up)
    G_down = compute_green_function(Phi_down)

    # Expectation values of number operators
    expectation_n_up = np.real(np.diag(G_up))
    expectation_n_down = np.real(np.diag(G_down))

     # Add interactions based on expectation values
    for i in range(num_sites):
        H_up[i, i] += U / 2 * expectation_n_down[i]
        H_down[i, i] += U / 2 * expectation_n_up[i]

   
    # # Apply chemical potential term to both spin-up and spin-down
    # for i in range(num_sites):
    #     H_up[i, i] += -mu  # Add chemical potential for spin-up states
    #     H_down[i, i] += -mu  # Add chemical potential for spin-down states

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
    print(f"Iteration {iteration + 1}: Energy = {E_curr}")

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
print("Final Slater Determinant (Slater Matrix):\n", slater_matrix)
print("Final expectation values (n_up):", np.real(np.diag(new_psi_up @ new_psi_up.T.conj())))
print("Final expectation values (n_down):", np.real(np.diag(new_psi_down @ new_psi_down.T.conj())))
print("Eigenvalues (spin-up):", eigenvalues_up)
print("Eigenvalues (spin-down):", eigenvalues_down)



##### electron-phonon coupling term

## using coherent state to calculate expectation value

import qutip as qt
import numpy as np

# Define parameters
N = 5  # Truncated Hilbert space per site (phonons)
L = 4  # Number of lattice sites
V = L  # Volume factor, assuming periodic boundary conditions

# Define discrete momentum indices
k_indices = range(L)
q_indices = range(L)

# Define electron creation and annihilation operators in real space
c_list = [qt.tensor([qt.destroy(2) if i == j else qt.qeye(2) for i in range(L)]) for j in range(L)]
c_dag_list = [c.dag() for c in c_list]

# Define phonon annihilation operators in real space
b_list = [qt.tensor([qt.destroy(N) if i == j else qt.qeye(N) for i in range(L)]) for j in range(L)]
b_dag_list = [b.dag() for b in b_list]

# Fourier transform to momentum space for electron operators
c_k = [(1/np.sqrt(L)) * sum(np.exp(-1j * (2 * np.pi * k_idx / L) * j) * c_list[j] for j in range(L)) for k_idx in k_indices]
c_k_dag = [c.dag() for c in c_k]

# Fourier transform to momentum space for phonon operators
b_q = [(1/np.sqrt(L)) * sum(np.exp(-1j * (2 * np.pi * q_idx / L) * j) * b_list[j] for j in range(L)) for q_idx in q_indices]
b_q_dag = [b.dag() for b in b_q]

# Define electron-phonon coupling function g_lambda_kq
g_lambda_kq = lambda k_idx, q_idx: np.exp(-((2 * np.pi * k_idx / L) - (2 * np.pi * q_idx / L))**2)  # Example function

# Create a coherent state instead of number states
alpha = 2.0  # Displacement parameter
coherent_state = qt.coherent(N, alpha)  # Coherent state
bosonic_state_real = qt.tensor([coherent_state for _ in range(L)])

# Compute expectation value in real space
boson_sum_real = sum(b_list[j] + b_dag_list[-j] for j in range(L))


#expectation_value_real = qt.expect(boson_sum_real, bosonic_state_real)
expectation_value_real = qt.expect(boson_sum_real, bosonic_state_real)

print("Expectation Value in Real Space:", expectation_value_real)


expectation_value_momentum = (1 / np.sqrt(L)) * sum(
        np.exp(-1j * 2 * np.pi * q_idx / L) * expectation_value_real
        for q_idx in range(L))

print("Expectation Value in momentum Space:", expectation_value_momentum)

# Define electron-phonon interaction Hamiltonian H_eph
H_eph = sum(
    (1/V) * g_lambda_kq(k_idx, q_idx) / 2 * (
        c_k_dag[(k_idx + q_idx) % L] * c_k[k_idx] * expectation_value_momentum
    )
    for k_idx in k_indices for q_idx in q_indices
)

# Print the corrected Hamiltonian
print("Corrected Electron-Phonon Hamiltonian H_eph:\n", H_eph)


