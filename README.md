# Quantum Bootstrap Embedding for Fermi-Bose Lattice Models
# NC State University

### Proposal for code directory organization structure:

- QBE/

    - lattice_hamiltonian_generation/

        - generate_bosonic_hamiltonians.py
        - generate_fermionic_hamiltonians.py
        - generate_hybrid_hamiltonians.py

    - mean_field_solvers/

        - boson_mean_field_solvers.py
        - fermion_mean_field_solvers.py
        - hybrid_mean_field_solvers.py

    - fragmentation/

        - generate_fermionic_fragments.py
        - generate_bosonic_fragments.py
        - generate_hybrid_fragments.py
        - generate_embedded_hamiltonians.py

    - circuits/

        - generate_fermionic_circuits.py
        - generate_boson_circuits.py
        - generate_hybrid_circuits.py

    - eigensolvers/

        - fermion_eigensolvers.py
        - boson_eigensolvers.py
        - hybrid_eigensolvers.py

    - be_optimization.py

        - boson_qbe_solvers.py
        - fermion_qbe_solvers.py
        - hybrid_qbe_solvers.py

- tests/

    - bosonic_tests.py
    - fermionic_tests.py
    - hybrid_tests.py

- examples/

    - bosonic_examples/
    - fermionic_tests/
    - hybrid_tests/