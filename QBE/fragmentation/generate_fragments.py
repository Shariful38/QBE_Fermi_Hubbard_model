from typing import Dict
from itertools import product
import numpy as np

class BaseLatticeFragmentsGenerator:

    def __init__(self,
                 fragmentation_lattice: np.ndarray | None = None,
                 lattice_generation_kwargs: Dict | None = None,
                 lattice_representation: str = "momentum"):
        
        self.lattice_generation_kwargs = lattice_generation_kwargs
        self.lattice_representation = lattice_representation

        if fragmentation_lattice is None:

            self.fragmentation_lattice = self.generate_uniform_fragmentation_lattice(lattice_generation_kwargs,
                                                                                     lattice_representation)
        else:

            self.fragmentation_lattice = fragmentation_lattice

        self.ints_to_lattice_map = self.get_lattice_points_int_map(self.fragmentation_lattice)

    def generate_uniform_fragmentation_lattice(self,
                                               lattice_generation_kwargs: Dict):
        
        k0 = lattice_generation_kwargs["k0"]
        num_dims = lattice_generation_kwargs["num_dims"]
        n_max = lattice_generation_kwargs["n_max"]

        # generate an array whose elements are a list of coordinates with dimension num_dims
        # that are in the cartesian product of an evenly-spaced list of points with spacing k0,
        # up to some max value n_max*k0.
        lattice_points = np.asarray([list(elem) for elem in product(k0*np.asarray(range(-n_max - 1, n_max + 1)), repeat= num_dims)])

        return lattice_points
    
    def get_lattice_points_int_map(self, lattice_points):

        lattice_to_ints_map = {n: coordinate for n, coordinate in enumerate(lattice_points)}

        return lattice_to_ints_map

class BosonicLatticeFragmentsGenerator(BaseLatticeFragmentsGenerator):

    def __init__(self,
                 fragmentation_lattice: np.ndarray | None = None,
                 lattice_generation_kwargs: Dict | None = None,
                 lattice_representation: str = "momentum"):
        
        super().__init__(self,
                         lattice_generation_kwargs=lattice_generation_kwargs,
                         fragmentation_lattice=fragmentation_lattice,
                         lattice_representation=lattice_representation)
    
    
    def generate_bosonic_fragments(self,
                                   bosonic_orbital_points: np.ndarray,
                                   threshold_integer: int):
        
        fragmentation_map = {}
        m = threshold_integer
        m_vector = np.full(shape=self.lattice_generation_kwargs["num_dims"],
                           fill_value=m)
        k0_vector = np.full(shape=self.lattice_generation_kwargs["num_dims"],
                            fill_value=self.lattice_generation_kwargs["k0"])

        if self.lattice_representation == "momentum":

            for n in self.ints_to_lattice_map:

                fragmentation_map[n] = []

                for orbital_coordinate in bosonic_orbital_points:

                    if np.abs(orbital_coordinate - self.ints_to_lattice_map[n]) <= np.abs(np.dot(m_vector, k0_vector)) or np.abs(orbital_coordinate - -1*self.ints_to_lattice_map[n]) <= np.abs(np.dot(m_vector, k0_vector)):

                        fragmentation_map[self.ints_to_lattice_map[n]].append(orbital_coordinate)

            #TODO: Implement way to check for and remove duplicate fragments.
            return fragmentation_map
        
        # Putting this here just in case we decide to implement real space fragmentation later.
        elif self.lattice_representation == "real_space":
            
            return None


class FermionicLatticeFragmentsGenerator:

    def __init__(self,
                 fermionic_orbital_points: np.ndarray,
                 fragmentation_lattice: np.ndarray | None = None,
                 lattice_generation_kwargs: Dict | None = None,
                 lattice_representation: str = "momentum"):
        
        super().__init__(self,
                         lattice_generation_kwargs=lattice_generation_kwargs,
                         fragmentation_lattice=fragmentation_lattice,
                         lattice_representation=lattice_representation)

        
        self.fermionic_orbital_points = fermionic_orbital_points


class HybridLatticeFragmentsGenerator(BosonicLatticeFragmentsGenerator,
                                      FermionicLatticeFragmentsGenerator):

    def __init__(self,
                 fermionic_orbital_points: np.ndarray,
                 bosonic_orbital_points: np.ndarray,
                 fragmentation_lattice: np.ndarray | None = None,
                 lattice_generation_kwargs: Dict | None = None,
                 lattice_representation: str = "momentum"):
        
        super().__init__(fermionic_orbital_points=fermionic_orbital_points,
                       bosonic_orbital_points=bosonic_orbital_points,
                       fragmentation_lattice=fragmentation_lattice,
                       lattice_generation_kwargs=lattice_generation_kwargs,
                       lattice_representation=lattice_representation)
        


