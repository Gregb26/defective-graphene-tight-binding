import numpy as np
from site_index import get_site_index

def fill_pristine_hamiltonian(Nx, Ny, t, neighbours=1):
    """
    Function that fills the tight binding hamiltonian matrix of pristine graphene in real space

    Inputs:
        Nx: integer: number of the unit cells in the x direction
        Ny: integer: number of the unit cells in the y direction
        t: float: hopping parameter
        neighbours: integer: order of neighbours electrons are allowed to hop to. Default is 1,
            i.e. nearest neighbour hopping only.

    Outputs:
        H_p: [Nx, Ny] array: tight binding hamiltonian matrix of pristine graphene in real space.
    """

    N = 2*Nx*Ny # number of atoms in the supercell. Nx*Ny unit cells and two carbon atoms per unit cell

    # nearest neighbours only
    if neighbours == 1:
        H_p = np.zeros((N, N))  # p subscript for pristine
        for x in range(Nx):
            for y in range(Ny):
                # A site index
                A = get_site_index(x, y, 0, Nx)

                # get id of nearest neighbours
                # there are three of them, all B-sites

                # neighbour 1 is in the same unit cell
                B1 = get_site_index(x, y, 1, Nx)
                # neighbour 2 is in unit cell at (x-1, y)
                B2 = get_site_index((x-1)%Nx, y, 1, Nx)
                # neighbour 3 is in unit cell at (x-1, y+1)
                B3 = get_site_index((x-1)%Nx, (y+1)%Ny, 1, Nx)

                # filling the  matrix
                for B in [B1, B2, B3]:
                    H_p[A, B] = -t
                    H_p[B, A] = -t # Hermitian

        return H_p

    # first and second-nearest neighbours only
    if neighbours == 2:

        # first-nearest neighbours
        H_p = fill_pristine_hamiltonian(Nx, Ny, t, neighbours=1)

        # second-nearest neighbours
        for x in range(Nx):
            for y in range(Ny):
                A = get_site_index(x, y, 0, Nx)
                # second neighbor 1 is in the unit cell (i-2, j)
                B1 = get_site_index((x-2)%Nx, y, 1, Nx)
                # second neighbour 2 is in unit cell at (i, j+1)
                B2 = get_site_index(x, (y+1)%Nx, 1, Nx)
                # second neighbour 3 is in unit cell at (i, j-1)
                B3 = get_site_index(x, (y-1)%Nx, 1, Nx)

                # filling matrix
                for B in [B1, B2, B3]:
                    H_p[A, B] = -t
                    H_p[B, A] = -t # Hermitian

        return H_p

def fill_vacancy_hamiltonian(Nx, Ny, vacancies, t, neighbours=1):
    """
    Function that that fills the tight binding hamiltonian matrix of graphene with vacancies in real space.
    Inputs:
        Nx: integer: number of unit cells in the x direction
        Ny: integer: number of unit cells in the y direction
        vacancies: list: list of the indices of atom to remove (vacancies)
        t: float: hopping parameter
        neighbours: integer: order of neighbours electrons are allowed to hop to. Default is 1 (nearest neighbour hopping only.)
    Returns:
        Tight binding hamiltonian matrix of graphene with vacancies in real space.
    """
    H_d = fill_pristine_hamiltonian(Nx, Ny, t, neighbours=neighbours) # d subscript for defective

    for id in vacancies:
        # no hopping allowed to vacant sites
        H_d[id, :] = 0
        H_d[:, id] = 0

    return H_d