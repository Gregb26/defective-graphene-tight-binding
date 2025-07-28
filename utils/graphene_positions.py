import numpy as np

def get_positions(Nx, Ny, a):
    """"
    Function that computes the atomic positions of atoms in a Nx x Ny supercell of graphene
    Inputs:
        Nx: integer: number of unit cells in the x direction
        Ny: integer: number of unit cells in the y direction
        a: lattice constant of graphene

    Outputs: [2*Nx*Ny, 2] array: the ith row is the (x, y) coordinates of the ith atom.
    """

    # lattice vectors of graphene
    a1 = (a / 2) * np.array([3, np.sqrt(3)])
    a2 = (a / 2) * np.array([3, -np.sqrt(3)])

    # positions of atoms in a unit cell (2 atoms per unit cell for graphene)
    tau1 = np.array([0, 0])
    tau2 = a * np.array([1, 0])

    positions = []

    for x in range(Nx):
        for y in range(Nx):
            R = x * a1 + y * a2
            positions.append(R + tau1)
            positions.append(R + tau2)

    return np.array(positions)