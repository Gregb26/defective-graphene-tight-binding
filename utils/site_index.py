def get_site_index(i, j, s, Nx):
    """
    Function that computes the site index of the atom of type A (s=0) or B (s=1) in the ith unit
    cell in the x direction and jth unit cell in the y direction.

    Inputs:
        i: integer: unit cell index in the x direction
        j: integer: unit cell index in the y direction
        s: 0 or 1: type of atom in the unit cell (A: s=0, B: s=1)
        Nx: integer: number of unit cells in the x direction

    Returns:
        id: integer: index of the atom in question
    """
    id = (j * Nx + i) * 2 + s
    return id

