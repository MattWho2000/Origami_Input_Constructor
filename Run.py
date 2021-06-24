import numpy as np
import sys
from Creases_to_Input import Fold_Pattern


#examples of use of the class. Uncomment and run accordingly :)


"""
#Sequence to set crease pattern from CUI (terminal) input
#suggested to start with this example to get an idea of how the code works
#see comments in "Fold_Pattern.crease_set()" for extra info on how to set a crease pattern
F = Fold_Pattern()
F.crease_set(2) #2 is the number of creases the user is going to be prompted to input
F.make_bonds_nodes_file()
F.make_bends_file()
F.make_folds_file()
"""





#Sequence to set crease pattern from file.
F = Fold_Pattern()
F.crease_set_file("input_creases.txt")
F.make_bonds_nodes_file()
F.make_bends_file()
F.make_folds_file()






"""
Sequence to set the Miura Ori fold pattern
    _______________________
    | \__/  \__/  \__/  \_|
    |_/  \__/  \__/  \__/ |
    | \__/  \__/  \__/  \_|
    |_/  \__/  \__/  \__/ |
    -----------------------
"""
"""
n = 5 #x_axis lattice points
m = 5 #y_axis lattice points
#Note: for valid folding, n and m must be odd

frac = 1/4
x_offset = frac/n
#ang = 10 #defining the angular offset of "vertical" creases from y-axis (in degrees)
#x_offset = (1/(2*m))*np.tan(np.radians(10))

F = Fold_Pattern()
lattice_y = np.linspace(0, 1, m+1)
lattice_x = [np.linspace(0, 1, n+1), np.linspace(0, 1, n+1)]
lattice_x[0][1:n] -= x_offset; lattice_x[1][1:n] += x_offset #establish zigzag pattern

for i in range(1, m):
    for j in range(1, n):
        p1 = (lattice_x[i%2][j], lattice_y[i])
        p21 = (lattice_x[i%2][j-1], lattice_y[i])
        p22 = (lattice_x[(i-1)%2][j], lattice_y[i-1])
        m_v1 = (-1)**(i+j)
        m_v2 = (-1)**(j)
        F.add_crease_pp(f"{p21[0]} {p21[1]}", f"{p1[0]} {p1[1]}", m_v1) #given point1 & point1, define crease
        F.add_crease_pp(f"{p22[0]} {p22[1]}", f"{p1[0]} {p1[1]}", m_v2) #given point1 & point1, define crease
        if i == m-1: #if you have reached topmost level
            p23 = (lattice_x[(i-1)%2][j], lattice_y[i+1])
            m_v3 = (-1)**(j)
            F.add_crease_pp(f"{p23[0]} {p23[1]}", f"{p1[0]} {p1[1]}", m_v3) #given point1 & point1, define crease

    p1 = (lattice_x[i%2][n], lattice_y[i])
    p2 = (lattice_x[i%2][n-1], lattice_y[i])
    m_v = (-1)**(i+n)
    F.add_crease_pp(f"{p2[0]} {p2[1]}", f"{p1[0]} {p1[1]}", m_v) #given point1 & point1, define crease

#print(F.crease_pattern)
F.make_bonds_nodes_file()
F.make_bends_file()
F.make_folds_file()
F.make_wcas()


"""
