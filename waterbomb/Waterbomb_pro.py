import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

"""
Node layout (with 1-5 as a mountain):

7     6     4

8     1     5

3     0     2
"""




def read_pos_file(fname):
    """
    Is given the directory of the waterbomb and then reads the file
    Returns a list of lists of positions for each node. [[T0], [n1, x1, y1, z1], [n2, x2, y2, z2] ...,
                                                         [T1], [n1, x1, y1, z1], [n2, x2, y2, z2] ...,
                                                                                                  ...]
    """

    file_pos = open(fname+"/pos.dat", 'r')
    lines = file_pos.readlines()

    pos_data = []

    for l in lines:
        if l[0] != '#' and l[0] != ' ' and l[0] != '\n':
            pos_data.append(l.split())

    P = []
    for i in range (len(pos_data)//10):
        temp = []
        for j in range(1, 10):
            p = pos_data[10*i + j][1:4]
            temp.append([float(p[0]), float(p[1]), float(p[2])])
        P.append(temp)

    return P, pos_data[10][0]



def vec_norm(v):
    return v/np.sqrt(v.dot(v))



def approx_theta(node_pos_at_t):
    """
    Accepts a list of lists of node positions
    Returns angle theta in radians. Averages out between angles of all mountain creases
    defined by 4 tuples of nodes:
    (1,5), (1,6), (1,8), (1,0).
    """
    npos = np.array(node_pos_at_t)
    a1 = npos[0] - npos[1]
    b1 = npos[6] - npos[1]
    a2 = npos[8] - npos[1]
    b2 = npos[5] - npos[1]

    cross1 = np.cross(a1,b1)
    cross2 = np.cross(a2,b2)
    normal = vec_norm(np.cross(cross2, cross1)) #defines "direction" of plane of paper

    #find angles of all mountain creases with the normal
    tu = [(1,5), (1,6), (1,8), (1,0)]
    angles = []
    for t in tu:
        direction = vec_norm(npos[t[1]] - npos[t[0]])
        angles.append(np.degrees(np.arccos(np.dot(direction, normal))))

    return angles


def make_plots(Hist, T):
    Hist = Hist/(len(P_data)-1)
    E = -T*np.log(Hist)

    plt.plot(Hist)
    plt.axvline(90, color='k', linestyle='dashed')
    plt.xlabel("Waterbomb state <---(θ)---> Preliminary state")
    plt.ylabel("Normalised PDF")
    plt.title("Probability Density Function vs state (θ)")
    plt.savefig(fname+"/graphs/pdf")
    plt.show()

    plt.plot(E)
    plt.axvline(90, color='k', linestyle='dashed')
    plt.xlabel("Waterbomb state <---(θ)---> Preliminary state)")
    plt.ylabel("Energy (arb. units)")
    plt.title("System energy vs state (θ)")
    plt.savefig(fname+"/graphs/energy")
    plt.show()


T = 30
fname = ["wb_2", "wb_tiny"]
Probs = []
Energies = []

for fn in fname:
    print(f"Considering file '{fn}'")

    waterbomb_thresh = 70
    w_cnt = 0
    preliminary_thresh = 110
    p_cnt = 0
    undef_cnt = 0

    No_bins = 180
    Hist = np.array([0]*No_bins)

    P_data, DT = read_pos_file(fn)
    for i in range (1, len(P_data)):
        p = P_data[i]
        angles = approx_theta(p)
        theta = sum(angles)/len(angles)

        if theta < waterbomb_thresh:
            w_cnt += 1
        elif theta > preliminary_thresh:
            p_cnt += 1
        else:
            undef_cnt += 1
        h_index = min(math.floor((theta/180)*No_bins), No_bins-1)
        Hist[h_index] += 1

    print('\n' + f"Percentage of time spent in each state:" + "\n" + f"{round(w_cnt/(len(P_data)-1)*100, 2)}% waterbomb, {round(p_cnt/(len(P_data)-1)*100, 2)}% preliminary" + '\n')
    Hist = Hist/(len(P_data)-1)
    Probs.append(Hist)
    Energies.append(-T*np.log(Hist))

for P in Probs:
    plt.plot(P)
plt.axvline(90, color = 'k', linestyle = 'dashed')
plt.xlabel("Waterbomb state <---(θ)---> Preliminary state")
plt.ylabel("Normalised PDF")
plt.title("Probability Density Function vs state (θ)")
plt.savefig("Complete_Probability")
plt.show()

Labels = ["Generic", "Diagonals off"]
col = ['b', 'g']
for i in range(len(Energies)):
    plt.plot(Energies[i], col[i], label = Labels[i])
plt.axvline(90, color = 'y', linestyle = 'dashed')
plt.xlabel("State No.1 <---(θ)---> State No.2")
plt.ylabel("Energy (arb. units)")
plt.ylim(110, 200)
plt.title("System energy vs state (θ)")
plt.legend(framealpha = 0.5)
plt.savefig("Complete_Energy")
plt.show()



"""
for p in P_data:
    print('\n')
    for npos in p:
        print(npos)
print(DT)
"""
