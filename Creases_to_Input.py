import numpy as np
import sys

class Fold_Pattern(object):
    """
    Base Class to store and manipulate a fold pattern. Saves it as a combination of
    initial paper dimensions and a list of folds that are represented as line segments.
    Functionality has been added to create the node, bend, bond etc. files from these line segments.
    Input can be provided from the terminal or given a file.
    """

    def __init__(self):
        """
        Simply initializes a rectangular piece of paper with no folds.

        User is prompted to input aspect ratio (long side (along y-axis) to short side (along x-axis)) and then the length of
        the short side. (e.g. for a rectangle of lengths of 2 along the x-direction and 3 along y, we input {3/2, 2}).

        You can also write the "A4" string for the ratio (with the tall side parallel to y_axis again).
        """

        ar = str(input("Please give aspect ratio as N/M (long/short): "))
        short = float(input("Please give the length of the short side: "))

        if ar == "A4":
            self.aspect_ratio = np.sqrt(2)
        else:
            n = float(eval(ar[(ar.find('/'))]))
            m = float(eval(ar[(ar.find('/')+1):]))
            self.aspect_ratio = float(n/m)

        self.paper_dims = (short, self.aspect_ratio*short)
        self.crease_pattern = []  #line segments ((point1), (point2)) in scaled coordinates





    def add_crease_pp(self, point1, point2, m_v): #given point point, define crease
        """
        A crease can be set, given the end-points of the line (in scaled or normal coordinates coordinates)
        point = string as string "x_coord y_coord"
        mv = +-1
        """
        point1 = point1.split()
        point1 = (eval(point1[0]), eval(point1[1]))
        point2 = point2.split()
        point2 = (eval(point2[0]), eval(point2[1]))
        self.crease_pattern.append(((point1, point2), m_v))





    def add_crease_pal(self, point, angle, length, m_v): #given point, angle with x-axis and length of crease in reduced coordinates
        """
        A crease can be set, given one point, angle with x-axis and length. Length can be given in both scaled and unscaled coordinates
        e.g. points (0, 0) and (1/2, 1/2) in scaled coordinates, define a scaled length of sqrt(2)/2.

        To define length in scaled coordinates write "S" in front of the number (e.g. Ssqrt(2)/2).
        Give angle always in unscaled coordinates.
        """
        point1 = point.split()
        point1 = (eval(point1[0]), eval(point1[1]))

        angle = np.radians(angle)
        g = np.tan(angle)/self.aspect_ratio #an intermediate function
        s = np.sqrt(g**2 / (1 + g**2)) #sine in scaled space
        c = np.sqrt(1 / (1 + g**2)) #cosine in sacled space

        if length[0] == "S":
            l = float(eval(length[1:]))
        else:
            l = length*np.sqrt((np.sin(angle)/np.paper_dims[1])**2 + (np.cos(angle)/np.paper_dims[0])**2)

        point2 = (point1[0]+l*c, point1[1]+l*s)
        if ((abs(point2[0]) < 0) or (abs(point2[0]) > 1) or (abs(point2[1]) < 0) or (abs(point2[1]) > 1)):
            raise Exception("coordinates of point 2 outside of allowed range")

        self.crease_pattern.append(((point1, point2), m_v))


    def crease_set_file(self, filename):
        """
        Same as crease set, but using a file
        """
        filec = open(filename, 'r')
        lines = filec.readlines()
        for i in range (1, len(lines), 5):
            print(f"Reading folds No. {i}")
            mode = int(lines[i])

            if mode == 1:
                p1 = lines[i+1]
                if p1 == "ps":
                    p1 = f"{self.crease_pattern[i-1][0][0][0]} {self.crease_pattern[i-1][0][0][1]}"
                elif p1 == "pe":
                    p1 = f"{self.crease_pattern[i-1][0][1][0]} {self.crease_pattern[i-1][0][1][1]}"

                p2 = lines[i+2]
                m_v = int(lines[i+3])
                print(f"{p1} {p2} {m_v}")
                self.add_crease_pp(p1, p2, m_v)



    def crease_set(self, creases):
        """
        Given the number of creases, it prompts the user to set what these creases are, 1 by 1. Order does not matter

        -A crease can be set, given the end-points of the line (in scaled coordinates) or given one point, angle
        with x-axis and length.

        -The second way of setting a crease, is more tricky, as the length has to be given in scaled coordinates
        e.g. points (0, 0) and (1/2, 1/2), define a length of sqrt(2)/2. Angle is always given in unscaled coordinates.
        """
        print('\n'+"There are two modes for crease input." + '\n' + "-> mode 1: By defining end-points in scaled coordinates" +"\n"+ "-> mode 2: By defining point, angle with x-axis and (scaled) length")

        for i in range (creases):
            print('\n' + f"Crease no. {i+1}")
            mode = int(input("Mode for this crease *(type 1 or 2): "))

            if mode == 1:
                if i >= 1:
                    print("Note: For point 1, you can input 'ps' or 'pe' which will use" + '\n' + "the previous crease's start or end (respectively) as point 1, for this crease")
                p1 = input("Please give point 1 as 'x_coord y_coord': ")
                if p1 == "ps":
                    p1 = f"{self.crease_pattern[i-1][0][0][0]} {self.crease_pattern[i-1][0][0][1]}"
                elif p1 == "pe":
                    p1 = f"{self.crease_pattern[i-1][0][1][0]} {self.crease_pattern[i-1][0][1][1]}"

                p2 = input("Please give point 2 as 'x_coord y_coord': ")
                m_v = int(input("Mountain (type '-1') or valley ('1')?: "))
                self.add_crease_pp(p1, p2, m_v)


            if mode == 2:
                if i >= 1:
                    print("Note: For point 1, you can input 'ps' or 'pe' which will use" + '\n' + "the previous crease's start or end (respectively) as point 1, for this crease")
                p1 = input("Please give point 1 as ('x_coord y_coord'): ")
                if p1 == "ps":
                    p1 = f"{self.crease_pattern[i-1][0][0][0]} {self.crease_pattern[i-1][0][0][1]}"
                elif p1 == "pe":
                    p1 = f"{self.crease_pattern[i-1][0][1][0]} {self.crease_pattern[i-1][0][1][1]}"

                angle = float(input("Please give angle of crease with x-axis: "))
                length = input("Please give the length of the crease: ")
                m_v = int(input("Mountain (type '-1') or valley ('1')?: "))
                self.add_crease_pal(p1, angle, length, m_v)

        self.crease_pattern = np.array(self.crease_pattern)





    def make_nodes_file(self):
        """
        Creates a node file. Implied that all nodes are initially inside the xy plane
        self.nodes = [(posx_i, posy_i)]
        """
        print("Constructing nodes file")

        node_file = open("nodes.inp", 'w')
        mass = 1
        drag = 100

        self.nodes_reduced_coord = []
        corners = np.array([(0,0), (0,1), (1,1), (1,0)]) #scaled coord
        self.nodes = []
        #Check if the nodes at the endpoints of each crease have been already appended. If not, append them
        for i in range(4):

            node_file.write(f"{self.paper_dims[0]*corners[i][0]} {self.paper_dims[1]*corners[i][1]} 0 {drag} {mass}" + '\n')
            self.nodes.append((self.paper_dims[0]*corners[i][0], self.paper_dims[1]*corners[i][1]))
            self.nodes_reduced_coord.append((corners[i][0], corners[i][1]))

        for cr in self.crease_pattern:

            if cr[0][0] not in self.nodes_reduced_coord:
                node_file.write(f"{self.paper_dims[0]*cr[0][0][0]} {self.paper_dims[1]*cr[0][0][1]} 0 {drag} {mass}" + '\n')
                self.nodes.append((self.paper_dims[0]*cr[0][0][0], self.paper_dims[1]*cr[0][0][1]))
                self.nodes_reduced_coord.append(cr[0][0])

            if cr[0][1] not in self.nodes_reduced_coord:
                node_file.write(f"{self.paper_dims[0]*cr[0][1][0]} {self.paper_dims[1]*cr[0][1][1]} 0 {drag} {mass}" + '\n')
                self.nodes.append((self.paper_dims[0]*cr[0][1][0], self.paper_dims[1]*cr[0][1][1]))
                self.nodes_reduced_coord.append(cr[0][1])

        node_file.close()
        #self.nodes = np.array(self.nodes)
        #self.nodes_reduced_coord = np.array(self.nodes_reduced_coord)
        #print(f"nodes in reduced coordinates: {self.nodes_reduced_coord}")
        #print(f"nodes in normal coordinates: {self.nodes}")





    def node_in_line(self, line, point):
        """
        Lines have the form (start_point, vector). If point inside (not including endpoints) (+ some error e)
        """
        e = 0.000001
        point_to_start = point - line[0]
        #if point in line's vector
        if (np.cross(point_to_start, line[1]) == 0) and (np.sqrt(line[1].dot(line[1])) > np.sqrt(point_to_start.dot(point_to_start))):
            if (point_to_start.dot(point_to_start) == 0):
                return False
            else:
                return True
        else:
            return False





    def make_bonds_nodes_file(self):
        """
        Creates a node and bond file. To be called after creases have been set
        It goes through the creases and "notes down" all lines that will have to become bonds, assuming a rectangular piece of paper
        It then goes over these lines and creates a list of nodes and bonds (in reduced coordinates).
        These are then used to write the nodes and bonds files.

        self.bonds_reduced = [(n1, n2, length, mv)]
        self.lines = [(start_point, vector, mv)]
        self.nodes_reduced_coord = [array([x_pos, y_pos])]
        """

        print("Constructing bonds file")
        print("Constructing nodes file")

        k = 1000 #strength of bonds
        R = 0.2 #maximum extension (in scaled coordinates)
        mass = 1
        drag = 10

        corners = np.array([(0,0), (0,1), (1,1), (1,0)]) #scaled coord

        self.lines = []

        """
        Go through all lines (including lines representing bonds that will become folds) that exist, and note them down.

        """

        for i in range (4):
            self.lines.append((corners[i], corners[(i+1)%4] - corners[i], 0))

        for cr in self.crease_pattern:
            self.lines.append((np.array(cr[0][0]), np.array(cr[0][1]) - np.array(cr[0][0]), cr[1]))

            for i in range (len(self.lines)):

                if self.node_in_line(self.lines[i], np.array(cr[0][0])):
                    st_p = self.lines[i][0]
                    vec = self.lines[i][1]
                    mv = self.lines[i][2]
                    self.lines.pop(i)
                    self.lines.append((st_p, np.array(cr[0][0]) - st_p, mv))
                    self.lines.append((np.array(cr[0][0]), st_p + vec - np.array(cr[0][0]), mv))

            for i in range (len(self.lines)):

                if self.node_in_line(self.lines[i], np.array(cr[0][1])):
                    st_p = self.lines[i][0]
                    vec = self.lines[i][1]
                    mv = self.lines[i][2]
                    self.lines.pop(i)
                    self.lines.append((st_p, np.array(cr[0][1]) - st_p, mv))
                    self.lines.append((np.array(cr[0][1]), st_p + vec - np.array(cr[0][1]), mv))

        """
        After you have the lines down, figure out the nodes and bonds between them. Nodes are sorted according to index
        """
        node_index = 0
        self.bonds_reduced = []
        self.nodes_reduced_coord = []

        #For all lines, note down their ends, and if they have been noted again, find and store their indices.
        for l in self.lines:
            if (self.arreqclose_in_list(l[0], self.nodes_reduced_coord)):
                for i in range(len(self.nodes_reduced_coord)):
                    if (self.nodes_reduced_coord[i].size == l[0].size and np.allclose(self.nodes_reduced_coord[i], l[0])):
                        id1 = i
                        break
            else:
                self.nodes_reduced_coord.append(l[0])
                id1 = node_index
                node_index += 1

            if (self.arreqclose_in_list(l[0]+l[1], self.nodes_reduced_coord)):
                for i in range(len(self.nodes_reduced_coord)):
                    if (self.nodes_reduced_coord[i].size == (l[0] + l[1]).size and np.allclose(self.nodes_reduced_coord[i], l[0]+l[1])):
                        id2 = i
                        break
            else:
                self.nodes_reduced_coord.append(l[0]+l[1])
                id2 = node_index
                node_index += 1

            self.bonds_reduced.append((id1, id2, np.sqrt(l[1].dot(l[1])), l[2]))

        #print('\n' + f"Bonds Reduced: {self.bonds_reduced}")
        #print(f"Nodes Reduced: {self.nodes_reduced_coord}")


        """Make the bond and node files"""
        bond_file = open("bonds.inp", 'w')
        node_file = open("nodes.inp", 'w')

        bond_file.write(f"{len(self.bonds_reduced)} NumBonds" + '\n')
        node_file.write(f"{len(self.nodes_reduced_coord)} NumNodes" + '\n')

        for n in self.nodes_reduced_coord:
            node_file.write(f"{self.paper_dims[0]*n[0]} {self.paper_dims[1]*n[1]} 0 {drag} {mass} x_y_z_drag_mass" + '\n')
        node_file.close()

        for b in self.bonds_reduced:
            vec_sc = self.nodes_reduced_coord[b[1]] - self.nodes_reduced_coord[b[0]]
            vec = np.array((vec_sc[0]*self.paper_dims[0], vec_sc[1]*self.paper_dims[1]))
            bond_file.write(f"{b[0]} {b[1]} {k} {R*np.sqrt(vec.dot(vec))} {np.sqrt(vec.dot(vec))} node1_node2_strength_R_D" + '\n')
        bond_file.close()





    #test for approximate equality (for floating point types)
    def arreqclose_in_list(self, myarr, list_arrays):
        return next((True for elem in list_arrays if elem.size == myarr.size and np.allclose(elem, myarr)), False)





    def find_angle_triplet(self, n_id1, n_id2, n_id3):
        """
        Returns angle in degrees between triplet (as seen from above). Right handed angles are assumed positive.
        """
        a = self.nodes_reduced_coord[n_id1] - self.nodes_reduced_coord[n_id2]
        b = self.nodes_reduced_coord[n_id3] - self.nodes_reduced_coord[n_id2]
        a = np.array([a[0]*self.paper_dims[0], a[1]*self.paper_dims[1]])
        b = np.array([b[0]*self.paper_dims[0], b[1]*self.paper_dims[1]])

        si = np.cross(a,b)/(np.sqrt(a.dot(a) * b.dot(b)))
        co = a.dot(b)/(np.sqrt(a.dot(a) * b.dot(b)))
        if co < -1: #rounding errors!!!!
            co = -1

        if si>=0:
            angle = np.degrees(np.arccos(co))
        else:
            angle = np.degrees(-1*np.arccos(co))

        return angle





    def find_adjacent_nodes(self, n_id):
        """
        Given the node id, find all adjacent nodes (closest neighbours). Returns a list with node ids.
        """
        list_neighbours = []
        for b in self.bonds_reduced:
            if b[0] == n_id:
                list_neighbours.append(b[1])
            elif b[1] == n_id:
                list_neighbours.append(b[0])

        return list_neighbours





    def find_angles(self, n_id1, n_id2):
        """
        Given the two of the three node ids, find angles for all the triplets that can be formed. Returns a list with node ids and respective angles.
        Assumes node 2 is the middle one! Node 1 is base, node 2 is the middle
        """

        list_angles = []
        nn = self.find_adjacent_nodes(n_id2)
        for n in nn:
            if n != n_id1:
                angle = self.find_angle_triplet(n_id1, n_id2, n)
                list_angles.append((angle, n)) #angle, node

        return list_angles





    def make_bends_file(self):
        """
        To create the bends file given self.bonds_reduced and self.nodes_reduced_coord
        """
        k = 100
        """
        Find all triplets of adjacent nodes. Examine all nodes one by one and make triplets using couples from their adjacent nodes.
        Note down the nodes participating and the angle defined for each triplet. Then for each node, delete the triplets
        that are not the simplest. By the term simplest, we mean that there is no triplet that defines an angle that is smaller.
        Only reason to exclude a triplet that is simple, is because it defines an angle of 180 degrees (and thus is tangent to some other crease).

        Note: There is no need to bother with the sign of the angle eventually as the energy is cosine dependent.
        """
        print("Constructing bends file")
        self.triplets = []
        for n1 in range (len(self.nodes_reduced_coord)): #for every middle
            used_node_pairs = []

            #print('\n' + f"Checking for node : {n1}")
            nn = self.find_adjacent_nodes(n1)
            #print(f"Found {len(nn)} nearest neighbours: {nn}")

            for n2 in nn: #for every edge node

                triplets_pos = []
                triplets_neg = []

                A = self.find_angles(n2, n1)
                #print(f"For {n2} as base and {n1} as middle, we have the following triplet(s): {A}")
                for a in A:
                    #print(f"Haven't covered triplet {n2} {n1} {a[1]} yet. Appending it")
                    if a[0] > 0:
                        triplets_pos.append([n2, n1, a[1], a[0]]) #node 1, node 2, node 3, angle
                    else:
                        triplets_neg.append([n2, n1, a[1], a[0]]) #node 1, node 2, node 3, angle

                triplets_pos.sort(key=lambda x: x[3])
                triplets_neg.sort(key=lambda x: x[3])

                if triplets_neg:
                    trneg = triplets_neg[-1] #node 1, node 2, node 3, angle
                    if ((trneg[2], trneg[0]) not in used_node_pairs):
                        used_node_pairs.append((trneg[0], trneg[2]))
                        self.triplets.append(trneg)


                if triplets_pos:
                    if triplets_pos[0][3] != 180:
                        trpos = triplets_pos[0] #node 1, node 2, node 3, angle
                        if ((trpos[2], trpos[0]) not in used_node_pairs):
                            used_node_pairs.append((trpos[0], trpos[2]))
                            self.triplets.append(trpos)


            #print('\n' + f"Here are the triplets: {self.triplets}" + '\n')

            bends_file = open("bends.inp", 'w')
            bends_file.write(f"{len(self.triplets)} NumBends" + '\n')
            for t in self.triplets:
                bends_file.write(f"{t[0]} {t[1]} {t[2]} {k} {180.0 - abs(t[3])} node1_node2_node3_strength_eqangle" + '\n')
            bends_file.close()


    def waterbomb_equilibrium(self, theta):
        """Given theta find equilibria for mountain and valley folds"""
        theta = np.radians(theta)
        c = np.arccos(np.cos(theta)**2)

        B = np.arcsin(np.sin(theta)/np.sin(c))
        E = np.cos(np.sin(c)/(np.cos(c) + 1))
        gamma_v = np.arccos(2*np.cos(c) - 1) - np.pi

        if theta <= np.pi/2:
            gamma_m = 2*(B + E - np.pi/2)
        else:
            gamma_m = 2*(np.arccos(1/np.sin(c) - 1/np.tan(c)) - B + np.pi/2)

        return abs(np.degrees(gamma_m)), abs(np.degrees(gamma_v))


    def make_folds_file(self):
        """
        To be called after bends are made
        Makes the folds file.
        Go over all creases (you can go over bonds_reduced, which contains info about whether the bond represents a crease or not
        in the variable mv, which is zero for non-creases, -1 for valleys and 1 for mountains), and then use the triplets to find suitable anchor points
        self.bonds_reduced = [(id1, id2, length, mv)]
        self.triplets = [[idi, idj, idk, angle]
        """
        print("Constructing folds file")
        kh = 1/2
        ke = 1/16
        phi_min = 115
        phi_max = 155

        #FOR Waterbomb only
        gamma_m, gamma_v = self.waterbomb_equilibrium(50)

        phi_eq = 135
        scale_potentials = False #if potential is scaled with size of crease
        self.folds = []

        for b in self.bonds_reduced:
            if b[3] != 0: #if it's supposed to be a crease
                ang0 = self.find_angles(b[0], b[1]) #returns list of [[angle, node]]
                ang1 = self.find_angles(b[1], b[0])
                ang0.sort(key=lambda x: x[0]) #Sort in ascending order according to angles
                ang1.sort(key=lambda x: x[0]) #Sort in ascending order according to angles

                for a0 in ang0:
                    if a0[0] > 0: #can add error message if angle is 180
                        n0 = a0[1]
                        break

                for a1 in ang1:
                    if a1[0] > 0:
                        n3 = a1[1]
                        break

                self.folds.append([n0, b[1], b[0], n3, b[3], b[2]]) #last is the length of the crease in scaled coordinates
        #print(f"The folds are the following: {self.folds}")

        folds_file = open("folds.inp", 'w')
        folds_file.write(f"{len(self.folds)} NumFolds" + '\n')
        for f in self.folds:
            theta = np.arccos((self.nodes_reduced_coord[f[2]][0] - self.nodes_reduced_coord[f[2]][0])/f[5])

            if scale_potentials == False:
                fact = 1
            else:
                fact = f[5]*np.sqrt((np.cos(theta) * self.paper_dims[0])**2 + (np.sin(theta) * self.paper_dims[1])**2)


            if f[4]>0:#if valley (this conditional is only for waterbomb)
                folds_file.write(f"{f[0]} {f[1]} {f[2]} {f[3]} {f[4]} {kh*fact} {ke*fact} {gamma_v-20} {gamma_v+20} {gamma_v} n1_n2_n3_n4_mv_khard_keasy_phimin_phimax_phieq" + '\n')
            else:
                folds_file.write(f"{f[0]} {f[1]} {f[2]} {f[3]} {f[4]} {kh*fact} {ke*fact} {gamma_m-20} {gamma_m+20} {gamma_m} n1_n2_n3_n4_mv_khard_keasy_phimin_phimax_phieq" + '\n')


            #folds_file.write(f"{f[0]} {f[1]} {f[2]} {f[3]} {f[4]} {kh/fact} {ke/fact} {phi_min} {phi_max} {phi_eq} n1_n2_n3_n4_mv_khard_keasy_phimin_phimax_phieq" + '\n')
        folds_file.close()



    def make_wcas(self):
        """
        Excluded volume interactions file. Naive version. Establishes repulsion
        in between all pairs of nodes. Improvement is possible (on paper). No time to do this during project
        """
        print("Constructing (naive) version of wcas file")

        en = 1
        sigma = 1
        EVI = []
        for i in range(len(self.nodes_reduced_coord)):
            for j in range(len(self.nodes_reduced_coord)):
                if i != j:
                    EVI.append([i, j, en, sigma])

        wcas_file = open("wcas.inp", 'w')
        for evi in EVI:
            wcas_file.write(f"{evi[0]} {evi[1]} {evi[2]} {evi[3]} node1_node2_energyscale_sigma" + '\n')

        wcas_file.close()
