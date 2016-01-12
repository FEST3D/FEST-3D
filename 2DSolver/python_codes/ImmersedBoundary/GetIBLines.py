from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# This is a script to generate the required IB lines in the following format
# Each line should contain the start point S, end point E and the outward unit normal vector N
# Hence, the ith line should contain
# Si_x Si_y Ei_x Ei_y Ni_x Ni_y

# Note that E_i = S_i+1

def put_to_file(start_points, end_points, normals, filename):
    # start_points, end_points and normals are vectors -> can be lists
    # or numpy arrays. Should be of the same length.
    # If 2D, then each entry of start_points should be a vector of length 2
    f = open(filename, 'w')
    for i in range(len(start_points)):
        sp = start_points[i]
        ep = end_points[i]
        nv = normals[i]
        line = list(sp) + list(ep) + list(nv)
        for l in line:
            f.write(l.__str__() + ' ')
        f.write('\n')
    f.close()


def get_IB_naca_4_dig(airfoil_str, chord_len= 1.0, num_points=51):
    # Input: airfoil_str: String of 4 numbers
    #        chord_len: length of chord (float or int)
    # Only implemented for symmetric 4 digit
    #TODO: Cambered 4 digit
    #TODO: Validation of input

    t = int(airfoil_str[-2:]) * 0.01
    c = chord_len
    x = np.linspace(0, c, num_points)
    y = 5 * t * c * \
        (0.2969*np.sqrt(x/c) - 0.126*(x/c) - 0.3516*(x/c)**2 + \
         0.2843*(x/c)**3 - 0.1036*(x/c)**4)

 #  plt.plot(x,y,'*')
 #  plt.show()

    # Now to create the start, end and normals list
    # Given n points, the number of line segments is n-1
    # But we want a closed airfoil with the trailing edge matching
    # Hence number of line segments = 2*(n-1)
    # But total number of points = 2*n - 2
    # We need outward normal
    start = []
    end = []
    normals = []

    # First storing the top surface line segments
    for i in range(len(x) - 1):
        start.append([x[i], y[i]])
        end.append([x[i+1], y[i+1]])
        # Given a line segment, find dx = x2 - x1 and dy = y2 - y1
        # The normal vector to the line segment is given by (dy, -dx)
        # We also need to normalise it
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
        ds = np.sqrt(dx**2 + dy**2)
        dx = dx / ds
        dy = dy / ds
        normals.append([-dy, dx])

    # Now we have the top surface line segments from LE to TE
    # Now the bottom surface line segments from TE to LE
    
    for i in range(len(x) - 1, 0, -1):
        start.append([x[i], -y[i]])
        end.append([x[i-1], -y[i-1]])
        # Given a line segment, find dx = x2 - x1 and dy = y2 - y1
        # The normal vector to the line segment is given by (dy, -dx)
        # We also need to normalise it
        dx = x[i-1] - x[i]
        dy = -y[i-1] + y[i]
        ds = np.sqrt(dx**2 + dy**2)
        dx = dx / ds
        dy = dy / ds
        normals.append([-dy, dx])

    return (start, end, normals)


def plot_start_end_normals(start, end, normals):
    # Works only for 2D!!
    #TODO: 3D extension??
    ax = plt.axes()
    scale_to_div = 20
    for i in range(len(start)):
        st = start[i]
        en = end[i]
        no = normals[i]
        no = [(no[0]/scale_to_div), (no[1]/scale_to_div)]

        plt.plot(st[0], st[1], 'ko')
        plt.plot(en[0], en[1], 'k*')
        plt.plot([st[0], en[0]], [st[1], en[1]], 'b')
        cen = [0.5*(st[0] + en[0]), 0.5*(st[1] + en[1])]
        ax.arrow(cen[0], cen[1], no[0], no[1], \
          head_width=0.003, head_length=0.003)
    
    # Not setting axes equal will make normals look like not perpendicular to line segment
    plt.axis('equal')
    plt.show()

