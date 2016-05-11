import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import pylab
from refinegrid import read_grid_3D

# NOTE: in numpy arrays, for an array of shape imx, jmx, kmx:
# a[0:imx, 1:jmx, 2:kmx-1] returns an array of shape (imx, jmx-1, kmx-3).
# Hence the index range of a is: ((0, imx-1), (0, jmx-2), (0, kmx-4))

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


p_list = [[0, 0, 0], [0, 0, 4], [-1, 5, 1], [-1, 5, 7], [5, 1, -2], [5, 1, 3], [3, 7, 2], [3, 7, 5]]
p = [[0,0,0], [2,0,0], [2,0,2], [0,0,2], [0,1,0], [2,1,0], [2,1,2], [0,1,2]]


def get_volume_tetrahedron(p1, p2, p3, p4):
    # Choosing p4 as common, the sides being pi-p4, i = 1, 2, 3
    # Converting pi - p4 to numpy arrays, if already not
    A = np.zeros((3,3))
    A[:, 0] = np.asarray(p1) - np.asarray(p4)
    A[:, 1] = np.asarray(p2) - np.asarray(p4)
    A[:, 2] = np.asarray(p3) - np.asarray(p4)
    
    volume_tet = np.linalg.det(A)/6.0
    return volume_tet


def get_volume_hexahedron(p_list):
    if len(p_list) != 8:
        print 'Need 8 points for hexahedra.. Duh! Exiting...'
        return None

    # Source: Split hex into 5 tetrahedron:
    # No assumptions about planarity seem to be made. All cuts were made
    # with a plane containing only 3 vertices at a time.
    # https://ieeexplore.ieee.org/ieee_pilot/articles/06/ttg2009061587/assets/img/article_1/fig_6/large.gif

    # Note: the list of points needs to be in the following order:
    # i, j, k
    # i+1, j, k
    # i+1, j+1, k
    # i, j+1, k
    # i, j, k+1
    # i+1, j, k+1
    # i+1, j+1, k+1
    # i, j+1, k+1

    # The indices of the 5 split tetrahedra can be visualised from the above
    # link. But since the volume of each tetrahedron depends on the determinant 
    # calculated, it is IMPERATIVE to ensure that a "correct" order is followed
    # for the 4 points. The logic to get the "correct" order is explained as 
    # follows (Refer wiki article on parallelepiped):
    # The determinant is taken of a matrix of pi - p4, i = 1, 2, 3. 
    # Graphically it denotes the sides with p4 as common vertex, with direction
    # outward from p4, i.e., directed from p4 to pi, i = 1, 2, 3
    # Hence, if you ensure that  cross(p1-p4, p2-p4) is along p3-p4, then the
    # determinant will be positive.
    # From the above link, a set of 5 tetrahedra was obtained. Each tetrahedra
    # has 4 points, and in the function calls below, care was taken to ensure
    # that the order is observed while passing parameters into the 
    # vol_tetrahedron function
    volume_hex = 0
    volume_hex += get_volume_tetrahedron(p_list[0], p_list[4], p_list[7], p_list[5])
    volume_hex += get_volume_tetrahedron(p_list[6], p_list[7], p_list[5], p_list[2])
    volume_hex += get_volume_tetrahedron(p_list[7], p_list[3], p_list[0], p_list[2])
    volume_hex += get_volume_tetrahedron(p_list[5], p_list[0], p_list[2], p_list[7])
    volume_hex += get_volume_tetrahedron(p_list[0], p_list[1], p_list[5], p_list[2])

    return volume_hex



def get_areas_and_normal_vectors(X,Y,Z):
    # This function gets all the normal vectors of the grid in the
    # format xnx, xny and so on.
    imx, jmx, kmx = X.shape
    print 'imx, jmx, kmx = ', imx, jmx, kmx 

    xnx = np.zeros((imx, jmx-1, kmx-1))
    xny = np.zeros((imx, jmx-1, kmx-1))
    xnz = np.zeros((imx, jmx-1, kmx-1))
    xA = np.zeros((imx, jmx-1, kmx-1))
    ynx = np.zeros((imx-1, jmx, kmx-1))
    yny = np.zeros((imx-1, jmx, kmx-1))
    ynz = np.zeros((imx-1, jmx, kmx-1))
    yA = np.zeros((imx-1, jmx, kmx-1))
    znx = np.zeros((imx-1, jmx-1, kmx))
    zny = np.zeros((imx-1, jmx-1, kmx))
    znz = np.zeros((imx-1, jmx-1, kmx))
    zA = np.zeros((imx-1, jmx-1, kmx))

    print 'shape of normals at first: ', xnx.shape
    # In order to simplify the code and reduce memory requirements
    # the equations have been simplified. The strategy is as follows:
    # Formulae for d1 and d2 as vectors in terms of X, Y, Z are derived
    # Then, we do cross(d1, d2) and then take the components
    # The components denote the respective normal components
    # The area is just the magnitude of the normals so calculated
    # Then the normal vectors are normalised as unit vectors
    d1_x = X[0:imx, 1:jmx, 1:kmx] - X[0:imx, 0:jmx-1, 0:kmx-1]
    d1_y = Y[0:imx, 1:jmx, 1:kmx] - Y[0:imx, 0:jmx-1, 0:kmx-1]
    d1_z = Z[0:imx, 1:jmx, 1:kmx] - Z[0:imx, 0:jmx-1, 0:kmx-1]
    d2_x = X[0:imx, 0:jmx-1, 1:kmx] - X[0:imx, 1:jmx, 0:kmx-1]
    d2_y = Y[0:imx, 0:jmx-1, 1:kmx] - Y[0:imx, 1:jmx, 0:kmx-1]
    d2_z = Z[0:imx, 0:jmx-1, 1:kmx] - Z[0:imx, 1:jmx, 0:kmx-1]
    xnx = d1_y*d2_z - d1_z*d2_y
    xny = d1_z*d2_x - d1_x*d2_z
    xnz = d1_x*d2_y - d1_y*d2_x
    xA = 0.5 * np.sqrt(xnx**2 + xny**2 + xnz**2)
    xnx = xnx/xA
    xny = xny/xA
    xnz = xnz/xA


    d1_x = X[1:imx, 0:jmx, 1:kmx] - X[0:imx-1, 0:jmx, 0:kmx-1]
    d1_y = Y[1:imx, 0:jmx, 1:kmx] - Y[0:imx-1, 0:jmx, 0:kmx-1]
    d1_z = Z[1:imx, 0:jmx, 1:kmx] - Z[0:imx-1, 0:jmx, 0:kmx-1]
    d2_x = X[1:imx, 0:jmx, 0:kmx-1] - X[0:imx-1, 0:jmx, 1:kmx]
    d2_y = Y[1:imx, 0:jmx, 0:kmx-1] - Y[0:imx-1, 0:jmx, 1:kmx]
    d2_z = Z[1:imx, 0:jmx, 0:kmx-1] - Z[0:imx-1, 0:jmx, 1:kmx]
    ynx = d1_y*d2_z - d1_z*d2_y
    yny = d1_z*d2_x - d1_x*d2_z
    ynz = d1_x*d2_y - d1_y*d2_x
    yA = 0.5 * np.sqrt(ynx**2 + yny**2 + ynz**2)
    ynx = ynx/yA
    yny = yny/yA
    ynz = ynz/yA


    d1_x = X[1:imx, 1:jmx, 0:kmx] - X[0:imx-1, 0:jmx-1, 0:kmx]
    d1_y = Y[1:imx, 1:jmx, 0:kmx] - Y[0:imx-1, 0:jmx-1, 0:kmx]
    d1_z = Z[1:imx, 1:jmx, 0:kmx] - Z[0:imx-1, 0:jmx-1, 0:kmx]
    d2_x = X[0:imx-1, 1:jmx, 0:kmx] - X[1:imx, 0:jmx-1, 0:kmx]
    d2_y = Y[0:imx-1, 1:jmx, 0:kmx] - Y[1:imx, 0:jmx-1, 0:kmx]
    d2_z = Z[0:imx-1, 1:jmx, 0:kmx] - Z[1:imx, 0:jmx-1, 0:kmx]
    znx = d1_y*d2_z - d1_z*d2_y
    zny = d1_z*d2_x - d1_x*d2_z
    znz = d1_x*d2_y - d1_y*d2_x
    zA = 0.5 * np.sqrt(znx**2 + zny**2 + znz**2)
    znx = znx/zA
    zny = zny/zA
    znz = znz/zA

  # print 'Random normals at 39,9,0 (python): '
  # print xnx[39][9][0], xny[39][9][0], xnz[39][9][0]
  # print ynx[39][9][0], yny[39][9][0], ynz[39][9][0]
  # print znx[39][9][0], zny[39][9][0], znz[39][9][0]

    return (xnx,xny,xnz,ynx,yny,ynz,znx,zny,znz,xA,yA,zA)


def create_cylindrical_meshgrid(r1, r2, t1, t2, z1, z2, dr=0.1, dt=0.5, dz=0.5):
    # t1, t2 in degrees
    r = np.linspace(r1, r2, round(abs((r2-r1))/dr + 1))
    t = np.linspace(t1, t2, round(abs((t2-t1))/dt + 1))
    z = np.linspace(z1, z2, round(abs((z2-z1))/dz + 1))

    R, T, Z = np.meshgrid(r,t,z,indexing='ij')
    X = R * np.cos(np.pi/180*T)
    Y = R * np.sin(np.pi/180*T)
    return (X,Y,Z)


def create_spherical_meshgrid(r1, r2, t1, t2, p1, p2, dr=0.1, dt=0.5, dp=0.5):
    # t1, t2 in degrees. Range: [-90, 90]
    # p1, p2 also in degrees. Range: [0, 180]
    r = np.linspace(r1, r2, round(abs((r2-r1))/dr + 1))
    t = np.linspace(t1, t2, round(abs((t2-t1))/dt + 1))
    p = np.linspace(p1, p2, round(abs((p2-p1))/dp + 1))

    R, P, T = np.meshgrid(r,p,t,indexing='ij')
    X = R * np.cos(np.pi/180*T) * np.cos(np.pi/180*P)
    Y = R * np.cos(np.pi/180*T) * np.sin(np.pi/180*P)
    Z = R * np.sin(np.pi/180*T)
    return (X,Y,Z)
    

def get_cylindrical_grid_volume(r1, r2, t1, t2, z1, z2, dr, dt, dz):
    X,Y,Z = create_cylindrical_meshgrid(r1,r2,t1,t2, z1, z2, dr, dt, dz)

    n1, n2, n3 = X.shape
    
    volume1 = 0
    volume2 = 0

    for i in range(n1-1):
        for j in range(n2-1):
            for k in range(n3-1):
                point_list = []
                point_list.append([X[i,j,k], Y[i,j,k], Z[i,j,k]])
                point_list.append([X[i+1,j,k], Y[i+1,j,k], Z[i+1,j,k]])
                point_list.append([X[i+1,j+1,k], Y[i+1,j+1,k], Z[i+1,j+1,k]])
                point_list.append([X[i,j+1,k], Y[i,j+1,k], Z[i,j+1,k]])
                point_list.append([X[i,j,k+1], Y[i,j,k+1], Z[i,j,k+1]])
                point_list.append([X[i+1,j,k+1], Y[i+1,j,k+1], Z[i+1,j,k+1]])
                point_list.append([X[i+1,j+1,k+1], Y[i+1,j+1,k+1], Z[i+1,j+1,k+1]])
                point_list.append([X[i,j+1,k+1], Y[i,j+1,k+1], Z[i,j+1,k+1]])
                volume1 += get_volume_hexahedron(point_list)

                # Getting the other diagonal
                point_list = get_other_diagonal(point_list)
                volume2 += get_volume_hexahedron(point_list)

    volume_analytical = 0.5*np.pi/180*(t2-t1)*(r2**2 - r1**2)*(z2 - z1)

    return (volume1, volume2, volume_analytical)


def get_spherical_grid_volume(r1, r2, t1, t2, p1, p2, dr, dt, dp):
    X,Y,Z = create_spherical_meshgrid(r1,r2,t1,t2, p1, p2, dr, dt, dp)

    n1, n2, n3 = X.shape
    
    volume1 = 0
    volume2 = 0

    for i in range(n1-1):
        for j in range(n2-1):
            for k in range(n3-1):
                point_list = []
                point_list.append([X[i,j,k], Y[i,j,k], Z[i,j,k]])
                point_list.append([X[i+1,j,k], Y[i+1,j,k], Z[i+1,j,k]])
                point_list.append([X[i+1,j+1,k], Y[i+1,j+1,k], Z[i+1,j+1,k]])
                point_list.append([X[i,j+1,k], Y[i,j+1,k], Z[i,j+1,k]])
                point_list.append([X[i,j,k+1], Y[i,j,k+1], Z[i,j,k+1]])
                point_list.append([X[i+1,j,k+1], Y[i+1,j,k+1], Z[i+1,j,k+1]])
                point_list.append([X[i+1,j+1,k+1], Y[i+1,j+1,k+1], Z[i+1,j+1,k+1]])
                point_list.append([X[i,j+1,k+1], Y[i,j+1,k+1], Z[i,j+1,k+1]])
                volume1 += get_volume_hexahedron(point_list)
                # Getting the other diagonal
                point_list = get_other_diagonal(point_list)
                volume2 += get_volume_hexahedron(point_list)

    volume_analytical = np.pi/180*(p2 - p1) * (r2**3 - r1**3)/3.0 * \
                          (np.sin(np.pi/180*t2) - np.sin(np.pi/180*t1))

    return (volume1, volume2, volume_analytical)


def get_cylindrical_grid_surface_area(r1, r2, t1, t2, z1, z2, dr, dt, dz):
    X,Y,Z = create_cylindrical_meshgrid(r1,r2,t1,t2, z1, z2, dr, dt, dz)

    n1, n2, n3 = X.shape

    t1_r = t1 * np.pi / 180
    t2_r = t2 * np.pi / 180
    
    print 'Mesh created'
    xnx,xny,xnz,ynx,yny,ynz,znx,zny,znz,xA,yA,zA = get_areas_and_normal_vectors(X,Y,Z)
    print 'Normals and areas obtained'

    SA = 0.0
    SA = xA[0, :, :].sum() + xA[-1, :, :].sum() + yA[:, 0, :].sum() + \
         yA[:, -1,  :].sum() + zA[:, :, 0].sum() + zA[:, :, -1].sum()

    SA_th = 2*(r2 - r1)*(z2 - z1) + (r2**2 - r1**2)*(t2_r - t1_r) + \
            (r1 + r2)*(t2_r - t1_r)*(z2 - z1)

    return (SA, SA_th)


def get_spherical_grid_surface_area(r1, r2, t1, t2, p1, p2, dr, dt, dp):
    X,Y,Z = create_spherical_meshgrid(r1,r2,t1,t2, p1, p2, dr, dt, dp)

    n1, n2, n3 = X.shape

    t1_r = t1 * np.pi / 180
    t2_r = t2 * np.pi / 180
    p1_r = p1 * np.pi / 180
    p2_r = p2 * np.pi / 180
    
    print 'Mesh created'
    xnx,xny,xnz,ynx,yny,ynz,znx,zny,znz,xA,yA,zA = get_areas_and_normal_vectors(X,Y,Z)
    print 'Normals and areas obtained'

    SA = 0.0
    SA = xA[0, :, :].sum() + xA[-1, :, :].sum() + yA[:, 0, :].sum() + \
         yA[:, -1,  :].sum() + zA[:, :, 0].sum() + zA[:, :, -1].sum()

    SA_th = (r2**2 - r1**2)*(t2_r - t1_r) + (r1**2 + r2**2)*(p2_r - p1_r)*(np.cos(t1_r) - np.cos(t2_r)) + \
            (np.sin(t1_r) + np.sin(t2_r))*0.5*(r2**2 - r1**2)*(p2_r - p1_r)

    return (SA, SA_th)


#################### 2D #######################################
def get_volume_quadrilateral(p_list):
    if len(p_list) != 4:
        print 'Need 4 points for quadrilateral.. Duh! Exiting...'
        return None

    p1, p2, p3, p4 = p_list
    # p1 = [x, y]. Hence, we have x = p1[0], y = p1[1]
    # p1 = [i, j]
    # p2 = [i+1, j]
    # p3 = [i+1, j+1]
    # p4 = [i, j+1]

    vol_quad = 0.0
    vol_quad = 0.5 * abs(\
               p1[0]*p4[1] - p1[1]*p4[0] + \
               p4[0]*p3[1] - p4[1]*p3[0] + \
               p3[0]*p2[1] - p3[1]*p2[0] + \
               p2[0]*p1[1] - p2[1]*p1[0]\
               )
    return vol_quad


def create_circular_meshgrid(r1, r2, t1, t2, dr=0.1, dt=0.5):
    # t1, t2 in degrees
    r = np.linspace(r1, r2, round(abs((r2-r1))/dr + 1))
    t = np.linspace(t1, t2, round(abs((t2-t1))/dt + 1))

    R, T = np.meshgrid(r,t,indexing='ij')
    X = R * np.cos(np.pi/180*T)
    Y = R * np.sin(np.pi/180*T)
    return (X,Y)


def get_circular_grid_volume(r1, r2, t1, t2, dr, dt):
    X,Y = create_circular_meshgrid(r1,r2,t1,t2,dr,dt)

    n1, n2 = X.shape
    
    volume1 = 0

    for i in range(n1-1):
        for j in range(n2-1):
            point_list = []
            point_list.append([X[i,j], Y[i,j]])
            point_list.append([X[i+1,j], Y[i+1,j]])
            point_list.append([X[i+1,j+1], Y[i+1,j+1]])
            point_list.append([X[i,j+1], Y[i,j+1]])
            volume1 += get_volume_quadrilateral(point_list)

    volume_analytical = 0.5*np.pi/180*(t2-t1)*(r2**2 - r1**2)

    return (volume1, volume_analytical)


def do_circular_mesh_convergence_study():
    # We are going to do a mesh convergence study, i.e., discretize the same
    # domain by different factors to see if approximate volume converges to
    # analytical volume
    k = np.linspace(1, 0.1, 10)
    dr = k * 0.2
    dt = k * 2
    r1, r2, t1, t2 = (2, 6.5, 0, 120)

    number_elements = [(round(abs((r2-r1))/dr[i] + 1) * \
                       round(abs((t2-t1))/dt[i] + 1)) \
                       for i in range(len(k))]
    vol_list1 = []
    vol_ana_list = []

    f = open('vol_list.csv', 'w')
    for i in range(len(k)):
        print 'Now doing ', i+1
        print 'Number elements: ', number_elements[i]
        print ' '
        v1, va = get_circular_grid_volume(r1, r2, t1, t2, dr[i], dt[i])
        vol_list1.append(v1)
        vol_ana_list.append(va)
        f.write(number_elements[i].__str__() + ', ' + va.__str__() + ', ' + \
            v1.__str__() + '\n')
    
    f.close()
   
   
    f = plt.figure(1)
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    
    plt.title('Circular Grid')
    plt.plot(number_elements, vol_list1, '.-', label='V_Computational_1')
    plt.plot(number_elements, vol_ana_list, '.-', label='V_Analytical')
    plt.xlabel('Number of elements')
    plt.ylabel('Volume of figure')
    plt.legend(loc='lower right')

    fig = plt.figure(2)
    ax = fig.add_subplot(111)
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    X,Y = create_circular_meshgrid(r1,r2,t1,t2,dr[-1],dt[-1])
    ax.scatter(X,Y) 
    pylab.ion()
    plt.show()

    plt.show()
    return (number_elements, vol_list1, vol_ana_list)


#################### 2D  ends###################################



def get_gridfile_volume(gridfile):
    imx, jmx, kmx, X, Y, Z = read_grid_3D(gridfile)

    n1, n2, n3 = X.shape
    
    volume1 = 0
    volume2 = 0

    for i in range(n1-1):
        for j in range(n2-1):
            for k in range(n3-1):
                point_list = []
                point_list.append([X[i,j,k], Y[i,j,k], Z[i,j,k]])
                point_list.append([X[i+1,j,k], Y[i+1,j,k], Z[i+1,j,k]])
                point_list.append([X[i+1,j+1,k], Y[i+1,j+1,k], Z[i+1,j+1,k]])
                point_list.append([X[i,j+1,k], Y[i,j+1,k], Z[i,j+1,k]])
                point_list.append([X[i,j,k+1], Y[i,j,k+1], Z[i,j,k+1]])
                point_list.append([X[i+1,j,k+1], Y[i+1,j,k+1], Z[i+1,j,k+1]])
                point_list.append([X[i+1,j+1,k+1], Y[i+1,j+1,k+1], Z[i+1,j+1,k+1]])
                point_list.append([X[i,j+1,k+1], Y[i,j+1,k+1], Z[i,j+1,k+1]])
                volume1 += get_volume_hexahedron(point_list)
                # Getting the other diagonal
                point_list = get_other_diagonal(point_list)
                volume2 += get_volume_hexahedron(point_list)

    return (volume1, volume2)



######################## Testing functions ###########################
def check_cylindrical_grid():
    X,Y,Z = create_cylindrical_meshgrid(2,4,0,30, 0, 5)

    # We are going to slow plot a random element to check order
    i,j,k = (8,8,8)

    point_list = []
    point_list.append([X[i,j,k], Y[i,j,k], Z[i,j,k]])
    point_list.append([X[i+1,j,k], Y[i+1,j,k], Z[i+1,j,k]])
    point_list.append([X[i+1,j+1,k], Y[i+1,j+1,k], Z[i+1,j+1,k]])
    point_list.append([X[i,j+1,k], Y[i,j+1,k], Z[i,j+1,k]])
    point_list.append([X[i,j,k+1], Y[i,j,k+1], Z[i,j,k+1]])
    point_list.append([X[i+1,j,k+1], Y[i+1,j,k+1], Z[i+1,j,k+1]])
    point_list.append([X[i+1,j+1,k+1], Y[i+1,j+1,k+1], Z[i+1,j+1,k+1]])
    point_list.append([X[i,j+1,k+1], Y[i,j+1,k+1], Z[i,j+1,k+1]])
    
    plot_hex_points_slow_mo(point_list)


def check_spherical_grid():
    X,Y,Z = create_spherical_meshgrid(2,4,0,90,0,180)

    # We are going to slow plot a random element to check order
    i,j,k = (5,4,4)

    point_list = []
    point_list.append([X[i,j,k], Y[i,j,k], Z[i,j,k]])
    point_list.append([X[i+1,j,k], Y[i+1,j,k], Z[i+1,j,k]])
    point_list.append([X[i+1,j+1,k], Y[i+1,j+1,k], Z[i+1,j+1,k]])
    point_list.append([X[i,j+1,k], Y[i,j+1,k], Z[i,j+1,k]])
    point_list.append([X[i,j,k+1], Y[i,j,k+1], Z[i,j,k+1]])
    point_list.append([X[i+1,j,k+1], Y[i+1,j,k+1], Z[i+1,j,k+1]])
    point_list.append([X[i+1,j+1,k+1], Y[i+1,j+1,k+1], Z[i+1,j+1,k+1]])
    point_list.append([X[i,j+1,k+1], Y[i,j+1,k+1], Z[i,j+1,k+1]])
    plot_hex_points_slow_mo(point_list)


def plot_mesh_grid(X, Y, Z):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    ax.view_init(32, -52)
    ax.scatter(X,Y,Z) 
    pylab.ion()
    plt.show()
    #time.sleep(2)
    #ax.scatter(X[0:10, 0:2, 0:2],Y[0:10, 0:2, 0:2],Z[0:10, 0:2, 0:2])
    #fig.canvas.draw()


def plot_hex_points_slow_mo(p_list):
    # Plots points one by one in order, in slow motion to determine
    # the ordering of points. This is to ensure that the order of points
    # as given in the volume_hexahedron is followed
    x = []
    y = []
    z = []

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # We wish to get the axes limits before hand. Otherwise, 
    # throughout the animation, the graph limits will change
    # We get the min and max of x,y,z and then set them as limits
    # but with an extra factor to see them clearly
    max_l = np.array(map(max, zip(*p_list)))
    min_l = np.array(map(min, zip(*p_list)))
    ax_range = max_l - min_l
    max_l = max_l + 0.5*ax_range
    min_l = min_l - 0.5*ax_range

    ax.set_xlim3d(min_l[0], max_l[0])
    ax.set_ylim3d(min_l[1], max_l[1])
    ax.set_zlim3d(min_l[2], max_l[2])
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    ax.view_init(32, -52)
    pylab.ion()
    plt.show()
    for p in p_list:
        x.append(p[0])
        y.append(p[1])
        z.append(p[2])
        ax.scatter(x,y,z)
        fig.canvas.draw()
        time.sleep(1)
    return ax    


def do_cylindrical_mesh_convergence_study():
    # We are going to do a mesh convergence study, i.e., discretize the same
    # domain by different factors to see if approximate volume converges to
    # analytical volume
    k = np.linspace(1, 0.1, 10)
    dr = k * 0.4
    dt = k * 4
    dz = k * 1
    r1, r2, t1, t2, z1, z2  = (2, 4, 0, 30, 0, 5)

    number_elements = [(round(abs((r2-r1))/dr[i] + 1) * \
                       round(abs((t2-t1))/dt[i] + 1) * \
                       round(abs((z2-z1))/dz[i] + 1)) for i in range(len(k))]
    vol_list1 = []
    vol_list2 = []
    vol_ana_list = []

    f = open('vol_list.csv', 'w')
    for i in range(len(k)):
        print 'Now doing ', i+1
        print 'Number elements: ', number_elements[i]
        print ' '
        v1, v2, va = get_cylindrical_grid_volume(r1, r2, t1, t2, z1, z2, dr[i], dt[i], dz[i])
        #v1 = round(v1*1000)/1000.0
        #v2 = round(v2*1000)/1000.0
        #va = round(va*1000)/1000.0
        vol_list1.append(v1)
        vol_list2.append(v2)
        vol_ana_list.append(va)
        f.write(number_elements[i].__str__() + ', ' + va.__str__() + ', ' + \
            v1.__str__() + ', ' + v2.__str__() + '\n')
    
    f.close()
   
   
    vol_list_avg = [0.5*(vol_list1[i] + vol_list2[i]) for i in range(len(vol_list1))]
    f = plt.figure(1)
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    
    plt.title('Cylindrical Grid')
    plt.plot(number_elements, vol_list1, '.-', label='V_Computational_1')
    plt.plot(number_elements, vol_list2, '.-', label='V_Computational_2')
    plt.plot(number_elements, vol_list_avg, '.-', label='V_Computational_avg')
    plt.plot(number_elements, vol_ana_list, '.-', label='V_Analytical')
    plt.xlabel('Number of elements')
    plt.ylabel('Volume of figure')
    plt.legend(loc='lower right')
    plt.show()
    return (number_elements, vol_list1, vol_list2, vol_list_avg, vol_ana_list)


def do_spherical_mesh_convergence_study():
    # We are going to do a mesh convergence study, i.e., discretize the same
    # domain by different factors to see if approximate volume converges to
    # analytical volume
    k = np.linspace(1, 0.1, 10)
    dr = k * 1
    dt = k * 4
    dp = k * 4
    r1, r2, t1, t2, p1, p2  = (2, 4, 0, 90, 0, 180)
    
    number_elements = [(round(abs((r2-r1))/dr[i] + 1) * \
                       round(abs((t2-t1))/dt[i] + 1) * \
                       round(abs((p2-p1))/dp[i] + 1)) for i in range(len(k))]

    vol_list1 = []
    vol_list2 = []
    vol_ana_list = []

    for i in range(len(k)):
        print 'Now doing ', i
        print 'Number elements: ', number_elements[i]
        print ' '
        v1, v2, va = get_spherical_grid_volume(r1, r2, t1, t2, p1, p2, dr[i], dt[i], dp[i])
        vol_list1.append(v1)
        vol_list2.append(v2)
        vol_ana_list.append(va)
    
    vol_list_avg = [0.5*(vol_list1[i] + vol_list2[i]) for i in range(len(vol_list1))]

    f = plt.figure(1)   
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)        
    plt.title('Spherical Grid')
    plt.plot(number_elements, vol_list1, label='V_Computational_1')
    plt.plot(number_elements, vol_list2, label='V_Computational_2')
    plt.plot(number_elements, vol_list_avg, label='V_Computational_avg')
    plt.plot(number_elements, vol_ana_list, label='V_Analytical')
    plt.xlabel('Number of elements')
    plt.ylabel('Volume of figure')
    plt.legend(loc='lower right')
    plt.show()
    return (number_elements, vol_list1, vol_list2, vol_list_avg, vol_ana_list)


def do_cylindrical_area_convergence_study():
    # We are going to do a mesh convergence study, i.e., discretize the same
    # domain by different factors to see if approximate volume converges to
    # analytical volume
    k = np.linspace(1, 0.1, 10)
    dr = k * 0.4
    dt = k * 4
    dz = k * 1
    r1, r2, t1, t2, z1, z2  = (2, 4, 0, 30, 0, 5)

    number_elements = [(round(abs((r2-r1))/dr[i] + 1) * \
                       round(abs((t2-t1))/dt[i] + 1) * \
                       round(abs((z2-z1))/dz[i] + 1)) for i in range(len(k))]
    sA_list = []
    sA_ana_list = []

    f = open('sA_list_cyl.csv', 'w')
    for i in range(len(k)):
        print 'Now doing ', i+1
        print 'Number elements: ', number_elements[i]
        print ' '
        sA, sA_th = get_cylindrical_grid_surface_area(r1, r2, t1, t2, z1, z2, dr[i], dt[i], dz[i])
        sA_list.append(sA)
        sA_ana_list.append(sA_th)
        f.write(number_elements[i].__str__() + ', ' + sA.__str__() + ', ' + \
            sA_th.__str__() + '\n')
    
    f.close()
   
    f = plt.figure(1)
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    
    plt.title('Cylindrical Grid - Area Convergence')
    plt.plot(number_elements, sA_list, '.-', label='sA_Computational')
    plt.plot(number_elements, sA_ana_list, '.-', label='sA_Analytical')
    plt.xlabel('Number of elements')
    plt.ylabel('Total Surface Area of figure')
    plt.legend(loc='lower right')
    plt.show()
    return (number_elements, sA_list, sA_ana_list)


def do_spherical_area_convergence_study():
    # We are going to do a mesh convergence study, i.e., discretize the same
    # domain by different factors to see if approximate volume converges to
    # analytical volume
    k = np.linspace(1, 0.1, 10)
    dr = k * 1
    dt = k * 4
    dp = k * 4
    r1, r2, t1, t2, p1, p2  = (2, 4, 0, 90, 0, 180)
    
    number_elements = [(round(abs((r2-r1))/dr[i] + 1) * \
                       round(abs((t2-t1))/dt[i] + 1) * \
                       round(abs((p2-p1))/dp[i] + 1)) for i in range(len(k))]

    sA_list = []
    sA_ana_list = []

    f = open('sA_list_sph.csv', 'w')
    for i in range(len(k)):
        print 'Now doing ', i
        print 'Number elements: ', number_elements[i]
        print ' '
        sA, sA_th = get_spherical_grid_surface_area(r1, r2, t1, t2, p1, p2, dr[i], dt[i], dp[i])
        sA_list.append(sA)
        sA_ana_list.append(sA_th)
        f.write(number_elements[i].__str__() + ', ' + sA.__str__() + ', ' + \
            sA_th.__str__() + '\n')
    
    f.close()
   
    f = plt.figure(1)
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    
    plt.title('Spherical Grid - Area Convergence')
    plt.plot(number_elements, sA_list, '.-', label='sA_Computational')
    plt.plot(number_elements, sA_ana_list, '.-', label='sA_Analytical')
    plt.xlabel('Number of elements')
    plt.ylabel('Total Surface Area of figure')
    plt.legend(loc='lower right')
    plt.show()
    return (number_elements, sA_list, sA_ana_list)


def get_centroid(point_list, index_list):
    # index list contains indices from 1 to 8
    centroid = np.zeros(np.shape(point_list[0]))
    for i in index_list:
        centroid += point_list[i-1]
    centroid = centroid/len(index_list)
    return centroid
       

def add_3D_line(ax, ps, pe):
    # Requires a figure already
    ax.plot([ps[0], pe[0]], [ps[1], pe[1]], zs=[ps[2], pe[2]])


global labels_and_points
labels_and_points = []

def update_position(e):
    global fig
    global ax
    global labels_and_points
    for label, x, y, z in labels_and_points:
        x2, y2, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
        label.xy = x2,y2
        label.update_positions(fig.canvas.renderer)
    fig.canvas.draw()


def test_areas_and_normals(X, Y, Z):
    print 'Mesh created'
    xnx,xny,xnz,ynx,yny,ynz,znx,zny,znz,xA,yA,zA = get_areas_and_normal_vectors(X,Y,Z)
    print 'Normals and areas obtained'

    # Choose random element
    i,j,k = (5,5,2)


    point_list = []
    point_list.append([X[i,j,k], Y[i,j,k], Z[i,j,k]])
    point_list.append([X[i+1,j,k], Y[i+1,j,k], Z[i+1,j,k]])
    point_list.append([X[i+1,j+1,k], Y[i+1,j+1,k], Z[i+1,j+1,k]])
    point_list.append([X[i,j+1,k], Y[i,j+1,k], Z[i,j+1,k]])
    point_list.append([X[i,j,k+1], Y[i,j,k+1], Z[i,j,k+1]])
    point_list.append([X[i+1,j,k+1], Y[i+1,j,k+1], Z[i+1,j,k+1]])
    point_list.append([X[i+1,j+1,k+1], Y[i+1,j+1,k+1], Z[i+1,j+1,k+1]])
    point_list.append([X[i,j+1,k+1], Y[i,j+1,k+1], Z[i,j+1,k+1]])
    
    global ax
    global fig
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    ax.view_init(32, -52)
    pylab.ion()
    x = []
    y = []
    z = []
    for p in point_list:
        x.append(p[0])
        y.append(p[1])
        z.append(p[2])
    ax.scatter(x,y,z, color='k', s=200)
    labels = ['i,j,k', 'i+1,j,k', 'i+1,j+1,k', 'i,j+1,k', 'i,j,k+1', 'i+1,j,k+1', 'i+1,j+1,k+1', 'i,j+1,k+1']
    
    for txt, xs, ys, zs in zip(labels, x, y, z):
        x2, y2, _ = proj3d.proj_transform(xs,ys,zs, ax.get_proj())
        label = plt.annotate(
            txt, xy = (x2, y2), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc3,rad=0'))
        labels_and_points.append((label, xs, ys, zs))

    # Plotting lines
    edge_list = [[1,2], [2,3], [3,4], [4,1], [5,6], [6,7], [7,8], [8,5], \
                 [1,5], [4,8], [2,6], [3,7]]
    
    for e in edge_list:
        add_3D_line(ax, point_list[e[0] - 1], point_list[e[1] - 1])
   
    labels2 = []
    x2 = []
    y2 = []
    z2 = []
    # Plotting normals
    # Starting point is the centroid of the face
    # Note: index given starting from 1, since it is easier to number
    # the vertices from 1 to 8
    centroid = get_centroid(point_list, [1,5,4,8])
    labels2.append('xn[i]')
    x2.append(centroid[0])
    y2.append(centroid[1])
    z2.append(centroid[2])
    normal = 0.4*abs(X[i+1,j,k]-X[i,j,k])*np.array([xnx[i,j,k], xny[i,j,k], xnz[i,j,k]])
    end_point = centroid + normal
    a = Arrow3D([centroid[0], end_point[0]], [centroid[1], end_point[1]], \
          [centroid[2], end_point[2]], mutation_scale=20, lw=2, arrowstyle="-|>", color = "b")
    ax.add_artist(a)
    ax.scatter(centroid[0], centroid[1], centroid[2])
    xnxmag = np.sqrt(xnx[i,j,k]**2 + xny[i,j,k]**2 + xnz[i,j,k]**2)
    print 'Mag of xn[i,j,k]: ', xnxmag
    print 'xn[i,j,k]: ', xnx[i,j,k], ', ', xny[i,j,k], ', ', xnz[i,j,k]
    print 'yn[i,j,k]: ', ynx[i,j,k], ', ', yny[i,j,k], ', ', ynz[i,j,k]
    print 'zn[i,j,k]: ', znx[i,j,k], ', ', zny[i,j,k], ', ', znz[i,j,k]

    
    centroid = get_centroid(point_list, [2,3,6,7])
    labels2.append('xn[i+1]')
    x2.append(centroid[0])
    y2.append(centroid[1])
    z2.append(centroid[2])
    normal = 0.4*abs(X[i+1,j,k]-X[i,j,k])*np.array([xnx[i+1,j,k], xny[i+1,j,k], xnz[i+1,j,k]])
    end_point = centroid + normal
    a = Arrow3D([centroid[0], end_point[0]], [centroid[1], end_point[1]], \
          [centroid[2], end_point[2]], mutation_scale=20, lw=2, arrowstyle="-|>", color = "b")
    ax.add_artist(a)
    ax.scatter(centroid[0], centroid[1], centroid[2])

    
    centroid = get_centroid(point_list, [1,2,3,4])
    labels2.append('zn[k]')
    x2.append(centroid[0])
    y2.append(centroid[1])
    z2.append(centroid[2])
    normal = 0.4*abs(Z[i,j,k+1]-Z[i,j,k])*np.array([znx[i,j,k], zny[i,j,k], znz[i,j,k]])
    end_point = centroid + normal
    a = Arrow3D([centroid[0], end_point[0]], [centroid[1], end_point[1]], \
          [centroid[2], end_point[2]], mutation_scale=20, lw=2, arrowstyle="-|>", color = "b")
    ax.add_artist(a)
    ax.scatter(centroid[0], centroid[1], centroid[2])
   

    centroid = get_centroid(point_list, [5,6,7,8])
    labels2.append('zn[k+1]')
    x2.append(centroid[0])
    y2.append(centroid[1])
    z2.append(centroid[2])
    normal = 0.4*abs(Z[i,j,k+1]-Z[i,j,k])*np.array([znx[i,j,k+1], zny[i,j,k+1], znz[i,j,k+1]])
    end_point = centroid + normal
    a = Arrow3D([centroid[0], end_point[0]], [centroid[1], end_point[1]], \
          [centroid[2], end_point[2]], mutation_scale=20, lw=2, arrowstyle="-|>", color = "b")
    ax.add_artist(a)
    ax.scatter(centroid[0], centroid[1], centroid[2])

    
    centroid = get_centroid(point_list, [1,2,5,6])
    labels2.append('yn[j]')
    x2.append(centroid[0])
    y2.append(centroid[1])
    z2.append(centroid[2])
    normal = 0.4*abs(Y[i,j+1,k]-Y[i,j,k])*np.array([ynx[i,j,k], yny[i,j,k], ynz[i,j,k]])
    end_point = centroid + normal
    a = Arrow3D([centroid[0], end_point[0]], [centroid[1], end_point[1]], \
          [centroid[2], end_point[2]], mutation_scale=20, lw=2, arrowstyle="-|>", color = "b")
    ax.add_artist(a)
    ax.scatter(centroid[0], centroid[1], centroid[2])

    
    centroid = get_centroid(point_list, [3,4,7,8])
    labels2.append('yn[j+1]')
    x2.append(centroid[0])
    y2.append(centroid[1])
    z2.append(centroid[2])
    normal = 0.4*abs(Y[i,j+1,k]-Y[i,j,k])*np.array([ynx[i,j+1,k], yny[i,j+1,k], ynz[i,j+1,k]])
    end_point = centroid + normal
    a = Arrow3D([centroid[0], end_point[0]], [centroid[1], end_point[1]], \
          [centroid[2], end_point[2]], mutation_scale=20, lw=2, arrowstyle="-|>", color = "b")
    ax.add_artist(a)
    ax.scatter(centroid[0], centroid[1], centroid[2])

    for txt, xs, ys, zs in zip(labels2, x2, y2, z2):
        x2, y2, _ = proj3d.proj_transform(xs,ys,zs, ax.get_proj())
        label = plt.annotate(
            txt, xy = (x2, y2), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'green', alpha = 0.1),
            arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc3,rad=0'))
        labels_and_points.append((label, xs, ys, zs))
    
    fig.canvas.mpl_connect('motion_notify_event', update_position)
    plt.show() 


def test_vec():
    # Just to check if the arrow module is working 
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    #ax.view_init(32, -52)
    pylab.ion()
    
    a = Arrow3D([0,1], [0,1], [0,1], mutation_scale=20, lw=1, arrowstyle="-|>", color = "k")
    ax.add_artist(a)
    plt.show()
    


def get_other_diagonal(p_list):
    # Just a function to change the order of points of p_list
    # p_list was obtained from a xkcd discussion on volume of hexahedron
    order_d = {0:1, 1:5, 2:6, 3:2, 4:0, 5:4, 6:7, 7:3}
    p_list_new = []
    for i in range(8):
        p_list_new.append(p_list[order_d[i]])
    return p_list_new


def change_p_list_order():
    # Just a function to change the order of points of p_list
    # p_list was obtained from a xkcd discussion on volume of hexahedron
    global p_list
    order_d = {0:0, 1:4, 2:5, 3:1, 4:2, 5:6, 6:7, 7:3}
    p_list_new = []
    for i in range(8):
        p_list_new.append(p_list[order_d[i]])
    p_list = p_list_new    
