import numpy as np
from matplotlib import pyplot as plt
from itertools import cycle
from pprint import pprint

density = []
x_speed = []
y_speed = []
pressure = []
x_density_left = {}
x_density_right = {}
y_density_left = {}
y_density_right = {}
x_x_speed_left = {}
x_x_speed_right = {}
y_x_speed_left = {}
y_x_speed_right = {}
x_y_speed_left = {}
x_y_speed_right = {}
y_y_speed_left = {}
y_y_speed_right = {}
x_pressure_left = {}
x_pressure_right = {}
y_pressure_left = {}
y_pressure_right = {}

imx = 97
jmx = 49

cellstyle = cycle(['ro', 'bo', 'go', 'co', 'mo', 'yo', ])
leftstyle = cycle(['b<', 'g<', 'c<', 'm<', 'y<', 'r<', ])
rightstyle = cycle(['g>', 'c>', 'm>', 'y>', 'r>', 'b>', ])


def read_cell_data(filename):
    global density
    global x_speed
    global y_speed
    global pressure

    f = open(filename)
    d = f.read()
    f.close()

    d = d.splitlines()

    d.pop(0)  # 'CELLDATA'
    v = d.pop(0)  # 'Density'
    v = d.pop(0)
    while True:
        try:
            density.append(float(v))
        except:
            break
        v = d.pop(0)

    v = d.pop(0)  # 'Velocity'
    v = d.pop(0)
    while True:
        s = v.split(' ')
        while True:
            try:
                s.remove('')
            except:
                break
        try:
            x_speed.append(float(s[0]))
        except:
            break
        y_speed.append(float(s[1]))
        v = d.pop(0)

    v = d.pop(0)  # 'Pressure'
    v = d.pop(0)
    while True:
        try:
            pressure.append(float(v))
        except:
            break
        if len(d) == 0:
            break
        v = d.pop(0)

    temp = {}
    for j in np.arange(1., jmx):
        temp[j] = {}
        for i in np.arange(1., imx):
            temp[j][i] = density.pop(0)
    density = temp

    temp = {}
    for j in np.arange(1., jmx):
        temp[j] = {}
        for i in np.arange(1., imx):
            temp[j][i] = x_speed.pop(0)
    x_speed = temp

    temp = {}
    for j in np.arange(1., jmx):
        temp[j] = {}
        for i in np.arange(1., imx):
            temp[j][i] = y_speed.pop(0)
    y_speed = temp

    temp = {}
    for j in np.arange(1., jmx):
        temp[j] = {}
        for i in np.arange(1., imx):
            temp[j][i] = pressure.pop(0)
    pressure = temp

def read_face_data(filename):

    global x_density_left
    global x_density_right
    global y_density_left
    global y_density_right
    global x_x_speed_left
    global x_x_speed_right
    global y_x_speed_left
    global y_x_speed_right
    global x_y_speed_left
    global x_y_speed_right
    global y_y_speed_left
    global y_y_speed_right
    global x_pressure_left
    global x_pressure_right
    global y_pressure_left
    global y_pressure_right

    f = open(filename)
    d = f.read()
    f.close()

    d = d.splitlines()

    d.pop(0)  # 'density'
    d.pop(0)  # 'xi_left'
    for j in np.arange(1., jmx-1.+1):
        x_density_left[j] = {}
        for i in np.arange(0-0.5, imx+1-0.5+1):
            x_density_left[j][i] = float(d.pop(0))
    d.pop(0)  # 'xi_right'
    for j in np.arange(1., jmx-1.+1):
        x_density_right[j] = {}
        for i in np.arange(0-0.5, imx+1-0.5+1):
            x_density_right[j][i] = float(d.pop(0))
    d.pop(0)  # 'eta_left'
    for j in np.arange(0-0.5, jmx+1-0.5+1):
        y_density_left[j] = {}
        for i in np.arange(1., imx-1.+1):
            y_density_left[j][i] = float(d.pop(0))
    d.pop(0)  # 'eta_right'
    for j in np.arange(0-0.5, jmx+1-0.5+1):
        y_density_right[j] = {}
        for i in np.arange(1., imx-1.+1):
            y_density_right[j][i] = float(d.pop(0))

    d.pop(0)  # 'x_speed'
    d.pop(0)  # 'xi_left'
    for j in np.arange(1., jmx-1.+1):
        x_x_speed_left[j] = {}
        for i in np.arange(0-0.5, imx+1-0.5+1):
            x_x_speed_left[j][i] = float(d.pop(0))
    d.pop(0)  # 'xi_right'
    for j in np.arange(1., jmx-1.+1):
        x_x_speed_right[j] = {}
        for i in np.arange(0-0.5, imx+1-0.5+1):
            x_x_speed_right[j][i] = float(d.pop(0))
    d.pop(0)  # 'eta_left'
    for j in np.arange(0-0.5, jmx+1-0.5+1):
        y_x_speed_left[j] = {}
        for i in np.arange(1., imx-1.+1):
            y_x_speed_left[j][i] = float(d.pop(0))
    d.pop(0)  # 'eta_right'
    for j in np.arange(0-0.5, jmx+1-0.5+1):
        y_x_speed_right[j] = {}
        for i in np.arange(1., imx-1.+1):
            y_x_speed_right[j][i] = float(d.pop(0))

    d.pop(0)  # 'y_speed'
    d.pop(0)  # 'xi_left'
    for j in np.arange(1., jmx-1.+1):
        x_y_speed_left[j] = {}
        for i in np.arange(0-0.5, imx+1-0.5+1):
            x_y_speed_left[j][i] = float(d.pop(0))
    d.pop(0)  # 'xi_right'
    for j in np.arange(1., jmx-1.+1):
        x_y_speed_right[j] = {}
        for i in np.arange(0-0.5, imx+1-0.5+1):
            x_y_speed_right[j][i] = float(d.pop(0))
    d.pop(0)  # 'eta_left'
    for j in np.arange(0-0.5, jmx+1-0.5+1):
        y_y_speed_left[j] = {}
        for i in np.arange(1., imx-1.+1):
            y_y_speed_left[j][i] = float(d.pop(0))
    d.pop(0)  # 'eta_right'
    for j in np.arange(0-0.5, jmx+1-0.5+1):
        y_y_speed_right[j] = {}
        for i in np.arange(1., imx-1.+1):
            y_y_speed_right[j][i] = float(d.pop(0))

    d.pop(0)  # 'pressure'
    d.pop(0)  # 'xi_left'
    for j in np.arange(1., jmx-1.+1):
        x_pressure_left[j] = {}
        for i in np.arange(0-0.5, imx+1-0.5+1):
            x_pressure_left[j][i] = float(d.pop(0))
    d.pop(0)  # 'xi_right'
    for j in np.arange(1., jmx-1.+1):
        x_pressure_right[j] = {}
        for i in np.arange(0-0.5, imx+1-0.5+1):
            x_pressure_right[j][i] = float(d.pop(0))
    d.pop(0)  # 'eta_left'
    for j in np.arange(0-0.5, jmx+1-0.5+1):
        y_pressure_left[j] = {}
        for i in np.arange(1., imx-1.+1):
            y_pressure_left[j][i] = float(d.pop(0))
    d.pop(0)  # 'eta_right'
    for j in np.arange(0-0.5, jmx+1-0.5+1):
        y_pressure_right[j] = {}
        for i in np.arange(1., imx-1.+1):
            y_pressure_right[j][i] = float(d.pop(0))

    #dirn_list = ['x', 'y']
    #var_list = ['density', 'x_speed', 'y_speed', 'pressure']
    #side_list = ['left', 'right']
    #return (eval(d+'_'+v+'_'+s for v in var_list for d in dirn_list for s in side_list))

def add(var, dirn, row, side='cell'):
    #var = 'density'
    #dirn = 'x'
    #side = 'left'
    #row = 1.

    if dirn == 'x':
        if side == 'cell':
            at = list(np.arange(1., imx))
            values = [eval(var+'[row][i]') for i in at]
        else:
            at = list(np.arange(0-0.5, imx+1-0.5+1))
            values = [eval(dirn+'_'+var+'_'+side+'[row][i]') for i in at]
    elif dirn == 'y':
        if side == 'cell':
            at = list(np.arange(1., jmx))
            values = [eval(var+'[j][row]') for j in at]
        else:
            at = list(np.arange(0-0.5, jmx+1-0.5+1))
            values = [eval(dirn+'_'+var+'_'+side+'[j][row]') for j in at]

    plt.plot(at, values, next(eval(side+'style')))

def init():
    plt.figure()

def show():
    plt.show(block=False)

def close():
    plt.close()
