#---------------------------------------------------------------------
# This program contains the module for plotting any form of output
# log as a function of no. of iterations, and saving it as an image.
# It is assumed that each log (example residue norm) is written to a
# single file.
# 
# Usage:
# 1. The name of the file and its description (What appears as the Y
#    axes label) are given as key, value pairs to the global plt_dict 
#    dictionary below
# 2. Program can be called as python plotandsavelogs.py
# 3. plt_dict = {'filename':'description'}
#
# The program does the following:
# 1. For the key, value pairs in the dictionary, it scans the current
#    directory for the filename. 
# 2. If present, it plots the value of parameter vs iterations and 
#    saves the png file of the same
#---------------------------------------------------------------------

import matplotlib.pyplot as plt
import os

#-------------------- Log file names and description -----------------
plt_dict = {\
'art_disspn':'Artificial Dissipation Parameter', \
'resnorms': 'Resnorm',\
}
#--------------------------------------------------------------------

def plot_and_save(filename, description):
    plotvar = []
    with open(filename, 'r') as f:
        for line in f:
            plotvar.append(float(line.split()[0]))
    plt.plot(plotvar)
    plt.xlabel('Iterations')
    plt.ylabel(description)
    plt.savefig(filename + '.png')

    plt.clf()


def get_all_plots(plt_dict):
    for key in plt_dict.keys():
      if os.stat(key).st_size != 1:
        plot_and_save(key, plt_dict[key])
    

if __name__ == '__main__':
    get_all_plots(plt_dict)
