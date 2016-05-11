# --------------------------------------------------------------------
# Usage: python visit_display_mesh.py <db_filename>
#
# This script will open visit command line interface and show the last
# state in the visit window
# --------------------------------------------------------------------


import visitpostprocessing as vs
import sys
import os


def visit_plot_mesh(db_name):
    vs.launch_visit()
    vs.draw_mesh(db_name)
    vs.draw_plots()
    t = raw_input('Any key to terminate')

if __name__ == '__main__':
    # Usage: python visit_display_mesh.py <db_filename>
    arg_len = len(sys.argv)
    
    if arg_len == 1:
        print 'Not enough arguments'
        print 'Usage: python visit_display_mesh.py <db_filename>'
        sys.exit(1)
    
    try:
        db_name = sys.argv[1]
    except IndexError:
        print 'Database file name necessary'
        print 'Usage: python visit_display_mesh.py <db_filename>'
        sys.exit(1)

    db_name = os.path.join(os.curdir, db_name)
    visit_plot_mesh(db_name)

