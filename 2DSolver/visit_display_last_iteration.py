# --------------------------------------------------------------------
# Usage:
# Run this script directly as 'python visit_display_last_iteration.py'
#
# This script will open visit command line interface and show the last
# state in the visit window
# --------------------------------------------------------------------

import visitpostprocessing as vs

plot_variable_data = vs.get_plot_variable_data()

db_name = './output*.vtk database'

vs.launch_visit()

vs.setup_visit(plot_variable_data, db_name)

vs.draw_plots()

vs.set_to_last_state()

t = raw_input('Any key to terminate')
