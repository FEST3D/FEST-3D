# --------------------------------------------------------------------
# Usage:
# Open visit by typing visit
# Go to controls->Launch CLI
# In the cli window, type 'execfile('visit_open_database_gui.py'
#
# This script will open the .vtk databse files and draw the required
# plots. This script is used to save the time of opening a database,
# adding a plot and displaying the plot.
# 
# The difference between this script and the 
# visit_display_last_iteration.py fil
# is that this is used in a visit
# gui instance. Hence a visit control GUI window is also present.
# --------------------------------------------------------------------
# import sys
# sys.path.append("/usr/local/visit/current/linux-x86_64/lib/site-packages")
# import visit

import visitpostprocessing as vs

plot_variable_data = vs.get_plot_variable_data()

db_name = './output*.vtk database'

vs.setup_visit(plot_variable_data, db_name)

vs.draw_plots()
