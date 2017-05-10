# This module contains all the required functions for post processing
# databases using visit
# All programs shall interface with this module to automate post
# processing in visit
import sys
sys.path.append("/usr/local/visit/current/linux-x86_64/lib/site-packages")
import visit

# Note: Follow the pattern shown as example in setup_visit() to read
# config files

# TODO: Adding subplots made generic
# TODO: Is AddPlot('Mesh', 'mesh') generic?

def get_plot_variable_data():
    plot_variable_data = {}
    
    config_file = open('config.md', 'r')
    # Extracting the value of gamma
    while config_file.readline() != '## gamma (ratio of specific heats)\n':
        pass
    # Read till the you get the line. The next line contains gamma value
    while config_file:
        gamma = config_file.readline()
        if not gamma.startswith('#'):
        	break
    gamma = gamma.strip()
    plot_variable_data['Gamma'] = gamma

    # Extracting the plot type: pseudocolor or vector or contour as
    # given in the config file 
    while config_file.readline() != '### Plot type\n':
        pass
    # The next line is the plot type
    while config_file:
        plottype = config_file.readline()
        if not plottype.startswith('#'):
            break
    plottype = plottype.strip()
    plot_variable_data['plot type'] = plottype
    
    # Leave two more blank lines as per config file
    config_file.readline()
    config_file.readline()
    
    # Exracting the plot variable name and its definition
    while config_file:
        plotvariable= config_file.readline()
        if not plotvariable.startswith('#'):
            break
    plotvariable = plotvariable.strip()
    config_file.close()
    plot_variable_data['plot variable'] = plotvariable
    
    # Extracting the variable definition from the output
    # variables list file
    visit_var_file = open('visit_output_variables.md', 'r')
    while visit_var_file.readline() != plottype + '\n':
        pass
    while visit_var_file:
        line = visit_var_file.readline()
        lst = line.replace(':',' ').replace('"',' ').replace('.',' ').split()
        if lst[1] == plotvariable:
            plotvariabledefinition = lst[2]
            break
    visit_var_file.close()
    # Replacing the value of gamma in the variable definition
    # If variable definition does not contain gamma, 
    # then nothing will happen
    plotvariabledefinition=plotvariabledefinition.replace('Gamma',gamma.__str__())    
    plot_variable_data['plot variable definition'] = plotvariabledefinition
    
    return plot_variable_data


def launch_visit():
    # This function launches the visit window. Used only via 
    # a cli version. Does not launch the entire GUI    
    visit.Launch(vdir='/usr/local/visit/bin/')


def setup_visit(plot_variable_data, db_name):
    print 'Opening database...'
    # Get the complete database name including directory path, db file
    # pattern and extension, followed by the word database
    # Example: db_filename = './output*.vtk database'
    visit.OpenDatabase(db_name)
    
    # Getting the list of current variables
    db_meta_data = visit.GetMetaData(db_name)
    scalar_list = [db_meta_data.GetScalars(i).name for i in range(db_meta_data.GetNumScalars())]
    vector_list = [db_meta_data.GetVectors(i).name for i in range(db_meta_data.GetNumVectors())]
    for vec in vector_list:
        scalar_list.append(vec + '_magnitude')

    # Check if current variable name exists in the database. If no, then
    # res becomes 0. Hence, we need to define a scalar expression for the
    # current variable
    if plot_variable_data['plot variable'] in scalar_list or \
       plot_variable_data['plot variable'] in vector_list:
        visit.AddPlot(plot_variable_data['plot type'], \
            plot_variable_data['plot variable'])
    else:
        print 'No existing variable of that name'
        visit.DefineScalarExpression(plot_variable_data['plot variable'], \
          plot_variable_data['plot variable definition'])
        visit.AddPlot(plot_variable_data['plot type'], \
          plot_variable_data['plot variable'])
    
    if plot_variable_data['plot type'] == 'Pseudocolor':
        p = visit.PseudocolorAttributes()
        p.SetCentering(1)
        p.colorTableName = 'hot_desaturated'
      # p.colorTableName = 'rainbow'
        visit.SetPlotOptions(p)
    # Making mesh transparent and overlaying it on the plot
 #  visit.AddPlot('Mesh', 'mesh')
 #  p = visit.MeshAttributes()
 #  p.SetOpacity(0.25)
 #  visit.SetPlotOptions(p)
    

def draw_plots():
    visit.DrawPlots()


def time_slider_next_state():
    visit.TimeSliderNextState()


def set_to_last_state():
    visit.SetTimeSliderState(visit.TimeSliderGetNStates()-1)


def check_and_advance_state(dbname):
    visit.CheckForNewStates(dbname)
    visit.TimeSliderNextState()


def close_visit():
    visit.Close()
