import os
import sys
import time
import datetime

from daemon import Daemon
from vtkwrite import translate_fortran_to_vtk as translate_fvtk_to_vtk
import visitpostprocessing as vs

class VisDaemon(Daemon):

    def __init__(self, pidfile, track_dir, stdin='/dev/null', 
            stdout='/dev/null', stderr='/dev/null'):
        self.track_dir = track_dir
        self.processed = []
        self.get_grid_file_name()
        self._database_open = False
        self.plot_variable_data = vs.get_plot_variable_data()
        super(VisDaemon, self).__init__(pidfile, stdin, stdout, stderr)

    def get_grid_file_name(self):
        # The grid file name is stored in the config file (config.md)
        config_file = open('config.md', 'r')
        # The grid file entry will have the appropriate header
        while config_file.readline() != '## Grid file\n':
            pass
        # The next line not beginning with a hash (#) is the grid file name
        while config_file:
            self._gridfile = config_file.readline()
            if not self._gridfile.startswith('#'):
                break
        self._gridfile = self._gridfile.strip()
        
    
    def run(self):
        print 'visdaemon listening for files.'
        vs.launch_visit()
        self._db = self.track_dir + '/' + 'output*.vtk database'
        while True:
            #TODO: Generalize filename (output*.fvtk)
            new_files = [f for f in os.listdir(self.track_dir) if \
                    os.path.isfile(self.track_dir + '/' + f) and \
                    f[:6] == 'output' and \
                    f[-5:] == '.fvtk' and \
                    f not in self.processed]
            new_files.sort()
            if new_files:
                next_file = new_files.pop(0)
                print 'Translating file ' + next_file
                output_file = self.track_dir + '/' + \
                        next_file.split('.')[0] + '.vtk'
                translate_fvtk_to_vtk(self.track_dir + '/' + self._gridfile, \
                    self.track_dir + '/' + next_file, \
                    output_file, \
                    'ccfd-iitm solver output')
                self.processed.append(next_file)

                if not self._database_open and len(self.processed) >= 2:
                    vs.setup_visit(self.plot_variable_data, self._db)
                    self._database_open = True                    
                    vs.draw_plots()
                    vs.time_slider_next_state()
                elif self._database_open:
                    print 'Advancing...'
                    vs.check_and_advance_state(self._db)
                
            time.sleep(0.1)

        vs.close_visit()

if __name__ == "__main__":
    pidfile = os.path.dirname(os.path.realpath(__file__)) + '/.visdaemon.pid'

    if len(sys.argv) >= 2:

        try:
            track_dir = os.path.abspath(sys.argv[2])
        except IndexError:
            if sys.argv[1] != 'stop':
                print 'Please specifiy a directory to track.'
                sys.exit(2)
            else:
                track_dir = None
        try:
            log = os.path.abspath(sys.argv[3])
        except IndexError:
            if sys.argv[1] != 'stop':
                print 'No log file set. Sending everything to /dev/null.'
                log = '/dev/null'
            else:
                log = None

        daemon = VisDaemon(pidfile, track_dir, stdout=log, stderr=log)

        if 'start' == sys.argv[1]:
            track_dir = None
            daemon.start()
        elif 'stop' == sys.argv[1]:
            daemon.stop()
        elif 'restart' == sys.argv[1]:
            daemon.restart()
        else:
            print "Unknown command"
            sys.exit(2)
        sys.exit(0)

    else:
        print "usage: %s start|stop|restart" % sys.argv[0]
        sys.exit(2)        
