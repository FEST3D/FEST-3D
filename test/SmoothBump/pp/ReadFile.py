""" Read tecplot files"""
import numpy as np


#-------- Tecplot class Start ----------#

class Tecplot:
    skiplist = [
                "Wall_distance",
                "u v w",
                "variables",
                "Varlocation",
                "STRANDID",
                "zone",
                "SOLUTIONTIME"
          ]
    Dimension = {}
    data = {}

    def __init__(self, FileNames, VarList):
        CheckList = self.skiplist + VarList
        for filename in FileNames:
            B = FileNames.index(filename)
            N = str(B)
            #print "Reading: ", B
            HeadLength, self.Dimension[B] = self.ReadHeader(filename, CheckList)
            self.data[B] = self.TecplotReader(filename,HeadLength, self.Dimension[B], VarList)



    def ReadHeader(self, file, CheckList):
        HeaderLength = 0
        dimension = []
        with open(file, 'r') as f:
            for line in f:
                if any(x in line for x in CheckList):
                    HeaderLength = HeaderLength + 1
                if 'zone' in line:
                    for str in line.split():
                        s = str.split("=")[-1]
                        if s.isdigit():
                            dimension.append(int(s))
        return HeaderLength, dimension


    def TecplotReader(self, file, HeadCount, Dim, VarList):
        arrays = []
        data = {}
        with open(file, 'r') as a:
            for idx, line in enumerate(a.readlines()):
                if idx < HeadCount:
                    continue
                else:
                    arrays.append([float(s) for s in line.split()])
        arrays = np.concatenate(arrays)
        NodeSize = Dim[0]*Dim[1]*Dim[2]
        CellSize = (Dim[0]-1)*(Dim[1]-1)*(Dim[2]-1)
        data["x"] = arrays[0*NodeSize:1*NodeSize].reshape(Dim, order="F")
        data["y"] = arrays[1*NodeSize:2*NodeSize].reshape(Dim, order="F")
        data["z"] = arrays[2*NodeSize:3*NodeSize].reshape(Dim, order="F")

        Start = 3*NodeSize
        count = 0
        ScaSize = [int(s) - 1 for s in Dim]
        for var in VarList:
            if var=="Velocity":
                data["u"] = arrays[Start+count*CellSize:Start+(count+1)*CellSize].reshape(ScaSize, order="F")
                count = count+1
                data["v"] = arrays[Start+count*CellSize:Start+(count+1)*CellSize].reshape(ScaSize, order="F")
                count = count+1
                data["w"] = arrays[Start+count*CellSize:Start+(count+1)*CellSize].reshape(ScaSize, order="F")
                count = count+1
                data["Velocity"] = np.sqrt(np.power(data["u"],2) + np.power(data["v"],2) + np.power(data["w"],2))
            else:
                data[var] = arrays[Start+count*CellSize:Start+(count+1)*CellSize].reshape(ScaSize, order="F")
                count = count + 1
        return data

#--- end Tecplot class ------------#


#------Control class  start---------#
class File:
    Control = "../system/control.md"
    Scheme  = "../system/fvscheme.md"
    Flow    = "../system/flow.md"
    Layout  = "../system/mesh/layout/layout.md"
    Output  = "../system/output_control.md"
    Resnorm = "../system/res_control.md"

    def __init__(self, FileInput):
        self.LoadFiles = self.GetFileNames(FileInput)

    def GetFileNames(self,FileInput):
        if not self.integer(FileInput[0]):
            #read name of the file to load
            try:
                with open(FileInput) as f:
                    LoadFiles = [line.strip() for line in f]
                    #print "Files to Load:", '\n'.join(LoadFiles)
            except IOError as e:
                print "unable to open file",FileInput
                print "Either file doesn't exist", 
                print " OR no read premissions"
        else:
            LoadFiles = ["../time_directories/"
                    +str(FileInput[0]).zfill(4)
                    +"/process_"+str(i).zfill(2)
                    +'.'+FileInput[2]
                    for i in range(FileInput[1])]
            #print "Files to Load:", '\n'.join(LoadFiles)
        return LoadFiles
            

    def integer(self,value):
        try:
            int(value)
            return True
        except:
            return False




class Key:
    Control =  ['CFL', 
                'RestartNumber', 
                'MaxIter', 
                'SaveIter', 
                'WriteFileFormat', 
                'WriteDataFormat', 
                'ReadFileFormat', 
                'ReadDataFormat', 
                'WritePrecision', 
                'ResnormInterval', 
                'PurgeNumber', 
                'Tolerance',
                'DebugLevel']

    FVScheme = ['FluxScheme',
                'FaceScheme',
                'Limiter',
                'TurbLimiter',
                'TurbulenceModel',
                'TransitionModel',
                'TimeStep',
                'TimeIntScheme',
                'HigherOrderBC']

    Flow =     ['VariableNumber',
                'DensityRef',
                'URef',
                'VRef',
                'WRef',
                'PressureRef',
                'TurbulenceIntensityRef',
                'MuRatioRef',
                'IntermittencyRef',
                'MURef',
                'MUVariation',
                'TRef',
                'SutherlandTemp',
                'Prandtl',
                'Gamma',
                'RGas']



class InputControl:

    def __init__(self,Input):
        self.ControlValueList  = self.ReadFile(File.Control)
        self.FVSchemeValueList = self.ReadFile(File.Scheme)
        self.FlowValueList     = self.ReadFile(File.Flow)
        self.ResnormValueList  = self.ReadFile(File.Resnorm)
        self.BCList            = self.GetBCList(File.Layout)
        self.VarList           = self.GetInputVarList(File.Output)

        self.ControlDict = dict(zip(Key.Control, self.ControlValueList))
        self.SchemeDict  = dict(zip(Key.FVScheme, self.FVSchemeValueList))
        self.FlowDict    = dict(zip(Key.Flow, self.FlowValueList))
        
        self.Files = File(Input)
        self.LoadFiles = self.Files.LoadFiles


    def Out(self):
        print self.ControlDict
        print self.SchemeDict
        print self.FlowDict
        print self.BCList
        print self.VarList
        print '\n'.join(self.ResnormValueList)


    def ReadFile(self, File):
        Remove = ["#", "=", "Control File", "CONFIG", "Flow File","{", "}", "Scheme File"]
        with open(File) as f:
            List = [entry.strip() for entry in f if not any(x in entry for x in Remove) ]
        return filter(None, List)


    def GetBCList(self, BCFile):
        with open(BCFile) as f:
            List = [l.strip().split() for l in f if l[0] is not "#"]
            del List[:2]
        BCList={}
        face = ["imin", "imax", "jmin", "jmax", "kmin", "kmax"]
        for item in List:
            BCList[List.index(item)] = {face[item[3:].index(i)]:int(i) for i in item[3:]}
        return BCList


    def GetInputVarList(self, VarFile):
        VarList=[]
        try:
            with open(VarFile) as f:
                next(f)
                for line in f:
                    if(line.strip()=='}'):
                        break
                    VarList.append(line.strip())
                #print "Variable to Load:", '\n'.join(VarList)
        except IOError as e:
            print "unable to open file", VarFile
            print "Either file doesn't exist", 
            print " OR no read premissions"
        return VarList
#------Control Class  end---------#
