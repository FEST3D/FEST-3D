# module setup which solution files to read
# and setup which variables to read from solution file
def GetFileNames(FileInput):
    if not integer(FileInput[0]):
        #read name of the file to load
        try:
            with open(FileInput) as f:
                LoadFiles = [line.strip() for line in f]
                print "Files to Load:", '\n'.join(LoadFiles)
        except IOError as e:
            print "unable to open file",FileInput
            print "Either file doesn't exist", 
            print " OR no read premissions"
    else:
        LoadFiles = ["../time_directories/"
                +str(FileInput[0]).zfill(4)
                +"/process_"+str(i).zfill(2)
                +'.vtk' 
                for i in range(FileInput[1])]
        print "Files to Load:", '\n'.join(LoadFiles)
    return LoadFiles
        

def integer(value):
    try:
        int(value)
        return True
    except:
        return False


def GetInputVarList(VarFile):
    VarList=[]
    try:
        with open(VarFile) as f:
            next(f)
            next(f)
            for line in f:
                if(line.strip()=='}'):
                    break
                VarList.append(line.strip())
            print "Variable to Load:", '\n'.join(VarList)
    except IOError as e:
        print "unable to open file", VarFile
        print "Either file doesn't exist", 
        print " OR no read premissions"
    return VarList


def GetBCList(BCFile):
    with open(BCFile) as f:
        List = [l.strip().split() for l in f if l[0] is not "#"]
        del List[:2]
    BCList={}
    face = ["imin", "imax", "jmin", "jmax", "kmin", "kmax"]
    for item in List:
        BCList[List.index(item)] = {face[item[3:].index(i)]:int(i) for i in item[3:]}
    return BCList


    
