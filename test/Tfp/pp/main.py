import numpy as np
import geometry
#from ReadFile import InputControl
#from Input import File
import ReadFile
import os
import Surface

if __name__=="__main__":
    LoadLevel = 1
    Blocks    = 2
    Type      = 'dat' #'vtk' or 'dat'
    nDrag     = [1.0, 0.0, 0.0]
    nLift     = [0.0, 1.0, 0.0]
    nSide     = [0.0, 0.0, 1.0]
    Input     = ReadFile.InputControl([LoadLevel, Blocks, Type])
    VarList   = Input.VarList
    BCList    = Input.BCList
    FileNames = Input.LoadFiles
    #GetDatabase(FileNames, VarList)
    assert Type in ['dat', 'vtk'], 'wrong file type: "dat", "vtk"'
    if Type is 'dat':
        db = ReadFile.Tecplot(FileNames, VarList) 
    else:
        db = ReadFile.Vtk(FileNames, VarList) 
    Database = db.data

    FaceData = Surface.Boundary(FileNames, Database).FaceData
    stress   = Surface.Stress(FaceData,Input,nDrag,nLift,nSide)
    #write data
    stress.WriteForces(Input)

