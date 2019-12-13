import numpy as np
import geometry
import ReadFile
import os
import entropy

if __name__=="__main__":
    LoadLevel = 1
    Blocks    = 2
    Type      = 'dat' #'vtk' or 'dat'
    Input     = ReadFile.InputControl([LoadLevel, Blocks, Type])
    VarList   = Input.VarList
    BCList    = Input.BCList
    FileNames = Input.LoadFiles
    assert Type in ['dat'], 'wrong file type; use: "dat"'
    db = ReadFile.Tecplot(FileNames, VarList) 
    Database = db.data
    entropy.Getentropy(FileNames, Database, Input)

