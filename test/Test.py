import os
import sys
from subprocess import call
import shutil

# How to run this script
# $python Test.py <arg1> <arg2> <arg3>
# <arg1> : Flux scheme
#        allowed option:
#            - ausm
#            - ldfss0
#            - ausmUP
#            - ausmP
#            - slau

# <arg2> : Higher order method
#       allowed options:
#           - muscl
#           - ppm
#           - weno

# <arg3> : Turbulence model
#       allowed options:
#           - sst
#           - sst2003
#           - kkl
#           - sa

# eg.:
# 1. $python Test.py ausm   muscl sst
# 2. $python Test.py slau   weno  sa
# 3. $python Test.py ausmUP ppm   kkl

def SetInput(Flux, HOM, Turbulence, Transition):
    Scheme={}
    Scheme['InviscidFlux'] = Flux
    Scheme['FaceState'] = HOM #Higher order method
    Scheme['Limiter'] = '0 0 0  0 0 0'
    Scheme['TurbulenceLimiter'] = '1 1 1'
    Scheme['TurbulenceModel']=Turbulence
    Scheme['TransitionModel']=Transition
    Scheme['TimeStep']='l'
    Scheme['TimeIntegration']='implicit'
    Scheme['HigherOrderBC']='0'
    return Scheme


def WriteFvSchemeFile(Scheme, RootDir):
    with open(RootDir+"/system/fvscheme.md", "w+") as file:
        file.write("\n")
        file.write("Scheme File\n")
        file.write("===========\n")
        file.write("## Inviscid Flux Scheme\n")
        file.write(str(Scheme['InviscidFlux'])+"\n\n")
        file.write("## Higher Order Method\n")
        file.write(str(Scheme['FaceState'])+"\n\n")
        file.write("## Switch: Limiter - Pressure based switch\n")
        file.write(str(Scheme['Limiter'])+"\n\n")
        file.write("## Turbulence Limiter Switch\n")
        file.write(str(Scheme['TurbulenceLimiter'])+"\n\n")
        file.write("## Turbulence model\n")
        file.write(str(Scheme['TurbulenceModel'])+"\n\n")
        file.write("## Transition model\n")
        file.write(str(Scheme['TransitionModel'])+"\n\n")
        file.write("## Time Step\n")
        file.write(str(Scheme['TimeStep'])+"\n\n")
        file.write("## Time Integration Method\n")
        file.write(str(Scheme['TimeIntegration'])+"\n\n")
        file.write("## Higher Order Boundary Conditions\n")
        file.write(str(Scheme['HigherOrderBC'])+"\n\n")


def run(case,f):
    os.chdir(case)
    call(['bash', "run.sh"], stdout=f, stderr=f)
    os.chdir("../")



if __name__=="__main__":
    Flux = sys.argv[1]
    HOM  = sys.argv[2] #Higher order method
    Turb = sys.argv[3]
    print " "
    print "  ----- Integrated Tests Started -----  "
    print "Total two processes will be used with MPICH library"
    with open("make.log", 'w+') as f:
        call(['make', 'all'], stdout=f, stderr=f)
    with open("Report.txt", "w+") as f:
        print "Running Test number 1  --->  Subsonic flow over a smooth bump"
        Scheme = SetInput(Flux, HOM, 'none', 'none')
        WriteFvSchemeFile(Scheme, 'SmoothBump')
        run('SmoothBump',f)
        print "Running Test number 2  --->  Laminar flow over a flat plate"
        Scheme = SetInput(Flux, HOM, 'none', 'none')
        WriteFvSchemeFile(Scheme, 'Lfp')
        run('Lfp', f)
        print "Running Test number 3  --->  Turbulent flow over a flat plate"
        print "This test might takes a few minutes ..."
        Scheme = SetInput(Flux, HOM, Turb, 'none')
        WriteFvSchemeFile(Scheme, 'Tfp')
        run('Tfp', f)
        print " ----- All tests completed -----"
        print " "

    with open("Report.txt", 'r') as f:
        data = f.read()
        print "Tests passed: ", str(data.count("Passed"))+" out of 3"
        print "Check test summary in 'Report.txt' file."
        print " "
