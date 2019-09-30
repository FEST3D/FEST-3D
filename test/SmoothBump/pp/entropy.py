import numpy as np
import sys
import geometry

def Getentropy(Filenames,Database, Input):

    Rgas=287.0
    gamma=1.4
    factor=Rgas/(gamma-1.)
    pInf=float(Input.FlowDict['PressureRef'])
    rInf=float(Input.FlowDict['DensityRef'])
    sInf = pInf/np.power(rInf, gamma)
    err=[]
    vol=[]
    for B in range(len(Filenames)):
        x = Database[B]["x"]
        y = Database[B]["y"]
        z = Database[B]["z"]
        # calculate cell volume
        Database[B]['Volume'] = geometry.compute_volume(x,y,z)
        Rho = Database[B]['Density']
        P   = Database[B]['Pressure']
        s = P/np.power(Rho, gamma)
        error = np.sqrt(np.sum(np.power((s-sInf)*Database[B]['Volume']/sInf,2)))
        err.append(error)
        vol.append(np.sum(Database[B]['Volume']))

        Ds = np.sqrt(sum(np.power(np.array(err),2))/sum(vol))
    print " ---------- Inviscid Test case: Smooth Bump ---------- "
    print " Flux Scheme        : "+ Input.SchemeDict['FluxScheme']
    print " Higher order method: "+ Input.SchemeDict['FaceScheme']
    print " Turbulence model   : "+ Input.SchemeDict['TurbulenceModel']
    print " Expected Change in entropy           : " + "{:.3E}".format(0.0)
    print " Calculated relative change in entropy: " + "{:.3E}".format(Ds)
    print " Difference                           : " + "{:.3E}".format(Ds*100) + " %"
    print " Allowed Tolerance                    : 0.1 %"
    if Ds < 1e-3:
        print "------------ >>> Test Passed  <<< --------------"
    else:
        print "------------ xxx Test Failed  xxx --------------"

