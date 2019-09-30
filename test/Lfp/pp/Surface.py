import numpy as np
import geometry
import ReadFile
import os

class Boundary:
    vol=[]
    FaceData={}
    def __init__(self,FileNames, Database):
        for B in range(len(FileNames)):
            N = str(B)
            x = Database[B]["x"]
            y = Database[B]["y"]
            z = Database[B]["z"]
            # calculate cell volume
            Database[B]['Volume'] = geometry.compute_volume(x,y,z)
            Database[B]['CellCentroid'] = geometry.compute_centroid(x,y,z)
            # Calculate face area
            Database[B]['in'] = geometry.compute_i_face_area(x,y,z)
            Database[B]['jn'] = geometry.compute_j_face_area(x,y,z)
            Database[B]['kn'] = geometry.compute_k_face_area(x,y,z)
            Database[B]['iArea'] = np.linalg.norm(Database[B]['in'], axis=3)
            Database[B]['jArea'] = np.linalg.norm(Database[B]['jn'], axis=3)
            Database[B]['kArea'] = np.linalg.norm(Database[B]['kn'], axis=3)
            Database[B]['iArea'] = Database[B]['iArea'].reshape(Database[B]['iArea'].shape+(1,), order='F')
            Database[B]['jArea'] = Database[B]['jArea'].reshape(Database[B]['jArea'].shape+(1,), order='F')
            Database[B]['kArea'] = Database[B]['kArea'].reshape(Database[B]['kArea'].shape+(1,), order='F')
            # Calculate face normal
            Database[B]['in'] = Database[B]['in']/Database[B]['iArea']
            Database[B]['jn'] = Database[B]['jn']/Database[B]['jArea']
            Database[B]['kn'] = Database[B]['kn']/Database[B]['kArea']

            # Initiallize face data base
            self.FaceData[B] = {}
            self.FaceData[B]["imin"] = {}
            self.FaceData[B]["jmin"] = {}
            self.FaceData[B]["kmin"] = {}
            self.FaceData[B]["imax"] = {}
            self.FaceData[B]["jmax"] = {}
            self.FaceData[B]["kmax"] = {}
        for key0, value0  in Database.items():
            for key, value in value0.items():
                self.FaceData[key0]["imin"][key] = value[0,:,:]
                self.FaceData[key0]["jmin"][key] = value[:,0,:]
                self.FaceData[key0]["kmin"][key] = value[:,:,0]
                self.FaceData[key0]["imax"][key] = value[-1,:,:]
                self.FaceData[key0]["jmax"][key] = value[:,-1,:]
                self.FaceData[key0]["kmax"][key] = value[:,:,-1]
        for key, value in self.FaceData.items():
            self.FaceData[key]['imax']['in'] = -self.FaceData[key]['imax']['in']
            self.FaceData[key]['jmax']['jn'] = -self.FaceData[key]['jmax']['jn']
            self.FaceData[key]['kmax']['kn'] = -self.FaceData[key]['kmax']['kn']
        self.vol.append(np.sum(Database[B]["Volume"]))
        #print sum(self.vol)


class Stress:

    def __init__(self, FaceData, Input, nD, nL, nS): 
        self.URef = float(Input.FlowDict['URef'])
        self.VRef = float(Input.FlowDict['VRef'])
        self.WRef = float(Input.FlowDict['WRef'])
        self.RhoRef = float(Input.FlowDict['DensityRef'])
        self.PRef = float(Input.FlowDict['PressureRef'])
        self.MuRef= float(Input.FlowDict['MURef'])
        self.Velocity = np.sqrt(self.URef*self.URef + self.VRef*self.VRef + self.WRef*self.WRef)
        self.DynamicVelocity = 0.5*self.RhoRef*self.Velocity*self.Velocity
        self.nDrag = nD
        self.nLift = nL
        self.nSide = nS
        self.BCList = Input.BCList
        self.FaceData = FaceData



    def WriteForces(self, Input):
        cd = []
        cl = []
        try:
            os.remove('surfaceData.dat')
        except OSError:
            pass
        try:
            os.remove('IntsurfaceData.txt')
        except OSError:
            pass
        ExpectedDragCoeff = 0.00133
        f = open('surfaceData.dat', 'ab')
        headerSurfaceData = "variables=X REX CF"
        f.write(headerSurfaceData+'\n')
        for block,  value in self.BCList.items():
            for face, BC in value.items():
                if BC==-5:
                    #print FaceData[block][face]['']
                    PDiff = (self.FaceData[block][face]['Pressure']- self.PRef)
                    x =0.5*(self.FaceData[block][face]['x'][1:,0] + self.FaceData[block][face]['x'][:-1,0])
                    nx   = self.FaceData[block][face]['jn'][:,:,0]
                    ny   = self.FaceData[block][face]['jn'][:,:,1]
                    nz   = self.FaceData[block][face]['jn'][:,:,2]
                    Area = self.FaceData[block][face]['jArea'][:,:,0]
                    dudx = self.FaceData[block][face]['Dudx']
                    dudy = self.FaceData[block][face]['Dudy']
                    dudz = self.FaceData[block][face]['Dudz']
                    dvdx = self.FaceData[block][face]['Dvdx']
                    dvdy = self.FaceData[block][face]['Dvdy']
                    dvdz = self.FaceData[block][face]['Dvdz']
                    dwdx = self.FaceData[block][face]['Dwdx']
                    dwdy = self.FaceData[block][face]['Dwdy']
                    dwdz = self.FaceData[block][face]['Dwdz']
                    delv = dudx + dvdy + dwdz
                    mu   = self.FaceData[block][face]['Mu']
                    txx = mu*((dudx + dudx) - 2.0*delv/3.0)
                    tyy = mu*((dvdy + dvdy) - 2.0*delv/3.0)
                    tzz = mu*((dwdz + dwdz) - 2.0*delv/3.0)
                    txy = mu*((dudy + dvdx))
                    tyz = mu*((dwdy + dvdz))
                    txz = mu*((dudz + dwdx))
                    tyx = txy
                    tzy = tyz
                    tzx = txz
                    Fx = (txx*nx + txy*ny + txz*nz)
                    Fy = (tyx*nx + tyy*ny + tyz*nz)
                    Fz = (tzx*nx + tzy*ny + tzz*nz)
                    Fn = Fx*nx + Fy*ny + Fz*nz
                    Fwallx = Fx - Fn*nx
                    Fwally = Fy - Fn*ny
                    Fwallz = Fz - Fn*nz
                    cf_x   = Fwallx/self.DynamicVelocity
                    cf_y   = Fwally/self.DynamicVelocity
                    cf_z   = Fwallz/self.DynamicVelocity
                    Fwall = np.sqrt(np.power(Fwallx,2) + np.power(Fwally,2) + np.power(Fwallz,2))
                    cp = PDiff/self.DynamicVelocity
                    cp_x = -cp*nx*Area
                    cp_y = -cp*ny*Area
                    cp_z = -cp*nz*Area
                    cd.append(np.sum((Fx*self.nDrag[0] + Fy*self.nDrag[1] + Fz*self.nDrag[2])*Area/self.DynamicVelocity))
                    cd.append(np.sum(cp_x*self.nDrag[0] + cp_y*self.nDrag[1] + cp_y*self.nDrag[2]))
                    cl.append(np.sum((Fx*self.nLift[0] + Fy*self.nLift[1] + Fz*self.nLift[2])*Area/self.DynamicVelocity))
                    cl.append(np.sum(cp_x*self.nLift[0] + cp_y*self.nLift[1] + cp_z*self.nLift[2]))
                    #print Cf_x.shape, Fwall.shape, nDrag[0].shape
                    Rex = self.RhoRef*self.Velocity*x/self.MuRef
                    #np.savetxt(f, np.c_[x, Rex, cp, cf_x, cf_y], delimiter="  ")
                    np.savetxt(f, np.c_[x, Rex, cf_x], delimiter="  ")
        CD = sum(cd)/0.05
        Difference = np.abs((0.00133-CD)*100/0.00133)
        print " ------ Laminar Test case: Flat plate ------ "
        print " Flux Scheme        : "+ Input.SchemeDict['FluxScheme']
        print " Higher order method: "+ Input.SchemeDict['FaceScheme']
        print " Turbulence model   : "+ Input.SchemeDict['TurbulenceModel']
        print " Expected drag coeffcient    : "+ "{:.3E}".format(0.00133)
        print " Calculated drag coefficient : "+ "{:.3E}".format(CD)
        print " Difference                  : "+ "{:.3E}".format(Difference) + " %"
        print " Allowed Tolerance           : 1 %"
        if Difference < 1:
            print "------------ >>> Test Passed  <<< --------------"
        else:
            print "------------ xxx Test Failed  xxx --------------"
