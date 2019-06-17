title: Mesh

<figure>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/Domain.png" alt="Domain" style="width:350px">
  </div>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <figcaption> Fig.1 Domain of interest divide into 4 sub-domain.</figcaption>
  </div>
</figure>


To solve a real-life problem, like flow over
an airfoil, first, a domain of computation has to be defined. This 
domain represents the system of interest which we
like to study. Now to apply a numerical method
to solve the governing equation (Navier-Stokes) on this domain, it needed
to be discretized. This process is not part of the FEST-3D
code and it assume a predefine discretized domain(grid) will
be provided as input. You can use softwares like Pointwise to generate the
good quality mesh. The FEST-3D expect structured grid in
 the following format
```
Imax  Jmax  Kmax
X(1,1,1)   Y(1,1,1)   Z(1,1,1)
X(2,1,1)   Y(2,1,1)   Z(2,1,1)
X(3,1,1)   Y(3,1,1)   Z(3,1,1)
...
X(Imax,1,1)   Y(Imax,1,1)   Z(Imax,1,1)
X(1,2,1)   Y(1,2,1)   Z(1,2,1)
X(2,2,1)   Y(2,2,1)   Z(2,2,1)
...
...
X(Imax,Jmax,1)   Y(Imax,Jmax,1)   Z(Imax,Jmax,1)
X(1,1,2)   Y(1,1,2)   Z(1,1,2)
X(2,1,2)   Y(2,1,2)   Z(2,1,2)
X(3,1,2)   Y(3,1,2)   Z(3,1,2)
...
...
...
X(Imax,Jmax,kmax)   Y(Imax,Jmax,kmax)   Z(Imax,Jmax,kmax)
```

In the case of multi-block, as shown in Fig. 1, each process requires a separate grid file
and block should be created such that every face of that
block employs only a single boundary condition. For now, two different types of boundary
conditions at a single face  of a block is not supported by the FEST-3D.
<figure>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/Mesh.png" alt="Mesh" style="width:350px">
  </div>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <figcaption> Fig.2 Example Mesh.</figcaption>
  </div>
</figure>

A very simple example of grid file is given below for 2D-smooth bump test case, shown in Fig. 3:<br>
`Blocks: 2` Number of blocks in I-direction are 2. So, two files will be written: `grid_00.txt` for first block and 
`grid__01.txt` for second block<br>
`Imax: 3` <br>
`Jmin: 2` <br>
`Kmin: 2` <br>
<figure>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/CoarseGrid.png" alt="CoarseGrid" style="width:350px">
  </div>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <figcaption> Fig.2 Simple coarse grid for 2D smooth bump test case.</figcaption>
  </div>
</figure>

### grid_00.txt
```
2 2 2
-1.5000e+00  2.3271e-26  0.0000e+00
 0.0000e+00  6.2500e-02  0.0000e+00
-1.5000e+00  8.0000e-01  0.0000e+00
 0.0000e+00  8.0000e-01  0.0000e+00
-1.5000e+00  2.3271e-26  1.0000e-01
 0.0000e+00  6.2500e-02  1.0000e-01
-1.5000e+00  8.0000e-01  1.0000e-01
 0.0000e+00  8.0000e-01  1.0000e-01
```

### grid_01.txt
```
2 2 2
0.0000e+00  6.2500e-02  0.0000e+00
1.5000e+00  2.3271e-26  0.0000e+00
0.0000e+00  8.0000e-01  0.0000e+00
1.5000e+00  8.0000e-01  0.0000e+00
0.0000e+00  6.2500e-02  1.0000e-01
1.5000e+00  2.3271e-26  1.0000e-01
0.0000e+00  8.0000e-01  1.0000e-01
1.5000e+00  8.0000e-01  1.0000e-01
```
