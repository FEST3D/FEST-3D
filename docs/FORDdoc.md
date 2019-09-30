---
project: FEST-3D
project_dir: ./
src_dir: ../src/
project_github: https://github.com/FEST3D/FEST-3D
project_download: https://github.com/FEST3D/FEST-3D/archive/master.zip
summary: <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/FEST3D.png" alt="FEST-3D" style="width:350px"></div>
media_dir: ./media
output_dir: ./html
page_dir: ./user_guide
tutorial_dir: ./Testcases
author: FEST-3D Team
author_description: Group of students from the Department of Aerospace Engineering at Indian Institute of Technology Madras (IITM), Chennai (600036), India.
github: https://github.com/FEST3D
email: fest3diitm[@]gmail.com
print_creation_date: true
year: 2019
parallel: 0
coloured_edges: true
graph_maxnodes: 50
fpp_extensions: fpp
docmark: <
display: public
         protected
         private
source: true
graph: true
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: by-nc
extra_filetypes: sh #
exclude: scheme.f90
exclude: face_interpolant.f90
md_extensions: markdown.extensions.toc
---

FEST-3D (Finite-volume Explicit STructured 3-Dimensional) is computational fluid dynamic solver written in Fortran 90 for solving the Navier-Stokes equations on structured grids using state-of-the-art finite-volume  methods. It is a modular, multiblock, finite-volume code developed to solve flow problems  in the field of aerodynamics.

## Installation
For installation instructions check out the [Documentation](./page/01_install.html) guide.

## Steps to run FEST-3D 
For detailed documentation on steps to run the FEST-3D after downlaod and installation, check out "[How to run FEST-3D](./page/04_Steps_to_run_FEST3D.html)" page.

## Tutorials
For tutorials check out the [test_case](./page/05_tutorials/index.html) descriptions.

@note
There are two separate GitHub repositories:<br>
1. <a href="https://github.com/FEST3D/FEST-3D" target="_blank">FEST-3D<a> directory, which contains the source code<br>
2. <a href="https://github.com/FEST3D/run" target="_blank">Run</a> directory, which contains the tutorials<br>
The run directory is a submodule of the FEST-3D code.
@endnote

## Modules and Subroutines
Although most of modules and subroutines are named such that its purpose is clear from its name, still more information is provided in [source File](./lists/files.html) and [Modules](./lists/modules.html) pages. Documentation of FEST-3D uses [Ford documentation generator](https://github.com/Fortran-FOSS-Programmers/ford) which generates documentation using the special comments written for each subroutine and module.

## License
FEST-3D is an open-source software available under GNU General Public License v3.0

## Reference and Citation
Details about the governing equations used in the FEST-3D code can be found in the publication by Jatinder Pal Singh Sandhu et al. (_Singh Sandhu, J. P., Girdhar, A., Ramakrishnan, R., Teja, R. D., & Ghosh, S., **A convergence study of solutions using two two-equation RANS turbulence models on a finite volume solver for structured grids**, AIAA 2018-3859_).

User are requested to cite this reference for any research publications made using the FEST-3D solver.


@Bug

  - Halt/Stop condition does not work on MPI
  - Can not run turbulent flow simulation on Windows Subsystem for Linux (WSL) due to outdated OS and libraries used by WSL.

@endbug

@todo
Add HLLC and Roe flux difference Schemes
