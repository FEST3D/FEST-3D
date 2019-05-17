---
project: FEST-3D
project_dir: ./
src_dir: ../src/
project_github: https://github.com/FEST3D/FEST-3D
project_download: https://fest3d.github.io/
media_dir: ./media
output_dir: ./html
page_dir: ./user_guide
tutorial_dir: ./Testcases
summary: ![FEST3D](|media|/FEST3D.png)
         {: style="text-align: center" }
author: FEST3D Team
author_description: Team consists of multiple students from the department of Aerospace engineering at Indian Institute of Technology Madras (IITM), Chennai (600036), India.
github: https://github.com/jayten
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

FEST3D (Finite-volume Explicit STructured 3-Dimensional) is computational fluid dynamic code written in Fortran 90 for solving Navier-Stokes equations on a structured grid using state of the art finite-volume numerical methods. It is a modular, multiblock finite-volume code developed to solve compressible flow problems encountered in the field of aerodynamics.

## Installation
For installation instruction check out [Documentation](./page/01_install.html) guide. 

## Tutorials
For tutorials check out [test_case](./page/05_tutorials/index.html) descriptions.

@Bug
 - Halt/Stop condition does not work on MPI
@endbug

@todo
Add HLLC and Roe flux difference Schemes



## FEST-3D Team:
### Developers:
    - Jatinder Pal Singh Sandhu
    - Rakesh Ramakrishnan
    - Anant Girdhar
    - R. D. Teja

### Advisor:
    - Dr. Santanu Ghosh
