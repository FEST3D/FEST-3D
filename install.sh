#!/usr/bin/env bash

#unset -xo pipefail

OS="`uname`"
dir="`pwd`"/bin/FEST3D
echo $dir
case ${OS} in
  'Linux')
    OS='Linux'
    PKG_MANAGER=$( command -v yum || command -v apt-get ) || echo "Neither yum nor apt-get found"
    echo "Following package manager found: ${PKG_MANAGER}"
    echo "Meeting dependencies"
    apt=`command -v apt-get`
    yum=`command -v yum`
    if [ -n "$apt" ]; then
      sudo apt-get update
      sudo apt-get -y install gfortran mpich python3-dev python3-numpy
    elif [ -n "$yum" ]; then
      sudo yum check-update
      sudo yum -y install gcc-gfortran mpich python3-devel
      sudo python3 -m pip install --upgrade pip
      python3 -m pip install --user numpy
      pip install --user numpy --upgrade
    else
      echo "Err: no path to apt-get or yum" >&2;
      exit 1;
    fi
    source ~/.bashrc
    if [[ -z "${FEST3D}" ]]; then
      echo 'FEST3D variable not set. Exporting FEST3D variable to .bashrc'
      echo 'export FEST3D="'${dir}'"' >> ~/.bashrc
      export FEST3D="${dir}"
    elif [[ "${FEST3D}" == "${dir}" ]]; then
      echo 'FEST3D env variable is set to: '${FEST3D}
    else
      echo 'export FEST3D="'${dir}'"' >> ~/.bashrc
      export FEST3D="${dir}"
      echo 'FEST3D variable is set, but to different binary path'
      echo 'FEST3D env variable is being set to: '${FEST3D}
      echo 'Please remove the old export line from ~/.bashrc'
    fi
    ;;
  'FreeBSD')
    OS='FreeBSD'
    echo "Install scirpt does not suport FreeBSD"
    ;;
  'WindowsNT')
    OS='Windows'
    echo "Install scirpt does not suport Windows"
    ;;
  'Darwin') 
    OS='Mac'
    #gfortran is part of gcc
    brew install gcc mpich numpy
    if [[ -z "${FEST3D}" ]]; then
      echo 'FEST3D variable not set. Exporting FEST3D variable to .bash_profile'
      echo 'export FEST3D="'${dir}'"' >> ~/.bash_profile
      export FEST3D="${dir}"
    elif [[ "${FEST3D}" == "${dir}" ]]; then
      echo 'FEST3D env variable is set to: '${FEST3D}
    else
      echo 'export FEST3D="'${dir}'"' >> ~/.bash_profile
      export FEST3D="${dir}"
      echo 'FEST3D variable is set, but to different binary path'
      echo 'FEST3D env variable is being set to: '${FEST3D}
      echo 'Please remove the old export line from ~/.bashrc'
    fi
    ;;
  'SunOS')
    OS='Solaris'
    ;;
  'AIX') ;;
  *) ;;
esac

# create a build folder for the out-of-source build
mkdir -p build
# switch to build directory
cd build
# run cmake; here we assume that the project's
# top-level CMakeLists.txt is located at '..'
cmake ..
# once CMake has done its job we just build using make as usual
 make
# if the project uses ctest we can run the tests like this
 make test
 cd ../tests/
 make
 python3 Test.py ausm muscl sst
 cd ../
