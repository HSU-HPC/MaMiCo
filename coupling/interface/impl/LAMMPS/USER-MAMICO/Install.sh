# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test ! -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# all package files
# only a few files have dependencies

for file in *.cpp *.h; do
  test -f ${file} && action $file
done

# add include directory for MAMICO to PKG_INC; currently, this path
# needs to exist in form of env.variable MAMICO_PATH, pointing to
# the directory that contains the folders coupling and tarch.
if (test $mode = 1) then
  
  if (test -e ../Makefile.package) then
    sed -i 's|-I${MAMICO_PATH} -I${LIB_EIGEN_PATH} -std=c++1z -DLAMMPS_MD -DMDDim3 -DMDCoupledParallel -DTarchParallel||g' ../Makefile.package
    sed -i 's|^PKG_INC =*|& -I${MAMICO_PATH} -I${LIB_EIGEN_PATH} -std=c++1z -DLAMMPS_MD -DMDDim3 -DMDCoupledParallel -DTarchParallel|' ../Makefile.package
    #sed -i -e 's|^PKG_INC =[ \t]*|&-I${MAMICO_PATH} -std=c++11 -DMDCoupledParallel |' ../Makefile.package
  fi
elif (test $mode = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's|[^ \t]*-I${MAMICO_PATH} -I${LIB_EIGEN_PATH} -std=c++1z -DLAMMPS_MD -DMDDim3 -DMDCoupledParallel -DTarchParallel[^ \t]* ||' ../Makefile.package
  fi

fi

