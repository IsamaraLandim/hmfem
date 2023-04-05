#!/bin/bash

nx=${1}
ny=${2}
nz=${3}

## size of domain 
Ix=${4} 
Iy=${5}
Iz=${6}

PERM=${7}


karst=${8}

type=${9}

WIFILE="./WI_pocoradial.wi"

##  To run name 
EXE="./run/exe"

echo "$layer"

PERMFILE="./permfields/spe_data.perm"

DIR="./AMG/AMG_Perm"$PERM"_"$nx"x"$ny"x"$nz"_"$Ix"x"$Iy"x"$Iz"_"$karst"_"$type""
mkdir -p $DIR

## Save terminal output
STDOUT=$DIR"/output.out"
STDERR=$DIR"/err.out"

## Run the programm
start=`date +%s`

# run the programm
$EXE $nx $ny $nz $Ix $Iy $Iz $PERM $karst $type $DIR $PERMFILE $WIFILE  2>$STDERR > $STDOUT 


end=`date +%s`
runtime=$((end-start))
echo "Terminating after "$runtime" seconds." 
eval "echo " $runtime "> "$DIR"/time.txt"


## plot flux
##IMG3=$DIR"/"$TNAME"_ufield3D.png"
##gnuplot -e "ufile='$DIR/flux.dat'; img1='$IMG3'; Lx='$Ix'; Ly='$Iy'; Lz='$Iz'; nsx='$nsx'; nsy='$nsy'; nsz='$nsz'; nx='$nx'; ny='$ny'; nz='$nz'" gnuplot/plotvector.plt 
##IMG3=$DIR"/"$TNAME"_ufieldPlane.png"
##gnuplot -e "ufile='$DIR/flux.dat'; img1='$IMG3'; Lx='$Ix'; Ly='$Iy'; Lz='$Iz'; nsx='1'; nsy='1'; nsz='1'; nx='$nx'; ny='$ny'; nz='$nz'" gnuplot/plotvector2D.plt

