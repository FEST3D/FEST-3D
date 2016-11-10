#!bin/bash 
WORKDIR=$(pwd)
SOURCE=$WORKDIR"/archive/setup/"
TARGET=$1
echo '---> Setting up' $WORKDIR"/"$1
mkdir $WORKDIR"/"$1
echo "=============================================================="
echo "--->copying content " $SOURCE " >>>------>>> " $TARGET
cp -r $SOURCE* $TARGET
echo "--->creating a softlink to executable"
ln -fs $WORKDIR/bin/FEST3D $TARGET
cd $TARGET
#mv sample.config.md config.md
echo "--->Changing work directory"
echo "============= Right now in " $TARGET " Directory ==============" 
