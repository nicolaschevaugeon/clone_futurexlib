#!/bin/bash
# check number of argument and silentely exit if incorect number
if [ $# -ne 2 ]
then
    exit 1;
fi
# first is the source to link from
SOURCE=$1
# second is the target to link to
TARGET=$2
# does the source exist
if [ -e $SOURCE ] 
then
    #does the link not allready done
    if  [ ! -h $TARGET  ]
    then
         ln -s $SOURCE $TARGET;
     fi
fi
exit 0;
