#!/bin/bash

#echo "tentative d'ecriture de ndiff_results.log dans :" $PWD
LOG_FILE="ndiff_results.log"
if [ -f $LOG_FILE ]  ; then     rm $LOG_FILE ; fi
ERR=0
bashver=$( bash --version | sed -Ee 's/GNU bash, version ([0-9.]+).*/\1/;q' )
echo $bashver

version="NONE"
ndiff=`(which ndiff)`>> $LOG_FILE
if [  "$ndiff" != "" ]; then
    echo  "1.">> $LOG_FILE
    version=`$ndiff -v 2>&1 ndiff_version.log ` 
   # vers=$(cat ndiff_version.log)
        echo "<"$version">"
else
    echo  "2.">> $LOG_FILE
    echo '---------------------------------------------------------------------------------' >> $LOG_FILE
    echo '|                      ERROR: ndiff is not installed.                           |' >> $LOG_FILE
    echo '|        The differences between execution results and reference results        |' >> $LOG_FILE
    echo '|                      can not be checked correctly                             |' >> $LOG_FILE
    echo '---------------------------------------------------------------------------------' >> $LOG_FILE
    exit 1
fi
echo $version >> $LOG_FILE
SUBSTRING=`(echo $version  | grep -oF 'extended-precision arithmetic' )`
echo "<"${SUBSTRING}">"  >> $LOG_FILE
reference="extended-precision arithmetic"
echo "<"${reference}">"  >> $LOG_FILE
if [ "${SUBSTRING}" = "${reference}" ] ;  then
    echo  "3.">> $LOG_FILE
else
    echo  "4.">> $LOG_FILE
    echo '---------------------------------------------------------------------------------' >> $LOG_FILE
    echo '|  WARNING: ndiff is not the `extended-precision arithmetic` diff tool needed.  |' >> $LOG_FILE
    echo '|        The differences between execution results and reference results        |' >> $LOG_FILE
    echo '|                      may not be checked correctly.                            |' >> $LOG_FILE
    echo '|         see:  https://www.math.utah.edu/~beebe/software/ndiff/	          |' >> $LOG_FILE
    echo '---------------------------------------------------------------------------------' >> $LOG_FILE
    exit 1
fi
echo "OK" >> $LOG_FILE


TESTCASE=`pwd`
if [ "$DEVROOT" != "" ] ;then
    echo ----------- treating $TESTCASE ------------    
    echo ----------- treating $TESTCASE ------------   >> $DEVROOT/$LOG_FILE
fi
q=0

for RESULT_FILE in `ls reference/`; do

    if [ -d reference/$RESULT_FILE ] ;     then  
	echo :::::::::: $RESULT_FILE is a directory not a reference file :::::::::: >> $LOG_FILE
    else

	echo :::::::::: treating $RESULT_FILE ::::::::::   >> $LOG_FILE
	if [ "$DEVROOT" != "" ] ;then
	    echo :::::::::: treating $RESULT_FILE ::::::::::   >> $DEVROOT/$LOG_FILE
	fi

	if [ -f $RESULT_FILE ] ;     then
	    if [ "$ndiff" != "" ] ;then
		$ndiff  -abserr 1.e-11 -separators '[ \t\(\),{}]'  -outfile temp.log $RESULT_FILE reference/$RESULT_FILE 
	    else 
		diff -q $RESULT_FILE reference/$RESULT_FILE > temp.log
	    fi

	    if [ $? = 0 ]
	    then

		q=0
	    else
		p=$(grep "absolute error" temp.log | wc -l)
		r=$(grep ">" temp.log | wc -l)
		s=$(grep "<" temp.log | wc -l)
		q=$(( ${m} + ${p}  + ${r}  + ${s} ))
	    fi
	    if [ $q -ne 0 ] 
	    then
		echo ERRORS : $q differences / errors messages : >> $LOG_FILE
		cat temp.log >> $LOG_FILE
		if [ "$DEVROOT" != "" ] ;then
		    cat temp.log >> $DEVROOT/$LOG_FILE
		fi
	    fi   
	else
	    echo    WARNING : ""$RESULT_FILE"" exists as reference but not as results ! >> $LOG_FILE
	    ERR=$(( ${ERR} + 1 ))
	fi
	ERR=$(( ${ERR} + ${q} ))
    fi
done

if [ $ERR -ne 0 ] 
then
    #   return 1 to signal a error to CTest
    exit 1
fi   

