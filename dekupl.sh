#!/bin/bash

printHelp() {
	echo "\n\t\t\tUser Commands\n
NAME
    DEKUPL   detect variant
OPTIONS
-test
    get data to test the result
-getconf
    get the configuration file to run dekulp
-nbthread
    fix the thread numbers when running, by default value=1
-run
    run dekupl
SYNOPSIS
    docker run --rm --privileged -v $(pwd):/data -w /data ebio/dekupl -nbthread NB -run"
}

while [ "$1" != "" ]; do
	case $1 in
		-test )				shift
								get_data_test=true
								;;
	    -getconf )			shift
 								get_conf_file=true
 								;;
		-run )				shift
								run_dekupl=true
								;;
		-nbthread )			shift
								nbthread=$1
								;;
		-h )				printHelp
								exit 0
								;;
		* )				echo "Use -h for help: " $1" unvalid"
						exit 1
	esac
 	shift
done

echo $nbthread

if [ "$get_data_test" = true ]; then
	echo "Get test data files"
	mkdir /data/data
	cp /bin/dekupl/data/* /data/data/
	chmod 666 /data/data/*
	exit 1
elif [ "$get_conf_file" = true ]; then
	echo "Get config file"
	cp /bin/dekupl/config.json /data
	chmod 666 /data/config.json
	exit 1
elif [ "$run_dekupl" = true ]; then
	echo "Run dekupl - HOME"
	FILE="/data/config.json"     
    if [ -f $FILE ]; then
       echo "File $FILE exists"
       cd /bin/dekupl
       if [ "$nbthread" = "" ]; then
           echo "run dekupl with 1 thread"
           wrapdocker docker pull itsjeffreyy/samtools
           docker pull ebio/jellyfish
           docker pull ebio/kallisto
           docker pull ebio/gsnap
           snakemake --timestamp --stats /data/log_file -s /data/Snakefile -j 1
        else
           echo "Run dekupl with $nbthread thread(s)"
           wrapdocker docker pull itsjeffreyy/samtools
           docker pull ebio/jellyfish
           docker pull ebio/kallisto
           docker pull ebio/gsnap
           snakemake --timestamp --stats /data/log_file -s /data/Snakefile -j $nbthread    
        fi
    else
       echo "File $FILE does not exist."
    fi
	exit 1
fi
