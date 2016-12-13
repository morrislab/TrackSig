#!/bin/bash

# Usage ./runRscript my_rscript.R [--extra-args "argument_for_r_script"; --job-id ...]
# Note that argument for r script should all be in quotes!
# R script names should be given with respect to $current_folder/src/ directory. Output or input directory names should be given with respect to $current_folder directory

rscript=$1

current_folder=/dupa-filer/yulia/mutational_signatures

if [ -z "$rscript" ]; then
    echo "Please provide the name of r script to run"
    exit 1
fi

if [[ $rscript != *.R ]]; then
        echo "The specified file $rscript is not an R script"
fi

extra_args=""
if [[ $2 = "--extra-args" ]]; then
        if [[ -z $3 ]]; then
                echo "Please provide extra argumets for R script"
                exit
        else
                extra_args=$3
        fi
fi

extra_job_id=""
if [[ $4 = "--job-id" ]]; then
        if [[ -z $5 ]]; then
                echo "Please provide extra argumets for R script"
                exit
        else
                extra_job_id=$5
        fi
fi


dir=${current_folder}/src/tmp
logdir=${current_folder}/log
if [[ !  -e "$dir" ]]; then
	mkdir $dir
fi

if [[ ! -e "$logdir" ]]; then
	mkdir $logdir
fi

for sc in ${rscript[@]}; do
	name=$(basename ${sc})
	name=${name%.R}
	if [[ ${extra_job_id} != "" ]];then
		name=${name}.${extra_job_id}
	fi
	script=${dir}/${name}.sh
 
	output_file=$logdir/${name}
	error_file=$logdir/${name}

	if [[ -e ${output_file} ]]; then
		rm ${output_file}
	fi
	if [[ -e ${error_file} ]]; then
		rm ${error_file}
	fi

#$ -l excl=true
	cat > $script <<EOF
#!/usr/bin/env bash

#$ -V
#$ -S /bin/bash
#$ -N "${name}"
#$ -e ${error_file}"
#$ -o ${output_file}"
#$ -l hostname="supa*"
#$ -l h_rt=48:00:00
#$ -l s_rt=40:00:11

set -eu

Rscript ${current_folder}/${sc} ${extra_args}

EOF
        qsub $script
done

#rm -r ${dir}
