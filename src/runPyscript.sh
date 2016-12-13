#!/bin/bash

# Usage ./runRscript my_rscript.R [--extra-args "argument_for_r_script"; --job-id ...]
# Note that argument for r script should all be in quotes!
# R script names should be given with respect to $current_folder/src/ directory. Output or input directory names should be given with respect to $current_folder directory

rscript=$1

current_folder=/home/q/qmorris/yulia
write_dir=/scratch/q/qmorris/yulia

if [ -z "$rscript" ]; then
    echo "Please provide the name of python script to run"
    exit 1
fi

if [[ $rscript != *.py ]]; then
        echo "The specified file $rscript is not an python script"
fi

extra_args=""
if [[ $2 = "--extra-args" ]]; then
        if [[ -z $3 ]]; then
                echo "Please provide extra argumets for python script"
                exit
        else
                extra_args=$3
        fi
fi

extra_job_id=""
if [[ $4 = "--job-id" ]]; then
        if [[ -z $5 ]]; then
                echo "Please provide extra argumets for python script"
                exit
        else
                extra_job_id=$5
        fi
fi

gravity=false
if [[ $6 = "--gravity" ]]; then
	gravity = true
fi

dir=${current_folder}/src/tmp
logdir=${write_dir}/log
if [[ !  -e "$dir" ]]; then
	mkdir $dir
fi

if [[ ! -e "$logdir" ]]; then
	mkdir $logdir
fi

for sc in ${rscript[@]}; do
	name=$(basename ${sc})
	name=${name%.py}
	if [[ ${extra_job_id} != "" ]];then
		name=${name}.${extra_job_id}
	fi
	script=${dir}/${name}.sh
 
	output_file=$logdir/${name}.o
	error_file=$logdir/${name}.e

	echo $name
	if [[ -e ${output_file} ]]; then
		rm ${output_file}
	fi
	if [[ -e ${error_file} ]]; then
		rm ${error_file}
	fi

	queue=""
	if $gravity; then
		queue="#PBS -q gravity"
	fi
	settings="nodes=1:ppn=8,walltime=47:00:00"
	if $gravity; then
		settings="nodes=1:ppn=12:gpus=2,walltime=12:00:00"
	fi

#$ -l excl=true
#PBS -q gravity
#nodes=1:ppn=12:gpus=2,walltime=12:00:00

	cat > $script <<EOF
#!/bin/bash

#PBS -V
#PBS -S /bin/bash
#PBS -N ${name}
#PBS -e ${error_file}
#PBS -o ${output_file}
#PBS -l ${settings}
#PBS -m abe
$queue

set -eu

module load gcc/4.9.0 intel/15.0.2 gsl/1.16-intel python/2.7.8
#pip2 install --user --upgrade https://storage.googleapis.com/tensorflow/linux/cpu/tensorflow-0.7.1-cp27-none-linux_x86_64.whl

python ${current_folder}/${sc} ${extra_args} 

EOF
        qsub $script
done

#rm -r ${dir}
