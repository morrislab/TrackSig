#!/bin/bash

sample_dir="/home/q/qmorris/yulia/prostate_cancer_200_fullresults/samples/"
phi_dir="/home/q/qmorris/yulia/prostate_cancer_200_fullresults/ssmphis/"
phi_script="/home/q/qmorris/yulia/src/exultant-pistachio/protocols/make-calibration-report/calc_ssm_phis.py"

for sample in $sample_dir/* ; do
	if [[ ! -d $sample ]]; then
		continue
	fi
	
	samplename=$(basename $sample)
	
	vcf=$sample/fh_$samplename.vcf 
	mutation_json=$sample/$samplename.muts.json
	# assignment_json=`ls $sample/[0-9]*.json`
	assignment_json=$sample/$samplename.mutass.zip
	tree_json=$sample/$samplename.summ.json

	gzip < $mutation_json > $mutation_json.gz

	gzip < $tree_json > $tree_json.gz

	#cd $sample/
	#zip $(basename $assignment_json).zip $(basename $assignment_json)
	#cd ..

	mutation_json=$mutation_json.gz
        # assignment_json=$assignment_json.zip
        tree_json=$tree_json.gz
	
	phi_file=$phi_dir/$samplename.txt

	echo $sample
	python $phi_script $tree_json $mutation_json $assignment_json > $phi_file
	
done

