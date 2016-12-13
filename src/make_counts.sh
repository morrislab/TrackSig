#!/bin/bash

vcf_path=$1 # Path to directory with vcf files
phi_path=$2 # Path to directory with phi.txt files produces by cals_ssm_phis.py
mutation_counts_path=$3 # Path where to write the mutation counts
mut_order_path=$4 # Path where to write mutation ordering (list of mutations sorted by phi). Needed to run get_clusters_at_timepoints.R. get_clusters_at_timepoints.R   composes list of tree node assignments for chunks of 100 mutations (prevalent tree node assignments at this time point)
mutation_types_path=$5 # Path where to write files listing mutation type (out of 96 trinucleotide-based types) for each mutation in vcf file sorted by phi
mutation_bootstrap_path=$6 # Path where to write bootstrapped mutation counts
group=$7 # Optional: the starting letter of the vcf files to be processed, for parallelization

do_bootstrap=false

N_BOOTSTRAPS=20

log_dir=/scratch/q/qmorris/yulia/

min_number() {
    printf "%s\n" "$@" | sort -g | head -n1
}


if [[ -z $vcf_path ]]; then
	echo "Please provide VCF path"
	exit
fi

if [[ -z $phi_path ]]; then
	echo "Please provide phi path"
	exit
fi

if [[ -z $mutation_counts_path ]]; then
	mutation_counts_path=.
fi

if [[ -z $mut_order_path ]]; then
	mut_order_path=.
fi

if [[ -z $mutation_types_path ]]; then
	mutation_types_path=.
fi

if [[ -z $mutation_bootstrap_path ]]; then
	mutation_bootstrap_path=$mutation_counts_path
fi

if [[ -z $ ]]; then
	group=""
fi  

for file in $vcf_path/$group*; do
    if [[ ! -f $file || $file != *.vcf ]]; then
      continue 
    fi

    echo $file 

    tumor_id=$(basename ${file})
    tumor_id=${tumor_id%.annotated.snv_mnv.vcf}
    tumor_id=${tumor_id%_filtered.vcf}

    tumor_part=""
    if [[ $tumor_id =~ ^.*tumor.*.vcf$ ]]; then
	tumor_part=${tumor_id:${#tumor_id}-5:1}
	tumor_part=tumor${tumor_part}.
    fi

    tumor_id=${tumor_id%.tumor*.vcf}
    tumor_id=${tumor_id%.vcf}
    tumor_id=${tumor_id#fh_}
	
    phi_file=$phi_path/$tumor_id.txt
    
    if [[ ! -f "$phi_file" ]]; then
	echo "$phi_file not found"
    	continue
    fi
   
   mutation_counts_file=$mutation_counts_path/$tumor_id.${tumor_part}phi.txt
   mutation_types_file=$mutation_types_path/$tumor_id.${tumor_part}mut_types.txt	
   mut_order_file=$mut_order_path/$tumor_id.${tumor_part}mut_order.txt

   echo $mutation_counts_file
   echo $mutation_types_file
   echo $mut_order_file

   if [ ! -f $mutation_types_file ] || [ ! -s  $mutation_types_file ]; then
	echo "Type file..."
	perl /home/q/qmorris/yulia/src/getMutationTypes.pl muse $file $phi_file  $mut_order_file >> $mutation_types_file 2>> $log_dir/log.txt
   fi

   num_mutations=$(min_number `cat $phi_file | wc -l` `cat $file | wc -l`)
   echo $num_mutations

   num_hundreds=$(($num_mutations/100 + ($num_mutations % 100 > 0))) 

   if [ ! -f $mutation_counts_file ]; then
	if [ $num_mutations -lt 100 ]; then
		echo "Less than 100 mutaions in a file $file or $phi_file" 
		touch $mutation_counts_file
	else
		echo "Count file..."
   		for i in `seq 1 $num_hundreds`; do 
  			python /home/q/qmorris/yulia/src/make_hundreds.py $mutation_types_file  $((i*100-100)) $((i*100-1)) >> $mutation_counts_file 2>>$log_dir/log.txt
  	 	done
	fi 
   fi

if [ "$do_bootstrap" = true ] ; then
   bootstrap_dir=$mutation_bootstrap_path/$tumor_id/
   if [ !  -z  $tumor_part  ]; then
	bootstrap_dir=${bootstrap_dir:0:${#bootstrap_dir}-1}
	bootstrap_dir=$bootstrap_dir.${tumor_part:0:${#tumor_part}-1}/
   fi

   echo $tumor_part
   echo $bootstrap_dir

   if [ ! -d $bootstrap_dir ]; then
	mkdir $bootstrap_dir
   fi
   
   for i in `seq 1 $N_BOOTSTRAPS`; do
	mutation_bootstrap_file=$bootstrap_dir/$tumor_id.${tumor_part}$i.mut.txt
	mutation_bootstrap_counts_file=$bootstrap_dir/$tumor_id.${tumor_part}$i.mut_counts.txt
	
	mutation_bootstrap_file_unsorted=$bootstrap_dir/$tumor_id.${tumor_part}$i.mut.unsorted.txt
        mutation_bootstrap_counts_file_unsorted=$bootstrap_dir/$tumor_id.${tumor_part}$i.mut_counts.unsorted.txt

	if [ ! -f $mutation_bootstrap_file_unsorted ] || [ ! -s  $mutation_bootstrap_file_unsorted ]; then

		if [ $num_mutations -lt 100 ]; then
                	touch $mutation_bootstrap_file
        	else
                	echo "Bootstrap..."
			python /home/q/qmorris/yulia/src/bootstrap_mutations.py $mutation_types_file  $num_mutations > $mutation_bootstrap_file_unsorted 2>>$log_dir/log.txt
		
			#sort -k3 -nr $mutation_bootstrap_file_unsorted > $mutation_bootstrap_file
		fi
	fi
	if [ ! -f $mutation_bootstrap_counts_file_unsorted ] || [ ! -s  $mutation_bootstrap_counts_file_unsorted ]; then
			for t in `seq 1 $num_hundreds`; do
				#python /home/q/qmorris/yulia/src/make_hundreds.py $mutation_bootstrap_file  $((t*100-100)) $((t*100-1)) >> $mutation_bootstrap_counts_file 2>>$log_dir/log.txt
				python /home/q/qmorris/yulia/src/make_hundreds.py $mutation_bootstrap_file_unsorted  $((t*100-100)) $((t*100-1)) >> $mutation_bootstrap_counts_file_unsorted 2>>$log_dir/log.txt
  			done
	fi
  
  done  
fi

done

#perl extractFeaturesPhi.pl muse pwgs/consprelim/0554ffe5-31f7-43f5-8372-2b73c9cf3527.annotated.snv_mnv.vcf pwgs/ssmphi/consprelim.sample-5000.cnvint/0554ffe5-31f7-43f5-8372-2b73c9cf3527.txt 0 99


