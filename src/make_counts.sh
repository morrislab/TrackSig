#!/bin/bash

# AUTHOR: Yulia Rubanova

vcf_file=$1 # Path to directory with vcf files
phi_file=$2 # Path to directory with phi.txt files produces by cals_ssm_phis.py

mutation_counts_path="data/counts/" # Path where to write the mutation counts
mut_order_path="data/mut_order/" # Path where to write mutation ordering (list of mutations sorted by phi). Needed to run get_clusters_at_timepoints.R. get_clusters_at_timepoints.R   composes list of tree node assignments for chunks of 100 mutations (prevalent tree node assignments at this time point)
mutation_types_path="data/mut_types/" # Path where to write files listing mutation type (out of 96 trinucleotide-based types) for each mutation in vcf file sorted by phi
mutation_bootstrap_path="data/bootstrap/" # Path where to write bootstrapped mutation counts

make_hundreds_script="src/make_hundreds.py"
get_mutation_type_script="src/getMutationTypes.pl"
bootstrap_mutations_script="src/bootstrap_mutations.py"

do_bootstrap=true

N_BOOTSTRAPS=30

log_dir="."

min_number() {
	printf "%s\n" "$@" | sort -g | head -n1
}


if [[ -z $vcf_file ]]; then
	echo "Please provide VCF path"
	exit
fi

if [[ -z $phi_file ]]; then
	echo "Please provide phi path"
	exit
fi

if [[ -z $mutation_counts_path ]]; then
	mutation_counts_path=.
else
   if [ ! -d "$mutation_counts_path" ]; then
	mkdir -p $mutation_counts_path
   fi
fi

if [[ -z $mut_order_path ]]; then
	mut_order_path=.
else
  if [ ! -d "$mut_order_path" ]; then
	mkdir -p $mut_order_path
  fi
fi

if [[ -z $mutation_types_path ]]; then
	mutation_types_path=.
else
   if [ ! -d "$mutation_types_path" ]; then
	mkdir -p $mutation_types_path
  fi
fi

if [[ -z $mutation_bootstrap_path ]]; then
	mutation_bootstrap_path=$mutation_counts_path
else
   if [ ! -d "$mutation_bootstrap_path" ]; then
	mkdir -p $mutation_bootstrap_path
  fi
fi


if [[ ! -f $vcf_file || $vcf_file != *.vcf ]]; then
	echo "Skipping $vcf_file" 
	continue 
fi

tumor_id=$(basename ${vcf_file})
tumor_id=${tumor_id%.vcf}  

if [[ ! -f "$phi_file" ]]; then
	echo "$phi_file not found. Aborting..."
exit
fi

mutation_counts_file=$mutation_counts_path/$tumor_id.${tumor_part}phi.txt
mutation_types_file=$mutation_types_path/$tumor_id.${tumor_part}mut_types.txt	
mut_order_file=$mut_order_path/$tumor_id.${tumor_part}mut_order.txt

if [ ! -f $mutation_types_file ] || [ ! -s  $mutation_types_file ]; then
	echo "Type file..."
	perl $get_mutation_type_script muse $vcf_file $phi_file  $mut_order_file >> $mutation_types_file #2>> $log_dir/log.txt
	rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
fi

if [ ! -f $mutation_types_file ] || [ ! -s  $mutation_types_file ]; then
	echo "ERROR: $mutation_types_file file not created. Aborting..."
	exit -1
fi


num_mutations=$(min_number `cat $phi_file | wc -l` `cat $vcf_file | wc -l`)

num_hundreds=$(($num_mutations/100 + ($num_mutations % 100 > 0))) 

if [ ! -f $mutation_counts_file ]; then
	if [ $num_mutations -lt 100 ]; then
		echo "Less than 100 mutaions in a file $vcf_file or $phi_file" 
		touch $mutation_counts_file
	else
		echo "Count file..."
			for i in `seq 1 $num_hundreds`; do
			if [ $num_mutations -ge $((i*100-1)) ]; then 
					###echo $i $((i*100-1))
				python $make_hundreds_script $mutation_types_file  $((i*100-100)) $((i*100-1)) >> $mutation_counts_file #2>>$log_dir/log.txt
				rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
				fi
		done
	fi 
fi

if [[ ! -f "$mutation_counts_file" ]]; then
	echo "ERROR: $mutation_counts_file file not created. Aborting..."
	exit -1
fi


if [ "$do_bootstrap" = true ] ; then
   bootstrap_dir=$mutation_bootstrap_path/$tumor_id/
   if [ !  -z  $tumor_part  ]; then
	bootstrap_dir=${bootstrap_dir:0:${#bootstrap_dir}-1}
	bootstrap_dir=$bootstrap_dir.${tumor_part:0:${#tumor_part}-1}/
   fi

   if [ ! -d $bootstrap_dir ]; then
	mkdir $bootstrap_dir
   fi
   
   for i in `seq 1 $N_BOOTSTRAPS`; do
	mutation_bootstrap_file=$bootstrap_dir/$tumor_id.${tumor_part}$i.mut.txt
	mutation_bootstrap_counts_file=$bootstrap_dir/$tumor_id.${tumor_part}$i.mut_counts.txt
	
	mutation_bootstrap_file_unsorted=$bootstrap_dir/$tumor_id.${tumor_part}$i.mut.unsorted.txt
	mutation_bootstrap_counts_file_unsorted=$bootstrap_dir/$tumor_id.${tumor_part}$i.mut_counts.unsorted.txt

	if [ ! -f $mutation_bootstrap_file_unsorted ] || [ ! -s  $mutation_bootstrap_file_unsorted ] || [ ! -f $mutation_bootstrap_file ] || [ ! -s  $mutation_bootstrap_file ]; then

		if [ $num_mutations -lt 100 ]; then
					touch $mutation_bootstrap_file
			else
			    # echo "Bootstrap mutations..."
				python $bootstrap_mutations_script $mutation_types_file  $num_mutations  > $mutation_bootstrap_file_unsorted #2>>$log_dir/log.txt
				rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
		
				sort -k3 -nr $mutation_bootstrap_file_unsorted > $mutation_bootstrap_file
		fi
	fi

	 if [ ! -f $mutation_bootstrap_counts_file ] || [ ! -s  $mutation_bootstrap_counts_file ]; then
			# echo "Bootstrap counts..."
			for t in `seq 1 $num_hundreds`; do
				if [ $num_mutations -ge $((t*100-1)) ]; then
					python $make_hundreds_script $mutation_bootstrap_file  $((t*100-100)) $((t*100-1)) >> $mutation_bootstrap_counts_file #2>>$log_dir/log.txt
					rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
				fi
			done
		fi

	# if [ ! -f $mutation_bootstrap_counts_file_unsorted ] || [ ! -s  $mutation_bootstrap_counts_file_unsorted ]; then
	# 		#echo "Bootstrap counts on unsorted mutations..."
	# 		for t in `seq 1 $num_hundreds`; do
	# 			if [ $num_mutations -ge $((t*100-1)) ]; then
	# 				python $make_hundreds_script $mutation_bootstrap_file_unsorted  $((t*100-100)) $((t*100-1)) >> $mutation_bootstrap_counts_file_unsorted # 2>>$log_dir/log.txt
	# 			fi
	# 		done
	# fi
  
  done 
fi

