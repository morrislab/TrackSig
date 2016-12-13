dir=$1
tmp_dir=$2

echo $(basename  $dir)

mkdir $tmp_dir/$(basename  $dir)

for d in `ls -d $dir/*`; do
	echo $(basename $d)
	archive_name=$tmp_dir/$(basename  $dir)/$(basename $d).tar
	tar -czf $archive_name $d
	n_files=`tar tzf $archive_name | grep -e "[^/]$" | wc -l`
	n_files_dir=`find $d -type f | wc -l`
	echo $n_files
	echo $n_files_dir
done

for f in `ls -p $dir | grep -v /`;do
	cp $f  $tmp_dir/$(basename  $dir)/ 
done

tar -czf $tmp_dir/$(basename  $dir).tar $tmp_dir/$(basename  $dir)
mv $tmp_dir/$(basename  $dir).tar $dir/..
