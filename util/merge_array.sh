#!/usr/bin/env bash
usage="
------------------------------
Merge Arrays if they shared elements
------------------------------
USAGE:
	$0 <arrays.csv>
eg:

save data below to test.csv:
		d,a,b,c
		b,c,a
		b,c,d
		c,d
		d,e
		z,x

and run:
  		$0 test.csv

Will get:
		a,b,c,d,e
		x,z

------------------------------
                   Songtao Gui
              songtaogui@sina.com
"
if [[ $# -ne 1 ]]; then 
	echo "$usage"
	exit 1
fi

in=$1

if [ ! -s "$in" ];then
	echo "No file: $in" >&2
	exit 1
fi
i=1
echo "Iterative run until md5sum values are same." >&2
cat $in | md5sum > ${in}.iter.md5sum
echo "Run round $i ..." >&2
cat $in | perl -F"," -lane '
@a=sort @F;
foreach $a (@a){
	push @{$hash{$a}},@a;
}
END{
	foreach $v (sort keys %hash){ 
		my %u;
		@u=grep {++$u{$_}==1} @{$hash{$v}};
		$,=",";
		print @u;
	}
}' | perl -ne 'print if ++$h{$_}==1' >${in}.round$i
if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> round$i: non-zero exit." >&2
	rm -f ${in}.round$i
	exit 1
fi
cat ${in}.round$i | md5sum >> ${in}.iter.md5sum &&\
echo "Done round $i ." >&2

while true
do
	md5_1=$(tac ${in}.iter.md5sum | sed -n '1p')
	md5_2=$(tac ${in}.iter.md5sum | sed -n '2p')
	if [ "$md5_1" == "$md5_2" ];then
		echo "Round ${i}: '$md5_1' == '$md5_2' Finished " >&2
		echo "Done. Final output in : ${in}.round$i" >&2
		exit;
	else 
		echo "Round ${i}: '$md5_1' != '$md5_2' Countinue ..." >&2
		cur_in=${in}.round$i
		let i++
		echo "Run round $i using $cur_in as input ..." >&2
		cat $cur_in | perl -F"," -lane '
		@a=sort @F;
		foreach $a (@a){
			push @{$hash{$a}},@a;
		}
		END{
			foreach $v (sort keys %hash){ 
				my %u;
				@u=grep {++$u{$_}==1} @{$hash{$v}};
				$,=",";
				print @u;
			}
		}' | perl -ne 'print if ++$h{$_}==1' >${in}.round$i
		if [[ $? -ne 0 ]] ; then 
			echo "[ERROR] --> round$i: non-zero exit." >&2
			rm -f ${in}.round$i
			exit 1
		fi
		cat ${in}.round$i | md5sum >> ${in}.iter.md5sum &&\
		echo "Done round $i ." >&2
	fi
done

