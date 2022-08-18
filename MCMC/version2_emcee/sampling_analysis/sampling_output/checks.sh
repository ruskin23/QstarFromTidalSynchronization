
> stored_steps.txt

for D in 1 2 3 4 5 6 7 8 10
do 
	for S in node_$D/system_*
       	do 
		spin_count=$(grep 'Final_Spin' $S/*.log | wc -l)
		N=64
		steps=$(($spin_count/$N))
		split1=$(echo $S | cut -d "/" -f 2)
		node=$(echo $S | cut -d "/" -f 1)
		split2=$(echo $split1 | cut -d "_" -f 2)
		echo "$node $split2 $spin_count $steps"	
       	done >>stored_steps.txt 
done
