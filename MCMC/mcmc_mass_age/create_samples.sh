for s in 121 123 124 125 126 128 129 130 131 132 133 134 137 138 139 140 141 142
do
    nohup python3 main.py -l $s -i 1 >output/output_1_$s.txt 2>&1 &
done

