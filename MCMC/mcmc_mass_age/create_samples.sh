for i in {11..15}
do
    nohup python3 main.py -l $i >output/output_$i.txt 2>&1 &
done

