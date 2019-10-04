for i in {113..130}
do
    nohup python3 main.py -l $i >output/output_$i.txt 2>&1 &
done

