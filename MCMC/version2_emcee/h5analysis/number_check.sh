 
while true
do
    echo "$(ps ax | grep 'python3 h5_analysis.py' | wc -l)"
    sleep 5
done
