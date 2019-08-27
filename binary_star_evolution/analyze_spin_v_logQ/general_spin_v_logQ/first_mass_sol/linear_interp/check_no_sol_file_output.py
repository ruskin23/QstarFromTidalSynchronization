with open('no_sol_files.csv','r') as f:
    for lines in f:
        with open(lines,'r') as f1:
            for lines1 in f1:
                print(lines1)
