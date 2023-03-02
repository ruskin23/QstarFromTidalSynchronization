import csv
import sys
import glob
import pandas as pd

def get_values(filename):

    random_list = []
    with open(filename, 'r') as f:
        for lines in f:
            lines = lines.strip('\n')

            if 'parameters' in lines:
                x = lines.split('=')
                if '[' not in x[1]:
                    random_list.append(x[1])
                if '[' in x[1]:
                    y = x[1][1:-2].split()
                    for vals in y:
                        random_list.append(vals)

            if 'Initial_Eccentricity' in lines or 'Final_Orbital_Period' in lines:
                x = lines.split(',')
                for vals in x:
                    random_list.append(vals.split('=')[1])


            if 'log_likelihood' in lines and '-inf' in lines:
                return

    return random_list


def write_on_file(filepath):

    #files = ['10031409_20220923105922_3405297.log']
    files = glob.glob(filepath + '/*.log')

    data_file = 'parsed.csv'
    if glob.glob(data_file): action = 'a'
    else: action = 'w'

    fnew = open(data_file, action)

    writer = csv.writer(fnew)

    for filename in files:
        values = get_values(filename)
        if values:
            writer.writerow(values)
    fnew.close()


def correction():

    data_file = 'raw_data_29.csv'
    fnew = open(data_file, 'w', newline='')
    writer = csv.writer(fnew)
    header = ['m1', 'm2', 'feh', 'age', 'wdisk', 'delta0', 'break', 'power', 'porb_i', 'e_i', 'porb_f', 'e_f']
    writer.writerow(header)

    with open('parsed.csv', 'r') as f:
        for i, lines in enumerate(f):
            lines = lines.strip('\n')
            x = lines.split(',')
            if len(x) == 13:
                del x[7]
            if len(x) == 15:
                indexes = [6, 8, 10]
                for idx in sorted(indexes, reverse=True):
                    del x[idx]
            if len(x) == 12: writer.writerow(x)
    fnew.close()


if __name__ == '__main__':

    filepath = sys.argv[1]
    write_on_file(filepath)
    correction()
