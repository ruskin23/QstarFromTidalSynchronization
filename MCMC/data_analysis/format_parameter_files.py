import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
import pickle
import argparse


class format_parameter_files():

    def __init__(self,
                 filename,
                 parameters,
                 instance,
                 format_file = None,
                 plot_file=None
                 ):

        self.instance = instance
        self.filename = filename
        self.parameters = parameters

        self.iter = {'iteration': []}
        self.data = {key:[] for key in self.parameters}

        self.value_array = {}

        self.new_fname = 'formatted_accepted_parameters_' + self.instance  + '.txt'
        print(self.new_fname)
        self.max_size = 0

        self.format_file = format_file
        self.plot_file = plot_file

    def store_file_data(self):

        with open(self.filename,'r') as f:
            reader =csv.reader(f, dialect='excel-tab')
            next(f)
            for line in reader:
                self.iter['iteration'].append(int(line[0]))
                for index,key in enumerate(self.parameters):
                    self.data[key].append(float(line[index+1]))


    def create_array(self,key):

        end_iteration = self.iter['iteration'][len(self.iter['iteration'])-1]
        self.max_size = end_iteration
        print(self.max_size)
        print(type(self.max_size))

        self.value_array[key] = np.zeros(self.max_size)

        for index,value in enumerate(self.iter['iteration']):
            self.value_array[key][value-1] = self.data[key][index]

    def fill_array(self,key):

        for index,value in enumerate(self.value_array[key]):
            if index>0 and self.value_array[key][index]==0:
                self.value_array[key][index] = self.value_array[key][index-1]

        return self.value_array

    def write_new_file(self):
        with open(self.new_fname,'w') as f:
            writer = csv.writer(f,delimiter='\t')
            writer.writerows(zip(*self.value_array.values()))


    def plot_output(self):
        logQ = []
        param = []
        index = self.parameters.index(self.plot_file)
        with open(self.new_fname,'r') as f:
            reader =csv.reader(f, dialect='excel-tab')
            for line in reader:
                param.append(float(line[index+1]))
                logQ.append(float(line[5]))

        plt.scatter(param,logQ,marker='.')
        plt.show()

    def format_files(self):
        self.store_file_data()

        for key in self.parameters:

            self.create_array(key)
            self.fill_array(key)

            for index,value in enumerate(self.value_array[key]):
                if self.value_array[key][index] > 0:
                    break
            self.value_array[key] = self.value_array[key][index:]

        self.iter['iteration'] = np.arange(1,self.max_size+1,1)
        self.iter['iteration'] = self.iter['iteration'][index:]
        self.value_array = {**self.iter,**self.value_array}

        self.write_new_file()

    def __call__(self):

        if self.format_file:
            print('formatting')
            self.format_files()
        if self.plot_file:
            print('plotting')
            self.plot_output()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action = 'store', dest = 'instance',
                        help='enter the instance of output file')
    parser.add_argument('-a', action = 'store_true', dest = 'all',
                        help='run for all files',
                        default =False)
    parser.add_argument('-p', action = 'store', dest = 'plot',
                        help='plot the given parameter')
    parser.add_argument('-f', action = 'store_true', dest = 'format_file',
                        help='format the parameter files and store in new file',
                        default =False)

    args = parser.parse_args()

    parameters = [  'age', 'teff', 'feh', 'Wdisk', 'logQ']

    if args.instance:
        fname = 'accepted_parameters_' + args.instance + '.txt'
        fill_values = format_parameter_files(fname,parameters,args.instance,args.format_file,args.plot)
        fill_values()

    elif args.all:
        for k in range(3,8):
            print(k)
            i = str(k)
            fname = 'accepted_parameters_' + i + '.txt'
            fill_values = format_parameter_files(fname,parameters,i,args.format_file,args.plot)
            fill_values()




