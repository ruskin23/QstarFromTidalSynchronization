import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
import pickle
import argparse


class format_parameter_files():

    def __init__(self,
                 parameters,
                 instance=None,
                 combine=None,
                 format_file=None,
                 plot_file=None
                 ):

        self.parameters = parameters
        self.instance = instance
        self.combine = combine
        self.format_file = format_file
        self.plot_file = plot_file

        self.iter = {'iteration': []}
        self.data = {key:[] for key in self.parameters}

        self.value_array = {}

        if self.instance:
            self.filename = 'accepted_parameters_' + self.instance + '.txt'
            self.new_fname = 'formatted_accepted_parameters_' + self.instance  + '.txt'
        if self.combine:
            self.fnames = ['formatted_accepted_parameters_'+str(k)+'.txt' for k in range(1,5)]
            self.new_fname = 'combined_accepted_parameters.txt'

        self.max_size = 0


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

    def combine_files(self):

        with open(self.new_fname,'w') as f_out:
            for name in self.fnames:
                with open(name,'r') as f_in:
                    for line in f_in:
                        f_out.write(line)

    def plot_output(self):
        logQ = []
        param = []
        index = self.parameters.index(self.plot_file)
        with open(self.new_fname,'r') as f:
            reader =csv.reader(f, dialect='excel-tab')
            for line in reader:
                param.append(float(line[index+1]))
                logQ.append(float(line[5]))

        plt.scatter(param,logQ,marker='x')
        #plt.hist2d(param,logQ,bins=50)
        plt.show()

    def __call__(self):

        if self.format_file:
            print('formatting')
            self.format_files()
        if self.combine:
            print('combing')
            self.combine_files()
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
                        help='plot a single file for given parameter')
    parser.add_argument('-f', action = 'store_true', dest = 'format_file',
                        help='format the parameter files and store in new file',
                        default =False)
    parser.add_argument('-c', action = 'store_true', dest = 'combine',
                        help='combine all files',
                        default =False)

    args = parser.parse_args()

    parameters = [  'teff', 'feh', 'Porb', 'logg', 'Wdisk', 'logQ']

    if args.instance:
        fill_values = format_parameter_files(parameters,args.instance,None,args.format_file,args.plot)
        fill_values()

    elif args.all:
        if args.combine: p = None
        else: p = args.plot
        for k in range(3,8):
            i = str(k)
            fill_values = format_parameter_files(parameters,i,None,args.format_file,p)
            fill_values()


    if args.combine:
        fill_values = format_parameter_files(parameters,None,True,None,args.plot)
        fill_values()



