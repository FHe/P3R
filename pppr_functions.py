"""
Functions used in pppr-gui

Authors/modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""
####################################################

import numpy as Num
import random
##############################  extract datasets and read parameters ############################
class DataSet:
    def __init__(self):
        self.x_label = ''
        self.y_label = ''
        self.n = 0
        self.weight = 1.0
        self.plot = 0
        self.lower = 0
        self.upper = 0

def eval_datastring(datastring):
    data = Num.ndarray((0,7),float)
    datasets = []
    datastring = str(datastring)
    tmp = datastring.split("\n")
    counter = 0
    for i in range(len(tmp)):
        if tmp[i] =='':
            pass
        else:
            tmp2 = str.rsplit(tmp[i])
            try:
                tmp3 = []
                for a in tmp2:
                    tmp3.append(float(a))
                for a in range(4):
                    tmp3.append(0.)
                data = Num.append(data, [tmp3], axis = 0)
                counter += 1
            except ValueError:
                if counter == 0: 
                    a = DataSet()
                    a.x_label = tmp2[0]
                    a.y_label = tmp2[1]
                    a.lower = counter
                    datasets.append(a)
                else:
                    a = DataSet()
                    a.x_label = tmp2[0]
                    a.y_label = tmp2[1]
                    a.lower = counter
                    datasets[-1].upper = counter
                    datasets[-1].n = datasets[-1].upper - datasets[-1].lower
                    datasets.append(a)
    datasets[-1].upper = counter
    datasets[-1].n = datasets[-1].upper - datasets[-1].lower
    data = Num.transpose(data)
    return data, datasets

def read_parameters(parameterfile):
    parameter={}
    param_labels = []
    f = open(parameterfile, 'r')
    data = f.readlines()
    f.close()
    for i in data:
        tmp = str.rsplit(i)
        if tmp[0] != '%':
            if len(tmp) == 5:
                parameter[tmp[0]]= [float(tmp[1]),float(tmp[2]),float(tmp[3]),False, 0., '']
                if tmp[4] == 'True':
                    parameter[tmp[0]][3] = True
            elif len(tmp) == 6:
                parameter[tmp[0]]= [float(tmp[1]),float(tmp[3]),float(tmp[4]),False, float(tmp[2]), '']
                if tmp[5] == 'True':
                    parameter[tmp[0]][3] = True
                
            param_labels.append(tmp[0])

    return parameter, param_labels
########################  writing parameter file  #######################################  
def write_par(parameter, param_labels, filename = 'parameters.new'):
    f = file(filename, 'w')
    f.write('% param_label              value     std-dev          min          max    '+\
            'refine_flag\n')
    for i in param_labels:
        line = "%13s %18.12f %12.8f %12.4f %12.4f %10s\n" % (i,parameter[i][0],
                                                      parameter[i][4],
                                                      parameter[i][1],
                                                      parameter[i][2],
                                                      parameter[i][3])
        f.write(line)
    f.close()
#

    
    
    



