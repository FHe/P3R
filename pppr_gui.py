"""
Functions and classes used to build the pppr GUI

Authors/modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""
###############################################################################

import wx
import wx.grid as gridlib
import os
import time
import random
import numpy as Num
from phreeqpy.iphreeqc.phreeqc_dll import IPhreeqc

from simplex import simplex, LM_fit
from pppr_functions import read_parameters, eval_datastring, write_par
from plotframe import createplotframe
###############################################################################
class wxPPPRFrame(wx.Frame):
    def __init__(self, parent, title, size):
        wx.Frame.__init__(self, parent, title= title, size=size)

        # A status bar
        self.dirname = ''
        self.statusbar = self.CreateStatusBar(2)

        # Setting up the ReadFilesmenu.
        filemenu= wx.Menu()
        menuReadModel = filemenu.Append(wx.ID_ANY,"&Read Modelfile"," Read in a Model file")
        menuReadData = filemenu.Append(wx.ID_ANY,"&Read Datafile"," Read in a Data file")
        menuReadDatabase = filemenu.Append(wx.ID_ANY,"&Load Database"," Load Database from a file")
        menuReadParameter = filemenu.Append(wx.ID_ANY,"&Read Parameterfile"," Read in a Parameter file")
        filemenu.AppendSeparator()
        menuExit = filemenu.Append(wx.ID_EXIT,"E&xit"," Terminate the program")

        #Setting up the WriteFiles
        writemenu = wx.Menu()
        menuWriteModel = writemenu.Append(wx.ID_ANY, "&Write .mod file", " Write PhreeqC Model with param. labels to a .mod file")
        menuWriteModel2 = writemenu.Append(wx.ID_ANY, "&Write .phrq file", " Write PhreeqC Model with best fit params to a .phrq file")
        menuWritePar = writemenu.Append(wx.ID_ANY, "&Write Parameter file", " Write parameters to a .par file")
        menuWriteData = writemenu.Append(wx.ID_ANY, "&Write Data input file", " Write data to a .dat file")
        menuWriteData2 = writemenu.Append(wx.ID_ANY, "&Write Results file", " Write data and fit results to a .dat file")

        
        # Creating the menubar.
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&Read Files") # Adding the "filemenu" to the MenuBar
        menuBar.Append(writemenu,"&Write Files")
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        
        # define events
        self.Bind(wx.EVT_MENU, self.OnReadData, menuReadData)
        self.Bind(wx.EVT_MENU, self.OnReadDatabase, menuReadDatabase)
        self.Bind(wx.EVT_MENU, self.OnReadModel, menuReadModel)
        self.Bind(wx.EVT_MENU, self.OnReadParameter, menuReadParameter)
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)

        self.Bind(wx.EVT_MENU, self.OnWriteModel, menuWriteModel)
        self.Bind(wx.EVT_MENU, self.OnWriteModel2, menuWriteModel2)
        self.Bind(wx.EVT_MENU, self.OnWritePar, menuWritePar)
        self.Bind(wx.EVT_MENU, self.OnWriteData, menuWriteData)
        self.Bind(wx.EVT_MENU, self.OnWriteData2, menuWriteData2)

        self.nb = PPPRNotebook(self)

        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.Show(True)
        

    def OnReadData(self,e):
        dlg = wx.FileDialog(self, "Choose a data file ", self.dirname, ".dat", "*.dat", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            f = open(self.dirname+'/'+filename, 'r')
            data = f.readlines()
            f.close()
            model = ''
            for line in data:
                model = model + line
            self.nb.datastring = model
            self.nb.DataPage.control.SetValue(self.nb.datastring)  
            self.nb.data, self.nb.datasets = eval_datastring(self.nb.datastring)
            self.nb.MainControlPage.fill_grid()
            self.nb.SetSelection(2)
        dlg.Destroy()

    def OnReadDatabase(self,e):
        dlg = wx.FileDialog(self, "Choose a database ", self.dirname, ".dat", "*.dat", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            self.nb.database = dirname+'/'+filename
            self.nb.database_loaded = True
        dlg.Destroy()

    def OnReadModel(self,e):
        """ Read in a Model file"""
        dlg = wx.FileDialog(self, "Choose a Model file", self.dirname, ".mod", "*.mod", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            f = open(self.dirname+'/'+filename, 'r')
            data = f.readlines()
            f.close()
            model = ''
            for line in data:
                model = model + line
            self.nb.batch_model = model
            self.nb.ModelPage.control.SetValue(self.nb.batch_model)  
            self.nb.SetSelection(1)
        dlg.Destroy()


    def OnReadParameter(self,e):
        """ Read in parameter file"""
        dlg = wx.FileDialog(self, "Choose a parameter file", self.dirname, ".par", "*.par", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.nb.parameter, self.nb.param_labels = read_parameters(self.dirname+'/'+filename)
           
            self.nb.DeletePage(3)
            self.nb.ParameterPage = ParameterPanel(self.nb)
            self.nb.AddPage(self.nb.ParameterPage, " Parameters ")
                        
            for i in range(len(self.nb.param_labels)):

                control0_tmp = wx.Button(self.nb.ParameterPage, 20000+i+7*len(self.nb.param_labels), label = self.nb.param_labels[i], pos=(20, 23*(i+1)+20), size=(120, 20))
                self.nb.ParameterPage.control0.append(control0_tmp)
                self.Bind(wx.EVT_BUTTON, self.nb.ParameterPage.clicklabel , self.nb.ParameterPage.control0[i])
                
                control1_tmp = wx.TextCtrl(self.nb.ParameterPage,20000+i, pos=(150, 23*(i+1)+20), size=(120,20))
                self.nb.ParameterPage.control1.append(control1_tmp)
                self.nb.ParameterPage.control1[i].SetValue(str(round(self.nb.parameter[self.nb.param_labels[i]][0], 12)))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparvalue, self.nb.ParameterPage.control1[i])
                
                control2_tmp = (wx.TextCtrl(self.nb.ParameterPage,20000+i+len(self.nb.param_labels), pos=(370, 23*(i+1)+20), size=(80,20)))
                self.nb.ParameterPage.control2.append(control2_tmp)
                self.nb.ParameterPage.control2[i].SetValue(str(self.nb.parameter[self.nb.param_labels[i]][1]))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparmin, self.nb.ParameterPage.control2[i])
                
                control3_tmp = (wx.TextCtrl(self.nb.ParameterPage,20000+i+2*len(self.nb.param_labels), pos=(460, 23*(i+1)+20), size=(80,20)))
                self.nb.ParameterPage.control3.append(control3_tmp)
                self.nb.ParameterPage.control3[i].SetValue(str(self.nb.parameter[self.nb.param_labels[i]][2]))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparmax, self.nb.ParameterPage.control3[i])
                                
                control4_tmp = (wx.CheckBox(self.nb.ParameterPage,20000+i+3*len(self.nb.param_labels), label = '', pos = (550, 23*(i+1)+25)))
                self.nb.ParameterPage.control4.append(control4_tmp)
                self.nb.ParameterPage.control4[i].SetValue(self.nb.parameter[self.nb.param_labels[i]][3])
                wx.EVT_CHECKBOX(self.nb.ParameterPage, self.nb.ParameterPage.control4[i].GetId(), self.nb.ParameterPage.editparstate)

                control5_tmp = wx.Button(self.nb.ParameterPage, 20000+i+4*len(self.nb.param_labels), label = '<', pos = (610, 23*(i+1)+20), size = (20,20))
                self.nb.ParameterPage.control5.append(control5_tmp)
                self.Bind(wx.EVT_BUTTON, self.nb.ParameterPage.toggleminus , self.nb.ParameterPage.control5[i])

                control6_tmp = (wx.TextCtrl(self.nb.ParameterPage,20000+i+5*len(self.nb.param_labels), pos=(640, 23*(i+1)+20), size=(50,20)))
                self.nb.ParameterPage.control6.append(control6_tmp)
                self.nb.ParameterPage.control6[i].SetValue('0')
                self.nb.ParameterPage.togglesteps.append(0)
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.togglestep, self.nb.ParameterPage.control6[i])

                control7_tmp = wx.Button(self.nb.ParameterPage,20000+i+6*len(self.nb.param_labels), label = '>', pos = (700, 23*(i+1)+20), size = (20,20))
                self.nb.ParameterPage.control7.append(control7_tmp)
                self.Bind(wx.EVT_BUTTON, self.nb.ParameterPage.toggleplus , self.nb.ParameterPage.control7[i])
                
                control8_tmp = wx.TextCtrl(self.nb.ParameterPage,20000+i+8*len(self.nb.param_labels), pos=(280, 23*(i+1)+20), size=(80,20))
                self.nb.ParameterPage.control8.append(control8_tmp)
                self.nb.ParameterPage.control8[i].SetValue(str(round(self.nb.parameter[self.nb.param_labels[i]][4], 8)))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparstddev, self.nb.ParameterPage.control8[i])

            self.nb.ParameterPage.SetScrollbars(0, 10, 0, int((len(self.nb.param_labels)+4)*2.3)+1)
            self.nb.SetSelection(3)
        dlg.Destroy()
            
    def OnExit(self,e):
        self.Close(True)  # Close the frame.

    def OnWriteModel(self,e): #write the model with parameter labels to a .mod file
        dlg = wx.FileDialog(self, "Write model to .mod file", self.dirname, ".mod", "*.mod", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            f = file(self.dirname+'/'+filename, 'w')
            f.write(self.nb.batch_model)
            f.close()
        dlg.Destroy()

    def OnWriteModel2(self,e): #write the model with best fit parameters inserted to a .phrq file - the output should be a vaild PhreeqC input file
        dlg = wx.FileDialog(self, "Write model to .phrq file", self.dirname,".phrq","*.phrq", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.nb.update_model()
            f = file(self.dirname+'/'+filename, 'w')
            f.write(self.nb.modelstring)
            f.close()
        dlg.Destroy()

    def OnWritePar(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Write parameters to .par file", dirname, ".par", "*.par", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            write_par(self.nb.parameter, self.nb.param_labels, filename)
        dlg.Destroy()

    def OnWriteData(self,e):
        dlg = wx.FileDialog(self, "Write Data .dat file", self.dirname, ".dat", "*.dat", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            f = file(self.dirname+'/'+filename, 'w')
            f.write(self.nb.datastring)
            f.close()
        dlg.Destroy()

    def OnWriteData2(self,e):
        dlg = wx.FileDialog(self, "Write Data and Fit Results to .dat file", self.dirname, ".dat", "*.dat", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            f = file(self.dirname+'/'+filename, 'w')
            f.write(' x_value        y_value        y_err          y_model        residual       chi**2    lower_bound    upper_bound   \n')
            dat = self.nb.data
            for i in range(len(dat[0])):
                line = " %12.6e   %12.6e   %12.6e   %12.6e   %12.6e   %12.6e   %12.6e   %12.6e\n" %\
                       (dat[0][i],dat[1][i],dat[2][i],dat[3][i],dat[4][i],dat[5][i], (dat[3][i]-dat[6][i]), (dat[3][i]+dat[6][i]))
                f.write(line)
            f.close()
        dlg.Destroy()
        
############################################################################################################
class PPPRNotebook(wx.Notebook):
    def __init__(self, parent):
        wx.Notebook.__init__(self,parent)

        self.Phreeq = None
        self.database_loaded = False
        self.database = ''
        self.Phreeq_age = 0

        self.batch_model = ''
        self.modelstring = ''
        self.datastring = ''
        self.data= []
        self.datasets = []
        self.chi2 = 0

        self.sensitivities = Num.array([])
        self.correl_matrix = Num.array([])
        self.cov_matrix = Num.array([])
        self.fpc = 0.0001
        self.used_params = ['None']

        self.mcn = 50
        self.mcmode = 0
        self.plotmode = 0

        self.parameter = {}
        self.param_labels = []
        self.frame = self.GetParent()

        self.FigureFrame = None
        self.Figure1 = None

        self.MainControlPage = MainControlPanel(self)
        self.ParameterPage = ParameterPanel(self)
        self.DataPage = DataPanel(self)
        self.ModelPage = ModelPanel(self)

        self.AddPage(self.MainControlPage, " Main Controls " )
        self.AddPage(self.ModelPage, " Model ")
        self.AddPage(self.DataPage, " Data ")
        self.AddPage(self.ParameterPage, " Parameters ")

    def phreeq_check(self):
        if self.Phreeq == None or self.Phreeq_age >= 20:
            self.Phreeq = None
            self.Phreeq = IPhreeqc()
            self.Phreeq.load_database(self.database)
            self.Phreeq.set_selected_output_file_off()
            self.Phreeq_age = 0
        else:
            self.Phreeq_age += 1

    def plot(self):
        if self.data != []:
            plot_dims = self.MainControlPage.plotdims
            if self.FigureFrame == None:
                self.FigureFrame = createplotframe(self, "PPPR Data and Model plot", (1200,900))
                self.FigureFrame.Show(True)
                self.Figure1 = self.FigureFrame.figure
            b=0
            for bl in self.datasets:
                if bl.plot > -1: b+=1
            if b > (plot_dims[0] * plot_dims[1]):
                print "\nPlot Dimensions are too small \nfor the number of Datasets to be displayed!"
                pass
            else:
                self.Figure1.clear()
                self.Figure1.suptitle('chi**2 = '+str(round(self.chi2,12)), fontsize=20)
                subplots = []
                other = []
                j = 0
                for i in range(len(self.datasets)):
                    if self.datasets[i].plot > -1:
                        subplots.append(self.Figure1.add_subplot(plot_dims[0],plot_dims[1],j+1))
                        tmp = self.datasets[i]
                        x = self.data[0][tmp.lower:tmp.upper]
                        subplots[j].plot(x,self.data[3][tmp.lower:tmp.upper],'k')
                        subplots[j].errorbar(x,self.data[1][tmp.lower:tmp.upper]\
                                         ,self.data[2][tmp.lower:tmp.upper], fmt = 'bo')
                        subplots[j].set_xlabel(tmp.x_label)
                        subplots[j].set_ylabel(tmp.y_label)
                        if self.plotmode == 1:
                            other.append(subplots[j].twinx())
                            other[-1].bar(x,self.data[4][tmp.lower:tmp.upper], width = 0.03, bottom = 0,\
                                          color = 'c', edgecolor = 'c', alpha = 0.25)
                        elif self.plotmode == 2:
                            y1 = self.data[3][tmp.lower:tmp.upper]+self.data[6][tmp.lower:tmp.upper]
                            y2 = self.data[3][tmp.lower:tmp.upper]-self.data[6][tmp.lower:tmp.upper]
                            subplots[j].plot(x, y1, 'g')
                            subplots[j].plot(x, y2, 'g')
                            subplots[j].fill_between(x, y1, y2, facecolor = 'green', alpha = 0.25)
                        if tmp.plot == 0:
                            subplots[j].set_yscale('linear')
                            subplots[j].set_xscale('linear')
                        elif tmp.plot == 1:
                            subplots[j].set_yscale('log')
                            subplots[j].set_xscale('linear')
                        elif tmp.plot == 2:
                            subplots[j].set_yscale('linear')
                            subplots[j].set_xscale('log')
                        elif tmp.plot == 3:
                            subplots[j].set_yscale('log')
                            subplots[j].set_xscale('log')
                        j +=1
                self.Figure1.canvas.draw()

    def plot_sens(self,dp,sens,param):
        if self.data != []:
            plot_dims = self.MainControlPage.plotdims
            if self.FigureFrame == None:
                self.FigureFrame = createplotframe(self, "PPPR Data and Model plot", (1200,900))
                self.FigureFrame.Show(True)
                self.Figure1 = self.FigureFrame.figure
            b=0
            for bl in self.datasets:
                if bl.plot > -1: b+=1
            if b > (plot_dims[0] * plot_dims[1]):
                print "\nPlot Dimensions are too small \nfor the number of Datasets to be displayed!"
                pass
            else:
                self.Figure1.clear()
                self.Figure1.suptitle('scaled sensitivities for parameter ' +param+ ', standard deviation = '+str(round(dp,12)), fontsize=16)
                subplots = []
                sensplots = []
                j = 0
                for i in range(len(self.datasets)):
                    if self.datasets[i].plot > -1:
                        subplots.append(self.Figure1.add_subplot(plot_dims[0],plot_dims[1],j+1))
                        sensplots.append(subplots[j].twinx())
                        tmp = self.datasets[i]
                        subplots[j].plot(self.data[0][tmp.lower:tmp.upper],self.data[3][tmp.lower:tmp.upper],'k')
                        subplots[j].errorbar(self.data[0][tmp.lower:tmp.upper],self.data[1][tmp.lower:tmp.upper]\
                                         ,self.data[2][tmp.lower:tmp.upper], fmt = 'bo')
                        subplots[j].set_xlabel(tmp.x_label)
                        subplots[j].set_ylabel(tmp.y_label)
                        sensplots[j].bar(self.data[0][tmp.lower:tmp.upper],sens[tmp.lower:tmp.upper], width = 0.03,\
                                         bottom = 0, color = 'r', edgecolor = 'r', alpha = 0.25)
                        if tmp.plot == 0:
                            subplots[j].set_yscale('linear')
                            subplots[j].set_xscale('linear')
                        elif tmp.plot == 1:
                            subplots[j].set_yscale('log')
                            subplots[j].set_xscale('linear')
                        elif tmp.plot == 2:
                            subplots[j].set_yscale('linear')
                            subplots[j].set_xscale('log')
                        elif tmp.plot == 3:
                            subplots[j].set_yscale('log')
                            subplots[j].set_xscale('log')
                        j +=1
                self.Figure1.canvas.draw()
    
    def update_model(self):
        keys = self.parameter.keys()
        for key in keys:
            if self.parameter[key][5] in keys:
                self.parameter[key][0] = self.parameter[self.parameter[key][5]][0]
        self.modelstring = self.batch_model[:]
        for key in keys:
            self.modelstring = self.modelstring.replace(key, str(self.parameter[key][0]))

    def get_used_params(self):
        self.used_params = []
        for i in self.parameter.keys():
            if self.parameter[i][3]:
                self.used_params.append(i)
        
    def model(self):
        self.phreeq_check()
        self.get_used_params()
        self.update_model()
        if self.modelstring != '':
            if self.database_loaded:
                self.Phreeq.run_string(self.modelstring)
                output = self.Phreeq.get_selected_output_array()
                output = output[(len(output)-len(self.data[0])):len(output)]
                if len(output) != len(self.data[0]):
                    print 'model output structure does not match data structure'
                    print output
                    print len(output), len(self.data[0])
                    pass
                else:
                    output = Num.transpose(output)
                    for i in range(len(output[0])):
                        if round(self.data[0][i],1) != round(output[0][i],1):
                            print 'Warning model conditions do not match data'
                            print round(self.data[0][i],1),  round(output[0][i],1)
                        self.data[3][i] = output[1][i]
                    self.data[4] = self.data[1]-self.data[3]
                    self.data[5] = (self.data[4]/self.data[2])**2
                    self.chi2 = 0
                    n = 0
                    for dat in self.datasets:
                        self.chi2 = self.chi2 + Num.sum(self.data[5][dat.lower:dat.upper])*dat.weight
                        n = n+ dat.n
                    self.chi2 = self.chi2/(n-len(self.used_params))
            else:
                print 'need to load a database before you can start modelling'
                pass
        else:
            pass

    def model_modified(self,param, h):
        self.phreeq_check()
        self.parameter[param][0] += h
        self.update_model()
        if self.modelstring != '':
            if self.database_loaded:
                self.Phreeq.run_string(self.modelstring)
                output = self.Phreeq.get_selected_output_array()
                output = output[(len(output)-len(self.data[0])):len(output)]
                output = Num.transpose(output)
        self.parameter[param][0] += -h
        self.update_model()
        return output[1]
        

    def statistics(self):
        #setup w matrix
        n = len(self.data[0])
        w = Num.zeros((n,n))
        for i in range(n):
            w[i][i] = (1/self.data[2][i])**2          
            
        #setup Sensitivity matrix X
        self.get_used_params()
        b = len(self.used_params)
        X = Num.zeros((b,n))
    
        for i in range(b):
            h = self.parameter[self.used_params[i]][0] * self.fpc
            if h == 0:
                h = self.fpc
            if self.parameter[self.used_params[i]][0] -0.5*h >= self.parameter[self.used_params[i]][1] \
               and self.parameter[self.used_params[i]][0] +0.5*h <= self.parameter[self.used_params[i]][2]:
            
                y1 = self.model_modified(self.used_params[i],-0.5*h)
                y2 = self.model_modified(self.used_params[i], 0.5*h)
                y1 = (y2 - y1)/h
            
            elif self.parameter[self.used_params[i]][0] -0.5*h < self.parameter[self.used_params[i]][1]:
                
                y1 = self.model_modified(self.used_params[i],0)
                y2 = self.model_modified(self.used_params[i],h)
                y1 = (y2 - y1)/h
            
            elif self.parameter[self.used_params[i]][0] +0.5*h > self.parameter[self.used_params[i]][2]:
                
                y1 = self.model_modified(self.used_params[i],-h)
                y2 = self.model_modified(self.used_params[i],0)
                y1 = (y2 - y1)/h
            for j in range(n):
                X[i][j] = y1[j] * self.parameter[self.used_params[i]][0] * Num.sqrt(w[j][j])

        self.model()
            
        V = Num.dot(X, Num.dot(w, Num.transpose(X)))
        C = Num.zeros((b,b))
        try:
            V = self.chi2 * Num.linalg.inv(V)
            for i in range(b):
                for j in range(b):
                    if i == j:
                        C[i][j] = Num.sqrt(V[i][j])
                        self.parameter[self.used_params[i]][4] = C[i][j]
                    else:
                        C[i][j] = V[i][j]/Num.sqrt(V[i][i]*V[j][j])
            keys = self.parameter.keys()
            for key in keys:
                if self.parameter[key][5] in keys:
                    self.parameter[key][4] = self.parameter[parameter[key][5]][4]

            self.sensitivities = X
            self.correl_matrix = C
            self.cov_matrix = V
        except:
            print "There's a hole in the matrix !!!\nOmit unused or insignificant parameters in the calculation.\n"
        
        return
    
    def single_param_sensitivity(self, param):
        #setup w matrix
        n = len(self.data[0])
        w = Num.zeros((n,n))
        for i in range(n):
            w[i][i] = (1/self.data[2][i])**2          
            
        #setup Sensitivity matrix X
        X = Num.zeros((n))
        h = self.parameter[param][0] * self.fpc
        if h == 0:
            h = self.fpc
        if self.parameter[param][0] -0.5*h >= self.parameter[param][1] \
            and self.parameter[param][0] +0.5*h <= self.parameter[param][2]:
            
            y1 = self.model_modified(param,-0.5*h)
            y2 = self.model_modified(param, 0.5*h)
            y1 = (y2 - y1)/h
            
        elif self.parameter[param][0] -0.5*h < self.parameter[param][1]:
                
            y1 = self.model_modified(param,0)
            y2 = self.model_modified(param,h)
            y1 = (y2 - y1)/h
            
        elif self.parameter[param][0] +0.5*h > self.parameter[param][2]:
                
            y1 = self.model_modified(param,-h)
            y2 = self.model_modified(param,0)
            y1 = (y2 - y1)/h
        for j in range(n):
            X[j] = y1[j] * self.parameter[param][0] * Num.sqrt(w[j][j])

        self.model()
            
        V = Num.dot(X, Num.dot(w, Num.transpose(X)))
        if V > 0:
            dp = Num.sqrt(self.chi2/V)
        else:
            dp = 0.

        return X, dp
############################################################################################################
class MainControlPanel(wx.Panel):
    def __init__(self,parent):
        wx.Panel.__init__(self,parent)

        self.fitalgo = 0
        self.plotdims = [3,3]
        self.simplex_params = [1.0,0.5,2.0,0.2,1e-6,10000, False]

        self.nb = self.GetParent()

        self.StopFit = False
        self.auto_fit = False
        self.auto_fit_n = 3
        self.FigureFrame = None
        self.Figure2 = None
        
        #Panel Headings
        wx.StaticText(self, label = 'Data Information: ', pos=(20, 12), size=(100, 20))
        wx.StaticText(self, label = 'Fit options: ', pos=(520, 12), size=(100, 20))
        wx.StaticLine(self, pos = (500,0), size = (5,650), style = wx.LI_VERTICAL)
        wx.StaticLine(self, pos = (500,115), size = (285,5), style = wx.LI_HORIZONTAL)
        wx.StaticLine(self, pos = (500,235), size = (285,5), style = wx.LI_HORIZONTAL)
        wx.StaticLine(self, pos = (500,345), size = (285,5), style = wx.LI_HORIZONTAL)
        wx.StaticLine(self, pos = (500,455), size = (285,5), style = wx.LI_HORIZONTAL)
        wx.StaticLine(self, pos = (500,565), size = (285,5), style = wx.LI_HORIZONTAL)

        #data set display in grid structure
        self.datagrid = gridlib.Grid(self, id = wx.ID_ANY, pos = (20,40), size = (460,600))
        self.datagrid.CreateGrid(25,5)
        for i in range(25):
            self.datagrid.SetRowSize(i,25)
            for j in range(5):
                self.datagrid.SetCellAlignment(i,j,wx.ALIGN_CENTRE,wx.ALIGN_CENTRE)
                if i == 0:
                    self.datagrid.SetColSize(j,70)
                    if j == 0:
                        self.datagrid.SetColLabelValue(j,'x label')
                    elif j == 1:
                        self.datagrid.SetColLabelValue(j,'y label')
                    elif j == 2:
                        self.datagrid.SetColLabelValue(j,'n')
                    elif j == 3:
                        self.datagrid.SetColLabelValue(j,'weight')
                    elif j == 4:
                        self.datagrid.SetColLabelValue(j,'plot type')
                if j < 3:
                    self.datagrid.SetReadOnly(i,j)
                elif j == 3:
                    self.datagrid.SetCellEditor(i,j,gridlib.GridCellFloatEditor())
                elif j == 4:
                    self.datagrid.SetCellEditor(i,j,gridlib.GridCellNumberEditor(-1,3))
        self.Bind(gridlib.EVT_GRID_CELL_CHANGE, self.OnCellChange)            

        #Fitting option use random start parameters
        self.getfitalgo = wx.ComboBox(self,-1, value="Levenberg-Marquardt",\
                                pos=(520, 35), size=(240,20),\
                                choices=["Levenberg-Marquardt","Downhill-Simplex"],style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.setfitalgo, self.getfitalgo)
        
        self.random_pars = wx.CheckBox(self, label = '  start Fit with random parameter values', pos = (520, 65))
        self.random_pars.SetValue(False)
        wx.EVT_CHECKBOX(self, self.random_pars.GetId(), self.setrandom_pars)

        #Fitting option auto_fit
        self.auto_fit_box = wx.CheckBox(self, label = '  auto fit,   number of params:', pos = (520, 90))
        self.auto_fit_box.SetValue(False)
        wx.EVT_CHECKBOX(self, self.auto_fit_box.GetId(), self.setauto_fit)

        self.auto_fit_n_box = wx.TextCtrl(self, pos=(700,90), size=(60,20))
        self.auto_fit_n_box.SetValue(str(self.auto_fit_n))
        self.Bind(wx.EVT_TEXT, self.setauto_fit_n, self.auto_fit_n_box)
        
        #### Downhill Simplex parameters and options #################################### 
        wx.StaticText(self, label = 'Simplex Parameters:  ', pos=(520, 130), size=(200, 20))

        wx.StaticText(self, label = 'alpha:    ', pos=(520, 157), size=(40, 20))
        self.alpha = wx.TextCtrl(self, pos=(570,155), size=(60,20))
        self.alpha.SetValue(str(self.simplex_params[0]))
        self.Bind(wx.EVT_TEXT, self.setalpha, self.alpha)

        wx.StaticText(self, label = 'beta:    ', pos=(640, 157), size=(50, 20))
        self.beta = wx.TextCtrl(self, pos=(700,155), size=(60,20))
        self.beta.SetValue(str(self.simplex_params[1]))
        self.Bind(wx.EVT_TEXT, self.setbeta, self.beta)

        wx.StaticText(self, label = 'gamma:    ', pos=(520, 182), size=(40, 20))
        self.gamma = wx.TextCtrl(self, pos=(570,180), size=(60,20))
        self.gamma.SetValue(str(self.simplex_params[2]))
        self.Bind(wx.EVT_TEXT, self.setgamma, self.gamma)

        wx.StaticText(self, label = 'delta:    ', pos=(640, 182), size=(50, 20))
        self.delta = wx.TextCtrl(self, pos=(700,180), size=(60,20))
        self.delta.SetValue(str(self.simplex_params[3]))
        self.Bind(wx.EVT_TEXT, self.setdelta, self.delta)

        wx.StaticText(self, label = 'Ftol:    ', pos=(520, 207), size=(40, 20))
        self.ftol = wx.TextCtrl(self, pos=(570,205), size=(60,20))
        self.ftol.SetValue(str(self.simplex_params[4]))
        self.Bind(wx.EVT_TEXT, self.setftol, self.ftol)

        wx.StaticText(self, label = 'maxiter:    ', pos=(640, 207), size=(50, 20))
        self.maxiter = wx.TextCtrl(self, pos=(700,205), size=(60,20))
        self.maxiter.SetValue(str(self.simplex_params[5]))
        self.Bind(wx.EVT_TEXT, self.setmaxiter, self.maxiter)
        
        ################### statistics input ###################################################
        wx.StaticText(self, label = 'Parameter Statistics:  ', pos=(520, 247), size=(200, 20))

        wx.StaticText(self, label = 'fract. param. change: ', pos=(520, 272), size=(110, 20))
        self.fpc_control = wx.TextCtrl(self, pos=(700,270), size=(60,20))
        self.fpc_control.SetValue(str(self.nb.fpc))
        self.Bind(wx.EVT_TEXT, self.setfpc, self.fpc_control)

        self.statisticsbutton = wx.Button(self, label = 'calculate parameter\nstatistics', pos =(640,295), size=(120,40))
        self.Bind(wx.EVT_BUTTON, self.OnClickStatistics, self.statisticsbutton)

        self.getparam = wx.ComboBox(self,-1, value=self.nb.used_params[-1],\
                                     pos=(520, 295), size=(110,20),\
                                     choices=self.nb.used_params,style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.setparam, self.getparam)

        #Data plotting
        wx.StaticText(self, label = 'Plotting Options (Fig. 1): ', pos=(520, 357), size=(200, 20))

        wx.StaticText(self, label = 'Plot dimensions', pos=(520, 382), size=(110, 20))
        self.getplotdims = wx.TextCtrl(self, pos=(700, 380), size=(60,20))
        self.getplotdims.SetValue(str(self.plotdims[0])+' '+str(self.plotdims[1]))
        self.Bind(wx.EVT_TEXT, self.setplotdims, self.getplotdims)

        self.getplotmode = wx.ComboBox(self,-1, value="Data+Model",\
                        pos=(520, 405), size=(110,20),\
                        choices=["Data+Model","Residuals","Monte-Carlo-range"],style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.setplotmode, self.getplotmode)

        self.button = wx.Button(self, label = 'Plot', pos =(640,405), size = (120,40))
        self.Bind(wx.EVT_BUTTON, self.OnClick, self.button)

        #Monte Carlo Perturbation
        wx.StaticText(self, label = 'Monte Carlo Parameter Perturbation: ', pos=(520, 467), size=(200, 20))
        wx.StaticText(self, label = 'number of samples', pos=(520, 492), size=(110, 20))
        self.getmcn = wx.TextCtrl(self, pos=(700, 490), size=(60,20))
        self.getmcn.SetValue(str(self.nb.mcn))
        self.Bind(wx.EVT_TEXT, self.setmcn, self.getmcn)

        self.getmcmode = wx.ComboBox(self,-1, value="Normal",\
                        pos=(520, 515), size=(110,20),\
                        choices=["Normal","Multivar.-Normal"],style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.setmcmode, self.getmcmode)
        self.mcbutton = wx.Button(self, label = 'Run Monte Carlo', pos =(640,515), size = (120,40))
        self.Bind(wx.EVT_BUTTON, self.OnClickMC, self.mcbutton)

        
        # Start and Stop Fit ###################################################################
        self.Startfitbutton = wx.Button(self, label = 'Start Fit', pos =(520,580), size=(170,55))
        self.Bind(wx.EVT_BUTTON, self.OnClickStartFit, self.Startfitbutton)
        self.Stopfitbutton = wx.Button(self, label = 'Stop Fit', pos =(700,580), size=(60,55))
        self.Bind(wx.EVT_BUTTON, self.OnClickStopFit, self.Stopfitbutton)
        
    #################### data grid functions ####################################################
    def fill_grid(self):
        for i in range(25):
            for j in range(5):
                self.datagrid.SetCellValue(i,j,'')
        for i in range(len(self.nb.datasets)):
            dat = self.nb.datasets[i]
            self.datagrid.SetCellValue(i,0,dat.x_label)
            self.datagrid.SetCellValue(i,1,dat.y_label)
            self.datagrid.SetCellValue(i,2,str(dat.n))
            self.datagrid.SetCellValue(i,3,str(dat.weight))
            if i < 9:
                self.datagrid.SetCellValue(i,4,str(0))
            else:
                self.datagrid.SetCellValue(i,4,str(-1))

    def OnCellChange(self, event):
        x = event.GetRow()
        y = event.GetCol()
        value = self.nb.MainControlPage.datagrid.GetCellValue(x,y)
        if y == 3:
            self.nb.datasets[x].weight = float(value)
        elif y == 4:
            self.nb.datasets[x].plot = int(value)
            self.nb.plot()
   ################################# Plotting options event functions #########################################
    def setplotdims(self, event):
        if event.GetString() == '':
            None
        else:
            dims = str.rsplit(str(event.GetString()))
            if len(dims) != 2:
                print 'need two integer numbers'
            else:
                self.plotdims = []
                for i in dims:
                    self.plotdims.append(int(i))

    def setplotmode(self, event):
        self.nb.plotmode = event.GetSelection()
        self.nb.plot()
                    
    def OnClick(self,e):
        self.nb.model()
        self.nb.plot()
        
    ################################# Fitting options event functions #########################################
    def setfitalgo(self,event):
        self.fitalgo = event.GetSelection()
        
    def setrandom_pars(self,e):
        self.simplex_params[6] = self.random_pars.GetValue()    

    def setauto_fit(self, e):
        self.auto_fit = self.auto_fit_box.GetValue()
        
    def setauto_fit_n(self, event):
        try:
            a = int(event.GetString())
            self.auto_fit_n = a
        except ValueError:
            pass
    ################################# Simplex options event functions #########################################
    def setalpha(self,event):
        try:
            a = float(event.GetString())
            if a <= 0 or a > 1:
                print 'alpha must be > 0 and <= 1'
            else:
                self.simplex_params[0] = a
        except ValueError:
            pass
    def setbeta(self,event):
        try:
            a = float(event.GetString())
            if a <= 0 or a > 1:
                print 'beta must be > 0 and <= 1'
            else:
                self.simplex_params[1] = a
        except ValueError:
            pass
    def setgamma(self,event):
        try:
            a = float(event.GetString())
            if a <= 0:
                print 'gamma must be > 0'
            else:
                self.simplex_params[2] = a
        except ValueError:
            pass
    def setdelta(self,event):
        try:
            a = float(event.GetString())
            if a <= 0 or a > 1:
                print 'delta must be > 0 and <= 1'
            else:
                self.simplex_params[3] = a
        except ValueError:
            pass
    def setftol(self,event):
        try:
            a = float(event.GetString())
            if a <= 0:
                print 'ftol must be > 0'
            else:
                self.simplex_params[4] = a
        except ValueError:
            pass
    def setmaxiter(self,event):
        try:
            a = int(event.GetString())
            if a <= 0:
                print 'maxiter must be > 0'
            else:
                self.simplex_params[5] = a
        except ValueError:
            pass
    ##############################################################################
    def setparam(self, event):
        param = event.GetSelection()
        sens = self.nb.sensitivities[param]
        parameter = self.nb.used_params[param]
        dp = self.nb.parameter[parameter][4]
        self.nb.plot_sens(dp, sens, parameter)
        pass
    
    def setfpc(self, event):
        try:
            a = float(event.GetString())
            if a <= 0 or a > 0.1:
                print 'fpc must be > 0 and <= 0.1'
            else:
                self.nb.fpc = a
        except ValueError:
            pass
    ################################################################################################################################
    def OnClickStatistics(self,e):
        self.nb.get_used_params()
        if len(self.nb.used_params) == 0:
            print "NO PARAMETERS SELECTED!"
            pass
        else:
            self.nb.frame.SetStatusText(' computing parameter statistics ', 0)
            self.nb.statistics()
        
            for i in range(len(self.nb.param_labels)):
                self.nb.ParameterPage.control8[i].SetValue(str(round(self.nb.parameter[self.nb.param_labels[i]][4], 12)))

            self.getparam.Clear()
            self.getparam.AppendItems(self.nb.used_params)
            self.getparam.SetSelection(0)

            if len(self.nb.correl_matrix) > 0:
                print '\n\n COVARIANCE MATRIX for: \n'
                print self.nb.used_params
                print '\n'
                print self.nb.cov_matrix
                print '\n\n CORRELATION MATRIX for: \n (numbers on the diagonal are std-deviations of the parameters)\n'
                print self.nb.used_params
                print '\n'
                print self.nb.correl_matrix
                b = len(self.nb.correl_matrix)
                print '\n\nParameter Correlations: \n'
                print '\nParameter Correlations between 0.95 and 1.00: \n'
                n = 0
                for i in range(b):
                    for j in range(b):
                        if j > i:
                            if (self.nb.correl_matrix[i][j] >= 0.95 and self.nb.correl_matrix[i][j] <= 1.0) or (self.nb.correl_matrix[i][j] <= -0.95 and self.nb.correl_matrix[i][j] >= -1.0):
                                n = n+1
                                print self.nb.used_params[i] + ' & ' + self.nb.used_params[j] + ': ' + str(round(self.nb.correl_matrix[i][j],5))
                if n == 0: print 'None \n'
                else: print '\n'
                n = 0
                print '\nParameter Correlations between 0.8 and 0.95: \n'
                for i in range(b):
                    for j in range(b):
                        if j > i:
                            if (self.nb.correl_matrix[i][j] >= 0.8 and self.nb.correl_matrix[i][j] < 0.95) or (self.nb.correl_matrix[i][j] <= -0.8 and self.nb.correl_matrix[i][j] > -0.95):
                                n = n+1
                                print self.nb.used_params[i] + ' & ' + self.nb.used_params[j] + ': ' + str(round(self.nb.correl_matrix[i][j],5))
                if n == 0: print 'None \n'
                else: print '\n'
                n = 0
                print '\nParameter Correlations between 0.5 and 0.8: \n'
                for i in range(b):
                    for j in range(b):
                        if j > i:
                            if (self.nb.correl_matrix[i][j] >= 0.5 and self.nb.correl_matrix[i][j] < 0.8) or (self.nb.correl_matrix[i][j] <= -0.5 and self.nb.correl_matrix[i][j] > -0.8):
                                n = n+1
                                print self.nb.used_params[i] + ' & ' + self.nb.used_params[j] + ': ' + str(round(self.nb.correl_matrix[i][j],5))
                if n == 0: print 'None \n'
                else: print '\n'
                print 'statistics calculation finished, chi**2 = '+str(round(self.nb.chi2,12))+'\n'
                print 'number of used variables = '+str(b)
            self.nb.frame.SetStatusText(' statistics calculation finished, chi**2 = '+str(round(self.nb.chi2,12)), 0)
            self.nb.plot()
            pass

    #########################################################################################################################
    def setmcn(self, event):
        try:
            self.nb.mcn = int(event.GetString())
        except ValueError:
            pass

    def setmcmode(self, event):
        self.nb.mcmode = event.GetSelection()
        
    def OnClickMC(self, event):
        flag = True
        self.nb.get_used_params()
        for param in self.nb.used_params:
            if param not in self.nb.batch_model:
                flag = False
                print 'Parameter "'+param+'" shell be perturbed by the Monte Carlo method,\n but is not used in the model!!'
            if self.nb.parameter[param][4] == 0:
                flag = False
                print 'Parameter "'+param+'" shell be perturbed by the Monte Carlo method,\n but has a Std-dev of zero!!'
        if flag:
            results = []
            means = []
            for param in self.nb.used_params:
                means.append(self.nb.parameter[param][0])
            if self.nb.mcmode == 0:
                self.nb.frame.SetStatusText('Performing Monte Carlo Parameter Perturbation',0)
                for i in range(self.nb.mcn):
                    while wx.GetApp().Pending():
                        wx.GetApp().Dispatch()
                        wx.GetApp().Yield(True)
                    self.nb.frame.SetStatusText(str(i),1)
                    for j in range(len(self.nb.used_params)):
                        self.nb.parameter[self.nb.used_params[j]][0] = Num.random.normal(means[j],self.nb.parameter[self.nb.used_params[j]][4])
                    model = self.nb.model_modified(self.nb.used_params[0], 0)
                    results.append(model)
                results = Num.transpose(results)
                for i in range(len(self.nb.data[0])):
                    self.nb.data[6][i] = Num.std(results[i])
                for i in range(len(self.nb.used_params)):
                    self.nb.parameter[self.nb.used_params[i]][0] = means[i]   
            elif self.nb.mcmode == 1:
                self.nb.frame.SetStatusText('Updating Parameter Statistics',0)
                self.nb.statistics()
                self.nb.frame.SetStatusText('Performing Monte Carlo Parameter Perturbation',0)
                for i in range(self.nb.mcn):
                    while wx.GetApp().Pending():
                        wx.GetApp().Dispatch()
                        wx.GetApp().Yield(True)
                    self.nb.frame.SetStatusText(str(i),1)
                    values = Num.random.multivariate_normal(means,self.nb.cov_matrix)
                    for j in range(len(self.nb.used_params)):
                        self.nb.parameter[self.nb.used_params[j]][0] = values[j]
                    model = self.nb.model_modified(self.nb.used_params[0], 0)
                    results.append(model)
                results = Num.transpose(results)
                for i in range(len(self.nb.data[0])):
                    self.nb.data[6][i] = Num.std(results[i])
                for i in range(len(self.nb.used_params)):
                    self.nb.parameter[self.nb.used_params[i]][0] = means[i]  

            self.nb.model()
            self.nb.frame.SetStatusText('Monte Carlo finished',0)
            self.nb.frame.SetStatusText('',1)
            self.nb.plotmode = 2
            self.getplotmode.SetStringSelection("Monte-Carlo-range")
            self.nb.plot()
        
    #########################################################################################################################
    def OnClickStartFit(self,e):
        flag = True
        self.nb.get_used_params()
        for param in self.nb.used_params:
            if param not in self.nb.batch_model:
                flag = False
                print 'Parameter "'+param+'" shell be adjusted in the fit,\n but is not used in the model!!'
        if flag == True and self.auto_fit == False:
            self.chi2 = -1
            self.StopFit = False
            if self.fitalgo == 0:
                LM_fit(self.nb.frame)
            elif self.fitalgo == 1:
                simplex(self.nb.frame)
            while wx.GetApp().Pending():
                wx.GetApp().Dispatch()
                wx.GetApp().Yield(True)
            for i in range(len(self.nb.param_labels)):
                if self.nb.parameter[self.nb.param_labels[i]][3] and self.nb.parameter[self.nb.param_labels[i]][5] == '':
                    self.nb.ParameterPage.control1[i].SetValue(str(round(self.nb.parameter[self.nb.param_labels[i]][0], 12)))
            self.nb.SetSelection(3)
        
        elif flag == True and self.auto_fit == True:
            print "STARTING AUTO FIT PROCEDURE"
            self.chi2 = -1
            self.StopFit = False
            param_list = []
            for key in self.nb.parameter.keys():
                if self.nb.parameter[key][3]:
                    param_list.append(key)
            print str(len(param_list))+" parameters are adjusted: " + str(param_list)
            if len(param_list) < self.auto_fit_n:
                self.auto_fit_n = len(param_list)
                self.auto_fit_n_box.SetValue(str(self.auto_fit_n))    
            j = 0
            while self.StopFit == False:
                j = j + 1
                parameter = random.sample(set(param_list), self.auto_fit_n)
                print "auto fit run " + str(j)
                print "parameter set : " + str(parameter)
                for key in self.nb.parameter.keys():
                    self.nb.parameter[key][3] = False
                for key in parameter:
                    self.nb.parameter[key][3] = True
                if self.fitalgo == 0:
                    LM_fit(self.nb.frame)
                elif self.fitalgo == 1:
                    simplex(self.nb.frame)
                while wx.GetApp().Pending():
                    wx.GetApp().Dispatch()
                    wx.GetApp().Yield(True)
                for i in range(len(self.nb.param_labels)):
                    if self.nb.parameter[self.nb.param_labels[i]][3] and self.nb.parameter[self.nb.param_labels[i]][5] == '':
                        self.nb.ParameterPage.control1[i].SetValue(str(round(self.nb.parameter[self.nb.param_labels[i]][0], 12)))
            
            for key in self.nb.parameter.keys():
                    self.nb.parameter[key][3] = False
            for key in param_list:
                self.nb.parameter[key][3] = True
            print "AUTO FIT STOPPED BY USER"
            self.nb.SetSelection(3)
        else:
            pass
        
    def OnClickStopFit(self,e):
        self.StopFit = True   
##########################################################################################################################
##########################################################################################################################
class ParameterPanel(wx.ScrolledWindow):
    def __init__(self,parent):
        wx.ScrolledWindow.__init__(self,parent)

        self.control0 = []
        self.control1 = []
        self.control2 = []
        self.control3 = []
        self.control4 = []
        self.control5 = []
        self.control6 = []
        self.control7 = []
        self.control8 = []
        self.togglesteps = []

        wx.StaticText(self, label = 'value', pos=(170, 10), size=(100, 20))
        wx.StaticText(self, label = 'std.-dev.', pos=(290, 10), size=(100, 20))
        wx.StaticText(self, label = 'min', pos=(390, 10), size=(70, 20))
        wx.StaticText(self, label = 'max', pos=(480, 10), size=(70, 20))
        wx.StaticText(self, label = 'refine', pos=(545, 10), size=(60, 15))
        wx.StaticText(self, label = ' - ', pos=(615, 10), size=(20, 20))
        wx.StaticText(self, label = ' step ', pos=(660, 10), size=(40, 20))
        wx.StaticText(self, label = ' + ', pos=(705, 10), size=(20, 20))

        self.uncheckbutton = wx.Button(self, label = 'u', pos =(549,25), size=(15,15))
        self.Bind(wx.EVT_BUTTON, self.uncheck, self.uncheckbutton)
        
        self.nb = self.GetParent()

    def uncheck(self,e):
        for i in range(len(self.control4)):
            self.nb.parameter[self.nb.param_labels[i]][3] = False
            self.control4[i].SetValue(False)

    def clicklabel(self, e):
        item = e.GetId()-7*len(self.nb.param_labels)-20000
        param = self.nb.param_labels[item]
        sens, dp = self.nb.single_param_sensitivity(param)
        self.nb.parameter[self.nb.param_labels[item]][4] = dp
        self.control8[item].SetValue(str(round(self.nb.parameter[self.nb.param_labels[item]][4],8)))
        self.nb.plot_sens(dp,sens,param)

    def editparvalue(self,event):
        item = event.GetId()-20000
        try:
            self.nb.parameter[self.nb.param_labels[item]][0] = float(event.GetString())
            self.nb.parameter[self.nb.param_labels[item]][5] = ''
            for key in self.nb.parameter.keys():
                if self.nb.parameter[key][5] == self.nb.param_labels[item]:
                    self.nb.parameter[key][0] = self.nb.parameter[self.nb.param_labels[item]][0]
        except ValueError:
            try:
                self.nb.parameter[self.nb.param_labels[item]][5] = event.GetString()
                self.nb.parameter[self.nb.param_labels[item]][3] = False
                self.control4[item].SetValue(False)
            except ValueError:
                pass
    def editparstddev(self,event):
        item = event.GetId()-8*len(self.nb.param_labels)-20000
        try:
            self.nb.parameter[self.nb.param_labels[item]][4] = float(event.GetString())
        except ValueError:
            pass
    def editparmin(self,event):
        item = event.GetId()-len(self.nb.param_labels)-20000
        try:
            self.nb.parameter[self.nb.param_labels[item]][1] = float(event.GetString())
        except ValueError:
            pass        
    def editparmax(self,event):
        item = event.GetId()-2*len(self.nb.param_labels)-20000
        try:
            self.nb.parameter[self.nb.param_labels[item]][2] = float(event.GetString())
        except ValueError:
            pass        
    def editparstate(self,event):
        item = event.GetId()- 3*len(self.nb.param_labels)-20000
        state = self.control4[item].GetValue()
        self.nb.parameter[self.nb.param_labels[item]][3] = state 
        if state and self.nb.parameter[self.nb.param_labels[item]][5] != '':
            self.nb.parameter[self.nb.param_labels[item]][5] = ''
            self.nb.ParameterPage.control1[item].SetValue(str(round(self.nb.parameter[self.nb.param_labels[item]][0], 12)))
    def toggleminus(self, event):
        item = event.GetId()-4*len(self.nb.param_labels)-20000
        step = self.togglesteps[item]
        self.nb.parameter[self.nb.param_labels[item]][0] = self.nb.parameter[self.nb.param_labels[item]][0] - step
        for key in self.nb.parameter.keys():
            if self.nb.parameter[key][5] == self.nb.param_labels[item]:
                self.nb.parameter[key][0] = self.nb.parameter[self.nb.param_labels[item]][0] 
        self.nb.parameter[self.nb.param_labels[item]][5] = ''    
        self.control1[item].SetValue(str(self.nb.parameter[self.nb.param_labels[item]][0]))
        self.nb.model()
        self.nb.plot()
    def togglestep(self,event):
        item = event.GetId()-5*len(self.nb.param_labels)-20000
        try:
            self.togglesteps[item] = float(event.GetString())
        except ValueError:
            pass
    def toggleplus(self, event):
        item = event.GetId()-6*len(self.nb.param_labels)-20000
        step = self.togglesteps[item]
        self.nb.parameter[self.nb.param_labels[item]][0] = self.nb.parameter[self.nb.param_labels[item]][0] + step
        for key in self.nb.parameter.keys():
            if self.nb.parameter[key][5] == self.nb.param_labels[item]:
                self.nb.parameter[key][0] = self.nb.parameter[self.nb.param_labels[item]][0]
        self.nb.parameter[self.nb.param_labels[item]][5] = '' 
        self.control1[item].SetValue(str(self.nb.parameter[self.nb.param_labels[item]][0]))
        self.nb.model()
        self.nb.plot()
############################################################################################################
############################################################################################################
class ModelPanel(wx.ScrolledWindow):
    def __init__(self,parent):
        wx.ScrolledWindow.__init__(self,parent)
        
        self.control = wx.TextCtrl(self, style=wx.TE_MULTILINE)
        self.Bind(wx.EVT_TEXT,self.editmodel, self.control)

        self.calcbutton = wx.Button(self, label = 'calculate', pos =(10,10), size=(120,25))
        self.Bind(wx.EVT_BUTTON, self.OnCalc, self.calcbutton)
        
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.calcbutton, 0, wx.SHAPED)
        self.sizer.Add(self.control, 1, wx.EXPAND)
        self.SetSizer(self.sizer)
        self.SetAutoLayout(1)
        self.sizer.Fit(self)
        self.nb = self.GetParent()

    def editmodel(self, event):
        self.nb.batch_model = event.GetString()

    def OnCalc(self, event):
        self.nb.model()
        self.nb.plot()
        
############################################################################################################
############################################################################################################
class DataPanel(wx.ScrolledWindow):
    def __init__(self,parent):
        
        wx.ScrolledWindow.__init__(self,parent)
        self.control = wx.TextCtrl(self, style=wx.TE_MULTILINE)
        self.Bind(wx.EVT_TEXT,self.editdata, self.control)
        
        self.refreshbutton = wx.Button(self, label = 'update datasets', pos =(10,10), size=(120,25))
        self.Bind(wx.EVT_BUTTON, self.OnRefresh, self.refreshbutton)
        
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.refreshbutton, 0, wx.SHAPED)
        self.sizer.Add(self.control, 1, wx.EXPAND)
        self.SetSizer(self.sizer)
        self.SetAutoLayout(1)
        self.sizer.Fit(self)
        
        self.nb = self.GetParent()

    def editdata(self, event):
        try:
            self.nb.datastring = event.GetString() 
        except:
            pass
        
    def OnRefresh(self, event):
        self.nb.data, self.nb.datasets = eval_datastring(self.nb.datastring)
        self.nb.MainControlPage.fill_grid()
        pass
    
############################################################################################################
############################################################################################################
app = wx.App(False)
frame = wxPPPRFrame(parent = None, title = "Python PhreeqC Parameter Refinement", size = (800,750))
frame.Show(True)
app.MainLoop()


