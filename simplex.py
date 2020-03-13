"""
Nelder Mead Simplex- and Levenberg-Marquardt optimization routines used in pppr

Authors/modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""
###############################################################################

import numpy as Num
from scipy.optimize import leastsq
import random
import wx
import time

############################### methods used by simplex ############################################################################################
def insert(used_params, point, parameter):
    for i in range(len(used_params)):
        key = used_params[i]
        parameter[key][0] = point[i]
    return parameter

def check_limits(used_params, point, parameter):
    for i in range(len(used_params)):
        key = used_params[i]
        if point[i] < parameter[key][1]: 
            point[i] = parameter[key][1]
        elif point[i] > parameter[key][2]: 
            point[i] = parameter[key][2]
    return point
    
def calc_average(points):
    av_point = Num.zeros((len(points[0])),float)
    for i in points:
        av_point = av_point + i
    av_point = av_point/(len(points))
    return av_point

def min_max(function_values):
    Ymin = function_values.min()
    mini = int(Num.where(function_values == Ymin)[0][0])
    Ymax = function_values.max()
    maxi = int(Num.where(function_values == Ymax)[0][0])
    return mini, maxi

def contraction(Xmax, Xav, beta, used_params, nb):
    #print 'contraction'
    Xcon = beta*Xmax+(1-beta)*Xav
    nb.parameter = insert(used_params, Xcon, nb.parameter)
    nb.model()
    Ycon = nb.chi2
    return Xcon,Ycon

def reflection(Xmax, Xav, alpha, used_params,nb):
    #print 'reflection'
    Xref = (1+alpha)*Xav - alpha*Xmax
    Xref = check_limits(used_params, Xref, nb.parameter)
    nb.parameter = insert(used_params, Xref, nb.parameter)
    nb.model()
    Yref = nb.chi2
    return Xref, Yref

def expansion(Xref, Xav, gamma, used_params, nb):
    Xexp = (1+gamma)*Xref - gamma*Xav
    Xexp = check_limits(used_params, Xexp, nb.parameter)
    nb.parameter = insert(used_params, Xexp, nb.parameter)
    nb.model()
    Yexp = nb.chi2
    return Xexp, Yexp

def compression(X, mini):
    #print 'compression'
    Y = Num.ndarray((0,len(X[0])),float)
    for x in X:
        x = (x + X[mini])/2
        Y = Num.append(Y,[x],axis = 0)
    return Y

def parameter_plot(fig,used_params,parameter,points,mini):
    fig.clear()
    fig.suptitle('Parameter Plot', fontsize = 20)
    plot = fig.add_subplot(111)
    plot.set_xticks(range(len(used_params)))
    plot.set_xticklabels(used_params,rotation = 90)
    low = []
    spread = []
    for i in used_params:
        low.append(parameter[i][1])
        spread.append(parameter[i][2]-parameter[i][1])
    for i in range(len(points)):
        if i == mini:
            pass
        else:
            plot.plot(range(len(used_params)),(points[i]-low)/spread,'bo')
        plot.plot(range(len(used_params)),(points[mini]-low)/spread,'ro')
    plot.set_xlim(-1,len(used_params)+1)
    plot.vlines(range(len(used_params)),ymin = 0,ymax = 1, color = 'k', linestyles = 'dashed')
    plot.set_ylim(0,1)
    fig.canvas.draw()

def calc_ftol(function_values):
    av = Num.sum(function_values)/len(function_values)
    sigma = 0
    for i in function_values:
        sigma = sigma + (i-av)**2
    sigma = Num.sqrt(sigma/len(function_values))
    return sigma
#################################Simplex main routine###################################################################################################
def simplex(frame):
    panel = frame.nb.MainControlPage
    alpha, beta, gamma, delta, ftol, maxiter, random_pars = panel.simplex_params
    frame.SetStatusText('Preparing Simplex',0)
    frame.nb.get_used_params()
    used_params = frame.nb.used_params
    used_params_values = []
    for key in used_params:
        used_params_values.append(frame.nb.parameter[key][0])        

    function_values = Num.ndarray((len(used_params)+1),float)
    points = Num.ndarray((len(used_params)+1,len(used_params)),float)
    for i in range(len(used_params)+1):
        wx.Yield()
        if i == 0 and not random_pars:
            points[i] = used_params_values
            frame.nb.parameter = insert(used_params, points[i], frame.nb.parameter)
            frame.nb.model()
            function_values[i] = frame.nb.chi2
        else:
            for j in range(len(points[i])):
                key = used_params[j]
                points[i][j] = used_params_values[j] + random.uniform(((frame.nb.parameter[key][1]-used_params_values[j])*delta), \
                                                                      ((frame.nb.parameter[key][2]-used_params_values[j])*delta))
            frame.nb.parameter = insert(used_params, points[i], frame.nb.parameter)
            frame.nb.model()
            function_values[i] = frame.nb.chi2
            
    not_converged = True
    z = 0
    mini, maxi = min_max(function_values)
    statusstring ='iteration '+str(z)+', best chi**2 = '+str(round(function_values[mini],8))
    frame.SetStatusText(statusstring,0)
    #if panel.FigureFrame == None:
    #    panel.FigureFrame = createplotframe(panel, "PPPR Parameter plot", (1000,500))
    #    panel.FigureFrame.Show(True)
    #    panel.Figure2 = panel.FigureFrame.figure
    #parameter_plot(panel.Figure2, used_params, frame.nb.parameter, points, mini)
    old_mini = function_values[mini]
    while not_converged:
        frame.SetStatusText(str(z),1)
        Xav = calc_average(points)
        wx.Yield()
        if panel.StopFit:
            frame.SetStatusText('Fit stopped after '+str(z)+' iterations',0)
            not_converged = False
            print('Fit aborted by user after '+str(z)+' iterations')
        Xref, Yref = reflection(points[maxi], Xav, alpha, used_params, frame.nb)
        if Yref < function_values[mini]:
            Xexp, Yexp = expansion(Xref, Xav, gamma, used_params, frame.nb)
            if Yexp < function_values[mini]:
                #print 'expansion'
                points[maxi] = Xexp
                function_values[maxi] = Yexp
            else:
                points[maxi] = Xref
                function_values[maxi] = Yref
        else:
            test = False
            for i in range(len(points)):
                if Yref < function_values[i]:
                    if i == maxi:
                        test = False
                    else:
                        test = True
            if test:
                points[maxi] = Xref
                function_values[maxi] = Yref
            else:
                if Yref < function_values[maxi]:
                    Xcon,Ycon = contraction(Xref, Xav, beta, used_params, frame.nb)
                else:
                    Xcon,Ycon = contraction(points[maxi], Xav, beta, used_params, frame.nb)

                if Ycon < function_values[maxi]:
                    points[maxi] = Xcon
                    function_values[maxi] = Ycon
                else:
                    points = compression(points, mini)
                    for i in range(len(points)):
                        frame.nb.parameter = insert(used_params, points[i], frame.nb.parameter)
                        frame.nb.model()
                        function_values[i] = frame.nb.chi2
        mini, maxi = min_max(function_values)
        act_ftol = calc_ftol(function_values)
        if function_values[mini]<old_mini:
            statusstring ='iteration '+str(z)+', best chi**2 = '+str(round(function_values[mini],8))+' ftol = '+str(round(act_ftol,8))
            frame.SetStatusText(statusstring,0)
            frame.nb.plot()
            #parameter_plot(panel.Figure2, used_params, frame.nb.parameter, points, mini)
            old_mini = function_values[mini]
            
        if act_ftol < ftol:
            not_converged = False
            print('\n CONVERGENCE REACHED DUE TO FTOL \n')
        if z >= maxiter:
            not_converged = False
            print('\n NO CONVERGENCE, STOP DUE TO MAXITER \n')
        z = z+1
    if not panel.StopFit: print(' Downhill Simplex stopped after '+str(z-1)+' iterations')
    print('best fit chi**2 = '+str(round(function_values[mini],8))+'\n')
    frame.SetStatusText('End of Downhill Simplex, best chi**2: '+str(round(function_values[mini],8)),0)
    frame.SetStatusText('',1)
    param_best = points[mini]
    frame.nb.parameter = insert(used_params, param_best, frame.nb.parameter)
    frame.nb.model()
    frame.nb.plot()
    return

                              
def LM_fit(frame):
    panel = frame.nb.MainControlPage
    frame.nb.get_used_params()
    used_params = frame.nb.used_params
    frame.SetStatusText('Starting Levenberg-Marquardt fit',0)
    vector = Num.array([])
    for key in used_params:
        vector = Num.append(vector, frame.nb.parameter[key][0])
        
    def target(vector, used_params, insert, nb):
        nb.parameter = insert(used_params, vector, nb.parameter)
        model = nb.model_modified(used_params[0],0)
        return nb.data[1] - model

    def Jacobi(vector, used_params, insert, nb):
        J = Num.ndarray((len(used_params),len(nb.data[0])), float)
        nb.parameter = insert(used_params, vector, nb.parameter)
        for i in range(len(used_params)):
            h = nb.parameter[used_params[i]][0] * nb.fpc
            if h == 0:
                h = nb.fpc
            if nb.parameter[used_params[i]][0] -0.5*h >= nb.parameter[used_params[i]][1] \
                and nb.parameter[used_params[i]][0] +0.5*h <= nb.parameter[used_params[i]][2]:
            
                y1 = nb.model_modified(used_params[i],-0.5*h)
                y2 = nb.model_modified(used_params[i], 0.5*h)
                J[i,:] = (y1 - y2)/h
            
            elif nb.parameter[used_params[i]][0] -0.5*h < nb.parameter[used_params[i]][1]:
                
                y1 = nb.model_modified(used_params[i],0)
                y2 = nb.model_modified(used_params[i],h)
                J[i,:] = (y1 - y2)/h
            
            elif nb.parameter[used_params[i]][0] +0.5*h > nb.parameter[used_params[i]][2]:
                
                y1 = nb.model_modified(used_params[i],-h)
                y2 = nb.model_modified(used_params[i],0)
                J[i,:] = (y1 - y2)/h
                
        return J

    result = leastsq(target, vector, args = (used_params, insert, frame.nb), Dfun = Jacobi, col_deriv = 1, full_output = 1, ftol = 1e-12)
    wx.Yield()
    print(result)
    frame.nb.parameter = insert(used_params, result[0], frame.nb.parameter)
    frame.nb.model()
    frame.SetStatusText('fitting finished, best chi**2 = '+str(round(frame.nb.chi2,12)),0)
    frame.nb.plot()
    return
                              
            
    
    
