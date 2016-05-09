"""
Functions and classes for plotting

Authors/modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""

################################################################################
import wx
import os
import numpy as Num
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

################################################################################
class wxPlotFrame(wx.Frame):
    def __init__(self, parent, title, size):
        wx.Frame.__init__(self, parent, title= title, size=size)
        
        self.menubar = wx.MenuBar()
        self.file = wx.Menu()
        menusavefig = self.file.Append(wx.ID_ANY, '&Save figure',\
                                       "Save the current figure")
        menuexit = self.file.Append(wx.ID_EXIT, \
                                    '&Exit', u"Close Fourier Frame")
        self.menubar.Append(self.file, '&File')
        self.SetMenuBar(self.menubar)
        
        self.figure = Figure()
        self.canvas = FigureCanvasWxAgg( self, -1, self.figure )
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas,1,wx.EXPAND)
        self.SetSizer(self.sizer)
        
        self.Bind(wx.EVT_MENU, self._OnSave,menusavefig)
        self.Bind(wx.EVT_MENU, self._OnClose,menuexit)
        self.Bind(wx.EVT_CLOSE, self._OnClose)

        self.toolbar = NavigationToolbar2Wx(self.canvas)
        self.toolbar.Realize()
        self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.toolbar.update()
        
        
    def _OnSave(self, event):
        dirname = ''
        dlg = wx.FileDialog(self, "Save the current Figure to a file",\
                            dirname, ".png", "*.png", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            self.figure.savefig(filename)
            
        dlg.Destroy()

    def _OnClose(self, event):
        parent = self.GetParent()
        parent.FigureFrame = None
        self.Destroy()
#######################################################################################
def createplotframe(parent, title, size):
    frame = wxPlotFrame(parent,title, size)
    return frame
