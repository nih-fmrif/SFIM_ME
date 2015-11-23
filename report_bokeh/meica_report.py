#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Javier Gonzalez Castillo 9/15/2015

"""
__version__="0.3"
import sys
import os

path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../melib'))
print("++ INFO [Main]: Using meBasic library located in: %s" % path)
if not path in sys.path:
    sys.path.insert(1, path)
del path
import meBasics as meb
import seaborn as sns
import numpy as np
import pandas as pd
import argparse
from scipy.stats import rankdata
from bokeh.models import ColumnDataSource, HoverTool, CustomJS, TapTool, Plot, Range1d
from bokeh.models.widgets import DataTable, TableColumn, NumberFormatter, BooleanFormatter, CheckboxEditor
from bokeh.models.widgets import Tabs, Panel, Select
from bokeh.document import Document
from bokeh.plotting import figure, gridplot, output_file, show, hplot, vplot
from bokeh.charts import BoxPlot, Histogram
from scipy.stats import scoreatpercentile
from glob import glob
from bokeh.models.glyphs import ImageURL

# READ INPUT PARAMETERS
# =====================
parser = argparse.ArgumentParser('Options')
parser.add_argument('-meicaDir',   dest = 'dir'  , help = 'directory where the output of meica was saved',default = None)
parser.add_argument('-runID',      dest = 'runID', help = 'data prefix',default = None)
parser.add_argument('-picDir',dest = 'pdir' ,      help = 'location of ICA maps (as photos)', default=None)
parser.add_argument('-nBins',dest = 'Nbins', help = 'Number of bins for the histograms. Default value = 200', default=200, type=int) 
options = parser.parse_args()

Nbins= int(options.Nbins)

# Check that all necessary options were provided
if options.pdir is None:
    print("++ ERROR: No path for ICA maps png files provided. Program will exit.")
    sys.exit() 
if options.dir is None:
    print("++ ERROR: No working directory provided. Program will exit.")
    sys.exit()
if options.runID is None:
    print("++ ERROR: No run ID provided. Program will exit.")
    sys.exit()
# Check if working directory exits
if not os.path.isdir(options.pdir):
    print("++ ERROR: ICA map directory [%s] does not exist." % options.pdir)
    print("++ ERROR: Exiting program now.")
    sys.exit()
    
if not os.path.isdir(options.dir):
    print("++ ERROR: Working directory [%s] does not exist." % options.dir)
    print("++ ERROR: Exiting program now.")
    sys.exit()

# If directory exits, check all necessary files are available.
if not os.path.exists(options.dir+'comp_table.txt'):
    print("++ ERROR: Component Table [%s] does not exists." % (options.dir+'comp_table.txt'))
    print("++ ERROR: Exiting program now.")
    sys.exit()
if not os.path.exists(options.dir+'accepted.txt'):
    print("++ ERROR: List of Accepted Components [%s] does not exists." % (options.dir+'accepted.txt'))
    print("++ ERROR: Exiting program now.")
    sys.exit()
if not os.path.exists(options.dir+'meica_mix.1D'):
    print("++ ERROR: Component Timeseries File [%s] does not exists." % (options.dir+'meica_mix.1D'))
    print("++ ERROR: Exiting program now.")
    sys.exit()
if not os.path.exists(options.dir+options.runID+'.chComp.FR2.nii'):
    print("++ ERROR: Per Component F_R2 Maps [%s] missing." % (options.dir+options.runID+'.chComp.FR2.nii'))
    print("++ ERROR: Exiting program now.")
    sys.exit()
if not os.path.exists(options.dir+options.runID+'.chComp.FS0.nii'):
    print("++ ERROR: Per Component F_S0 Maps [%s] missing." % (options.dir+options.runID+'.chComp.FS0.nii'))
    print("++ ERROR: Exiting program now.")
    sys.exit()
if not os.path.exists(options.dir+options.runID+'.chComp.Kappa_mask.nii'):
    print("++ ERROR: Per Component Kappa Masks [%s] missing." % (options.dir+options.runID+'.chComp.Kappa_mask.nii'))
    print("++ ERROR: Exiting program now.")
    sys.exit()


meicaDir    = os.path.abspath(options.dir)
ICAmap_dir  = os.path.abspath(options.pdir)
Program_dir = os.path.dirname(__file__)
print("++ INFO [Main]: ME-ICA Output Directory is %s" % meicaDir)
print("++ INFO [Main]: Path to PNG ICA files is %s" % ICAmap_dir)
print("++ INFO [Main]: Number of bins for histograms %d" % Nbins)
# Load Inputs into memory
# =======================

accepted_components = np.loadtxt(os.path.join(meicaDir,'accepted.txt')).astype('int')
comp_table          = np.loadtxt(os.path.join(meicaDir,'comp_table.txt'))
comp_timeseries     = np.loadtxt(os.path.join(meicaDir,'meica_mix.1D'))
comp_ffts           = np.loadtxt(os.path.join(meicaDir,'meica_mix_fft.1D'))
freq_axis           = np.loadtxt(os.path.join(meicaDir,'meica_mix_fft_freqs.1D'))

FR2_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.chComp.FR2.nii')
print("++ INFO [Main]: FR2 Maps [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.chComp.FR2.nii',FR2_maps.shape))
cR2_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.chComp.cR2.nii')
print("++ INFO [Main]: FR2 Maps [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.chComp.cR2.nii',cR2_maps.shape))
FS0_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.chComp.FS0.nii')
print("++ INFO [Main]: FS0 Maps [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.chComp.FS0.nii',FS0_maps.shape))
cS0_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.chComp.cS0.nii')
print("++ INFO [Main]: FS0 Maps [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.chComp.cS0.nii',cS0_maps.shape))
KM_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.chComp.Kappa_mask.nii'); KM_maps=(KM_maps>0);
print("++ INFO [Main]: Kappa Masks [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.chComp.Kappa_mask.nii',KM_maps.shape))
RM_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.chComp.Rho_mask.nii'); RM_maps=(RM_maps>0);
print("++ INFO [Main]: Rho Masks [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.chComp.Rho_mask.nii',RM_maps.shape))
Kappa_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.chComp.Kappa.nii')
print("++ INFO [Main]: Kappa Maps [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.chComp.Kappa.nii',Kappa_maps.shape))
Rho_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.chComp.Rho.nii')
print("++ INFO [Main]: Rho Maps [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.chComp.Rho.nii',Rho_maps.shape))
ICA_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.ICA.Zmaps.nii')
print("++ INFO [Main]: ICA Maps [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.ICA.Zmaps.nii',ICA_maps.shape))
ICAM_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.ICA.Zmaps.mask.nii'); ICAM_maps=(ICAM_maps>0);
print("++ INFO [Main]: ICA Masks [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.ICA.Zmaps.mask.nii',ICAM_maps.shape))
mask_orig,_,_ = meb.niiLoad(options.dir+options.runID+'.mask.orig.nii'); mask_orig = (mask_orig>0);
print("++ INFO [Main]: Original Mask [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.mask.orig.nii',mask_orig.shape))
Wgth_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.chComp.weightMaps.nii')
print("++ INFO [Main]: Weight Maps [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.chComp.weightMaps.nii',Wgth_maps.shape))

Nt, Nc              = comp_timeseries.shape
fica_psel           = np.zeros((Nc,))
fica_psel[accepted_components] = 1

ICAmap_default_path = os.path.join(Program_dir,'AxSag_Component_XX.png')
print(ICAmap_default_path)
ICAmap_paths        = glob(str(os.path.join(ICAmap_dir,'AxSag*.png')))

print( "++ INFO: Number of ICA Maps: %d" % len(ICAmap_paths))

# Setting up output file
# ======================
output_file(options.runID+".html", title="Static Bokeh Report for " +options.runID)

# Reading the different features of interest
# ==========================================
Nc,Nf  = comp_table.shape
kappa  = comp_table[:,1]
rho    = comp_table[:,2]
var    = comp_table[:,3]
ratio  = comp_table[:,6]
cID    = comp_table[:,0]
maxFR2 = comp_table[:,4]
maxFS0 = comp_table[:,5]
maxZICA= comp_table[:,7]
NvZmask= comp_table[:,8]
NvFR2mask = comp_table[:,9]
NvFS0mask = comp_table[:,10]
NvKapmask = comp_table[:,11]
NvRhomask = comp_table[:,12]

kappa_ranked = np.sort(kappa)[::-1]
rho_ranked   = np.sort(rho)[::-1]
var_ranked   = np.sort(var)[::-1]
ratio_ranked = np.sort(ratio)[::-1]

loc_by_kappa = Nc - rankdata(kappa)
loc_by_rho   = Nc - rankdata(rho)
loc_by_var   = Nc - rankdata(var)
loc_by_ratio = Nc - rankdata(ratio)

component_colormap      = { "1.0" : "#00ff00", "0.0" : "#ff0000"}
component_color         = [component_colormap[str(x)] for x in fica_psel]
component_status_labels = { "1.0" : "Accepted", "0.0" : "Rejected"}
component_status        = [component_status_labels[str(x)] for x in fica_psel]

Source = ColumnDataSource(data=dict(cID = cID, 
                                    kappa = kappa, loc_by_kappa = loc_by_kappa, 
                                    rho = rho, loc_by_rho = loc_by_rho, 
                                    var = var, loc_by_var = loc_by_var, 
                                    ratio = ratio, loc_by_ratio = loc_by_ratio,
                                    comp_color  = component_color,
                                    comp_status = component_status,
                                    maxFR2      = maxFR2,
                                    maxFS0      = maxFS0,
                                    maxZICA     = maxZICA,
                                    NvZmask     = NvZmask,
                                    NvFR2mask     = NvFR2mask,
                                    NvFS0mask     = NvFS0mask,
                                    NvKapmask     = NvKapmask,
                                    NvRhomask     = NvRhomask))

# ==============================================================================
#                                 FEATURE TABLE

comp_table_columns = [
    TableColumn(field="cID",  title="ID", formatter=NumberFormatter(format='0o')),
    TableColumn(field="comp_status",title="Status", editor=CheckboxEditor()),
    TableColumn(field="var", title="Variance", formatter=NumberFormatter(format='0.00')),
    TableColumn(field="kappa",title="Kappa", formatter=NumberFormatter(format='0.00')),
    TableColumn(field="rho",  title="Rho", formatter=NumberFormatter(format='0.00')),
    TableColumn(field="ratio",title="Ka/Rh", formatter=NumberFormatter(format='0.00')),
    TableColumn(field="maxFR2",title="maxFR2", formatter=NumberFormatter(format='0.000')),
    TableColumn(field="maxFS0",title="maxFS0", formatter=NumberFormatter(format='0.000')),
    TableColumn(field="maxZICA",title="maxZICA", formatter=NumberFormatter(format='0.000')),
    TableColumn(field="NvZmask",title="Nv(Z mask)"),
    TableColumn(field="NvFR2mask",title="Nv(FR2 mask)"),
    TableColumn(field="NvFS0mask",title="Nv(FS0 mask)"),
    TableColumn(field="NvKapmask",title="Nv(Kappa mask)"),
    TableColumn(field="NvRhomask",title="Nv(Rho mask)")
]
comp_table_DTABLE = DataTable(source=Source,columns=comp_table_columns,width=1350, height=250, editable=True, selectable=True, sortable=False)

# ==============================================================================
#                                 FEATURE PLOTS
# Feaute Plots Tools
# ==================
TOOLS = "tap,box_zoom,reset"
HoverKappa = HoverTool(tooltips=[("Component", "@cID"),("Kappa",     "@kappa"),("Rho",       "@rho"),
                       ("Variance",  "@var"),("Ratio",     "@ratio"),("Status", "$color[swatch]:comp_color")])
HoverRho   = HoverTool(tooltips=[("Component", "@cID"),("Kappa",     "@kappa"),("Rho",       "@rho"),
                       ("Variance",  "@var"),("Ratio",     "@ratio"),("Status", "$color[swatch]:comp_color")])
HoverVar   = HoverTool(tooltips=[("Component", "@cID"),("Kappa",     "@kappa"),("Rho",       "@rho"),
                       ("Variance",  "@var"),("Ratio",     "@ratio"),("Status", "$color[swatch]:comp_color")])
HoverRatio = HoverTool(tooltips=[("Component", "@cID"),("Kappa",     "@kappa"),("Rho",       "@rho"),
                       ("Variance",  "@var"),("Ratio",     "@ratio"),("Status", "$color[swatch]:comp_color")])
HoverKvsR  = HoverTool(tooltips=[("Component", "@cID"),("Kappa",     "@kappa"),("Rho",       "@rho"),
                       ("Variance",  "@var"),("Ratio",     "@ratio"),("Status", "$color[swatch]:comp_color")])
                                      
# Feature Plots
# =============
sp_kappa = figure(tools=[TOOLS, HoverKappa],width=325, height=250, y_axis_label='Kappa', toolbar_location='left')
sp_kappa.circle('loc_by_kappa','kappa',size=5,color='comp_color',source=Source)
sp_kappa.yaxis.axis_label_text_font_size = "12pt"
sp_tab_kappa  = Panel(child=sp_kappa, title='Sorted by Kappa')
sp_tabs_kappa = Tabs(tabs=[sp_tab_kappa])
 
sp_rho   = figure(tools=[TOOLS, HoverRho],width=325, height=250, y_axis_label='Rho', toolbar_location=None)
sp_rho.circle('loc_by_rho','rho',size=5,color='comp_color',source=Source)
sp_rho.yaxis.axis_label_text_font_size = "12pt"
sp_tab_rho  = Panel(child=sp_rho, title='Sorted by Rho')
sp_tabs_rho = Tabs(tabs=[sp_tab_rho])

sp_var   = figure(tools=[TOOLS,HoverVar],width=325, height=250, y_axis_label='Variance', toolbar_location=None)
sp_var.circle('loc_by_var','var',size=5,color='comp_color',source=Source)
sp_var.yaxis.axis_label_text_font_size = "12pt"
sp_tab_var  = Panel(child=sp_var, title='Sorted by Variance')
sp_tabs_var = Tabs(tabs=[sp_tab_var])

sp_ratio = figure(tools=[TOOLS,HoverRatio],width=325, height=250, y_axis_label='K/R Ratio', toolbar_location=None)
sp_ratio.circle('loc_by_ratio','ratio',size=5,color='comp_color',source=Source)
sp_ratio.yaxis.axis_label_text_font_size = "12pt"

sp_kvr = figure(tools=[TOOLS,HoverKvsR],width=325, height=250, y_axis_label='rho', x_axis_label='kappa', toolbar_location=None)
sp_kvr.circle('kappa','rho',size=5,color='comp_color',source=Source)
sp_kvr.xaxis.axis_label_text_font_size = "12pt"
sp_kvr.yaxis.axis_label_text_font_size = "12pt"

sp_tab_ratio = Panel(child=sp_ratio, title='Kappa/Rho Ratio')
sp_tab_kvr   = Panel(child=sp_kvr, title='Kappa vs. Rho')
sp_tabs_left = Tabs(tabs=[sp_tab_ratio,sp_tab_kvr])

# ==============================================================================

# ==============================================================================
#                          TIME SERIES PLOTS

# Load Default Data on Plots
# ==========================
default_ts_x  = range(Nt)
default_ts_y  = np.zeros((Nt,))
default_fft_x = freq_axis
default_fft_y = np.zeros((Nt,))

# Generate Plots
# ==============
sp_ts = figure(tools=[],toolbar_location=None, width=680,height=200, x_axis_label='Time [TR]', 
               title='Component Timeseries',x_range=(0,Nt), title_text_font_size='12pt')
sp_ts.xaxis.axis_label_text_font_size = "12pt"

sp_fft = figure(tools=[],toolbar_location=None, width=680,height=200, x_axis_label='Frequency', 
                title='Component Spectrum',x_range=(min(default_fft_x), max(default_fft_x)), title_text_font_size='12pt')
sp_fft.xaxis.axis_label_text_font_size = "12pt"

# Generate Data Sources for interactivity
# =======================================
timeseries_to_display = ColumnDataSource(data=dict(cID=cID, x=default_ts_x, y=default_ts_y))
available_timeseries  = ColumnDataSource(data=dict(cID=cID, x=default_ts_x, y=comp_timeseries.T,comp_color  = component_color))
sp_ts.line('x','y', source=timeseries_to_display, line_width=3)

ffts_to_display = ColumnDataSource(data=dict(cID=cID, x=default_fft_x, y=default_fft_y))
available_ffts  = ColumnDataSource(data=dict(cID=cID, x=default_fft_x, y=comp_ffts.T,comp_color  = component_color))
sp_fft.line('x','y',source=ffts_to_display, line_width=3)
# ==============================================================================

# ==============================================================================
#                    BRAIN MAP PLOTS
# Convert Input maps into usable mosaics
# ======================================
available_ICAmaps = ColumnDataSource(data=dict(cID=cID,urls=ICAmap_paths))
ICAmap_to_display = ColumnDataSource(data=dict(x=[0], y=[0], w=[1246], h=[86], url=[ICAmap_default_path]))
xdr               = Range1d(start=-630, end=630)
ydr               = Range1d(start=-45, end=40)
ICAmapFigure      = figure(tools=[],title="ICA maps", x_range=xdr, y_range=ydr,width=1350, height=400, x_axis_type=None, y_axis_type=None, toolbar_location=None, title_text_font_size='12pt')
ICAmapImg         = ImageURL(url="url", x="x", y="y", w="w", h="h", anchor="center")
ICAmapFigure.add_glyph(ICAmap_to_display,ICAmapImg)
ICAmapFigure.outline_line_color='#ffffff'

map_tab_ICA = Panel(child=ICAmapFigure, title='ICA Z-maps')
map_tabs = Tabs(tabs=[map_tab_ICA])
# ==============================================================================

# ==============================================================================
#                     TABS FOR HISTOGRAMS AND BOX PLOTS

FR2_hist=np.zeros((Nc,Nbins))
FR2_Ledges=np.zeros((Nc,Nbins))
FR2_Redges=np.zeros((Nc,Nbins))
FS0_hist=np.zeros((Nc,Nbins))
FS0_Ledges=np.zeros((Nc,Nbins))
FS0_Redges=np.zeros((Nc,Nbins))

cR2_hist=np.zeros((Nc,Nbins))
cR2_Ledges=np.zeros((Nc,Nbins))
cR2_Redges=np.zeros((Nc,Nbins))
cS0_hist=np.zeros((Nc,Nbins))
cS0_Ledges=np.zeros((Nc,Nbins))
cS0_Redges=np.zeros((Nc,Nbins))

Kappa_hist=np.zeros((Nc,Nbins))
Kappa_Ledges=np.zeros((Nc,Nbins))
Kappa_Redges=np.zeros((Nc,Nbins))
Rho_hist=np.zeros((Nc,Nbins))
Rho_Ledges=np.zeros((Nc,Nbins))
Rho_Redges=np.zeros((Nc,Nbins))

ICA_hist=np.zeros((Nc,Nbins))
ICA_Ledges=np.zeros((Nc,Nbins))
ICA_Redges=np.zeros((Nc,Nbins))
uICA_hist=np.zeros((Nc,Nbins))
uICA_Ledges=np.zeros((Nc,Nbins))
uICA_Redges=np.zeros((Nc,Nbins))

Wgth_hist=np.zeros((Nc,Nbins))
Wgth_Ledges=np.zeros((Nc,Nbins))
Wgth_Redges=np.zeros((Nc,Nbins))
uWgth_hist=np.zeros((Nc,Nbins))
uWgth_Ledges=np.zeros((Nc,Nbins))
uWgth_Redges=np.zeros((Nc,Nbins))

for c in range(Nc):
    histUseDensity = False
    # BOLD Model
    aux_mask_Kappa = np.squeeze(KM_maps[:,:,:,c])
    aux_FR2        = np.squeeze(FR2_maps[:,:,:,c])
    aux_input      = aux_FR2[aux_mask_Kappa]
    FR2_hist[c,:], aux = np.histogram(aux_input,Nbins, density=histUseDensity)
    FR2_Ledges[c,:]    = aux[:-1]
    FR2_Redges[c,:]    = aux[1:]
    _,_,Nz         = aux_FR2.shape
    print(Nz)
    
    aux_cR2       = np.squeeze(cR2_maps[:,:,:,c])
    aux_input     = aux_cR2[aux_mask_Kappa]
    cR2_hist[c,:], aux = np.histogram(aux_input,Nbins,density=histUseDensity)
    cR2_Ledges[c,:]    = aux[:-1]
    cR2_Redges[c,:]    = aux[1:]
    
    aux_Kappa     = np.squeeze(Kappa_maps[:,:,:,c])
    aux_input     = aux_Kappa[aux_mask_Kappa]
    Kappa_hist[c,:], aux = np.histogram(aux_input, Nbins,density=histUseDensity);
    Kappa_Ledges[c,:] = aux[:-1]
    Kappa_Redges[c,:] = aux[1:]
    
    # Non-BOLD Model
    aux_mask_Rho = np.squeeze(RM_maps[:,:,:,c])
    aux_FS0      = np.squeeze(FS0_maps[:,:,:,c])
    aux_input    = aux_FS0[aux_mask_Rho]
    FS0_hist[c,:], aux = np.histogram(aux_input,Nbins,density=histUseDensity)
    FS0_Ledges[c,:]    = aux[:-1]
    FS0_Redges[c,:]    = aux[1:]
    
    aux_cS0   = np.squeeze(cS0_maps[:,:,:,c])
    aux_input = aux_cS0[aux_mask_Rho]
    cS0_hist[c,:], aux = np.histogram(aux_input,Nbins,density=histUseDensity)
    cS0_Ledges[c,:]    = aux[:-1]
    cS0_Redges[c,:]    = aux[1:]
    
    aux_Rho   = np.squeeze(Rho_maps[:,:,:,c])
    aux_input = aux_Rho[aux_mask_Rho]
    Rho_hist[c,:], aux = np.histogram(aux_input, Nbins, density=histUseDensity)
    Rho_Ledges[c,:]    = aux[:-1]
    Rho_Redges[c,:]    = aux[1:]
    
    #ICA Maps
    aux_mask_ICA = np.squeeze(ICAM_maps[:,:,:,c])
    aux_ICA      = np.squeeze(ICA_maps[:,:,:,c])
    aux_input    = aux_ICA[aux_mask_ICA] ##### <--------------- Need to think of appropriate masking here
    ICA_hist[c,:], aux = np.histogram(aux_input, Nbins)
    ICA_Ledges[c,:]    = aux[:-1]
    ICA_Redges[c,:]    = aux[1:]
    
    aux_ICA   = np.squeeze(ICA_maps[:,:,:,c])
    aux_input = aux_ICA[mask_orig]
    uICA_hist[c,:], aux = np.histogram(aux_input, Nbins)
    uICA_Ledges[c,:]    = aux[:-1]
    uICA_Redges[c,:]    = aux[1:]
    
    #Weight Maps
    _,_,Nzz = aux_mask_ICA.shape
    print(Nzz)
    if ~(Nzz==Nz):
        aux_mask_ICA = aux_mask_ICA[:,:,np.arange(Nz)]
    aux_Wgth   = np.squeeze(Wgth_maps[:,:,:,c])
    aux_input  = aux_Wgth[aux_mask_ICA] ##### <--------------- Need to think of appropriate masking here
    Wgth_hist[c,:], aux = np.histogram(aux_input, Nbins)
    Wgth_Ledges[c,:]    = aux[:-1]
    Wgth_Redges[c,:]    = aux[1:]
    
    aux_Wgth   = np.squeeze(Wgth_maps[:,:,:,c])
    aux_input  = aux_Wgth[mask_orig]
    uWgth_hist[c,:], aux = np.histogram(aux_input, Nbins)
    uWgth_Ledges[c,:]    = aux[:-1]
    uWgth_Redges[c,:]    = aux[1:]
    
    
    
Hist_dis_cs = ColumnDataSource(data=dict(FR2_hist=np.ones((Nbins,)),   FR2_Redges=range(Nbins),     FR2_Ledges=range(Nbins),
                                         FS0_hist=np.ones((Nbins,)),   FS0_Redges=range(Nbins),     FS0_Ledges=range(Nbins),
                                         cR2_hist=np.ones((Nbins,)),   cR2_Redges=range(Nbins),     cR2_Ledges=range(Nbins),
                                         cS0_hist=np.ones((Nbins,)),   cS0_Redges=range(Nbins),     cS0_Ledges=range(Nbins),
                                         Kappa_hist=np.ones((Nbins,)), Kappa_Redges=range(Nbins),   Kappa_Ledges=range(Nbins),
                                         Rho_hist=np.ones((Nbins,)),   Rho_Redges=range(Nbins),     Rho_Ledges=range(Nbins),
                                         ICA_hist=np.ones((Nbins,)),   ICA_Redges=range(Nbins),     ICA_Ledges=range(Nbins),
                                         uICA_hist=np.ones((Nbins,)),  uICA_Redges=range(Nbins),    uICA_Ledges=range(Nbins),
                                         Wgth_hist=np.ones((Nbins,)),   Wgth_Redges=range(Nbins),   Wgth_Ledges=range(Nbins),
                                         uWgth_hist=np.ones((Nbins,)),  uWgth_Redges=range(Nbins),  uWgth_Ledges=range(Nbins)))
                                             
Hist_avl_cs  = ColumnDataSource(data=dict(FR2_hist=FR2_hist,     FR2_Redges=FR2_Redges,     FR2_Ledges=FR2_Ledges,
                                          FS0_hist=FS0_hist,     FS0_Redges=FS0_Redges,     FS0_Ledges=FS0_Ledges,
                                          cR2_hist=cR2_hist,     cR2_Redges=cR2_Redges,     cR2_Ledges=cR2_Ledges,
                                          cS0_hist=cS0_hist,     cS0_Redges=cS0_Redges,     cS0_Ledges=cS0_Ledges,
                                          Kappa_hist=Kappa_hist, Kappa_Redges=Kappa_Redges, Kappa_Ledges=Kappa_Ledges,
                                          Rho_hist=Rho_hist,     Rho_Redges=Rho_Redges,     Rho_Ledges=Rho_Ledges,
                                          ICA_hist=ICA_hist,     ICA_Redges=ICA_Redges,     ICA_Ledges=ICA_Ledges,
                                          uICA_hist=uICA_hist,    uICA_Redges=uICA_Redges,    uICA_Ledges=uICA_Ledges,
                                          Wgth_hist=Wgth_hist,     Wgth_Redges=Wgth_Redges,     Wgth_Ledges=Wgth_Ledges,
                                          uWgth_hist=uWgth_hist,    uWgth_Redges=uWgth_Redges,    uWgth_Ledges=uWgth_Ledges
                                          ))

F_Hists = figure(width=500, height=300, title='F-stat Histograms',title_text_font_size='12pt',x_axis_label='FS0 or FR2',y_axis_label='# Voxels')
F_Hists.xaxis.axis_label_text_font_size = "12pt"
F_Hists.yaxis.axis_label_text_font_size = "12pt"
F_Hists.quad(top='FR2_hist',bottom=0,left='FR2_Ledges',right='FR2_Redges', source=Hist_dis_cs, line_color='green', fill_color='green', alpha=0.5, legend='F R2')
F_Hists.quad(top='FS0_hist',bottom=0,left='FS0_Ledges',right='FS0_Redges', source=Hist_dis_cs, line_color='red', fill_color='red', alpha=0.5, legend='F S0')
#F_Hists.line('FR2_Ledges','FR2_hist',source=Hist_dis_cs, line_width=3, color='green')
#F_Hists.line('FS0_Ledges','FS0_hist',source=Hist_dis_cs, line_width=3, color='red')
C_Hists = figure(width=500, height=300, title='TE-Fits Histograms',title_text_font_size='12pt',x_axis_label='cS0 or cR2',y_axis_label='# Voxels')
C_Hists.xaxis.axis_label_text_font_size = "12pt"
C_Hists.yaxis.axis_label_text_font_size = "12pt"
C_Hists.quad(top='cR2_hist',bottom=0,left='cR2_Ledges',right='cR2_Redges', source=Hist_dis_cs, line_color='green', fill_color='green', alpha=0.5, legend='c R2')
C_Hists.quad(top='cS0_hist',bottom=0,left='cS0_Ledges',right='cS0_Redges', source=Hist_dis_cs, line_color='red', fill_color='red', alpha=0.5, legend='c S0')
FC_Tab01_C = Panel(child=C_Hists, title='Fit Coefficients')
FC_Tab02_F = Panel(child=F_Hists, title='F-stat')
FC_Tabs    = Tabs(tabs=[FC_Tab01_C, FC_Tab02_F])

KR_Hists = figure(width=500, height=300, title='Kappa/Rho Histograms',title_text_font_size='12pt',x_axis_label='Kappa or Rho',y_axis_label='# Voxels')
KR_Hists.xaxis.axis_label_text_font_size = "12pt"
KR_Hists.yaxis.axis_label_text_font_size = "12pt"
KR_Hists.quad(top='Kappa_hist',bottom=0, left='Kappa_Ledges',right='Kappa_Redges', source=Hist_dis_cs, line_color='green', fill_color='green', alpha=0.5, legend='Kappa')
KR_Hists.quad(top='Rho_hist',  bottom=0, left='Rho_Ledges',  right='Rho_Redges',   source=Hist_dis_cs, line_color='red', fill_color='red', alpha=0.5, legend='Rho')

ICA_Hists = figure(width=500, height=300, title='ICA-Z Histograms',title_text_font_size='12pt',x_axis_label='Weight',y_axis_label='# Voxels')
ICA_Hists.xaxis.axis_label_text_font_size = "12pt"
ICA_Hists.yaxis.axis_label_text_font_size = "12pt"
ICA_Hists.quad(top='uICA_hist',bottom=0, left='uICA_Ledges',right='uICA_Redges', source=Hist_dis_cs, line_color='black', fill_color='black', alpha=0.3, legend='ICA-Z - No Threshold')
ICA_Hists.quad(top='ICA_hist',bottom=0, left='ICA_Ledges',right='ICA_Redges', source=Hist_dis_cs, line_color='blue', fill_color='blue', alpha=0.5, legend='ICA-Z - Inside mask')
Wgth_Hists = figure(width=500, height=300, title='Weight Histograms',title_text_font_size='12pt',x_axis_label='Weight',y_axis_label='# Voxels')
Wgth_Hists.xaxis.axis_label_text_font_size = "12pt"
Wgth_Hists.yaxis.axis_label_text_font_size = "12pt"
Wgth_Hists.quad(top='uWgth_hist',bottom=0, left='uWgth_Ledges',right='uWgth_Redges', source=Hist_dis_cs, line_color='black', fill_color='black', alpha=0.3, legend='Weights - No Threshold')
Wgth_Hists.quad(top='Wgth_hist',bottom=0, left='Wgth_Ledges',right='Wgth_Redges', source=Hist_dis_cs, line_color='blue', fill_color='blue', alpha=0.5, legend='Weights - Inside mask')
IW_Tab01_I = Panel(child=ICA_Hists, title='ICA Maps')
IW_Tab02_W = Panel(child=Wgth_Hists, title='Weight Maps')
IW_Tabs    = Tabs(tabs=[IW_Tab01_I, IW_Tab02_W])

# ==============================================================================

# ==============================================================================
#                   JAVA SCRIPT INTERACTIVITY
 
update_ts = CustomJS(args=dict(timeseries_to_display=timeseries_to_display, 
                               comp_ts=available_timeseries, 
                               ffts_to_display=ffts_to_display, 
                               comp_fft=available_ffts,
                               ICApaths=available_ICAmaps,
                               ICAmap_to_display=ICAmap_to_display,
                               Hist_dis_cs=Hist_dis_cs,Hist_avl_cs=Hist_avl_cs
                               ), 
       code="""
         var c            = cb_obj.get('selected')['1d'].indices
         
         var data2disp_ts = timeseries_to_display.get('data')
         x2disp_ts        = data2disp_ts['x']
         y2disp_ts        = data2disp_ts['y']
         var comp_ts      = comp_ts.get('data')
         ts_x             = comp_ts['x']
         ts_y             = comp_ts['y'];
         for (i = 0; i < x2disp_ts.length; i++) {
            y2disp_ts[i]  = ts_y[c][i];
         }
         
         var data2disp_fft = ffts_to_display.get('data')
         x2disp_fft       = data2disp_fft['x']
         y2disp_fft       = data2disp_fft['y']
         var comp_fft     = comp_fft.get('data')
         fft_x            = comp_fft['x']
         fft_y            = comp_fft['y']
         for (i=0; i < x2disp_fft.length; i++) {
             y2disp_fft[i] = fft_y[c][i]
         }
         
         var ICA2display  = ICAmap_to_display.get('data')
         url2display      = ICA2display['url']
         var availICAurls = ICApaths.get('data')
         allICAurls       = availICAurls['urls']
         url2display[0]   = allICAurls[c]
         
         var Hist2disp = Hist_dis_cs.get('data')
         var HistAvl   = Hist_avl_cs.get('data')
         FR2hist2disp_y   = Hist2disp['FR2_hist']
         FR2hist2disp_re  = Hist2disp['FR2_Redges']
         FR2hist2disp_le  = Hist2disp['FR2_Ledges']
         FR2histAvl_y     = HistAvl['FR2_hist']
         FR2histAvl_re    = HistAvl['FR2_Redges']
         FR2histAvl_le    = HistAvl['FR2_Ledges']
         FS0hist2disp_y   = Hist2disp['FS0_hist']
         FS0hist2disp_re  = Hist2disp['FS0_Redges']
         FS0hist2disp_le  = Hist2disp['FS0_Ledges']
         FS0histAvl_y     = HistAvl['FS0_hist']
         FS0histAvl_re    = HistAvl['FS0_Redges']
         FS0histAvl_le    = HistAvl['FS0_Ledges']
         cR2hist2disp_y   = Hist2disp['cR2_hist']
         cR2hist2disp_re  = Hist2disp['cR2_Redges']
         cR2hist2disp_le  = Hist2disp['cR2_Ledges']
         cR2histAvl_y     = HistAvl['cR2_hist']
         cR2histAvl_re    = HistAvl['cR2_Redges']
         cR2histAvl_le    = HistAvl['cR2_Ledges']
         cS0hist2disp_y   = Hist2disp['cS0_hist']
         cS0hist2disp_re  = Hist2disp['cS0_Redges']
         cS0hist2disp_le  = Hist2disp['cS0_Ledges']
         cS0histAvl_y     = HistAvl['cS0_hist']
         cS0histAvl_re    = HistAvl['cS0_Redges']
         cS0histAvl_le    = HistAvl['cS0_Ledges']
         Kappahist2disp_y   = Hist2disp['Kappa_hist']
         Kappahist2disp_re  = Hist2disp['Kappa_Redges']
         Kappahist2disp_le  = Hist2disp['Kappa_Ledges']
         KappahistAvl_y     = HistAvl['Kappa_hist']
         KappahistAvl_re    = HistAvl['Kappa_Redges']
         KappahistAvl_le    = HistAvl['Kappa_Ledges']
         Rhohist2disp_y   = Hist2disp['Rho_hist']
         Rhohist2disp_re  = Hist2disp['Rho_Redges']
         Rhohist2disp_le  = Hist2disp['Rho_Ledges']
         RhohistAvl_y     = HistAvl['Rho_hist']
         RhohistAvl_re    = HistAvl['Rho_Redges']
         RhohistAvl_le    = HistAvl['Rho_Ledges']
         ICAhist2disp_y   = Hist2disp['ICA_hist']
         ICAhist2disp_re  = Hist2disp['ICA_Redges']
         ICAhist2disp_le  = Hist2disp['ICA_Ledges']
         ICAhistAvl_y     = HistAvl['ICA_hist']
         ICAhistAvl_re    = HistAvl['ICA_Redges']
         ICAhistAvl_le    = HistAvl['ICA_Ledges']
         uICAhist2disp_y   = Hist2disp['uICA_hist']
         uICAhist2disp_re  = Hist2disp['uICA_Redges']
         uICAhist2disp_le  = Hist2disp['uICA_Ledges']
         uICAhistAvl_y     = HistAvl['uICA_hist']
         uICAhistAvl_re    = HistAvl['uICA_Redges']
         uICAhistAvl_le    = HistAvl['uICA_Ledges']
         Wgthhist2disp_y   = Hist2disp['Wgth_hist']
         Wgthhist2disp_re  = Hist2disp['Wgth_Redges']
         Wgthhist2disp_le  = Hist2disp['Wgth_Ledges']
         WgthhistAvl_y     = HistAvl['Wgth_hist']
         WgthhistAvl_re    = HistAvl['Wgth_Redges']
         WgthhistAvl_le    = HistAvl['Wgth_Ledges']
         uWgthhist2disp_y   = Hist2disp['uWgth_hist']
         uWgthhist2disp_re  = Hist2disp['uWgth_Redges']
         uWgthhist2disp_le  = Hist2disp['uWgth_Ledges']
         uWgthhistAvl_y     = HistAvl['uWgth_hist']
         uWgthhistAvl_re    = HistAvl['uWgth_Redges']
         uWgthhistAvl_le    = HistAvl['uWgth_Ledges']
         for (i = 0; i < FR2hist2disp_y.length; i++) {
            FR2hist2disp_y[i]     = FR2histAvl_y[c][i];
            FR2hist2disp_re[i]    = FR2histAvl_re[c][i];
            FR2hist2disp_le[i]    = FR2histAvl_le[c][i];
            FS0hist2disp_y[i]     = FS0histAvl_y[c][i];
            FS0hist2disp_re[i]    = FS0histAvl_re[c][i];
            FS0hist2disp_le[i]    = FS0histAvl_le[c][i];
            cR2hist2disp_y[i]     = cR2histAvl_y[c][i];
            cR2hist2disp_re[i]    = cR2histAvl_re[c][i];
            cR2hist2disp_le[i]    = cR2histAvl_le[c][i];
            cS0hist2disp_y[i]     = cS0histAvl_y[c][i];
            cS0hist2disp_re[i]    = cS0histAvl_re[c][i];
            cS0hist2disp_le[i]    = cS0histAvl_le[c][i];
            Kappahist2disp_y[i]   = KappahistAvl_y[c][i];
            Kappahist2disp_re[i]  = KappahistAvl_re[c][i];
            Kappahist2disp_le[i]  = KappahistAvl_le[c][i];
            Rhohist2disp_y[i]     = RhohistAvl_y[c][i];
            Rhohist2disp_re[i]    = RhohistAvl_re[c][i];
            Rhohist2disp_le[i]    = RhohistAvl_le[c][i];
            ICAhist2disp_y[i]     = ICAhistAvl_y[c][i];
            ICAhist2disp_re[i]    = ICAhistAvl_re[c][i];
            ICAhist2disp_le[i]    = ICAhistAvl_le[c][i];
            uICAhist2disp_y[i]     = uICAhistAvl_y[c][i];
            uICAhist2disp_re[i]    = uICAhistAvl_re[c][i];
            uICAhist2disp_le[i]    = uICAhistAvl_le[c][i];
            Wgthhist2disp_y[i]     = WgthhistAvl_y[c][i];
            Wgthhist2disp_re[i]    = WgthhistAvl_re[c][i];
            Wgthhist2disp_le[i]    = WgthhistAvl_le[c][i];
            uWgthhist2disp_y[i]     = uWgthhistAvl_y[c][i];
            uWgthhist2disp_re[i]    = uWgthhistAvl_re[c][i];
            uWgthhist2disp_le[i]    = uWgthhistAvl_le[c][i];
        }
         Hist_dis_cs.trigger('change'); 
         ICAmap_to_display.trigger('change');
         timeseries_to_display.trigger('change');
         ffts_to_display.trigger('change');
    """)

# Additional Code to link with time series
kappa_taptool          = sp_kappa.select(type=TapTool)
kappa_taptool.callback = update_ts 
rho_taptool          = sp_rho.select(type=TapTool)
rho_taptool.callback = update_ts 
var_taptool          = sp_var.select(type=TapTool)
var_taptool.callback = update_ts 
ratio_taptool          = sp_ratio.select(type=TapTool)
ratio_taptool.callback = update_ts 
kvr_taptool          = sp_kvr.select(type=TapTool)
kvr_taptool.callback = update_ts

# ==============================================================================

# ==============================================================================
#                       GRAPH   LAYOUT
top_left  = hplot(sp_tabs_kappa, sp_tabs_rho)
top_right = hplot(sp_tabs_var, sp_tabs_left)
top       = hplot(top_left, top_right)
middle    = hplot(sp_ts, sp_fft)
pl        = vplot(comp_table_DTABLE,top, middle, map_tabs)
h         = vplot(FC_Tabs, IW_Tabs, KR_Hists)
p         = hplot(pl,h)
show(p)
# ==============================================================================
