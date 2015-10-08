#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Benjamin E. Gutierrez 9/15/2015

Do some prelimary work for Javier on bokeh.

Current big issue:  The linked brushing (being able to click on a component and see the same component on other plots) is difficult to work with
		    I am still not enitely sure how to make linked brushing related to component # while being able to freely adjust the order of
		    elements in the plot.  Will need to do more testing to fully understand what is going on.
"""
__version__="0.3"
import sys
import os

path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../melib'))
print "++ INFO [Main]: Using meBasic library located in: %s" % path
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

options = parser.parse_args()

# Check that all necessary options were provided
if options.pdir is None:
        print "++ ERROR: No path for ICA maps png files provided. Program will exit."
        sys.exit() 
if options.dir is None:
	print "++ ERROR: No working directory provided. Program will exit."
	sys.exit()
if options.runID is None:
	print "++ ERROR: No run ID provided. Program will exit."
	sys.exit()
# Check if working directory exits
if not os.path.isdir(options.pdir):
	print "++ ERROR: ICA map directory [%s] does not exist." % options.pdir
	print "++ ERROR: Exiting program now."
	sys.exit()
	
if not os.path.isdir(options.dir):
	print "++ ERROR: Working directory [%s] does not exist." % options.dir
	print "++ ERROR: Exiting program now."
	sys.exit()

# If directory exits, check all necessary files are available.
if not os.path.exists(options.dir+'comp_table.txt'):
	print "++ ERROR: Component Table [%s] does not exists." % (options.dir+'comp_table.txt')
	print "++ ERROR: Exiting program now."
	sys.exit()
if not os.path.exists(options.dir+'accepted.txt'):
	print "++ ERROR: List of Accepted Components [%s] does not exists." % (options.dir+'accepted.txt')
	print "++ ERROR: Exiting program now."
	sys.exit()
if not os.path.exists(options.dir+'meica_mix.1D'):
	print "++ ERROR: Component Timeseries File [%s] does not exists." % (options.dir+'meica_mix.1D')
	print "++ ERROR: Exiting program now."
	sys.exit()
if not os.path.exists(options.dir+options.runID+'.chComp.FR2.nii'):
	print "++ ERROR: Per Component F_R2 Maps [%s] missing." % (options.dir+options.runID+'.chComp.FR2.nii')
	print "++ ERROR: Exiting program now."
	sys.exit()
if not os.path.exists(options.dir+options.runID+'.chComp.FS0.nii'):
	print "++ ERROR: Per Component F_S0 Maps [%s] missing." % (options.dir+options.runID+'.chComp.FS0.nii')
	print "++ ERROR: Exiting program now."
	sys.exit()
if not os.path.exists(options.dir+options.runID+'.chComp.Kappa_mask.nii'):
	print "++ ERROR: Per Component Kappa Masks [%s] missing." % (options.dir+options.runID+'.chComp.Kappa_mask.nii')
	print "++ ERROR: Exiting program now."
	sys.exit()


meicaDir    = os.path.abspath(options.dir)
ICAmap_dir  = os.path.abspath(options.pdir)
Program_dir = os.path.dirname(__file__)
print "++ INFO [Main]: ME-ICA Output Directory is %s" % meicaDir
print "++ INFO [Main]: Path to PNG ICA files is %s" % ICAmap_dir

# Load Inputs into memory
# =======================

accepted_components = np.loadtxt(os.path.join(meicaDir,'accepted.txt')).astype('int')
comp_table          = np.loadtxt(os.path.join(meicaDir,'comp_table.txt'))
comp_timeseries     = np.loadtxt(os.path.join(meicaDir,'meica_mix.1D'))
comp_ffts           = np.loadtxt(os.path.join(meicaDir,'meica_mix_fft.1D'))
freq_axis           = np.loadtxt(os.path.join(meicaDir,'meica_mix_fft_freqs.1D'))

FR2_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.chComp.FR2.nii')
print "++ INFO [Main]: FR2 Maps [%s] loaded successfully. [%s]" % (options.dir+options.runID+'.chComp.FR2.nii',FR2_maps.shape)
FS0_maps,_,_ = meb.niiLoad(options.dir+options.runID+'.chComp.FS0.nii')
print "++ INFO [Main]: FS0 Maps [%s] loaded successfully." % (options.dir+options.runID+'.chComp.FS0.nii')
Kappa_masks,_,_ = meb.niiLoad(options.dir+options.runID+'.chComp.Kappa_mask.nii')
print "++ INFO [Main]: Kappa Masks [%s] loaded successfully." % (options.dir+options.runID+'.chComp.Kappa_mask.nii')

Nt, Nc              = comp_timeseries.shape
fica_psel           = np.zeros((Nc,))
fica_psel[accepted_components] = 1

ICAmap_default_path = os.path.join(Program_dir,'AxSag_Component_XX.png')
print ICAmap_default_path
ICAmap_paths        = glob(str(os.path.join(ICAmap_dir,'AxSag*.png')))

print  "++ INFO: Number of ICA Maps: %d" % len(ICAmap_paths)

# Setting up output file
# ======================
output_file(options.runID+".html", title="Static Bokeh Report for " +options.runID)

# Reading the different features of interest
# ==========================================
Nc,Nf = comp_table.shape
kappa = comp_table[:,1]
rho   = comp_table[:,2]
var   = comp_table[:,3]
ratio = comp_table[:,7]
cID   = comp_table[:,0]
maxFR2= comp_table[:,5]
maxFS0= comp_table[:,6]

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
                                    maxFS0      = maxFS0))

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
]
comp_table_DTABLE = DataTable(source=Source,columns=comp_table_columns,width=1350, height=150, editable=True, selectable=True, sortable=False)

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
                                    
# Feature Plots
# =============
sp_kappa = figure(tools=[TOOLS, HoverKappa],width=325, height=250, y_axis_label='Kappa', toolbar_location='left')
sp_kappa.circle('loc_by_kappa','kappa',size=5,color='comp_color',source=Source)
sp_kappa.yaxis.axis_label_text_font_size = "12pt"

sp_rho   = figure(tools=[TOOLS, HoverRho],width=325, height=250, y_axis_label='Rho', toolbar_location=None)
sp_rho.circle('loc_by_rho','rho',size=5,color='comp_color',source=Source)
sp_rho.yaxis.axis_label_text_font_size = "12pt"

sp_var   = figure(tools=[TOOLS,HoverVar],width=325, height=250, y_axis_label='Variance', toolbar_location=None)
sp_var.circle('loc_by_var','var',size=5,color='comp_color',source=Source)
sp_var.yaxis.axis_label_text_font_size = "12pt"

sp_ratio = figure(tools=[TOOLS,HoverRatio],width=325, height=250, y_axis_label='K/R Ratio', toolbar_location=None)
sp_ratio.circle('loc_by_ratio','ratio',size=5,color='comp_color',source=Source)
sp_ratio.yaxis.axis_label_text_font_size = "12pt"

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
# ==============================================================================

# ==============================================================================
#                     TABS FOR HISTOGRAMS AND BOX PLOTS

## Let's compute the histograms
#for c in range(Nc)
#    aux_fr2_data  = FR2_maps[:,:,:,c]
#    aux_fr2_mask  = (Kappa_masks[:,:,:,c]==1)
#    aux_fr2_input = aux_fr2_data[aux_fr2_mask]
#    aux_fr2_
#fr2_data =[] 
#fr2_idx  =[]
#for c in range(Nc):
#    aux_fr2_data  = FR2_maps[:,:,:,c]
#    aux_fr2_mask  = (Kappa_masks[:,:,:,c]==1)
#    aux_fr2_input = aux_fr2_data[aux_fr2_mask]
#    fr2_data.append(aux_fr2_input)
#    aux_fr2_Nv    = aux_fr2_input.shape
#    fr2_idx       = fr2_idx.append(np.repeat(c,aux_fr2_Nv))
#fr2_na  = np.vstack((fr2_data,fr2_idx)).T
#fr2_df  = pd.DataFrame(fr2_na,columns=['data','cid'])
#fr2_his = Histogram(fr2_df,values='data')
#fr2_bpl = BoxPlot(fr2_df,values='data',label='cid')
#fr2_dot = figure(plot_width=600,plot_height=300)


#aux_data  = FR2_maps[:,:,:,10]
#aux_mask  = (Kappa_masks[:,:,:,10]==1)
#aux_input = aux_data[aux_mask]
#aux_Nv    = aux_input.shape
#print "[%f,%f]" % (aux_input.min(),aux_input.max())
#aux_hist, aux_edges = np.histogram(aux_data,bins=200)
#aux_PD = pd.DataFrame(np.vstack((aux_input,np.repeat(1,aux_Nv))).T,columns=['data','id'])
#fr2_his = figure(plot_width=600,plot_height=300)
#fr2_his.quad(top=aux_hist, bottom=0, left=aux_edges[:-1], right=aux_edges[1:], line_color="#033649")
#fr2_bpl = BoxPlot(aux_PD,values='data',label='id')
#fr2_sel = Select(title="Component ID:", value="000", options=[str(c).zfill(3) for c in range(Nc)])
#fr2_dot = figure(plot_width=600,plot_height=300)
#fr2_dot.circle(x=aux_input,y=aux_input)

#fr2_all = vplot(fr2_sel,fr2_his,fr2_bpl)
#fr2_tab = Panel(child=fr2_all, title='FR2') 

#fs0_fig = figure(plot_width=300,plot_height=800)
#fs0_sel = Select(title="Component ID:", value="010", options=[str(c).zfill(3) for c in range(Nc)])
#fs0_all = vplot(fs0_sel,fs0_fig)
#fs0_tab = Panel(child=fs0_all, title='FS0') 

#tabs = Tabs(tabs=[fr2_tab, fs0_tab])
 
# ==============================================================================

# ==============================================================================
#                   JAVA SCRIPT INTERACTIVITY
 
update_ts = CustomJS(args=dict(timeseries_to_display=timeseries_to_display, 
                               comp_ts=available_timeseries, 
                               ffts_to_display=ffts_to_display, 
                               comp_fft=available_ffts,
                               ICApaths=available_ICAmaps,
                               ICAmap_to_display=ICAmap_to_display), 
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

# ==============================================================================

# ==============================================================================
#                       GRAPH   LAYOUT
top_left  = hplot(sp_kappa, sp_rho)
top_right = hplot(sp_var, sp_ratio)
top       = hplot(top_left, top_right)
middle    = hplot(sp_ts, sp_fft)
pl        = vplot(comp_table_DTABLE,top, middle, ICAmapFigure)
#p         = hplot(pl,tabs)
show(pl)
# ==============================================================================
