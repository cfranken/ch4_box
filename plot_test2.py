### python 3

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy as np
from matplotlib.figure import Figure
import matplotlib.pylab as plt

import scipy.io as spio

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict        

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

case1 = loadmat('case1_test.mat')
case2 = loadmat('case2_test.mat')
case3 = loadmat('case3_test.mat')
case4 = loadmat('case4_test.mat')
case5 = loadmat('case5_test.mat')

case1_nh_ch4 = case1['case1']['nh_ch4_ems'] 
case1_sh_ch4 = case1['case1']['sh_ch4_ems'] 
case1_nh_oh = case1['case1']['nh_oh_ems'] 
case1_sh_oh = case1['case1']['sh_oh_ems'] 
case1_nh_co = case1['case1']['nh_co_ems'] 
case1_sh_co = case1['case1']['sh_co_ems'] 



case2_nh_ch4 = case2['case2']['nh_ch4_ems'] 
case2_sh_ch4 = case2['case2']['sh_ch4_ems'] 
case2_nh_oh = case2['case2']['nh_oh_ems'] 
case2_sh_oh = case2['case2']['sh_oh_ems'] 
case2_nh_co = case2['case2']['nh_co_ems'] 
case2_sh_co = case2['case2']['sh_co_ems'] 



case3_nh_ch4 = case3['case3']['nh_ch4_ems'] 
case3_sh_ch4 = case3['case3']['sh_ch4_ems'] 
case3_nh_oh = case3['case3']['nh_oh_ems'] 
case3_sh_oh = case3['case3']['sh_oh_ems'] 
case3_nh_co = case3['case3']['nh_co_ems'] 
case3_sh_co = case3['case3']['sh_co_ems'] 

case4_nh_ch4 = case4['case4']['nh_ch4_ems'] 
case4_sh_ch4 = case4['case4']['sh_ch4_ems'] 
case4_nh_oh = case4['case4']['nh_oh_ems'] 
case4_sh_oh = case4['case4']['sh_oh_ems'] 
case4_nh_co = case4['case4']['nh_co_ems'] 
case4_sh_co = case4['case4']['sh_co_ems'] 
"""
# extract kx 
nh_kx = case5['case5']['nh_kx']; 
sh_kx = case5['case5']['sh_kx']; 
"""

case5_nh_ch4 = case5['case5']['nh_ch4_ems'] 
case5_sh_ch4 = case5['case5']['sh_ch4_ems'] 
case5_nh_oh = case5['case5']['nh_oh_ems'] 
case5_sh_oh = case5['case5']['sh_oh_ems'] 
case5_nh_co = case5['case5']['nh_co_ems'] 
case5_sh_co = case5['case5']['sh_co_ems'] 



# add up nh and sh emissions 
# Methane CH4 
case1_ch4_ems = case1_nh_ch4 + case1_sh_ch4
case2_ch4_ems = case2_nh_ch4 + case2_sh_ch4
case3_ch4_ems = case3_nh_ch4 + case3_sh_ch4
case4_ch4_ems = case4_nh_ch4 + case4_sh_ch4
case5_ch4_ems = case5_nh_ch4 + case5_sh_ch4

# Hydroxyl OH 
case1_oh_ems = case1_nh_oh + case1_sh_oh
case2_oh_ems = case2_nh_oh + case2_sh_oh
case3_oh_ems = case3_nh_oh + case3_sh_oh
case4_oh_ems = case4_nh_oh + case4_sh_oh
case5_oh_ems = case5_nh_oh + case5_sh_oh

# Carbon monoxide CO 
case1_co_ems = case1_nh_co + case1_sh_co
case2_co_ems = case2_nh_co + case2_sh_co 
case3_co_ems = case3_nh_co + case3_sh_co 
case4_co_ems = case4_nh_co + case4_sh_co 
case5_co_ems = case5_nh_co + case5_sh_co 

### Extract observations for plotting fits
nh_ch4_obs = case3['obs']['nh_ch4']
sh_ch4_obs = case3['obs']['sh_ch4']

nh_co_obs = case3['obs']['nh_co']
sh_co_obs = case3['obs']['sh_co']

time = np.linspace(1980, 2016,37)

# plot sources from inversion 
fig = Figure()
FigureCanvas(fig)
ch4_subplot = fig.add_subplot(2,2,1)
ch4_subplot.plot(time , case1_ch4_ems, color = 'black')
ch4_subplot.plot(time, case2_ch4_ems, color='green')
ch4_subplot.plot(time, case3_ch4_ems, color = 'red')
ch4_subplot.plot(time, case4_ch4_ems, color = 'blue')
ch4_subplot.plot(time, case5_ch4_ems, color = 'purple')
ch4_subplot.set_title(r'\text{$CH_4$} Emissions from Inversion')
ch4_subplot.set_xlabel('Years')
ch4_subplot.set_ylabel('Terragrams')

oh_subplot = fig.add_subplot(2,2,2)
oh_subplot.plot(time , case1_oh_ems, color = 'black')
oh_subplot.plot(time, case2_oh_ems, color = 'green')
oh_subplot.plot(time, case3_oh_ems, color = 'red')
oh_subplot.plot(time, case4_oh_ems, color = 'blue')
oh_subplot.plot(time, case5_oh_ems, color = 'purple')
oh_subplot.set_title('OH Source from Inversion')
oh_subplot.set_xlabel('Years')
oh_subplot.set_ylabel(r'$\frac{molecules}{cm^3}$')

# plot CO emissions 
co_subplot = fig.add_subplot(2,2,3)
co_subplot.plot(time , case1_co_ems, color = 'black')
co_subplot.plot(time, case2_co_ems, color = 'green')
co_subplot.plot(time, case3_co_ems, color = 'red')
co_subplot.plot(time, case4_co_ems, color = 'blue')
co_subplot.plot(time, case5_co_ems, color = 'purple')
co_subplot.set_title('CO Source from Inversion')
co_subplot.set_xlabel('Years')
co_subplot.set_ylabel('Teragrams')
"""
kx_subplot = fig.add_subplot(2,2,4)
kx_subplot.plot(time, nh_kx, time, sh_kx)
kx_subplot.set_title('Arbitrary Reaction with OH')
kx_subplot.set_xlabel('years')
kx_subplot.set_ylabel('1/seconds')
"""

fig.legend(['Case 1' , 'Case 2', 'Case 3', 'Case 4', 'case 5'], loc='center right')
fig.subplots_adjust(wspace=0.4)
fig.tight_layout()
fig.savefig('case_emissions')

### plot the concentrations
# case 1
case1_nh_ch4_con = case1 ['case1']['concentrations']['nh_ch4']
case1_nh_oh_con = case1 ['case1']['concentrations']['nh_oh']
case1_nh_co_con = case1 ['case1']['concentrations']['nh_co']

# case 2
case2_nh_ch4_con = case2['case2']['concentrations']['nh_ch4']
case2_nh_oh_con = case2['case2']['concentrations']['nh_oh']
case2_nh_co_con = case2['case2']['concentrations']['nh_co']


# Case 3
case3_nh_ch4_con = case3 ['case3']['concentrations']['nh_ch4']
case3_nh_oh_con = case3 ['case3']['concentrations']['nh_oh']
case3_nh_co_con = case3 ['case3']['concentrations']['nh_co']
case3_sh_co_con = case3 ['case3']['concentrations']['sh_co']

# Case 3
case4_nh_ch4_con = case4 ['case4']['concentrations']['nh_ch4']
case4_nh_oh_con = case4 ['case4']['concentrations']['nh_oh']
case4_nh_co_con = case4 ['case4']['concentrations']['nh_co']
case4_sh_co_con = case4 ['case4']['concentrations']['sh_co']

# case 5
case5_nh_ch4_con = case5 ['case5']['concentrations']['nh_ch4']
case5_nh_oh_con = case5 ['case5']['concentrations']['nh_oh']
case5_nh_co_con = case5 ['case5']['concentrations']['nh_co']
case5_sh_co_con = case5 ['case5']['concentrations']['sh_co']




fig2 = plt.figure()
ch4_con = fig2.add_subplot(2,2,1)
ch4_con.scatter(time, nh_ch4_obs)
ch4_con.plot(time , case1_nh_ch4_con, color = 'black')
ch4_con.plot(time , case2_nh_ch4_con, color = 'green')
ch4_con.plot(time , case3_nh_ch4_con, color = 'red')
ch4_con.plot(time , case4_nh_ch4_con, color = 'blue')
ch4_con.plot(time , case5_nh_ch4_con, color = 'purple')
ch4_con.set_xlabel('time')
ch4_con.set_ylabel('ppb')
ch4_con.set_title(r'\text{$CH_4$} Concentrations')

# Plot OH

oh_con = fig2.add_subplot(2,2,2)
oh_con.plot(time , case1_nh_oh_con, color = 'black')
oh_con.plot(time , case2_nh_oh_con, color = 'green')
oh_con.plot(time , case3_nh_oh_con, color = 'red')
oh_con.plot(time , case4_nh_oh_con, color = 'blue')
oh_con.plot(time , case5_nh_oh_con, color = 'purple')
oh_con.set_xlabel('time')
oh_con.set_ylabel('ppb')
oh_con.set_title('OH Concentrations')

# Plot CO concentrations
co_con = fig2.add_subplot(2,2,3)
co_con.scatter(time, nh_co_obs)
co_con.plot(time , case1_nh_co_con, color = 'black')
co_con.plot(time , case2_nh_co_con, color = 'green')
co_con.plot(time , case3_nh_co_con, color = 'red')
co_con.plot(time , case4_nh_co_con, color = 'blue')
co_con.plot(time , case5_nh_co_con, color = 'purple')
co_con.set_xlabel('time')
co_con.set_ylabel('ppb')
co_con.set_title('CO Concentrations')

fig2.legend(['Case 1', 'case 2', ' case 3'], loc='lower right')
fig2.subplots_adjust(wspace=0.4)
fig2.tight_layout()
plt.savefig('case_concentrations.png')
