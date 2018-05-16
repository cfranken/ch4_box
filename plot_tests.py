### python 3

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy as np
from matplotlib.figure import Figure


import scipy.io as spio
### Defining helper functions to translate Matlab structs to python dictionaries 

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


### call in matfiles and translate to python dictionaries 
case1 = loadmat('case1_test.mat')
case2 = loadmat('case2_test.mat')
case3 = loadmat('case3_test.mat')


### import emissions for case 1
case1_nh_ch4 = case1['case1']['nh_ch4_ems'] 
case1_sh_ch4 = case1['case1']['sh_ch4_ems'] 
case1_nh_oh = case1['case1']['nh_oh_ems'] 
case1_sh_oh = case1['case1']['sh_oh_ems'] 

### Import emissions from case 2
case2_nh_ch4 = case2['case2']['nh_ch4_ems'] 
case2_sh_ch4 = case2['case2']['sh_ch4_ems'] 
case2_nh_oh = case2['case2']['nh_oh_ems'] 
case2_sh_oh = case2['case2']['sh_oh_ems'] 


### import emissions for case 3
case3_nh_ch4 = case3['case3']['nh_ch4_ems'] 
case3_sh_ch4 = case3['case3']['sh_ch4_ems'] 
case3_nh_oh = case3['case3']['nh_oh_ems'] 
case3_sh_oh = case3['case3']['sh_oh_ems'] 

# Defining time vector (1980 to 2016)
time = np.linspace(1980, 2016,37)

# plot sources from inversion 
fig = Figure()
FigureCanvas(fig)
ch4_subplot = fig.add_subplot(1,2,1)
ch4_subplot.plot(time , case1_nh_ch4, color = 'black')

ch4_subplot.plot(time, case2_nh_ch4, color='green')
ch4_subplot.plot(time, case3_nh_ch4, color = 'red')
ch4_subplot.set_title('Methane Emissions from Inversion')
ch4_subplot.set_xlabel('Years')
ch4_subplot.set_ylabel('Terragrams')

oh_subplot = fig.add_subplot(1,2,2)
oh_subplot.plot(time , case1_nh_oh, color = 'black')
oh_subplot.plot(time, case2_nh_oh, color = 'green')
oh_subplot.plot(time, case3_nh_oh, color = 'red')
oh_subplot.set_title('OH Source from Inversion')
oh_subplot.set_xlabel('Years')
oh_subplot.set_ylabel('Teragrams')

fig.legend(['Case 1' , 'Case 2', 'Case 3'], loc='center right')
fig.subplots_adjust(wspace=0.4)
fig.savefig('case_timeseries')

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

fig2 = Figure()
ch4_con = fig2.add_subplot(2,2,1)
ch4_con.plot(time , case1_nh_ch4_con, color = 'black')
ch4_con.plot(time , case2_nh_ch4_con, color = 'green')
ch4_con.plot(time , case3_nh_ch4_con, color = 'red')
ch4_con.set_xlabel('time')
ch4_con.set_ylabel('ppb')
ch4_con.set_title('CH4 Concentrations')

# Plot OH
oh_con = fig2.add_subplot(2,2,2)
oh_con.plot(time , case1_nh_oh_con, color = 'black')
oh_con.plot(time , case2_nh_oh_con, color = 'green')
oh_con.plot(time , case3_nh_oh_con, color = 'red')
oh_con.set_xlabel('time')
oh_con.set_ylabel('ppb')
oh_con.set_title('OH Concentrations')

# Plot CO concentrations
co_con = fig2.add_subplot(2,2,3)
co_con.plot(time , case1_nh_co_con, color = 'black')
co_con.plot(time , case2_nh_co_con, color = 'green')
co_con.plot(time , case3_nh_co_con, color = 'red')
co_con.set_xlabel('time')
co_con.set_ylabel('ppb')
co_con.set_title('CO Concentrations')

fig2.legend(['Case 1', 'case 2', ' case 3'], loc='lower right')
fig2.subplots_adjust(wspace=0.4)
fig.savefig('case_concentrations')
