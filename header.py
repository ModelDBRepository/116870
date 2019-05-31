
# target values of the file 'vm_trace.txt'
# ge=0.032057
# gi=0.096171
# se=0.008014
# si=0.024043


vmFile = 'vm_trace.txt'         # input file
resfile = 'g_dist_parms.txt'    # output file

Iext = 0.                   # constant injected current in nA
gtot = 0.146928             # total input conductance in uS
C = 0.33                    # capacitance in nF
gl = 0.0187                 # leak conductance in uS
Vl = -99.                   # leak reversal potential in mV
Ve = 0.                     # rev. potential of exc. in mV
Vi = -75.                   # rev. potential of inh. in mV
te = 2.728                  # exc. corr. time constant in ms 
ti = 10.49                  # inh. corr. time constant in ms 


vt=-20.                     # threshold for spike detection
dt=0.048                    # time step in ms
t_pre = 5.                  # excluded time preceding spike
t_post = 10.                # excluded time after spike

n_smooth = 3
n_ival = 100                # nb of intervals analysed

n_minISI = 1000             # min nb of datapoints in interval
n_maxISI = 10000            # max nb of datapoints in interval

g_start = [0.03, 0.001, 0.001]

pre=int(t_pre/dt)           # excluded pre spike steps
ahp=int(t_post/dt)          # excluded post spike steps

he1 = 1.-dt/te
he2 = dt/te
hi1 = 1.-dt/ti
hi2 = dt/ti
