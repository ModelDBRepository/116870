
from header import *
from Numeric import *
from pylab import *
from simplex import *
from os import system
from methods import *


h = {}
h['lin'] = zeros(1, Float)
h['cross'] = zeros(1, Float)
h['sq'] = zeros(1, Float)
h['rest'] = zeros(1, Float)

hp = {}
hp['lin'] = zeros(1, Float)
hp['sq'] = zeros(1, Float)


a = zeros(1, Float)
b = zeros(1, Float)


def norm(ge, gi, se, si):
    """computes the norm"""
    hlin_e = -ge*(he1-1.)*he2*te/(2.*dt*se**2)
    hcross_e = he1*te/(2.*dt*se**2)
    hsq_e = -(1. + he1**2)*te/(4.*dt*se**2)
    hrest_e = -ge**2*he2**2*te/(4.*dt*se**2)
    hlin_i = -gi*(hi1-1.)*hi2*ti/(2.*dt*si**2)
    hcross_i = hi1*ti/(2.*dt*si**2)
    hsq_i = -(1. + hi1**2)*ti/(4.*dt*si**2)
    hrest_i = -gi**2*hi2**2*ti/(4.*dt*si**2)
#
    hlin0_e = ge*(2.*dt - he1*he2*te)/(2.*dt*se**2)
    hsq0_e = -(2.*dt + he1**2*te)/(4.*dt*se**2)
    hrest0_e = -ge**2/(2.*se**2)
    hlinf_e = ge*he2*te/(2.*dt*se**2)
    hsqf_e = -te/(4.*dt*se**2)
    hplin_e = hlin0_e
    hpsq_e = hsq0_e
    hlin0_i = gi*(2.*dt - hi1*hi2*ti)/(2.*dt*si**2)
    hsq0_i = -(2.*dt + hi1**2*ti)/(4.*dt*si**2)
    hrest0_i = -gi**2/(2.*si**2)
    hlinf_i = gi*hi2*ti/(2.*dt*si**2)
    hsqf_i = -ti/(4.*dt*si**2)
    hplin_i = hlin0_i
    hpsq_i = hsq0_i
#
#   logarithm of norm
    nrml = log(-pi/hsq0_e)/2. - hlin0_e**2/4./hsq0_e + hrest0_e
    nrml += log(-pi/hsq0_i)/2. - hlin0_i**2/4./hsq0_i + hrest0_i
#
    for k in range(1,n-1):
        hplin_e = hlin_e - hplin_e*hcross_e/2./hpsq_e
        hpsq_e = hsq_e - hcross_e**2/4./hpsq_e
        hplin_i = hlin_i - hplin_i*hcross_i/2./hpsq_i
        hpsq_i = hsq_i - hcross_i**2/4./hpsq_i
#
        nrml += log(-pi/hpsq_e)/2. - hplin_e**2/4./hpsq_e + hrest_e
        nrml += log(-pi/hpsq_i)/2. - hplin_i**2/4./hpsq_i + hrest_i

#
    hplin_e = hlinf_e-hplin_e*hcross_e/2./hpsq_e
    hpsq_e = hsqf_e-hcross_e**2/4./hpsq_e
    hplin_i = hlinf_i-hplin_i*hcross_i/2./hpsq_i
    hpsq_i = hsqf_i-hcross_i**2/4./hpsq_i
#
    nrml += log(-pi/hpsq_e)/2. - hplin_e**2/4./hpsq_e + hrest_e
    nrml += log(-pi/hpsq_i)/2. - hplin_i**2/4./hpsq_i + hrest_i
#
#
    return nrml




def int(ge, gi, se, si):
    """computes the integrated probability"""
    a0 = a[1:size(a)-1]
    ap = a[2:]
    b0 = b[1:size(b)-1]
    bp = b[2:]
    bm = b[:size(b)-2]
#
    # arrays shifted: h[''][k] refers to g_e^{k+1}
    h['lin'][1:n-2] = (-ge*(he1-1.)*he2*te/se**2 + a0*(-b0*(1.+hi1**2)+gi*hi2+hi1*(bm+bp-gi*hi2))*ti/si**2)/2./dt
    h['cross'][1:n-2] = (he1*te/se**2+a0*ap*hi1*ti/si**2)/2./dt
    h['sq'][1:n-2] = -((1.+he1**2)*te/se**2 + a0**2*(1.+hi1**2)*ti/si**2)/4./dt
    h['rest'][1:n-2] = -te*(ge*he2/se)**2/4./dt - (b0-bm*hi1-gi*hi2)**2*ti/(4.*dt*si**2)
#
#
    h['lin'][0] = -(ge*(he1-1.)*he2*te/se**2 + a[0]*(b[0]-b[1]*hi1+b[0]*hi1**2+gi*(hi1-1.)*hi2)*ti/si**2)/2./dt
    h['cross'][0] = (he1*te/se**2+a[0]*a[1]*hi1*ti/si**2)/2./dt
    h['sq'][0] = -((1.+he1**2)*te/se**2 + a[0]**2*(1.+hi1**2)*ti/si**2)/4./dt
    h['rest'][0] = -te*(ge*he2/se)**2/4./dt - (b[0]-gi*hi2)**2*ti/4./dt/si**2
    h['lin'][n-2] = (ge*he2*te/se**2 + a[n-2]*(b[n-3]*hi1-b[n-2]+gi*hi2)*ti/si**2)/2./dt
    h['sq'][n-2] = -(te/se**2 + a[n-2]**2*ti/si**2)/4./dt
    h['rest'][n-2] = -te*(ge*he2/se)**2/4./dt - (b[n-2]-b[n-3]*hi1-gi*hi2)**2*ti/(4.*dt*si**2)
#
    hlin0_e = ge*(2.*dt-he1*he2*te)/(2.*dt*se**2)
    hlin0_i = (2.*dt*gi+hi1*(b[0]-gi*hi2)*ti)/(2.*dt*si**2) 
    hcross0_e = he1*te/(2.*dt*se**2)
    hcross0_i = a[0]*hi1*ti/(2.*dt*si**2)
    hsq0_e = -(2.*dt+he1**2*te)/(4.*dt*se**2)
    hsq0_i = -(2.*dt + hi1**2*ti)/(4.*dt*si**2)
    hrest0 = -(ge/se)**2/2. - (gi/si)**2/2.
    hp['lin'][0] = h['lin'][0] - hcross0_e*hlin0_e/2./hsq0_e - hcross0_i*hlin0_i/2./hsq0_i
    hp['sq'][0] = h['sq'][0] - hcross0_e**2/4./hsq0_e - hcross0_i**2/4./hsq0_i

    for k in range(1, n-1):
        hp['lin'][k] = h['lin'][k] - hp['lin'][k-1]*h['cross'][k-1]/2./hp['sq'][k-1]
        hp['sq'][k] = h['sq'][k] - h['cross'][k-1]**2/4./hp['sq'][k-1]

#
# prl: logarithm of probability
    prlv = -hp['lin']**2/4./hp['sq'] + log(-pi/hp['sq'])/2. + h['rest']
    prl = -(hlin0_e**2/hsq0_e+hlin0_i**2/hsq0_i)/4. + log(pi**2/hsq0_e/hsq0_i)/2. + hrest0         
    prl += sum(prlv)
#
    return prl




def prob(para):
    """returns the logarithm of the inverse probability (to be minimised)"""
    [ge, se, si] = [para[0],para[1],para[2]]
    gi = gtot-gl-ge
    return -(int(ge, gi, se, si)-norm(ge, gi, se, si))



sf = open(resfile, 'w')         # open file and save result
sf.write('#\t\tcurrent ival\t\t\ttemporary average\n')
sf.write('# npts    ge\t  gi\t  se\t  si\t  ge\t  gi\t  se\t  si\n\n')
sf.close()

n_ival = ivals(vmFile,d)          # get intervals from file
res = zeros((n_ival,9,), Float)

for iInt in range(n_ival):
    vm=smooth_g(d['%d'%(iInt+1)],n_smooth)     # smooth voltage trace
    n = size(vm)
    
    for name in h.keys():
        h[name] = resize(h[name],(n-1,))
		
    for name in hp.keys():
	hp[name] = resize(hp[name],(n-1,))
		
    vm0 = vm[:size(vm)-1]
    vmp = vm[1:]
    a = -(vm0-Ve)/(vm0-Vi)
    b = -(gl*(vm0-Vl)+C*(vmp-vm0)/dt-Iext)/(vm0-Vi)
    print "processing ISI", iInt+1, "with", n, "datapoints"
    bp = minimise(prob, g_start)    # maximise probability
#
    res[iInt,0] = n                 # nb of dataponts
    res[iInt,1] = bp[0]             # exc. mean
    res[iInt,2] = gtot-gl-bp[0]     # inh. mean
    res[iInt,3:5] = bp[1:3]         # exc. and inh. SDs
    # momentary averages
    res[iInt,5:9]=mean(res[:iInt+1,:])[1:5]
        
    sf = open(resfile, 'a')         # open file and save result
    line='%d\t'%(res[iInt,0])
    for i in range(1,9):
        line+='%2.5f\t'%(res[iInt,i])
            
    sf.write(line+'\n')
    sf.close()











