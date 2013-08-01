#! /usr/bin/env python
from matplotlib import use
from glob import glob
import os,sys
use('tkagg')
PROVOKE=False
PAPERMWA=False
SQRT=True
import numpy as n
from pylab import *

dlogk = 0.2
myk = 10**n.arange(n.log10(0.06),n.log10(3),dlogk)

def logrebin(ki,yi,dlogk,binsizes=None,ks=None):
    #upweight by the input bin size
    #if the bin sizes are given, use them, otherwise, do diff and throw out the highest k point
    if not binsizes is None:
        yi *= binsizes
    else:
        yi = yi[:-1]*n.diff(ki)
        ki = ki[:-1]
    if ks is None:
        kmin = ki.min()
        kmax = ki.max()
        ks = 10**n.arange(n.log10(kmin)-dlogk/2,n.log10(kmax)+dlogk/2,dlogk)


    dmodel = lambda k: 10**n.poly1d(n.polyfit(n.log10(ki),n.log10(yi),1))(n.log10(k))
    rebinned = dmodel(ks)/(dlogk*ks)
    return ks,rebinned



h = 0.71

"""
Plot the sensitivity of various telescopes
"""

def stepwise(A):
    #put two values for every value input. useful for plotting something with beginning and ending step points
    return n.vstack([A,A]).T.ravel()

"""
The PAPER sensitivity curve is in units of mK^2  and k is in units of h Mpc^-1

"""
PAPER_noise = n.array([[0.2,0.5,20],
            [.5,1.5,200],
            [1.5    ,4,2000]])
#PAPER_noise_midpoints = n.array([[0.35,20],[0.75,200],[2.75,2000]]) #upweight by the dk, divide by the new dk later
#PAPER_noise_model = lambda k: 10**n.poly1d(n.polyfit(n.log10(PAPER_noise_midpoints[:,0]),
#                     n.log10(PAPER_noise_midpoints[:,1]),1))(n.log10(k))
#PAPER_noise_logks = 10**n.arange(n.log10(0.2),n.log10(4),dlogk)
#PAPER_noise_dk = n.diff(PAPER_noise_logks)
#PAPER_noise_logks = PAPER_noise_logks[:-1]
#PAPER_noise_rebinned = PAPER_noise_model(PAPER_noise_logks)/PAPER_noise_dk
#PAPER_noise_logks,PAPER_noise_rebinned = logrebin(PAPER_noise_midpoints[:,0],PAPER_noise_midpoints[:,1],
#    dlogk,binsizes=PAPER_bin_sizes)
"""
The MWA sensitivity curve is in units of densitity ratio (mK^2/<mK^2>) and k is in units of Mpc^-1
"""

MWA_noise = n.array([
[0.03,0.05,0.4e-3],
[0.05,0.08,2.5e-3],
[0.08,0.15,3e-3  ],
[0.15,0.2,7e-3   ],
[0.2,0.3,20e-3   ],
[0.3,0.5,60e-3   ],
])
#make a model of the noise for rebinning
MWA_noise_midpoints = n.array([[(l[1]+l[0])/2,l[2]*(l[1]-l[0])] for l in MWA_noise])#the bin weighted values
MWA_noise_model = lambda k: 10**n.poly1d(n.polyfit(n.log10(MWA_noise_midpoints[:,0]),
                            n.log10(MWA_noise_midpoints[:,1]),1))(n.log10(k))
#make a new set of bins
MWA_noise_logks = 10**n.arange(n.log10(0.03),n.log10(0.5),dlogk)
MWA_noise_dks = n.diff(MWA_noise_logks)
MWA_noise_logks = MWA_noise_logks[:-1]
MWA_noise_rebinned = MWA_noise_model(MWA_noise_logks)/MWA_noise_dks
MWA_noise_logks,MWA_noise_rebinned = logrebin(n.array([(l[1]+l[0])/2 for l in MWA_noise]),
    n.array([l[2] for l in MWA_noise]),dlogk,binsizes=n.array([(l[1]-l[0])/2 for l in MWA_noise]),ks=myk*h)

#MWA_noise_horizon = n.loadtxt('MWA_noise1200_horizonWedge.txt')
#MWA_noise_fov = n.loadtxt('MWA_noise1200_fovWedge.txt')
"""
Jonnies Sensitivity
The HERA sensitivity curve is in units of mK^2 and k is in units of h Mpc^-1
"""

F = n.load('hera24x24_14m_sense.npz')
HERA_noise = n.ma.masked_invalid(F['T_errs'])
HERA_ks = F['ks']
kmin,kmax = HERA_ks.min(),HERA_ks.max()
HERA_logks = 10**n.arange(n.log10(kmin/2),n.log10(kmax),dlogk) #bins
digitized = n.digitize(HERA_ks, HERA_logks)
HERA_noise_logk = n.ma.array([HERA_noise[digitized == i].mean() for i in range(0, len(HERA_logks))])

"""
Josh's sensitivity
"""
dfiles = glob('dillon/*csv')
HERAstages = [os.path.basename(f)[:-4] for f in dfiles]
dsense_conservative = {}
dsense_fiducial = {}
dsense_beardsley = {}
HERA_ks = {}
nlim = 10**4
for HERAstage,dfile in zip(HERAstages,dfiles):
    #the first two k-bins are bad in the conservative case
    D = n.loadtxt(dfile,skiprows=7,delimiter=',',usecols=(0,1))
    dsense_conservative[HERAstage] = D[D[:,1]<nlim,:] #throw out points with really high noise (nlim=10^4mK^2)
    D = n.loadtxt(dfile,skiprows=5,delimiter=',',usecols=(3,4))
    dsense_fiducial[HERAstage] = D[D[:,1]<nlim,:]
    D = n.loadtxt(dfile,skiprows=5,delimiter=',',usecols=(6,7))
    dsense_beardsley[HERAstage] = D[D[:,1]<nlim,:]
    kmin =n.min([dsense_conservative[HERAstage][:,0].min(),
                dsense_fiducial[HERAstage][:,0].min(),
                dsense_beardsley[HERAstage][:,0].min()])
    kmax = n.max([dsense_conservative[HERAstage][:,0].max(),
                dsense_fiducial[HERAstage][:,0].max(),
                dsense_beardsley[HERAstage][:,0].max()])
    if HERAstage.endswith('547'):print kmin
    HERA_ks[HERAstage] = 10**(n.arange(n.log10(kmin),n.log10(kmax),dlogk))


"""
Build the PAPER sensitivity to the same k bins as Josh's HERA
"""
PAPER_noise_midpoints = n.array([[0.06,1],[1,100]])
PAPER_dlogk = 1.0
PAPER_noise_model = lambda k: 10**(n.poly1d(n.polyfit(
                    n.log10(PAPER_noise_midpoints[:,0]),n.log10(PAPER_noise_midpoints[:,1]),
                            1))(n.log10(k)))
PAPER_noise_logks = myk #use the ks from the biggest scope
PAPER_noise_rebinned = PAPER_noise_model(PAPER_noise_logks)*PAPER_dlogk/dlogk
print 'PAPER_noise_logks.min()=',PAPER_noise_logks.min()

"""
Load in a model of the power spectrum
"""

A = n.loadtxt('/Users/danny/Work/radio_astronomy/Models/pritchard/powerTb_lowtau.dat')
z = A[:,0]
k = A[0,:]*h


zs = [7,8,9]
ks = [0.1]

ps = []
deltas = []
pritchard = {}
for Z in zs:
    zi = n.argwhere(n.abs(z-Z)==n.min(n.abs(z-Z))).squeeze()+1
    print zi
    pritchard[Z] = k**3 * A[zi,:] * 1.e6/ (2* n.pi**2)
    if SQRT:
        pritchard[Z] = n.sqrt(pritchard[Z])
    #loglog(k,k**3 * A[zi,:] * 1.e6/ (2* n.pi**2) ,'k')


if SQRT:
    """
    Square root the sensitivities and model
    """
    PAPER_noise_rebinned = n.sqrt(PAPER_noise_rebinned)
    MWA_noise_rebinned = n.sqrt(MWA_noise_rebinned)
    for HERAstage in HERAstages:
        dsense_conservative[HERAstage] = n.sqrt(dsense_conservative[HERAstage])
        dsense_fiducial[HERAstage] = n.sqrt(dsense_fiducial[HERAstage])
        dsense_beardsley[HERAstage] = n.sqrt(dsense_beardsley[HERAstage])
    ytext = '$\Delta [mK]$'
else:
    ytext = '$\Delta^2$ [mK$^2$]'

"""
End of power spectrum loading
"""

z0=8.
T0 = 28.*((1+z0)/10)**(1/2)

#print "T0 = ",T0

"""
Check that we are doing the binning correctly
"""
figure()
ax = subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
#plot the pritchard model
for Z in zs:
    loglog(k,pritchard[Z] ,'k')
#for array in ['HERA37','HERA127','HERA547'][::-1]:
array = 'HERA37'
hlk,hls =logrebin(dsense_beardsley[array][:,0],dsense_beardsley[array][:,1],
                dlogk,ks=myk)
step(dsense_beardsley[array][:,0],dsense_beardsley[array][:,1],lw=1,color='#FF8C00',where='post',label='HERA37')
step(hlk,hls,lw=3,color='#FF8C00',where='post',label='HERA37-binned')
step(PAPER_noise_logks,PAPER_noise_rebinned,color='r',label='PAPER 1200hrs',lw=3,where='post')
#pks = 10**n.arange(n.log10(0.5),n.log10(1),PAPER_dlogk)
pks = n.array([0.1,1,10])
PAPER_noise = n.vstack([10**(n.log10(pks)-PAPER_dlogk/2),10**(n.log10(pks)+PAPER_dlogk/2),PAPER_noise_model(pks)]).T
print PAPER_noise.shape
plot(PAPER_noise[:,:2].ravel(),stepwise(PAPER_noise[:,2]))
print PAPER_noise[:,:2].ravel()
#step(pks,PAPER_noise_model(pks),where='mid',label='PAPER raw',color='r',ls=':')
#plot(pks,PAPER_noise_model(pks),'.')
#plot(myk,PAPER_noise_model(myk))
print pks,PAPER_noise_model(pks)
ylim([10**-1,10**5])
xlim([10**-2,10**1]) 
xlabel('$k$ [h Mpc$^{-1}$]')
ylabel(ytext)
legend(loc='upper left')
subplots_adjust(left=0.15,bottom=0.175)
print 'HERA_binning.png'
savefig('HERA_binning.png')
#annotate(array,[hlk[-1],hls[-1]],xytext=[0.85,0.6+i*.05],textcoords='axes fraction',arrowprops=dict(arrowstyle="->",
#                                connectionstyle="arc3,rad=.2"))


"""
Plot the sensitivities.
"""
#plot the sensitivities
figure()
ax =subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
#loglog(MWA_noise[:,:2].ravel()/h,stepwise(MWA_noise[:,2])*T0**2/h**3,label='MWA 1600hrs',lw=3)
#loglog(PAPER_noise[:,:2].ravel(),stepwise(PAPER_noise[:,2]),'r',label='PAPER 1200hrs',lw=3)
#plot(n.linspace(0.1,4),PAPER_noise_model(n.linspace(0.1,4)))
#plot(PAPER_noise_midpoints[:,0],PAPER_noise_midpoints[:,1],'.')
#plot(n.linspace(0.01,1)/h**3,MWA_noise_model(n.linspace(0.01,1))*T0**2/h**3)
#plot(MWA_noise_midpoints[:,0]/h**3,MWA_noise_midpoints[:,1]*T0**2/h**3,'.')
step(MWA_noise_logks/h,MWA_noise_rebinned*T0**2/h**3,color='b',label='MWA 1600hrs',lw=3,where='post')
step(PAPER_noise_logks,PAPER_noise_rebinned,color='r',label='PAPER 1200hrs',lw=3,where='post')
#plot the pritchard model
for Z in zs:
    loglog(k,pritchard[Z] ,'k')
ylim([10**-1,10**5])
xlim([10**-2,10**1])    

xlabel('$k$ [h Mpc$^{-1}$]')
ylabel(ytext)
legend(loc='upper left')
subplots_adjust(left=0.15,bottom=0.175)
savefig('PAPER_MWA_sensitivity.png')
print 'PAPER_MWA_sensitivity.png'
#step(HERA_ks,HERA_noise,label='HERA 720hrs',lw=3,color='#FF8C00',where='mid')
#step(HERA_logks,HERA_noise_logk,label='HERA 720hrs',lw=3,color='#FF8C00',where='pre')
i=0
for array in ['HERA37','HERA127','HERA547'][::-1]:
    hlk,hls =logrebin(dsense_beardsley[array][:,0],dsense_beardsley[array][:,1],
                dlogk,ks=myk)
    #print "hlk.min()=",hlk.min()

    if array.endswith('547'):
        step(dsense_beardsley[array][:,0],dsense_beardsley[array][:,1],lw=1,color='#FF8C00',where='post',label='HERA')
        step(hlk,hls,lw=3,color='#FF8C00',where='post',label='HERA')


    else:
 #       step(dsense_beardsley[array][:,0],dsense_beardsley[array][:,1],lw=1,color='#FF8C00',where='post')
        step(hlk,hls,lw=3,color='#FF8C00',where='post')
    annotate(array,[hlk[-1],hls[-1]],xytext=[0.85,0.6+i*.05],textcoords='axes fraction',arrowprops=dict(arrowstyle="-",
                                connectionstyle="arc3,rad=.2"))
    i +=1


legend(loc='upper left')
savefig('HERA_sensitivity_compare.png')
print 'HERA_sensitivity_compare.png'
"""
MAke a plot of the estimated fg window
"""

figure()
ax =subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')    
#loglog(MWA_noise[:,:2].ravel()/h,stepwise(MWA_noise[:,2])*T0**2/h**3,label='MWA 1600hrs',lw=3)
#loglog(PAPER_noise[:,:2].ravel(),stepwise(PAPER_noise[:,2]),'r',label='PAPER 1200hrs',lw=3)
#plot(n.linspace(0.1,4),PAPER_noise_model(n.linspace(0.1,4)))
#plot(PAPER_noise_midpoints[:,0],PAPER_noise_midpoints[:,1],'.')
#plot(n.linspace(0.01,1)/h**3,MWA_noise_model(n.linspace(0.01,1))*T0**2/h**3)
#plot(MWA_noise_midpoints[:,0]/h**3,MWA_noise_midpoints[:,1]*T0**2/h**3,'.')
step(MWA_noise_logks/h,MWA_noise_rebinned*T0**2/h**3,color='b',label='MWA 1600hrs',lw=3,where='post')
step(PAPER_noise_logks,PAPER_noise_rebinned,color='r',label='PAPER 1200hrs',lw=3,where='post')
#plot the pritchard model
for Z in zs:
    loglog(k,pritchard[Z] ,'k')
ylim([10**-1,10**5])
xlim([10**-2,10**1])    

xlabel('$k$ [h Mpc$^{-1}$]')
ylabel(ytext)
legend(loc='upper left')
subplots_adjust(left=0.15,bottom=0.175)


#step(dsense_fiducial['HERA547'][:,0],dsense_fiducial['HERA547'][:,1],ls='.-',lw=3,color='g',where='pre',label='fiducial')
#step(dsense_conservative['HERA547'][:,0],dsense_conservative['HERA547'][:,1],ls=':',lw=3,color='g',where='pre',label='conservative')
kmin = n.max((dsense_beardsley['HERA547'][:,0].min(),dsense_conservative['HERA547'][:,0].min()))
kmax = n.min((dsense_beardsley['HERA547'][:,0].max(),dsense_conservative['HERA547'][:,0].max()))

wlogks,woptimistic = logrebin(dsense_beardsley['HERA37'][:,0],dsense_beardsley['HERA37'][:,1],0.05,ks=myk)
print wlogks
wlogks,wpessmistic = logrebin(dsense_conservative['HERA37'][:,0],dsense_conservative['HERA37'][:,1],0.05,ks=myk)
print wlogks
plot(wlogks,woptimistic)
plot(wlogks,wpessmistic)
fill_between(wlogks,woptimistic,y2=wpessmistic,color='#FF8C00',alpha=0.2)
print (wpessmistic-woptimistic)/wpessmistic
print 'HERA_sensitivity_fg_window.png'
savefig('HERA_sensitivity_fg_window.png')





#Add current best sensitivity limits
figure()
subplot(111)
arrowlen = 10.

"""
RESULTS
"""


GMRT_paciga_2013 = n.array([[0.1,2e5],[0.13,4e5],[0.16,1e5],[0.19,1.9e5],[0.275,2e5],[0.31,4e5],[0.4,6e5],[0.5,8e4]])
#plot the MWA and PAPER upper limits
MWA_X13 = n.array([[0.3,(5e3)**2],[0.4,(9e3)**2],[0.9,(15e3)**2]])#the three lowest k bins that are noise dominated in Fig 9 of Dillon et al
PAPER_psa32mr = n.array([[0.25,2e3],[0.27,2.5e3],[0.3,8e3],[0.33,4e3]])#the four lowest k measurements
#load the actual results from a file
lines = open('PSA32_pspec.txt').readlines()
PAPER_psa32mr = []
for line in lines:
    psak = float(line.strip().split(',')[0])
    psa_ul = float(line.strip().split()[1])
    psaval = float(line.strip().split()[2][1:])
    psa2sig = float(line.split('+/-')[1].split(')')[0])
#    if (psaval-psa2sig)<0: 
#    if psaval<10**8:
    if True:
        PAPER_psa32mr.append([psak,psa_ul])
PAPER_psa32mr = n.array(PAPER_psa32mr)

if SQRT:
    """
    sqrt the results
    """
    GMRT_paciga_2013[:,1] = n.sqrt(GMRT_paciga_2013[:,1])
    PAPER_psa32mr[:,1] = n.sqrt(PAPER_psa32mr[:,1])
    MWA_X13[:,1] = n.sqrt(MWA_X13[:,1])


for l in GMRT_paciga_2013:
    plot(l[0],l[1],'vy')
    annotate('',xy=[l[0],l[1]/arrowlen],xycoords='data',xytext=[l[0],l[1]],textcoords='data',
        arrowprops={'arrowstyle':'-','color':'y','lw':1})
#text(GMRT_paciga_2013[-1,0]*1.1,GMRT_paciga_2013[-1,1],'GMRT-2007 \n(2013 revised)\n50 hours')
text(GMRT_paciga_2013[0,0]/4,GMRT_paciga_2013[0,1],'GMRT\n Paciga 2010\n2013(revised)')
if PROVOKE:
    #plot the old GMRT result
    plot(0.15,2e5,'^y',zorder=-10)
    plot(0.65,4e3,'^y',zorder=-10)
    plot(0.15,2e5,'xk',ms=20)
    plot(0.65,4e3,'xk',ms=20)


#plot the model from Pritchard
for Z in zs:
    loglog(k,pritchard[Z] ,'k')


#ylim([10**-1,10**3])
ylim([1e-2,PAPER_psa32mr[:,1].max()*1.1])

xlim([10**-2,10**1])    

if PAPERMWA:
    #plot the sensitivity curves againb
    #loglog(MWA_noise[:,:2].ravel()/h,stepwise(MWA_noise[:,2])*T0**2/h**3,'b',label='MWA 1600hrs',lw=3)
    #loglog(PAPER_noise[:,:2].ravel(),stepwise(PAPER_noise[:,2]),'r',label='PAPER 1200hrs',lw=3)
    step(MWA_noise_logks/h,MWA_noise_rebinned*T0**2/h**3,color='b',label='MWA 1600hrs',lw=3,where='post')
    annotate('MWA 2013',[MWA_noise_logks[0]/h,MWA_noise_rebinned[0]*T0**2/h**3],xytext=[-100,10],color='b',
                    textcoords='offset points',arrowprops=dict(arrowstyle="-",
                    connectionstyle="arc3,rad=.2",color='b'))
    
    step(PAPER_noise_logks,PAPER_noise_rebinned,color='r',label='PAPER 1200hrs',lw=3,where='post')
    annotate('PAPER 2013',[PAPER_noise_logks[0],PAPER_noise_rebinned[0]],xytext=[-100,10],color='r',
                    textcoords='offset points',arrowprops=dict(arrowstyle="-",
                    connectionstyle="arc3,rad=.2",color='red'))
    xlabel('$k$ [h Mpc$^{-1}$]')
    ylabel(ytext)
    subplots_adjust(left=0.15,bottom=0.175)
    
    savefig('GMRT_limits_with_PAPER_MWA_sensitivity.png')
    print 'GMRT_limits_with_PAPER_MWA_sensitivity.png'

for l in MWA_X13:
    plot(l[0]/h,l[1]*h**3,'xb')
    annotate('',xy=[l[0]/h,l[1]*h**3/arrowlen],xycoords='data',xytext=[l[0]/h,l[1]*h**3],textcoords='data',
        arrowprops={'arrowstyle':'-','color':'b'})
#text(MWA_X13[-1,0]+0.5,MWA_X13[-1,1]/20,'MWA-2011\n5 hours\n32ants',color='b')
text(MWA_X13[-1,0]+0.5,MWA_X13[-1,1]/20,'MWA-2011\nDillon 2013',color='b')

#text(0.02,0.3,'MWA 2013',color='b')

for l in PAPER_psa32mr:
    plot(l[0],l[1],'+r')
    annotate('',xy=[l[0],l[1]/arrowlen],xycoords='data',xytext=[l[0],l[1]],textcoords='data',
        arrowprops={'arrowstyle':'-','color':'r','lw':1})
#text(PAPER_psa32mr[0,0]/2.5,PAPER_psa32mr[-1,1]/10,'PAPER-2011\n55 hours\n32ants',color='r')
text(PAPER_psa32mr[0,0]*1.5,PAPER_psa32mr[-1,1]/10,'PAPER-2011\nParsons 2013',color='r')
#print "PAPER label",PAPER_psa32mr[0,1],PAPER_psa32mr[-1,1]
#text(2,3e3,'PAPER 2013',color='r')






xlabel('$k$ [h Mpc$^{-1}$]')
ylabel(ytext)
subplots_adjust(left=0.15,bottom=0.175)

ylim([1e-2,PAPER_psa32mr[:,1].max()*1.1])
print "MWA_PAPER_with_limits_2013.png"
savefig('MWA_PAPER_with_limits_2013.png')



#step(HERA_ks,HERA_noise,label='HERA 720hrs',lw=3,color='#FF8C00')
#step(HERA_logks,HERA_noise_logk,label='HERA 720hrs',lw=3,color='#FF8C00',where='pre')
i=0
for array in ['HERA37','HERA127','HERA547'][::-1]:
    hlk,hls =logrebin(dsense_beardsley[array][:,0],dsense_beardsley[array][:,1],
                dlogk,ks=myk)
    #print "hlk.min()=",hlk.min()

    if array.endswith('547'):
#        step(dsense_beardsley[array][:,0],dsense_beardsley[array][:,1],lw=1,color='#FF8C00',where='post',label='HERA')
        step(hlk,hls,lw=3,color='#FF8C00',where='post',label='HERA')


    else:
 #       step(dsense_beardsley[array][:,0],dsense_beardsley[array][:,1],lw=1,color='#FF8C00',where='post')
        step(hlk,hls,lw=3,color='#FF8C00',where='post')
    annotate(array,[hlk[-1],hls[-1]],xytext=[0.85,0.5+i*.05],textcoords='axes fraction',arrowprops=dict(arrowstyle="-",
                                connectionstyle="arc3,rad=.2"))
    i +=1

text(2.0,0.9,'HERA 2016',color='#FF8C00')
#text(0.02,10**7,'z=8')
print "HERA_with_limits.png"
savefig("HERA_with_limits.png")
