#! /usr/bin/env python
import numpy as n, pylab as p, aipy as a

def aconv2(pos):
    ipos = n.fft.fft2(pos)
    uv = n.abs(n.fft.ifft2(ipos * n.conj(ipos)))
    return uv

def grid(xy,dim):
    pos = n.zeros((dim,dim), dtype=n.float)
    xy_g = n.round(xy).astype(n.int)
    a.utils.add2array(pos,xy_g,n.ones((xy.shape[0],), dtype=n.float))
    return pos

def beam(uv):
    bm = n.abs(n.fft.ifft2(uv))
    return bm / bm.max()

DIM = 2048
CEN = (DIM/2,DIM/2)
#grid_ax = n.arange(-11.5,12,1)
grid_ax = n.arange(-12,12,1)
xy = []
# add in the HERA grid
xy = [[i,j] for i in grid_ax for j in grid_ax]
# add dithered elements at end of each row/col (doesn't do much)
#for i in grid_ax[:-1]:
#    xy.append([i+.5,12])
#    xy.append([12,i+.5])

## add in some outliers
xyo0,xyo1,xyo2 = [],[],[]

xyo0.append([ 12.5,10.5])
xyo0.append([-11.5,12.5])

#xyo.append([  0,-24])
#xyo.append([-24,-24])
#xyo.append([-24,  0])
#xyo.append([-24, 23])

xyo1.append([-12, 35])
xyo1.append([ 11, 35])
xyo1.append([ 35, 35])
xyo1.append([ 35, 11])
xyo1.append([ 35,-12])
xyo1.append([ 35,-36])

xyo2.append([ 35,-60])
xyo2.append([ 11,-60])
xyo2.append([-12,-60])
xyo2.append([-36,-60])
xyo2.append([-60,-60])
xyo2.append([-60,-36])
xyo2.append([-60,-12])
xyo2.append([-60, 11])
xyo2.append([-60, 35])
xyo2.append([-60, 59])

xy += xyo0 + xyo1 + xyo2 + [[-i-.5,-j-.5] for i,j in xyo1] + [[-i-1.5,-j-1.5] for i,j in xyo2]

#xy = 14 * n.array(xy)
SCALAR = 6.
xy = SCALAR * n.array(xy).astype(n.float)
pos = grid(xy, DIM)
inv_bm = a.img.gaussian_beam(SCALAR/2/n.sqrt(2), shape=(DIM,DIM))
_inv_bm = n.fft.fft2(inv_bm)
inv_bm2 = a.img.gaussian_beam(SCALAR/2, shape=(DIM,DIM))
_inv_bm2 = n.fft.fft2(inv_bm2)
pos /= pos.max()
uv = aconv2(pos)
uv[0,0] = 0
pos = n.fft.fft2(n.fft.ifft2(pos) * _inv_bm).astype(n.float)
pos = n.ma.masked_less(pos, 1e-4)
uv_clip = uv.clip(0,1)
uv = n.fft.fft2(n.fft.ifft2(uv) * _inv_bm2).astype(n.float)
bm = beam(uv)
uv_clip = n.fft.fft2(n.fft.ifft2(uv_clip) * _inv_bm2).astype(n.float)
bm_clip = beam(uv_clip)
uv = n.ma.masked_less(uv, 1e-4)
uv_clip = n.ma.masked_less(uv_clip, 1e-4)
EXT = DIM/2/SCALAR * 14 # m
RES = 1. / EXT / a.const.arcmin
RES_EXT = DIM/2 * RES
p.subplot(231)
p.imshow(a.img.recenter(pos, CEN), vmin=0, vmax=1, 
    origin='lower', interpolation='nearest',
    extent=(-EXT,EXT,-EXT,EXT))
p.xlim(-1000,1000); p.ylim(-1000,1000)
p.xlabel('Position (m)')
p.subplot(232)
p.imshow(a.img.recenter(uv, CEN), vmin=0, 
    origin='lower', interpolation='nearest',
    extent=(-EXT/2,EXT/2,-EXT/2,EXT/2))
p.xlim(-1000,1000); p.ylim(-1000,1000)
p.xlabel('u ($\\lambda$)')
p.subplot(233)
p.imshow(a.img.recenter(uv_clip, CEN), vmin=0, 
    origin='lower', interpolation='nearest',
    extent=(-EXT/2,EXT/2,-EXT/2,EXT/2))
p.xlim(-1000,1000); p.ylim(-1000,1000)
p.xlabel('u ($\\lambda$)')
p.subplot(235)
p.imshow(a.img.recenter(10*n.log10(bm), CEN), vmax=0, vmin=-40, 
    origin='lower', interpolation='nearest',
    extent=(-RES_EXT,RES_EXT,-RES_EXT,RES_EXT))
p.xlim(-600,600); p.ylim(-600,600)
p.xlabel('$\\theta$ (arcmin)')
p.subplot(236)
p.imshow(a.img.recenter(10*n.log10(bm_clip), CEN), vmax=0, vmin=-40, 
    origin='lower', interpolation='nearest',
    extent=(-RES_EXT,RES_EXT,-RES_EXT,RES_EXT))
p.xlim(-600,600); p.ylim(-600,600)
p.xlabel('$\\theta$ (arcmin)')
p.show()

