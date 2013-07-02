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

DIM = 512
CEN = (DIM/2,DIM/2)
#grid_ax = n.arange(-11.5,12,1)
grid_x = n.arange(-3,4,1)
grid_y = n.arange(-8,8,1) * 4
xy = []
# add in the HERA grid
xy = [[i,j] for i in grid_x for j in grid_y]
# add dithered elements at end of each row/col (doesn't do much)

## add in some outliers
xyo = []

xyo.append([ 2.0, -27.00])
xyo.append([ 2.0,  30.00])
xyo.append([ 2.5,  31.50])
xyo.append([ 2.5,   2.50])

xyo.append([ 9.0, -33.25])
xyo.append([ 8.5,  29.08])
xyo.append([15.5, -32.50])
xyo.append([15.0,  28.42])
xyo.append([22.0, -31.75])
xyo.append([21.5,  28.25])
xyo.append([28.5, -31.00])
xyo.append([28.0,  27.75])
xyo.append([35.0, -30.25])
xyo.append([34.5,  27.08])
xyo.append([41.5, -18.00])
xyo.append([41.0,  13.00])
xyo.append([48.0,   0.00])

xy += xyo

SCALAR = 4
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
p.xlim(-500,500); p.ylim(-200,800)
p.xlabel('Position (m)')
p.subplot(232)
p.imshow(a.img.recenter(uv, CEN), vmin=0, 
    origin='lower', interpolation='nearest',
    extent=(-EXT/2,EXT/2,-EXT/2,EXT/2))
p.xlim(-500,500); p.ylim(-500,500)
p.xlabel('u ($\\lambda$)')
p.subplot(233)
p.imshow(a.img.recenter(uv_clip, CEN), vmin=0, 
    origin='lower', interpolation='nearest',
    extent=(-EXT/2,EXT/2,-EXT/2,EXT/2))
p.xlim(-500,500); p.ylim(-500,500)
p.xlabel('u ($\\lambda$)')
p.subplot(235)
p.imshow(a.img.recenter(10*n.log10(bm), CEN), vmax=0, vmin=-40, 
    origin='lower', interpolation='nearest',
    extent=(-RES_EXT,RES_EXT,-RES_EXT,RES_EXT))
p.xlim(-600,600); p.ylim(-600,600)
p.xlabel('$\\theta$ (arcmin)')
p.subplot(236)
p.imshow(a.img.recenter(10*n.log10(bm_clip), CEN), vmax=0, vmin=-20, 
    origin='lower', interpolation='nearest',
    extent=(-RES_EXT,RES_EXT,-RES_EXT,RES_EXT))
p.colorbar()
p.xlim(-600,600); p.ylim(-600,600)
p.xlabel('$\\theta$ (arcmin)')
p.show()

