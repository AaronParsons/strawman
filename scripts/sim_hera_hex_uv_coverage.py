#! /usr/bin/env python
import numpy as n, pylab as p, aipy as a

def gen_hex_pos(N):
    pos = []
    for i in xrange(N,2*N-1):
        y = n.sqrt(3)/2 * (2*N-1-i)
        for j in xrange(-i/2,i/2):
            x = j + float(i%2)/2 + 0.5
            pos.append((x, y))
            pos.append((x,-y))
    i = 2*N-1
    for j in xrange(-i/2,i/2):
        x = j + float(i%2)/2 + 0.5
        pos.append((x,0))
    return pos

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
#DIM = 1024
CEN = (DIM/2,DIM/2)
#grid_hex = []
#for w in [16,20,24,28,24,20,16]:
#    for i in xrange(4): grid_hex.append(n.arange(-w/2,w/2))
#
#xy = []
#for i,row in enumerate(grid_hex):
#    for j in row:
#        xy.append([i-len(grid_hex)/2,j])

N = 14.

### add in some outliers
xyo0,xyo1,xyo2 = [],[],[]
##
xyo0.append([ 0.75*N+0.5, n.sqrt(3)/2*(N/2-1)+1./n.sqrt(3)])
xyo0.append([ 0.0*N+0.5, n.sqrt(3)/2*(N+1)-1./n.sqrt(3)])
xyo0.append([-0.75*N+0.0, n.sqrt(3)/2*(N/2-0)+1./n.sqrt(3)])
##
xyo1.append([ 1.5*N-0.5, n.sqrt(3)/2*(N-1)])
xyo1.append([ 0.0*N+0.5, n.sqrt(3)/2*(2*N-1)])
xyo1.append([-1.5*N+1.0, n.sqrt(3)/2*(N-0)])
 
xyo2.append([ 3.0*N-1.5, n.sqrt(3)/2*(-1)])
xyo2.append([ 3.0*N-1.0, n.sqrt(3)/2*(2*N-2)])
xyo2.append([ 1.5*N-0.0, n.sqrt(3)/2*(3*N-2)])
xyo2.append([ 0.0*N+1.0, n.sqrt(3)/2*(4*N-2)])
xyo2.append([-1.5*N+1.5, n.sqrt(3)/2*(3*N-1)])
xyo2.append([-3.0*N+2.0, n.sqrt(3)/2*(2*N-0)])
#
SCALAR = 14.
xy = gen_hex_pos(int(N))
print len(xy)
#p.plot(pos[:,0], pos[:,1], '.')
#p.show()
#import sys; sys.exit(0)

pos0 = grid(SCALAR*n.array(xy).astype(n.float), DIM)
#xy += xyo0 + xyo1 + xyo2 + [[-i+0.5,-j-1+1./n.sqrt(3)] for i,j in xyo1] #+ [[-i,-j+1./n.sqrt(3)-0.5] for i,j in xyo2]
xy += xyo0 + xyo1 + xyo2 + [[-i+0.5,-j+n.sqrt(3)/6] for i,j in xyo1] + [[-i+0.5,-j+n.sqrt(3)/6] for i,j in xyo2]
#xy += xyo1 + xyo2

xy = SCALAR * n.array(xy).astype(n.float)
pos = grid(xy, DIM)
inv_bm = a.img.gaussian_beam(SCALAR/2/n.sqrt(2), shape=(DIM,DIM))
_inv_bm = n.fft.fft2(inv_bm)
inv_bm2 = a.img.gaussian_beam(SCALAR/2, shape=(DIM,DIM))
_inv_bm2 = n.fft.fft2(inv_bm2)
pos /= pos.max()

uv = aconv2(pos) #- aconv2(pos0)
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

