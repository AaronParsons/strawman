#! /usr/bin/env python
import numpy as n, pylab as p, aipy as a
import random, time, matplotlib

cdict = {'red':  ((0.00, 1.0, 1.0), #  0
                  (0.25, 0.0, 0.0), #  1
                  (0.50, 1.0, 1.0), #  2
                  (0.75, 0.0, 0.0), #  3
                  (1.00, 0.0, 0.0)),#  4
         'green':((0.00, 1.0, 1.0),
                  (0.25, 0.0, 0.0),
                  (0.50, 0.0, 0.0),
                  (0.75, 1.0, 1.0),
                  (1.00, 1.0, 1.0)),
         'blue': ((0.00, 1.0, 1.0),
                  (0.25, 0.0, 0.0),
                  (0.50, 1.0, 1.0),
                  (0.75, 1.0, 1.0),
                  (1.00, 1.0, 1.0)),
}
vmin,vmax = 0, 100
cmap1 = matplotlib.colors.LinearSegmentedColormap('steps',cdict,256)

cmap2 = p.get_cmap('gist_yarg')
random.seed()

def aconv2(pos):
    ipos = n.fft.fft2(pos)
    uv = n.abs(n.fft.ifft2(ipos * n.conj(ipos)))
    return uv

def _grid(xy,dim):
    dim *= 4
    pos = n.zeros((dim,dim), dtype=n.float)
    #x = rp[:,0] * n.cos(2*n.pi*rp[:,1])
    #y = rp[:,0] * n.sin(2*n.pi*rp[:,1])
    #xy = n.array([x,y]).transpose()
    xy_g = n.round(xy * dim/8).astype(n.int)
    a.utils.add2array(pos,xy_g,n.ones((xy.shape[0],), dtype=n.float))
    return pos

def _gen_r(dim):
    dim *= 4
    rx,ry = n.indices((dim,dim))
    r = (rx - dim/2)**2 + (ry - dim/2)**2
    r = n.sqrt(r.astype(n.float)).clip(.25,dim) / dim
    return a.img.recenter(r, (dim/2,dim/2))

def _energy(r,pos,uv,uv2,N):
    N2 = N*(N+1)/2
    return n.sum(uv2/r)/(N2*(N2+1)/2) + n.sum(uv/r) / N2

def comp_mr(N, start_xy=None, others=None, start_phs=3, maxiter=1000, doit=True):
    """Compute optimal minimum redundancy configuration.  Starts antennas at
    specified positions (start_xy) on a grid with resolution specified by
    start_phs.  Each iteration, antenna positions are dithered by adding noise.
    After churning for a while on a given resolution, if no
    improvement is made, the resolution of the grid is increased to fine-tune
    the positions."""
    try:
        if not start_xy is None: best_xy = start_xy
        else: best_xy = n.random.uniform(-1,1,size=(N,2))
        cnt,phs = n.Inf, start_phs - 1
        while phs <= 7:
            cnt += 1
            if cnt > maxiter:
                cnt = 0
                phs += 1
                T = 2.**-phs
                dim = 2**(phs)
                print 'Doing DIM=%d, T=%f' % (dim, T),
                r = _gen_r(dim)
                if not others == None: other_pos = _grid(others,dim)
                else: other_pos = 0
                pos = _grid(best_xy,dim) + other_pos
                uv = aconv2(pos); uv2 = aconv2(uv)
                best_E = _energy(r,pos,uv,uv2,N)
                print 'E=', best_E
            xy = best_xy + n.random.normal(scale=T, size=best_xy.shape)
            xy[:,0] = xy[:,0].clip(-1,1)
            xy[:,1] = xy[:,1].clip(-1,1)
            if not doit: xy = best_xy
            pos = _grid(xy,dim) + other_pos
            uv = aconv2(pos); uv2 = aconv2(uv)
            E = _energy(r,pos,uv,uv2,N)
            if E < best_E:
                cnt = 0
                print dim, phs, E, time.time()
                best_E, best_xy = E, xy
            if not doit: cnt = n.Inf
    except(KeyboardInterrupt): print 'Ending'
    return best_E, best_xy

DIM = 2**6
N = 40
others = []
for i in xrange(-12,12):
    for j in xrange(-12,12):
        #r = n.sqrt(i**2 + j**2) / 24 / 0.66
        #ph = n.arctan2(j,i) / (2*n.pi)
        #x = r * n.cos(2*n.pi*ph)
        #y = r * n.sin(2*n.pi*ph)
        #others.append([r,ph])
        others.append([i/16.,j/16.])
others=n.array(others)
#E, rp = comp_mr(N, start_rp=None, doit=False);
E, rp = comp_mr(N, start_phs=3, others=others, maxiter=2000)
print E
print rp

# Show the antpos, UV, and UV*UV matrices
pos = _grid(rp,DIM) + _grid(others,DIM)
uv = aconv2(pos)
uv2 = aconv2(uv)
uv[0,0] = 0
uv2[0,0] = 0
DIM = pos.shape[0]
pos = a.img.recenter(pos, (DIM/2,DIM/2))[3*DIM/8:5*DIM/8,3*DIM/8:5*DIM/8]
uv = a.img.recenter(uv, (DIM/2,DIM/2))[DIM/4:3*DIM/4,DIM/4:3*DIM/4]
uv2 = a.img.recenter(uv2, (DIM/2,DIM/2))
p.subplot(131);
p.imshow(pos,origin='lower',extent=(-75,75,-75,75),
    vmin=vmin, vmax=vmax, cmap=cmap1, interpolation='nearest')
p.subplot(132);
p.imshow(uv, origin='lower',extent=(-150,150,-150,150),
    vmin=vmin, vmax=vmax, cmap=cmap1, interpolation='nearest')
p.subplot(133);
p.imshow(uv2,origin='lower',extent=(-300,300,-300,300),
    cmap=cmap2, interpolation='nearest')
p.show()

