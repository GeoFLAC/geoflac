#!/usr/bin/env python
'''Interpolate flac fields to a high resolution regular grid suitable for plotting

Usage: flac_interpolate.py frame field

frame (starts from 1)
field can be either 'phase', 'temperature', 'aps', 'density', 'eII', 'sII',
                    'sxx', 'szz', 'sxz', 'srII', 'pres', 'diss', 'visc'
'''
import sys, os
import numpy as np
import flac


# domain bounds (in km)
xmin = 300
xmax = 700
zmin = -150
zmax = 10

# grid resolution (in km)
dx = .8
dz = .2

# extra padding space for interpolation
padding = 12


def excluding(x, z, f, xmin, xmax, zmin, zmax):
    '''Excluding nodes of (x, z, f) that are outside the domain bounds'''
    # bool * bool is element-wise logical AND
    ind = (xmin <= x) * (x <= xmax) * (zmin <= z) * (z <= zmax)
    return x[ind], z[ind], f[ind]


def clip_topo(x, z, f, x0, z0):
    '''Setting f=NaN for nodes above (x0, z0).
    x, z, f x0, z0 are 2d arrays.'''
    xx0 = x0[:,0]
    zz0 = z0[:,0]
    xx = x[:,0]
    zz = np.interp(xx, xx0, zz0)

    for i in range(len(xx)):
        ind = z[i,:] > zz[i]
        f[i,ind] = np.nan
    return np.ma.masked_array(f, mask=np.isnan(f))


def interpolate(frame, field):
    # read coordinate of nodes
    fl = flac.Flac()
    xx, zz = fl.read_mesh(frame)
    x, z = flac.make_uniform_grid(xmin, xmax, zmin, zmax, dx, dz)

    if field == 'phase':
        ## phase
        # read marker location, age and phase
        mx, mz, mage, mphase, mid = fl.read_markers(frame)
        mx, mz, mphase = excluding(mx, mz, mphase, xmin-padding, xmax+padding, zmin-padding, zmax+padding)
        ph = flac.nearest_neighbor_interpolation2d(mx, mz, mphase, x, z)
        f = ph.astype(np.float32)
    elif field in ('aps', 'density', 'eII', 'sII',
                   'sxx', 'szz', 'sxz', 'srII', 'pres', 'diss', 'visc'):
        # read field on elements
        cf = getattr(fl, 'read_'+field)(frame)
        cx, cz = flac.elem_coord(xx, zz)
        cx, cz, cf = excluding(cx, cz, cf, xmin-padding, xmax+padding, zmin-padding, zmax+padding)
        f = flac.gaussian_interpolation2d(cx, cz, cf, x, z)
    elif field in ('temperature', ):
        # nodal field is interpolated by scipy bilinear interpolation
        # read vector field on nodes
        cf = getattr(fl, 'read_'+field)(frame)
        cx, cz, cf = excluding(xx, zz, cf, xmin-padding, xmax+padding, zmin-padding, zmax+padding)
        f = flac.bilinear_interpolation2d(cx, cz, cf, x, z)
    elif field in ('velo'):
        # read vector field on nodes
        cf = getattr(fl, 'read_'+field)(frame)
        cfx = cf[0]
        cx, cz, cfx = excluding(x, z, cfx, xmin-padding, xmax+padding, zmin-padding, zmax+padding)
        fx = flac.bilinear_interpolation2d(cx, cz, cfx, x, z)
        cfz = cf[1]
        cx, cz, cfz = excluding(x, z, cfz, xmin-padding, xmax+padding, zmin-padding, zmax+padding)
        fz = flac.bilinear_interpolation2d(cx, cz, cfz, x, z)
    else:
        raise RuntimeError('unknown field %s' % field)

    f = clip_topo(x, z, f, xx, zz)
    return x, z, f


if __name__ == '__main__':

    if len(sys.argv) < 3:
        print __doc__
        exit(1)

    frame = int(sys.argv[1])
    field = sys.argv[2].lower()
    xx, zz, f = interpolate(frame, field)
    flac.printing(xx, zz, f)
