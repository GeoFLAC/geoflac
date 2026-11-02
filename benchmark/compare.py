#!/usr/bin/python3
from __future__ import print_function
import sys, os
import numpy as np
import flac

class Stuff():
    pass


def read_data(fl, frame):
    stuff = Stuff()
    stuff.T = fl.read_temperature(frame)
    stuff.x, stuff.z = fl.read_mesh(frame)
    stuff.vx, stuff.vz = fl.read_vel(frame)

    stuff.aps = fl.read_aps(frame)
    stuff.sxx = fl.read_sxx(frame)
    stuff.szz = fl.read_szz(frame)
    stuff.sxz = fl.read_sxz(frame)

    stuff.m_x, stuff.m_z, stuff.m_age, stuff.m_phase, stuff.m_ID, *unused = \
        fl.read_markers(frame)
    return stuff


def reldiff(oldf, newf):
    m = np.abs(oldf).max()
    diff = np.abs(newf - oldf)
    return diff.max()/m, diff.std()/m


def compare(old, new):
    max, sigma = reldiff(old.T, new.T)
    print('  Temperature:\t', max, sigma)

    max, sigma = reldiff(old.x, new.x)
    print('  X coordinate:\t', max, sigma)

    max, sigma = reldiff(old.z, new.z)
    print('  Z coordinate:\t', max, sigma)

    max, sigma = reldiff(old.vx, new.vx)
    print('  X velocity:\t', max, sigma)

    max, sigma = reldiff(old.vz, new.vz)
    print('  Z velocity:\t', max, sigma)

    max, sigma = reldiff(old.aps, new.aps)
    print('  Pl. strain:\t', max, sigma)

    max, sigma = reldiff(old.sxx, new.sxx)
    print('  Stress xx:\t', max, sigma)

    max, sigma = reldiff(old.szz, new.szz)
    print('  Stress zz:\t', max, sigma)

    max, sigma = reldiff(old.sxz, new.sxz)
    print('  Stress xz:\t', max, sigma)

    max, sigma = reldiff(old.m_x, new.m_x)
    print('  Marker X:\t', max, sigma)

    max, sigma = reldiff(old.m_z, new.m_z)
    print('  Marker Z:\t', max, sigma)

    max, sigma = reldiff(old.m_phase, new.m_phase)
    print('  Marker Phase:\t', max, sigma)

    return


olddir = sys.argv[1]
frame = int(sys.argv[2])

curdir = os.getcwd()
newdir = curdir

# name holder
old = 0
new = 0

try:
    # read old and new results
    flo = flac.Flac(wdir=olddir)
    old = read_data(flo, frame)

    fln = flac.Flac(wdir=newdir)
    new = read_data(fln, frame)

    # compare results
    print()
    print('Relative difference (max, stddev) of frame =', frame,
          ' step =', int(flo.steps[frame]))
    compare(old, new)

finally:
    pass
