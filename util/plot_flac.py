#!/usr/bin/env python3
'''Plot FLAC model outputs.
'''

import sys
import os
import glob
import argparse
import numpy as np
import flac
import flac_interpolate as fi
import flac_gravity

def find_trench_index(z):
    '''Returns the i index of trench location.'''
    zz = z[:,0]
    # the highest point defines the forearc
    imax = zz.argmax()
    # the trench is the lowest point west of forearc
    i = zz[:imax].argmin()
    return i


def plot_topo(frame, output_file=None):
    if output_file:
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    fl = flac.Flac()
    x, z = fl.read_mesh(frame)
    nx, nz = fl.nx, fl.nz

    xx = x[:,0]
    topo = z[:,0]

    sizeoffloat = 4
    offset = (frame-1) * sizeoffloat * nx
    if not os.path.exists('dtopo.0'):
        print("Error: dtopo.0 not found in the current directory.")
        sys.exit(1)
        
    with open('dtopo.0', 'rb') as f:
        f.seek(offset)
        erosion_rate = -np.fromfile(f, dtype=np.single, count=nx)

    fig, ax1 = plt.subplots()
    ax1.plot(xx, topo, 'b-')
    ax1.set_xlabel('x (km)')
    ax1.set_ylabel('topo (km)', color='b')
    ax1.tick_params('y', colors='b')

    ax2 = ax1.twinx()
    ax2.plot(xx, erosion_rate, 'r--')
    ax2.set_ylabel('erosion rate (mm/yr)', color='r')
    ax2.tick_params('y', colors='r')

    if output_file:
        plt.savefig(output_file, dpi=150)
        print(f"Saved topography plot to {output_file}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description='Plot FLAC model outputs.')
    parser.add_argument('frame', type=int, help='Output frame number')
    parser.add_argument('-m', '--mode', choices=['phase', 'strain_rate', 'stress', 'all', 'topo'], default='phase',
                        help="Plot mode: 'phase' (plot3), 'strain_rate' (plot4), 'stress' (plot5), 'all' stacked (plot6), or 'topo' (plot_topo)")
    parser.add_argument('--left', type=float, default=-200, help='Left domain bound in km (relative to trench)')
    parser.add_argument('--right', type=float, default=250, help='Right domain bound in km (relative to trench)')
    parser.add_argument('--up', type=float, default=10, help='Upper domain bound in km')
    parser.add_argument('--down', type=float, default=-200, help='Lower domain bound in km')
    parser.add_argument('--dx', type=float, default=None, help='Grid spacing x (km)')
    parser.add_argument('--dz', type=float, default=None, help='Grid spacing z (km)')
    parser.add_argument('-o', '--output', type=str, default=None, help='Custom output filename')

    args = parser.parse_args()

    if args.mode == 'topo':
        plot_topo(args.frame, args.output)
        return

    # GMT-based plotting modes
    frame = args.frame
    left = args.left
    right = args.right
    up = args.up
    down = args.down

    # Resolve default resolutions
    dx = args.dx
    dz = args.dz
    if dx is None:
        if args.mode in ('phase', 'all'):
            dx = 1.5
        else:
            dx = 4.0
    if dz is None:
        if args.mode in ('phase', 'all'):
            dz = 0.6
        elif args.mode == 'strain_rate':
            dz = 2.5
        elif args.mode == 'stress':
            dz = 1.5

    # Locate CPT files relative to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    phcpt = os.path.join(script_dir, 'phase15.cpt')
    strain_ratecpt = os.path.join(script_dir, 'strain_rate.cpt')
    stresscpt = os.path.join(script_dir, 'stress.cpt')

    fl = flac.Flac()
    x, z = fl.read_mesh(frame)
    itrench = find_trench_index(z)
    xtrench = x[itrench,0]

    xmin = xtrench + left
    xmax = xtrench + right
    zmin = down
    zmax = up

    # Read topography and gravity at uniform spacing
    px, topo, topomod, gravity = flac_gravity.compute_gravity(frame, fill_basin=True, erode_forearc=True, ref_grav="right")
    px *= 1e-3
    topo *= 1e-3
    topomod *= 1e-3
    gravity *= 1e5

    gfile = 'topo-grav.%d' % frame
    with open(gfile, 'w') as f:
        flac.printing(px, topo, gravity, topomod, stream=f)

    # temperature data file
    tfile = 'intp3-T.%d' % frame
    if not os.path.exists(tfile):
        T = fl.read_temperature(frame)
        with open(tfile, 'w') as f:
            f.write('%d %d\n' % x.shape)
            flac.printing(x, z, T, stream=f)

    tgrd = 'temperature3.%d.grd' % frame
    if not os.path.exists(tgrd):
        cmd = 'tail -n +2 %(tfile)s | gmt surface -G%(tgrd)s -Ll0 -I%(dx)f/%(dz)f -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f' % locals()
        os.system(cmd)

    model = os.path.split(os.getcwd())[-1]
    
    # Setup plotting grids
    if args.mode in ('phase', 'all'):
        phasefile = 'intp3-phase.%d' % frame
        if not os.path.exists(phasefile):
            fi.xmin = xtrench + left
            fi.xmax = xtrench + right
            fi.zmin = down
            fi.zmax = up
            fi.dx = dx
            fi.dz = dz
            xx, zz, ph = fi.interpolate(frame, 'phase')
            with open(phasefile, 'w') as f:
                f.write('%d %d\n' % xx.shape)
                flac.printing(xx, zz, ph, stream=f)
        phgrd = 'phase3.%d.grd' % frame
        if not os.path.exists(phgrd):
            cmd = 'tail -n +2 %(phasefile)s | gmt xyz2grd -G%(phgrd)s -I%(dx)f/%(dz)f -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f' % locals()
            os.system(cmd)

    if args.mode in ('strain_rate', 'all'):
        strain_ratefile = 'intp3-srII.%d' % frame
        if not os.path.exists(strain_ratefile):
            fi.xmin = xtrench + left
            fi.xmax = xtrench + right
            fi.zmin = down
            fi.zmax = up
            fi.dx = dx
            fi.dz = dz
            xx, zz, sr = fi.interpolate(frame, 'srII')
            with open(strain_ratefile, 'w') as f:
                f.write('%d %d\n' % xx.shape)
                flac.printing(xx, zz, sr, stream=f)
        strrgrd = 'strain_rate3.%d.grd' % frame
        if not os.path.exists(strrgrd):
            cmd = 'tail -n +2 %(strain_ratefile)s | gmt xyz2grd -G%(strrgrd)s -I%(dx)f/%(dz)f -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f' % locals()
            os.system(cmd)

    if args.mode in ('stress', 'all'):
        stressfile = 'intp3-sII.%d' % frame
        if not os.path.exists(stressfile):
            fi.xmin = xtrench + left
            fi.xmax = xtrench + right
            fi.zmin = down
            fi.zmax = up
            fi.dx = dx
            fi.dz = dz
            xx, zz, s = fi.interpolate(frame, 'sII')
            with open(stressfile, 'w') as f:
                f.write('%d %d\n' % xx.shape)
                flac.printing(xx, zz, s, stream=f)
        strgrd = 'stress3.%d.grd' % frame
        if not os.path.exists(strgrd):
            cmd = 'tail -n +2 %(stressfile)s | gmt xyz2grd -G%(strgrd)s -I%(dx)f/%(dz)f -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f' % locals()
            os.system(cmd)

    # GMT script generation and execution
    aspect_ratio = float(up - down) / (right - left)
    width = 6.5
    height = width * aspect_ratio
    height2 = 1.0

    gravgridsize = 50
    gmin = int(gravity.min() / gravgridsize - 1) * gravgridsize
    gmax = int(gravity.max() / gravgridsize + 1) * gravgridsize
    gravann = max(abs(gmin), abs(gmax))

    topogridsize = 2
    tpmin = int(topo.min() / topogridsize - 1) * topogridsize
    tpmax = int(topo.max() / topogridsize + 1) * topogridsize
    topoann = max(abs(tpmin), abs(tpmax))
    
    cint = 200

    if args.mode == 'all':
        psfile = args.output or 'result6.%d.ps' % frame
        pngfile = os.path.splitext(psfile)[0] + '.png'
        shiftz = height + 1

        cmd = '''
rm -f .gmtcommands* .gmtdefaults*
gmt gmtset MEASURE_UNIT = inch
gmt gmtset LABEL_FONT_SIZE=14 ANNOT_FONT_SIZE_PRIMARY=10
gmt gmtset PAPER_MEDIA A3 

gmt psbasemap -JX%(width)f/%(height)f -Ba100f10/a20f5::WSne:."stress": -R%(left)f/%(right)f/%(zmin)f/%(zmax)f -X0.9 -P -K > %(psfile)s
gmt grdimage %(strgrd)s -C%(stresscpt)s -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f -J -P -O -K >> %(psfile)s
gmt grdcontour %(tgrd)s -C%(cint)f -A200 -W1p -J -R -P -K -O >> %(psfile)s

gmt psbasemap -JX%(width)f/%(height)f -Ba100f10/a20f5::WSne:."strain rate": -R%(left)f/%(right)f/%(zmin)f/%(zmax)f -Y4.5 -P -O -K >> %(psfile)s
gmt grdimage %(strrgrd)s -C%(strain_ratecpt)s -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f -J -P -O -K >> %(psfile)s
gmt grdcontour %(tgrd)s -C%(cint)f -A200 -W1p -J -R -P -K -O >> %(psfile)s

gmt psbasemap -JX%(width)f/%(height)f -Ba100f10/a20f5::WSne:."phase": -R%(left)f/%(right)f/%(zmin)f/%(zmax)f -Y4.5 -P -O -K >> %(psfile)s
gmt grdimage %(phgrd)s -C%(phcpt)s -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f -J -P -O -K >> %(psfile)s
gmt grdcontour %(tgrd)s -C%(cint)f -A200 -W1p -J -R -P -K -O >> %(psfile)s

gmt psxy -JX%(width)f/%(height2)f -R%(xmin)f/%(xmax)f/-1/1 -A -Wthin,0,. -Y%(shiftz)f -P -K -O <<END>> %(psfile)s
%(xmin)f 0
%(xmax)f 0
END

cut -f1,3 %(gfile)s | gmt psxy -JX%(width)f/%(height2)f -R%(xmin)f/%(xmax)f/-%(gravann)f/%(gravann)f -B/a%(gravann)ff%(gravgridsize)f:"@;red;mGal@;;":E -A -Wthin,red -P -K -O >> %(psfile)s
cut -f1,2 %(gfile)s | gmt psxy -JX%(width)f/%(height2)f -R%(xmin)f/%(xmax)f/-%(topoann)f/%(topoann)f -B/a%(topoann)ff%(topogridsize)f:"@;blue;km@;;":W -A -Wthin,blue -P -K -O >> %(psfile)s
cut -f1,4 %(gfile)s | gmt psxy -J -R -A -Wthin,blue,-- -P -K -O >> %(psfile)s

gmt psscale -C%(phcpt)s -D6.75/-2.5/2.5/0.3 -Ba1f0.5 -V -O -K >> %(psfile)s
gmt psscale -C%(strain_ratecpt)s -D6.75/-7/2.5/0.3 -Ba1f0.5 -V -O -K >> %(psfile)s
gmt psscale -C%(stresscpt)s -D6.75/-11.5/2.5/0.3 -Ba2e8f5e7 -V -O -K >> %(psfile)s

echo %(xmin)f %(topogridsize)d 14 0 1 LB "Model=%(model)s   Frame=%(frame)d   Trench location=%(xtrench).3f km" | gmt pstext -D0/1 -N -J -R -P -O >> %(psfile)s
convert -density 150 %(psfile)s %(pngfile)s
''' % locals()
        os.system(cmd)
        print(f"Generated 3-panel plot: {pngfile}")

    else:
        # Single panel modes: phase, strain_rate, stress
        if args.mode == 'phase':
            cpt_file = phcpt
            grd_file = phgrd
            label = 'phase'
            psfile = args.output or 'result3.%d.ps' % frame
        elif args.mode == 'strain_rate':
            cpt_file = strain_ratecpt
            grd_file = strrgrd
            label = 'strain_rate'
            psfile = args.output or 'result4.%d.ps' % frame
        elif args.mode == 'stress':
            cpt_file = stresscpt
            grd_file = strgrd
            label = 'stress'
            psfile = args.output or 'result5.%d.ps' % frame

        pngfile = os.path.splitext(psfile)[0] + '.png'
        shiftz = height + 0.3

        cmd = '''
rm -f .gmtcommands* .gmtdefaults*
gmt gmtset MEASURE_UNIT = inch
gmt gmtset LABEL_FONT_SIZE=14 ANNOT_FONT_SIZE_PRIMARY=10

gmt psbasemap -JX%(width)f/%(height)f -Ba100f10/a20f5::WSne -R%(left)f/%(right)f/%(zmin)f/%(zmax)f -X0.9 -Y4 -P -K > %(psfile)s
gmt grdimage %(grd_file)s -C%(cpt_file)s -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f -J -P -O -K >> %(psfile)s
gmt grdcontour %(tgrd)s -C%(cint)f -A200 -W1p -J -R -P -K -O >> %(psfile)s

gmt psxy -JX%(width)f/%(height2)f -R%(xmin)f/%(xmax)f/-1/1 -A -Wthin,0,. -Y%(shiftz)f -P -K -O <<END>> %(psfile)s
%(xmin)f 0
%(xmax)f 0
END

cut -f1,3 %(gfile)s | gmt psxy -JX%(width)f/%(height2)f -R%(xmin)f/%(xmax)f/-%(gravann)f/%(gravann)f -B/a%(gravann)ff%(gravgridsize)f:"@;red;mGal@;;":E -A -Wthin,red -P -K -O >> %(psfile)s
cut -f1,2 %(gfile)s | gmt psxy -JX%(width)f/%(height2)f -R%(xmin)f/%(xmax)f/-%(topoann)f/%(topoann)f -B/a%(topoann)ff%(topogridsize)f:"@;blue;km@;;":W -A -Wthin,blue -P -K -O >> %(psfile)s
cut -f1,4 %(gfile)s | gmt psxy -J -R -A -Wthin,blue,-- -P -K -O >> %(psfile)s

echo %(xmin)f %(topogridsize)d 14 0 1 LB "Model=%(model)s   Frame=%(frame)d   Trench location=%(xtrench).3f km" | gmt pstext -D0/1 -N -J -R -P -O >> %(psfile)s
convert -density 150 %(psfile)s %(pngfile)s
''' % locals()
        os.system(cmd)
        print(f"Generated single-panel plot: {pngfile}")


if __name__ == '__main__':
    main()
