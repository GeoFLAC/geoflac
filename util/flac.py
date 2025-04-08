#!/usr/bin/python3

from __future__ import print_function
import sys, os, zlib, base64, glob

try:
    import numpy as np
except ImportError:
    print('Error: Failed importing "numpy" module.')
    print('Please install the module and add it to PYTHONPATH environment variable.')
    sys.exit(1)

try:
    # for interpolation
    from scipy import interpolate
except ImportError:
    pass


# the precision of saved file
doubleprecision = False
sizeofint = 4
if doubleprecision:
    sizeoffloat = 8
    default_dtype = np.double
else:
    sizeoffloat = 4
    default_dtype = np.single



class Flac(object):
    '''Read Flac data file. Most data are 2D arrays with shape (nx, nz) or (nex, nez).
    '''

    def __init__(self, swap_endian=False):
        self.swap_endian = swap_endian
        self.reload()
        return


    def reload(self):
        # read record number
        tmp = np.fromfile('_contents.0', sep=' ')
        tmp.shape = (-1, 3)
        self.frames = tmp[:,0]
        self.steps = tmp[:,1]
        self.time = tmp[:,2]
        self.nrec = len(self.time)

        # number of elements in x and z
        nex, nez = np.fromfile('nxnz.0', sep=' ', dtype=int)
        # number of nodes in x and z
        self.nx, self.nz = nex+1, nez+1
        self.nnodes = self.nx * self.nz
        self.nelements = nex * nez

        return


    def read_mesh(self, frame):
        columns = 2
        f = open('mesh.0')
        offset = (frame-1) * columns * self.nnodes * sizeoffloat
        f.seek(offset)
        x, z = self._read_data(f, columns)
        self._reshape_nodal_fields(x, z)
        return x, z


    def read_vel(self, frame):
        columns = 2
        f = open('vel.0')
        offset = (frame-1) * columns * self.nnodes * sizeoffloat
        f.seek(offset)
        vx, vz = self._read_data(f, columns)
        self._reshape_nodal_fields(vx, vz)
        return vx, vz


    def read_temperature(self, frame):
        columns = 1
        f = open('temperature.0')
        offset = (frame-1) * columns * self.nnodes * sizeoffloat
        f.seek(offset)
        T = self._read_data(f, columns)
        self._reshape_nodal_fields(T)
        return T


    def read_aps(self, frame):
        columns = 1
        f = open('aps.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        aps = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(aps)
        return aps


    def read_density(self, frame):
        columns = 1
        f = open('density.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        density = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(density)
        return density


    def read_area(self, frame):
        columns = 1
        f = open('area.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        area = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(area)
        return area


    def read_area(self, frame):
        columns = 1
        f = open('area.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        area = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(area)
        return area


    def read_strain(self, frame):
        columns = 1
        f = open('exx.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        exx = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(exx)

        f = open('ezz.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        ezz = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(ezz)

        f = open('exz.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        exz = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(exz)
        return exx, ezz, exz


    def read_eII(self, frame):
        columns = 1
        f = open('eII.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        eII = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(eII)
        return eII


    def read_sII(self, frame):
        columns = 1
        f = open('sII.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        sII = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(sII)
        return sII


    def read_sxx(self, frame):
        columns = 1
        f = open('sxx.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        sxx = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(sxx)
        return sxx


    def read_sxz(self, frame):
        columns = 1
        f = open('sxz.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        sxz = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(sxz)
        return sxz


    def read_syy(self, frame):
        columns = 1
        f = open('syy.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        syy = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(syy)
        return syy


    def read_szz(self, frame):
        columns = 1
        f = open('szz.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        szz = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(szz)
        return szz


    def read_srII(self, frame):
        columns = 1
        f = open('srII.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        srII = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(srII)
        return srII


    def read_pres(self, frame):
        columns = 1
        f = open('pres.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        pres = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(pres)
        return pres


    def read_fmelt(self, frame):
        columns = 1
        f = open('fmelt.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        fmelt = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(fmelt)
        return fmelt


    def read_fmagma(self, frame):
        columns = 1
        f = open('fmagma.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        fmagma = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(fmagma)
        return fmagma


    def read_eff_melt(self, frame):
        columns = 1
        f = open('eff_melt.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        eff_melt = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(eff_melt)
        return eff_melt


    def read_diss(self, frame):
        columns = 1
        f = open('diss.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        diss = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(diss)
        return diss


    def read_visc(self, frame):
        columns = 1
        f = open('visc.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        visc = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(visc)
        return visc


    def read_phase(self, frame):
        columns = 1
        f = open('phase.0')
        offset = (frame-1) * columns * self.nelements * sizeofint
        f.seek(offset)
        phase = self._read_data(f, columns, count=self.nelements, dtype=np.int32)
        self._reshape_elemental_fields(phase)
        return phase


    def read_markers(self, frame):
        # read tracer size
        tmp = np.fromfile('_markers.0', sep=' ')
        tmp.shape = (-1, 4)
        n = int(tmp[frame-1,2])

        suffix = '.%06d.0' % frame
        f2 = open('marker2' + suffix)
        dead = self._read_data(f2, count=n, dtype=np.int32)
        tmp = self._read_data(f2, count=n, dtype=np.int32)
        phase = self._remove_dead_markers(tmp, dead)
        tmp = self._read_data(f2, count=n, dtype=np.int32)
        ntriag = self._remove_dead_markers(tmp, dead)
        f2.close()

        f1 = open('marker1' + suffix)
        tmp = self._read_data(f1, count=n)
        x = self._remove_dead_markers(tmp, dead)

        tmp = self._read_data(f1, count=n)
        z = self._remove_dead_markers(tmp, dead)

        tmp = self._read_data(f1, count=n)
        age = self._remove_dead_markers(tmp, dead)

        tmp = self._read_data(f1, count=n)
        a1 = self._remove_dead_markers(tmp, dead)

        tmp = self._read_data(f1, count=n)
        a2 = self._remove_dead_markers(tmp, dead)

        tmp = np.arange(1, n+1)
        ID = self._remove_dead_markers(tmp, dead)
        f1.close()
        return x, z, age, phase, ID, a1, a2, ntriag


    def read_tracers(self):
        # read tracer size
        tmp = np.fromfile('_tracers.0', sep=' ')
        tmp.shape = (-1, 4)
        ntracerrec = tmp.shape[0]
        ntracers = int(tmp[0,1])
        n = ntracerrec * ntracers

        time = tmp[:,3]

        x = self._read_data('tracerx.0', count=n)
        x.shape = (ntracerrec, ntracers)

        z = self._read_data('tracerz.0', count=n)
        z.shape = (ntracerrec, ntracers)

        T = self._read_data('tracert.0', count=n)
        T.shape = (ntracerrec, ntracers)

        p = self._read_data('tracerp.0', count=n)
        p.shape = (ntracerrec, ntracers)

        sxx = self._read_data('tracersxx.0', count=n)
        sxx.shape = (ntracerrec, ntracers)

        szz = self._read_data('tracerszz.0', count=n)
        szz.shape = (ntracerrec, ntracers)

        sxz = self._read_data('tracersxz.0', count=n)
        sxz.shape = (ntracerrec, ntracers)

        e = self._read_data('tracere.0', count=n)
        e.shape = (ntracerrec, ntracers)

        edot = self._read_data('traceredot.0', count=n)
        edot.shape = (ntracerrec, ntracers)

        phase = self._read_data('tracerph.0', count=n, dtype=np.int32)
        phase.shape = (ntracerrec, ntracers)

        return x, z, T, p, sxx, szz, sxz, e, edot, phase


    def _read_data(self, fileobj, columns=1,
                   count=None, dtype=None):
        '''Read data from a file-like object 'fileobj'.

        The 'dtype' specifies the storage type, default to single precision
        float.
        '''

        # number of nodes
        if count is None:
            count = self.nnodes

        # total number of items
        n = columns * count

        if dtype is None:
            dtype = default_dtype

        result = np.fromfile(fileobj, dtype, n)
        if self.swap_endian:
            result.byteswap(True)

        if columns == 1:
            return result
        else:
            result.shape = (columns, -1)
            return tuple(result[i,:] for i in range(columns))


    def _reshape_nodal_fields(self, *argv):
        for x in argv:
            x.shape = (self.nx, self.nz)
        return


    def _reshape_elemental_fields(self, *argv):
        for x in argv:
            x.shape = (self.nx-1, self.nz-1)
        return


    def _remove_dead_markers(self, a, dead):
        b = a[dead==1]
        return b


class FlacFromVTK(object):
    '''Read flac data from vts or vtp files.
       Try to maintain the same API as Flac()
    '''

    def __init__(self, swap_endian=False):
        self.last_frame_read = None
        self.cached_vts = None

        allvts = glob.glob('flac.*.vts')
        allvts.sort()
        if allvts[0] != 'flac.000001.vts':
            print('Error: missing first vts frame!')
            sys.exit(1)
        self.nrec = len(allvts)
        self.frames = list(range(1, self.nrec+1))
        self.steps = np.zeros(self.nrec)
        self.time = np.zeros(self.nrec)
        self._read_vtk(1) # get grid size
        return


    def _read_vtk(self, frame):
        filename = 'flac.%06d.vts' % self.frames[frame]
        print('Reading from', filename)
        f = open(filename, 'r')
        d = f.readlines()
        for n, line in enumerate(d):
            d[n] = line.strip()
        self.last_frame_read = frame
        self.cached_vts = d

        # parse element size
        s = '<StructuredGrid WholeExtent="'
        if d[2][:len(s)] != s:
            print('Error: invalid vtk data')
            sys.exit(2)
        a = d[2][len(s):].split(' ')
        nex = int(a[1])
        nez = int(a[3])
        self.nx, self.nz = nex+1, nez+1
        self.nnodes = self.nx * self.nz
        self.nelements = nex * nez

        # parse time and steps
        s = '<DataArray type="Float32" Name="TIME" '
        i = 4
        if d[i][:len(s)] != s:
            print('Error: invalid vtk data')
            sys.exit(2)
        self.time[frame] = float(d[i+1])

        s = '<DataArray type="Float32" Name="CYCLE" '
        i = 7
        if d[i][:len(s)] != s:
            print('Error: invalid vtk data')
            sys.exit(2)
        self.steps[frame] = float(d[i+1])
        return d


    def _get_vtk_data(self, frame):
        if frame == self.last_frame_read:
            data = self.cached_vts
        else:
            data = self._read_vtk(frame)
        return data


    def _unpack_vtk(self, line):
        # 4 int32 (16 bytes) were base64-encoded, became 16*3/2 = 24 bytes
        header = base64.standard_b64decode(line[:24].encode('ascii'))
        data = line[24:].encode('ascii')
        if int.from_bytes(header[:4], 'little', signed=True) != 1:
            print('Error: invalid compressed vtk data')
            sys.exit(2)

        # size of original data
        n = int.from_bytes(header[4:8], 'little', signed=True)
        if int.from_bytes(header[8:12], 'little', signed=True) != n:
            print('Error: invalid compressed vtk data')
            sys.exit(2)

        # size of zlib compressed data
        m = int.from_bytes(header[12:], 'little', signed=True)

        a = zlib.decompress(base64.standard_b64decode(data), bufsize=n)
        return a


    def _locate_line(self, data, name, dtype="Float32"):
        s = '<DataArray type="%s" Name="%s"' % (dtype, name)
        for n, line in enumerate(data):
            if line[:len(s)] == s: break
        else:
            print('Error: reading data', s)
            sys.exit(1)

        line = self._unpack_vtk(data[n+1])
        return line


    def read_mesh(self, frame):
        data = self._get_vtk_data(frame)
        s = '<Points>'
        for n, line in enumerate(data):
            if line[:len(s)] == s: break
        else:
            print('Error: reading data', s)
            sys.exit(1)

        s = '<DataArray type="Float32"  NumberOfComponents="3" format="binary">'
        if data[n+1] != s:
            print('Error: reading data', s)
            sys.exit(1)

        a = self._unpack_vtk(data[n+2])
        a = np.frombuffer(a, dtype=np.float32).reshape(-1,3)
        x, z = a[:,0], a[:,1]
        x.shape = (self.nx, self.nz)
        z.shape = (self.nx, self.nz)
        return x.transpose(), z.transpose()


    def read_vel(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Velocity")
        a = np.frombuffer(a, dtype=np.float32).reshape(-1,3)
        vx, vz = a[:,0], a[:,1]
        vx.shape = (self.nx, self.nz)
        vz.shape = (self.nx, self.nz)
        return vx.transpose(), vz.transpose()


    def read_temperature(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Temperature")
        T = np.frombuffer(a, dtype=np.float32)
        T.shape = (self.nx-1, self.nz-1)
        return T.transpose()


    def read_aps(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Plastic strain")
        aps = np.frombuffer(a, dtype=np.float32)
        aps.shape = (self.nx-1, self.nz-1)
        return aps.transpose()


    def read_density(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Density")
        density = np.frombuffer(a, dtype=np.float32)
        density.shape = (self.nx-1, self.nz-1)
        return density.transpose()


    def read_area(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Area")
        area = np.frombuffer(a, dtype=np.float32)
        area.shape = (self.nx-1, self.nz-1)
        return area.transpose()


    def read_strain(self, frame):
        return NotImplemented


    def read_eII(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "eII")
        eII = np.frombuffer(a, dtype=np.float32)
        eII.shape = (self.nx-1, self.nz-1)
        return eII.transpose()


    def read_sII(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "sII")
        sII = np.frombuffer(a, dtype=np.float32)
        sII.shape = (self.nx-1, self.nz-1)
        return sII.transpose()


    def read_sxx(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Sxx")
        sxx = np.frombuffer(a, dtype=np.float32)
        sxx.shape = (self.nx-1, self.nz-1)
        return sxx.transpose()


    def read_sxz(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Sxz")
        sxz = np.frombuffer(a, dtype=np.float32)
        sxz.shape = (self.nx-1, self.nz-1)
        return sxz.transpose()


    def read_szz(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Szz")
        szz = np.frombuffer(a, dtype=np.float32)
        szz.shape = (self.nx-1, self.nz-1)
        return szz.transpose()


    def read_srII(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Strain rate")
        srII = np.frombuffer(a, dtype=np.float32)
        srII.shape = (self.nx-1, self.nz-1)
        return srII.transpose()


    def read_pres(self, frame):
        return NotImplemented


    def read_fmelt(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Melt fraction")
        fmelt = np.frombuffer(a, dtype=np.float32)
        fmelt.shape = (self.nx-1, self.nz-1)
        return fmelt.transpose()


    def read_fmagma(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Magma fraction")
        fmagma = np.frombuffer(a, dtype=np.float32)
        fmagma.shape = (self.nx-1, self.nz-1)
        return fmagma.transpose()


    def read_eff_melt(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Eff melt fraction")
        eff_melt = np.frombuffer(a, dtype=np.float32)
        eff_melt.shape = (self.nx-1, self.nz-1)
        return eff_melt.transpose()


    def read_diss(self, frame):
        return NotImplemented


    def read_visc(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Viscosity")
        visc = np.frombuffer(a, dtype=np.float32)
        visc.shape = (self.nx-1, self.nz-1)
        return visc.transpose()


    def read_phase(self, frame):
        data = self._get_vtk_data(frame)
        a = self._locate_line(data, "Phase", "Int32")
        phase = np.frombuffer(a, dtype=np.int32)
        phase.shape = (self.nx-1, self.nz-1)
        return phase.transpose()


    def read_markers(self, frame):
        filename = 'flacmarker.%06d.vtp' % self.frames[frame]
        print('Reading from', filename)
        f = open(filename, 'r')
        d = f.readlines()
        for n, line in enumerate(d):
            d[n] = line.strip()
        data = d

        # # parse # of markers (not dead)
        # s = '<Piece NumberOfPoints="'
        # if d[3][:len(s)] != s:
        #     print('Error: invalid vtk data')
        #     sys.exit(2)
        # n = int(d[3][len(s):-2])

        # read point data
        a = self._locate_line(data, "age")
        age = np.frombuffer(a, dtype=np.float32)

        a = self._locate_line(data, "phase", dtype="Int32")
        phase = np.frombuffer(a, dtype=np.int32)

        a = self._locate_line(data, "ID", dtype="Int32")
        ID = np.frombuffer(a, dtype=np.int32)

        a = self._locate_line(data, "a1")
        a1 = np.frombuffer(a, dtype=np.float32)

        a = self._locate_line(data, "a2")
        a2 = np.frombuffer(a, dtype=np.float32)

        a = self._locate_line(data, "ntriag", dtype="Int32")
        ntriag = np.frombuffer(a, dtype=np.int32)

        # read point coordinate
        s = '<Points>'
        for n, line in enumerate(data):
            if line[:len(s)] == s: break
        else:
            print('Error: reading data', s)
            sys.exit(1)

        s = '<DataArray type="Float32"  NumberOfComponents="3" format="binary">'
        if data[n+1] != s:
            print('Error: reading data', s)
            sys.exit(1)

        a = self._unpack_vtk(data[n+2])
        a = np.frombuffer(a, dtype=np.float32).reshape(-1,3)
        x, z = a[:,0], a[:,1]

        # True if not using thermochron
        if (True): return x, z, age, phase, ID, a1, a2, ntriag

        # Thermochronology from Chase Shyu's work
        '''
        a = self._locate_line(data, "temp")
        temp = np.frombuffer(a, dtype=np.float32)

        a = self._locate_line(data, "tempmax")
        tempmax = np.frombuffer(a, dtype=np.float32)

        a = self._locate_line(data, "cooling rate")
        coolingrate = np.frombuffer(a, dtype=np.float32)

        a = self._locate_line(data, "zft age")
        zft_age = np.frombuffer(a, dtype=np.float32)

        a = self._locate_line(data, "zft if", dtype="Int32")
        zft_if = np.frombuffer(a, dtype=np.int32)

        a = self._locate_line(data, "zft Temp.")
        zft_temp = np.frombuffer(a, dtype=np.float32)

        a = self._locate_line(data, "zhe age")
        zhe_age = np.frombuffer(a, dtype=np.float32)

        a = self._locate_line(data, "zhe if", dtype="Int32")
        zhe_if = np.frombuffer(a, dtype=np.int32)

        a = self._locate_line(data, "zhe Temp.")
        zhe_temp = np.frombuffer(a, dtype=np.float32)

        chronage = [zft_age, zhe_age]
        chronif = [zft_if, zhe_if]
        chrontemp = [zft_temp, zhe_temp]
        return x, z, age, temp, tempmax, coolingrate, phase, ID, \
               chronage, chronif, chrontemp
        '''

################################################################
### Replace Flac by FlacFromVTK if only vts files are available
################################################################
#
#if (not os.path.exists('_contents.0')) and os.path.exists('flac.000001.vts'):
#    Flac = FlacFromVTK
#
################################################################



def marker_ntriag2elem(ntriag, nz):
    '''Convert markers' ntriag to element number (i,j) and triangle number (k)
    nz is the # of nodes in z direction.'''
    k = ((ntriag - 1) % 2) + 1
    j = ((ntriag - k) // 2) % (nz - 1)
    i = (ntriag - k) // 2 // (nz - 1)
    return i, j, (k-1)


def marker_interpolate_elem(ntriag, nz, efield):
    '''Interpolation elemental field (e.g. stress) onto markers
    '''
    if efield.shape[1] != nz-1:
        raise ValueError('The array length in 2nd dimension is %d, excepting %d' % (nfield.shape[1], nz-1))

    i, j, k = marker_ntriag2elem(ntriag, nz)
    f = np.zeros(ntriag.shape, dtype=efield.dtype)
    f = efield[i[:], j[:]]
    return f


def marker_interpolate_node(ntriag, a1, a2, nz, nfield):
    '''Interpolation nodal field (e.g. temperature) onto markers
    '''
    if nfield.shape[1] != nz:
        raise ValueError('The array length in 2nd dimension is %d, excepting %d' % (nfield.shape[1], nz))

    f = np.zeros_like(a1)
    a3 = 1.0 - a1 - a2

    # Upper triangle:
    # 1 --- 3
    # |   /
    # |  /
    # | /
    # 2
    #
    # Lower triangle:
    #       1
    #     / |
    #    /  |
    #   /   |
    # 2 --- 3
    i, j, k = marker_ntriag2elem(ntriag, nz)
    u = (k == 0)  # uppper triangles
    l = (k == 1)  # lower triangles

    f[u] = (nfield[i[u]  , j[u]  ] * a1[u] +
            nfield[i[u]  , j[u]+1] * a2[u] +
            nfield[i[u]+1, j[u]  ] * a3[u])
    f[l] = (nfield[i[l]+1, j[l]  ] * a1[l] +
            nfield[i[l]  , j[l]+1] * a2[l] +
            nfield[i[l]+1, j[l]+1] * a3[l])
    return f


def elem_coord(x, z):
    '''Turning nodal coordinates to element coordinates'''
    cx = (x[:-1, :-1] + x[:-1, 1:] + x[1:, :-1] + x[1:, 1:]) / 4
    cz = (z[:-1, :-1] + z[:-1, 1:] + z[1:, :-1] + z[1:, 1:]) / 4
    return cx, cz


def make_uniform_grid(xmin, xmax, zmin, zmax, dx, dz):
    # grid size
    nx = (xmax - xmin) / dx + 1
    nz = (zmax - zmin) / dz + 1

    # generate uniform grid
    xx = np.linspace(xmin, xmax, nx)
    zz = np.linspace(zmin, zmax, nz)

    # the order of argument ensures the shape of arrays is (nx, nz)
    z, x = np.meshgrid(zz, xx)
    return x, z


def nearest_neighbor_interpolation2d(x0, z0, f0, x, z):
    '''Interpolating field f0, which is defined on (x0, z0)
    to a new grid (x, z) using nearest neighbor method'''

    if x0.shape != z0.shape:
        raise Exception('x0 and z0 arrays have different shape')

    if x0.shape != f0.shape:
        raise Exception('x0 and f0 arrays have different shape')

    if x.shape != z.shape:
        raise Exception('x and z arrays have different shape')

    points = np.vstack((x0.flat, z0.flat)).T
    f = interpolate.griddata(points, f0.flat, (x, z), method='nearest')

    return f


def bilinear_interpolation2d(x0, z0, f0, x, z):
    '''Interpolating field f0, which is defined on (x0, z0)
    to a new grid (x, z) using bilinear method'''

    if x0.shape != z0.shape:
        raise Exception('x0 and z0 arrays have different shape')

    if x0.shape != f0.shape:
        raise Exception('x0 and f0 arrays have different shape')

    if x.shape != z.shape:
        raise Exception('x and z arrays have different shape')

    points = np.vstack((x0.flat, z0.flat)).T
    f = interpolate.griddata(points, f0.flat, (x, z), method='linear')

    return f


def gaussian_interpolation2d(x0, z0, f0, x, z):
    '''Interpolating field f0, which is defined on (x0, z0)
    to a new grid (x, z) using weighted gaussian method'''

    if x0.shape != z0.shape:
        raise Exception('x0 and z0 arrays have different shape')

    if x0.shape != f0.shape:
        raise Exception('x0 and f0 arrays have different shape')

    if x.shape != z.shape:
        raise Exception('x and z arrays have different shape')

    # using 1d index for x0, z0, f0
    x0 = x0.flat
    z0 = z0.flat
    f0 = f0.flat

    dx = 1.5 * (x[1,0] - x[0,0])
    dz = 1.5 * (z[0,1] - z[0,0])

    f = np.zeros(x.shape)
    g = np.zeros(x.shape)
    for i in range(len(x0)):
        weight = np.exp(-((x - x0[i]) / dx)**2 - ((z - z0[i]) / dz)**2)
        f += weight * f0[i]
        g += weight

    return f / g


def printing(*args, **kwd):
    '''printing(arg1, arg2, ..., stream=None)

    stream: None -- output to stdout
            filename -- output to file
            file-like object -- output to the object
    '''

    stream = kwd.get('stream')
    if stream == None:
        stream = sys.stdout
    elif isinstance(stream, (str, unicode)):
        filename = stream
        stream = open(filename, 'w')

    # flatten numpy arrays
    try:
        args = tuple(x.flat for x in args)
    except AttributeError:
        pass

    narg = len(args)
    fmt = '%.15e' + '\t%.15e'*(narg-1)

    for items in zip(*args):
        print(fmt % tuple(items), file=stream)
    return



## Before running this module, 'cd' to a directory containing the flac data.
if __name__ == '__main__':
    import sys

    if (len(sys.argv) < 3):
        print('''Usage: flac.py frame fieldname''')
        sys.exit(1)

    fl = Flac()

    frame = int(sys.argv[1])
    if frame == -1: frame = fl.nrec

    #results = []
    fieldname = sys.argv[2]
    fn = fl.__getattribute__('read_' + fieldname)
    field = fn(frame)
    #results.append(field)

    # print to screen
    printing(field)

    #print('# time =', fl.time[fl.nrec-1], 'Myrs')

