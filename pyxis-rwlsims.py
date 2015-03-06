## simulator 

# import Pyxis essentials
import Pyxis 
import ms
import im
import mqt
import lsm
from Pyxis.ModSupport import *

# import python essentials
import os
import sys
import numpy
import math

# some useful constants
PI = math.pi
FWHM = math.sqrt( math.log(256) )

# other python packages
import pyfits

# use simms to create empty MSs (https://github.com/SpheMakh/simms)
from simms import simms

import im.lwimager 
import im.argo


def make_empty_ms(msname='$MS', observatory='$OBSERVATORY', antennas='$ANTENNAS',
                  synthesis='$SYNTHESIS', dtime='$DTIME', freq0='$FREQ0',
                  dfreq='$DFREQ',nchan='$NCHAN', **kw):
    """ creates an empty MS """
    
    msname, observatory, antennas, synthesis, dtime, freq0, dfreq, nchan =\
            interpolate_locals('msname observatory antennas synthesis dtime freq0 dfreq nchan')

    if not exists(msname) or MS_REDO:
        simms(tel=observatory, pos=antennas, msname=msname, synthesis=float(synthesis),
              dtime=float(dtime), freq0=freq0, dfreq=dfreq, nchan=int(nchan), **kw)
    
    v.MS = msname
    return msname


def simsky(msname='$MS', lsmname='$LSM', tdlsec='$TDLSEC', tdlconf='$TDLCONF',
           column='$COLUMN', noise=0,recenter_lsm=True, args=[],**kw):
    """ Simulates visibilities into an MS """

    msname, lsmname, column, tdlsec, tdlconf = \
        interpolate_locals('msname lsmname column tdlsec tdlconf')
    
    fits = True if verify_sky(lsmname) is 'FITS' else False

    v.MS = msname

    # use LWIMAGER to predict visibilities is skymodel is a FITS file
    if fits:
        lsmname = conform(fitsname=lsmname)
        v.LSM = lsmname
        _column = 'MODEL_DATA' if noise else column
        im.lwimager.predict_vis(image=lsmname, wprojplanes=128, column=_column, padding=1.5,**kw)

        if noise:
            simnoise(noise=noise,addToCol=_column,column=column)
    # use MeqTrees to predict visibilities if skymodel is Tigger Model or an ASCII file
    else:
        if recenter_lsm:
            # save temp lsm in temporary file
            tlsm = tempfile.NamedTemporaryFile(suffix='.lsm.html')
            tlsm.flush
            tlsmname = tlsm.name
            x.sh('tigger-convert --recenter=$DIRECTION $lsmname $tlsmname -f')
            v.LSM = lsmname
        else:
            v.LSM = lsmname

        args = ["${ms.MS_TDL} ${lsm.LSM_TDL}"] + list(args)

        options = {}
        options['ms_sel.output_column'] = column

        if noise:
            options['noise_stddev'] = noise

        options.update(kw)
        mqt.run(SIMCRIPT, job='_tdl_job_1_simulate_MS',
                config=tdlconf, section=tdlsec, options=options, args=args)
        
        tlsm.close()

document_globals(simsky,"MS LSM COLUMN TDLSEC TDLCONF")    


def driver():
    make_empty_ms()
    simsky()
    im.lwimager.make_image()


def conform(msname='$MS',fitsname='$LSM',outfile=None):
    """ conforms FITS file to a structure acceptable to LWIMAGER.
        This assunes a 2D FITS file. It will also work with a 3D FITS file, 
        with the 3rd axis being frequency.
    """

    msname, fitsname = interpolate_locals('msname fitsname')

    ## Insure conformance by letting LWIMAGER create the header
    # Lets use a temp file
    tf = tempfile.NamedTemporaryFile(suffix='.fits',dir='.')
    im.argo.make_empty_image(msname=msname,image=tf.name)
    
    hdu = pyfits.open(fitsname)
    hdr = pyfits.open(tf.name)[0].header # header from the lwimager empty image
    data = hdu[0].data
    shape = list(data.shape)

    if len(shape) == 2:
        shape + [1,1]
    elif len(shape) == 3:
        shape = shape[:2] + [1] + [shape[2]]
    else:
        warn('FITS file has more than 3 axes, not touching it!')
        return

    hdu[0].data = numpy.reshape(data,shape)
    hdu[0].header = hdr

    outfile = outfile or fitsname.replace('.fits','4d.fits')
    hdu.writeto(outfile,clobber=True)
    
    return outfile


def verify_sky(fname):
    """ verify if skymodel is compitable with simulator """

    ext = fname.split('.')[-1]
    if ext.lower() == 'fits':
        return 'FITS'
    elif ext.lower() == 'txt' or fname.endswith('.lsm.html'):
        return 'TIGGER'
    else:
        raise TypeError('Sky model "%s" has to be either one of FITS,ASCII,Tigger Model (lsm.html) '%fname)


def compute_vis_noise (noise=0,sefd=None):
    """Computes nominal per-visibility noise"""

    tab = ms.ms()
    spwtab = ms.ms(subtable="SPECTRAL_WINDOW")
    sefd = sefd or SEFD 
   
    freq0 = spwtab.getcol("CHAN_FREQ")[ms.SPWID,0]
    wavelength = 300e+6/freq0
    bw = DFREQ or spwtab.getcol("CHAN_WIDTH")[ms.SPWID,0]
    dt = DTIME or tab.getcol("EXPOSURE",0,1)[0]
    dtf = (tab.getcol("TIME",tab.nrows()-1,1)-tab.getcol("TIME",0,1))[0]
    # close tables properly, else the calls below will hang waiting for a lock...
    tab.close()
    spwtab.close()

    info(">>> $MS freq %.2f MHz (lambda=%.2fm), bandwidth %.2g kHz, %.2fs integrations, %.2fh synthesis"%(freq0*1e-6,wavelength,bw*1e-3,dt,dtf/3600))
    if not noise:
        noise = sefd/math.sqrt(2*bw*dt)
        info(">>> SEFD of %.2f Jy gives per-visibility noise of %.2f mJy"%(sefd,noise*1000))
    else:
        info(">>> using per-visibility noise of %.2f mJy"%(noise*1000))
    return noise


def simnoise (noise=0,rowchunk=100000,skipnoise=False,addToCol=None,scale_noise=1.0,column='MODEL_DATA'):
    """ generate noise and add to a given column in the MS """

    spwtab = ms.ms(subtable="SPECTRAL_WINDOW")
    freq0 = spwtab.getcol("CHAN_FREQ")[ms.SPWID,0]/1e6
    tab = ms.msw()
    dshape = list(tab.getcol('DATA').shape)
    nrows = dshape[0]

    noise = noise or compute_vis_noise()
    if addToCol: colData = tab.getcol(addToCol)

    for row0 in range(0,nrows,rowchunk):
        nr = min(rowchunk,nrows-row0)
        dshape[0] = nr
        data = noise*(numpy.random.randn(*dshape) + 1j*numpy.random.randn(*dshape)) * scale_noise
        if addToCol:
            data+=colData[row0:(row0+nr)]
            info(" $addToCol + noise --> $column (rows $row0 to %d)"%(row0+nr-1))
        else : info("Adding noise to $column (rows $row0 to %d)"%(row0+nr-1))
        tab.putcol(column,data,row0,nr);
    tab.close()
