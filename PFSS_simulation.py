import pickle
import os
import astropy
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sunpy.io.fits
import sunpy.map
from sunpy.net import Fido, attrs as a
from sunpy.image import coalignment
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, quantity_support
import eispac
import pfsspy
import pfsspy.tracing as tracing

def save(filename, obj):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

def restore(filename):
    with open(filename, "rb") as f: # "rb" because we want to read in binary mode
        state = pickle.load(f)
    return state


aia_submap = restore('20140202/20140202_aia_submap')
eis_fixed = restore('20140202/20140202_eis_map')
aia = restore('20140202/20140202_aia_map')
pfss_in = restore('20140202/20140202_pfss_in')

hp_lon = np.linspace(eis_fixed.bottom_left_coord.Tx/u.arcsec, eis_fixed.top_right_coord.Tx/u.arcsec, len(aia_submap.data[0])) * u.arcsec
hp_lat = np.linspace(eis_fixed.bottom_left_coord.Ty/u.arcsec, eis_fixed.top_right_coord.Ty/u.arcsec, len(aia_submap.data[0:])) * u.arcsec
lon, lat = np.meshgrid(hp_lon, hp_lat)
seeds = SkyCoord(lon.ravel(), lat.ravel(),
                 frame=aia.coordinate_frame)

m = pfss_in.map
pfss_out = pfsspy.pfss(pfss_in)

print('Currently Tracing')
tracer = tracing.FortranTracer(max_steps=80000, step_size=0.65)
flines = tracer.trace(seeds, pfss_out)

save('20140202/20140202_flines', flines)
fline_list = [i.coords for i in flines]
save('20140202/20140202_fline_list', fline_list)
