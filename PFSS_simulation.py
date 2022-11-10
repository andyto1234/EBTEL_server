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
import json
from functions_pickle import *
from tqdm import tqdm
from ebtel_adapt import adapt2pfsspy
from glob import glob


date = '20110415/'
aia_submap = restore(date+'aia_submap')
eis_fixed = restore(date+'eis_map')
aia = restore(date+'aia_map')
pfss_in = restore(date+'pfss_in')
adapt_map = glob(date+'adapt*')[0]


hp_lon = np.linspace(eis_fixed.bottom_left_coord.Tx/u.arcsec, eis_fixed.top_right_coord.Tx/u.arcsec, len(eis_fixed.data[0])) * u.arcsec
hp_lat = np.linspace(eis_fixed.bottom_left_coord.Ty/u.arcsec, eis_fixed.top_right_coord.Ty/u.arcsec, len(eis_fixed.data[0:])) * u.arcsec
lon, lat = np.meshgrid(hp_lon, hp_lat)
seeds = SkyCoord(lon.ravel(), lat.ravel(),
                 frame=aia.coordinate_frame)

pfss_model = adapt2pfsspy(adapt_map,rss=2.5)

# nrho = 45
# rss = 2.5
# pfss_input = pfsspy.Input(pfss_in, nrho, rss)

# pfss_out = pfsspy.pfss(pfss_input)

print('Currently Tracing')
tracer = pfsspy.tracing.FortranTracer(max_steps=80000, step_size=0.8)
flines = tracer.trace(seeds, pfss_model)

save(date+'flines', flines)
fline_list = [i.coords for i in flines]
save(date+'fline_list', fline_list)

print('Tracing Finished')

def get_loop_length(line):
    c = line.coords.cartesian.xyz
    s = np.append(0., np.linalg.norm(np.diff(c.value, axis=1), axis=0).cumsum()) * c.unit
    return np.diff(s).sum()

gauss = []
length = []

print(f'Compiling list for {len(flines)} B-fields and loop lengths')

for i in tqdm(range(len(flines))):
    gauss.append(np.mean([sum(b**2)**0.5 for b in flines[i].b_along_fline]))
    length.append(get_loop_length(flines[i]).to_value(u.m))

print(f'Average loop b-field: {np.mean(gauss)}; length: {np.mean(length)/1e6}')
heating = (0.0492*((29e6/np.array(length))*(np.array(gauss)/76)))
for i in range(len(gauss)):
    with open(date+'mean_field.txt', 'a+') as outfile:  
        outfile.write(f'{gauss[i]}\n')
    with open(date+'length.txt', 'a+') as outfile:  
        outfile.write(f'{length[i]}\n')
    with open(date+'heating.txt', 'a+') as outfile:  
        outfile.write(f'{heating[i]}\n')

print('finished writing into txt file')
    

