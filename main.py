from astropy.io import fits
import numpy as np
import astropy.coordinates
import astropy.units as u
from functions_pickle import *
from functions_data import data_download_from_date,hmi_daily_download, aia_correction
import matplotlib.pyplot as plt
import pfsspy
import sunpy
from tqdm import tqdm
from aiapy.calibrate import correct_degradation, update_pointing
from astropy.coordinates import SkyCoord
from sunpy.net import Fido,attrs
# from sys import setrecursionlimit
import faulthandler
faulthandler.enable()

# setrecursionlimit(10**6) 
change_obstime = lambda x,y: SkyCoord(x.replicate(observer=x.observer.replicate(obstime=y), obstime=y))
change_obstime_frame = lambda x,y: x.replicate_without_data(observer=x.observer.replicate(obstime=y), obstime=y)

# Restoring maps pickled
date = '20160520/'
m_hmi = restore(date+'m_hmi_full_res.pickle')
m_aia = restore(date+'m_aia.pickle')
m_cutout = restore(date+'m_aia_cutout.pickle')

# Locating the center of the active region for cutout
ar_center = SkyCoord(lon=219*u.deg, lat=-8*u.deg, frame=m_hmi.coordinate_frame)
ar_width = 600*u.arcsec
ar_height = 600*u.arcsec

# Resample PFSS
m_hmi_resample = m_hmi.resample((1800, 1080 )*u.pix)

# Setting parameters for PFSSpy
print('Setting parameters for PFSSpy')
nrho = 70
rss = 2

# PFSS input map and output PFSS format
pfss_input = pfsspy.Input(m_hmi_resample, nrho, rss)

print('Setting PFSS output')
pfss_output = pfsspy.pfss(pfss_input)
save(f"{date}pfss_output.pickle", pfss_output)

# Locating seeds using AIA cutout
new_frame = change_obstime_frame(m_hmi.coordinate_frame, m_cutout.date)
blc_ar_synop = change_obstime(m_cutout.bottom_left_coord.transform_to(new_frame),
                              m_hmi.date)
trc_ar_synop = change_obstime(m_cutout.top_right_coord.transform_to(new_frame),
                              m_hmi.date)

# Pinning the seeds onto synoptic map using the coordinates above
masked_pix_y, masked_pix_x = np.where(m_hmi_resample.data > 10)
seeds = m_hmi_resample.pixel_to_world(masked_pix_x*u.pix, masked_pix_y*u.pix,).make_3d()
in_lon = np.logical_and(seeds.lon > blc_ar_synop.lon, seeds.lon < trc_ar_synop.lon)
in_lat = np.logical_and(seeds.lat > blc_ar_synop.lat, seeds.lat < trc_ar_synop.lat)
seeds = seeds[np.where(np.logical_and(in_lon, in_lat))]

# Setting parameters for the PFSS tracer
ds = 0.01
max_steps = int(np.ceil(10 * nrho / ds))

# Tracing PFSS fieldlines
print('Currently Tracing')
tracer = pfsspy.tracing.FortranTracer(step_size=ds, max_steps=max_steps)
# tracer = pfsspy.tracing.FortranTracer(max_steps=max_steps)
fieldlines = tracer.trace(SkyCoord(seeds), pfss_output,)

save(date+'flines.pickles', fieldlines)
fline_list_closed = [i.coords for i in fieldlines.closed_field_lines]
fline_list_open = [j.coords for j in fieldlines.open_field_lines]
save(date+'fline_list_closed.pickles', fline_list_closed)
save(date+'fline_list_open.pickles', fline_list_open)

print('Tracing Finished')

def get_loop_length(line):
    c = line.cartesian.xyz
    s = np.append(0., np.linalg.norm(np.diff(c.value, axis=1), axis=0).cumsum())
    return np.diff(s).sum()

gauss = []
length = []

print(f'Compiling list for {len(fieldlines)} B-fields and loop lengths')

def get_ebtel_param(fline_list,description):
    print(f'Working on {description}')
    for i in tqdm(fline_list):
        if len(i) > 0:
            data = np.sqrt(np.nansum(np.square(pfss_output.get_bvec(i)),axis=1))
            data[data == 0] = np.nan
            gauss.append(np.nanmean(data))
            length.append(get_loop_length(i))

    print(f'{description} : Average loop b-field: {np.mean(gauss)}; length: {np.mean(length)/1e6}')
    heating = (0.0492*((29e6/np.array(length))*(np.array(gauss)/76)))
    for i in tqdm(range(len(gauss))):
        with open(date+f'mean_field_{description}.txt', 'a+') as outfile:  
            outfile.write(f'{gauss[i]}\n')
        with open(date+f'length_{description}.txt', 'a+') as outfile:  
            outfile.write(f'{length[i]}\n')
        with open(date+f'heating_{description}.txt', 'a+') as outfile:  
            outfile.write(f'{heating[i]}\n')

get_ebtel_param(fline_list_closed,'closed')
get_ebtel_param(fline_list_open,'open')

print('finished writing into txt file')
