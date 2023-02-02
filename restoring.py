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

date = '20160520/'

fieldlines = restore(date+'flines.pickles')
pfss_output = restore(date+'pfss_output.pickle')
fline_list_closed = [i.coords for i in fieldlines.closed_field_lines]
fline_list_open = [j.coords for j in fieldlines.open_field_lines]

print(f'Closed fields: {len(fline_list_closed)}; Open fields: {len(fline_list_open)}')


def get_loop_length(line):
    c = line.cartesian.xyz
    s = np.append(0., np.linalg.norm(np.diff(c.value, axis=1), axis=0).cumsum())
    return np.diff(s).sum()

print(f'Compiling list for {len(fieldlines)} B-fields and loop lengths')

def get_ebtel_param(fline_list,description):
    from os import remove
    from glob import glob
    
    print(f'Working on {description}')
    gauss = []
    length = []
    for i in tqdm(fline_list, desc="Curating Loop Length & Field Strength"):
        if len(i) > 0:
            data = np.sqrt(np.nansum(np.square(pfss_output.get_bvec(i)),axis=1))
            data[data == 0] = np.nan
            gauss.append(np.nanmean(data))
            length.append(get_loop_length(i))

    heating = (0.0492*((29e6/np.array(length))*(np.array(gauss)/76)))
    print(f'Sanity Check: Gauss array length:{len(gauss)}; Length array length:{len(length)}; Heating array length:{len(heating)}')
    try:
        for f in glob(f"*{description}*.txt"):
            remove(f)
    except OSError:
        print('No txt file found')

    for i in tqdm(range(len(gauss)), desc="Writing EBTEL parameters"):
        with open(date+f'mean_field_{description}.txt', 'a+') as outfile:  
            outfile.write(f'{gauss[i]}\n')
        with open(date+f'length_{description}.txt', 'a+') as outfile:  
            outfile.write(f'{length[i]}\n')
        with open(date+f'heating_{description}.txt', 'a+') as outfile:  
            outfile.write(f'{heating[i]}\n')
    print(f'{description} : Average loop b-field: {np.mean(gauss)}; length: {np.mean(length)/1e6}')
    return gauss, length, heating

gauss_closed, length_closed, heating_closed = get_ebtel_param(fline_list_closed,'closed')
gauss_open, length_open, heating_open = get_ebtel_param(fline_list_open,'open')

gauss_total = gauss_closed + gauss_open
length_total = length_closed + length_open
heating_total = list(heating_closed) + list(heating_open)

for i in tqdm(range(len(gauss_total)), desc="combining two lists of EBTEL parameters"):
    with open(date+f'mean_field_total.txt', 'a+') as outfile:  
        outfile.write(f'{gauss_total[i]}\n')
    with open(date+f'length_total.txt', 'a+') as outfile:  
        outfile.write(f'{length_total[i]}\n')
    with open(date+f'heating_total.txt', 'a+') as outfile:  
        outfile.write(f'{heating_total[i]}\n')
