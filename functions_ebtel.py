import numpy as np
from functions_pixels import *
from tqdm import tqdm
from functions_pickle import *
from scipy.io import readsav
import glob
import re
import sunpy
import sys
import multiprocessing as mp
from functools import partial
import json
from pathlib import Path
from datetime import datetime


def get_var(date_string, description):
    with open(date_string+f'length_{description}.txt', "r") as f:
        length = f.readlines()
        length = [float(l.replace('\n','')) for l in length]

    with open(date_string+f'mean_field_{description}.txt', "r") as f:
        gauss = f.readlines()
        gauss = [float(l.replace('\n','')) for l in gauss]

    with open(date_string+f'heating_{description}.txt', "r") as f:
        heating = f.readlines()
        heating = [float(l.replace('\n','')) for l in heating]
    flines = restore(date_string+'flines.pickles')
    print('Order of variables: length, gauss, heating')
    print(f'List length of loop length: {len(length)}; Magnetic field: {len(gauss)}; Heating: {len(heating)}')
    
    return length, gauss, heating, flines

def get_new_var(date_string):
    with open(date_string+'intensity.txt', "r") as f:
        intensity = f.readlines()
        intensity = [list(l.replace('\n','')) for l in intensity]
    print('Order of variables: length, gauss, heating')
    print(f'List length of intensities: {len(intensity)}')
    
    return intensity


def check_file(python_index, file_list):
    file_heating = [float('.'.join(('0',re.findall(r'\d+', i)[-2]))) for i in file_list]
    file_length = [int(re.findall(r'\d+', i)[-1]) for i in file_list]
    synth_file_from_idl = np.array(file_length)+(np.array(file_heating)*1e4)
    synth_file_real_data = np.array(length)+(np.array(heating)*1e4)

    idl_index = list(synth_file_from_idl).index(min(synth_file_from_idl, key=lambda x:abs(x-synth_file_real_data[python_index])))
    # print(f'real heating: {heating_a[python_index]}, idl heating:{file_heating[idl_index]}')
    # print(f'real length: {length_a[python_index]}, idl length:{file_length[idl_index]}')
    return idl_index


def fline_multi(date_string, length, flines, aia_submap, file_list, i):
    try:
        if length[i] > 5e6:
            idl_index = check_file(i, file_list)
            pix_x, pix_y = filter_pix(flines[i].coords, aia_submap)
            intensity = readsav(file_list[idl_index])['int']
            # print(i)
            dict = {'pix_x':pix_x, 'pix_y':pix_y, 'int':intensity}
            # print(dict)
            # save(f'{date_string}simulated_intensities/{i}', dict)
            with open(f'{date_string}simulated_intensities/{i}', 'wb') as outp:  # Overwrites any existing file.
                pickle.dump(dict, outp, pickle.HIGHEST_PROTOCOL)
            return pix_x, pix_y, intensity
        else:
            pass
    except:
        pass

def synthetic_map(blank_array, date):
    files_intensity = glob.glob(date+"simulated_intensities/*")
    for file in tqdm(files_intensity):
        dict = restore(file)
        for num, time in enumerate(dict['int']):
            if dict['int'] != np.inf:
                blank_array[num][dict['pix_y'], dict['pix_x']] += time
            else:
                blank_array[num][dict['pix_y'], dict['pix_x']] += 0
    return blank_array

def inf_check(date):
    files_intensity = glob.glob(date+"simulated_intensities/*")
    int_list=[]
    for file in tqdm(files_intensity):
        dict = restore(file)
        int_list.append(dict['int'])
    return int_list

if __name__ == "__main__":
    '''
    calling statement: python functions_ebtel.py '20160520/'
    '''
    # date = '20110415/'
    date = sys.argv[1]
    aia_submap = restore(date+'m_aia_cutout.pickle')
    # eis_fixed = restore(date+'eis_map')
    aia = restore(date+'m_aia.pickle')
    # pfss_in = restore(date+'pfss_in')
    """
    change this later - aia_submap to eis_fixed
    """
    blank_data = np.zeros(len(aia_submap.data[0:])*len(aia_submap.data[0])).reshape(len(aia_submap.data[0:]),len(aia_submap.data[0]))
    blank_data = np.zeros((1800, blank_data.shape[0], blank_data.shape[1]))
    length, gauss, heating, flines = get_var(date, 'total')
    file_list_multi = glob.glob(date+"simulation_results/*.sav")
    Path(f'{date}simulated_intensities/').mkdir(parents=True, exist_ok=True)

    start=datetime.now()
    fline_partial = partial(fline_multi, date, length, flines, aia_submap, file_list_multi)
    print(f'Spreading process into multiple cores; Processing {len(flines)} fieldlines')
    for i in tqdm(range(len(flines))):
        fline_partial(i)

    # with mp.Pool(processes = 40) as p:
    #     p.map(fline_partial, range(len(flines)))

    print('Creating synthetic map')
    blank_data = synthetic_map(blank_data, date)

    # for result in results:
    #     blank_data[result[1], result[0]] += result[2]

    # for i in tqdm(range(len(flines))):
    #     if length[i] > 5e6:
    #         idl_index = check_file(i, file_list)
    #         try:
    #             pix_x, pix_y = filter_pix(flines[i].coords, eis_fixed)
    #             blank_data[pix_y, pix_x] += readsav(file_list[idl_index])['int'][400]
    #         except:
    #             failed_list.append(i)
    print('Saving syntheic data')
    # synth_map_multi = sunpy.map.Map(blank_data, aia_submap.meta)
    save(date+'synth_map_data.pickle', blank_data)
    print('Saved syntheic data')
    print(f'Total time used: {datetime.now()-start}')
    print(f'Synthetic map directory: {date+"synth_map.pickle"}')
