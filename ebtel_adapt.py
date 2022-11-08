import pfsspy
import sunpy
import numpy as np
from functions_pickle import *

def adapt2pfsspy(filepath, #must already exist on your computer
                 rss=2.5, # Source surface height
                 nr=60, # number of radial gridpoiints for model
                 realization="mean", #which slice of the adapt ensemble to choose
                 return_magnetogram = False # switch to true for function to return the input magnetogram
                ):

    # Load the FITS file into memory
    # ADAPT includes 12 "realizations" - model ensembles
    # pfsspy.utils.load_adapt is a specific function that knows
    # how to handel adapt maps
    adaptMapSequence = pfsspy.utils.load_adapt(filepath)
    # If realization = mean, just average them all together
    if realization == "mean" : 
        br_adapt_ = np.mean([m.data for m in adaptMapSequence],axis=0)
        adapt_map = sunpy.map.Map(br_adapt_,adaptMapSequence[0].meta)
        
    # If you enter an integer between 0 and 11, the corresponding
    # realization is selected
    elif isinstance(realization,int) : adapt_map = adaptMapSequence[realization]
    else : raise ValueError("realization should either be 'mean' or type int ") 
    # pfsspy requires that the y-axis be in sin(degrees) not degrees
    # pfsspy.utils.car_to_cea does this conversion
    adapt_map_strumfric = pfsspy.utils.car_to_cea(adapt_map)
    # Option to return the magnetogram
    if return_magnetogram : 
        return adapt_map_strumfric
    # Otherwise run the PFSS Model and return
    else :
        adapt_map_input = sunpy.map.Map(adapt_map_strumfric.data,
                                        adapt_map_strumfric.meta)
        peri_input = pfsspy.Input(adapt_map_input, nr, rss)
        peri_output = pfsspy.pfss(peri_input)
        return peri_output

