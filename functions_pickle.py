import pickle

def save(filename, obj):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)
        
def restore(filename):
    with open(filename, "rb") as f: # "rb" because we want to read in binary mode
        state = pickle.load(f)
    return state