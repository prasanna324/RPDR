
def load_RPDR_mrn_multiple(dir_data, path_labs, delimiter='|'):
    ''' load_RPDR_mrn_multiple(dir_data, path_labs, delimiter='|'):
        Sequentially loads all files from RPDR data dump when multiple files have the same name. 
        
        1. Starts in dir_data (should have trailing slash), grabs all sub-folders' names automatically, then sequentially loads: dir_data/[sub-folders]/path_labs (where path_labs is the name of the file)
        * Note for whatever reason, on a multiple-split file dump from RPDR the labs, demographics, etc files are all named the exact same, just in different zips
        2. Calls the traditional load MRN function on each file
        3. Concatenates all results and returns 1 DF
        
        '''
    import os
    import pandas as pd
    
    # get list of subdirectories
    subdirectories = [x[0] for x in os.walk(dir_data)][1:]
    
    first=True
    # for each subdir, use the traditional load function to load data and concat
    for subdir in subdirectories:
        path_to_labs=subdir+'/'+path_labs
        mrn = load_RPDR_mrn(path=path_to_labs, delimiter=delimiter)
        
        if first==True:
            concat_pd = mrn
            first=False
        else:
            concat_pd=pd.concat([concat_pd, mrn],ignore_index=True)
    
    return concat_pd

def load_RPDR_mrn(path, delimiter='|'):
    import pandas as pd
    
    mrn = pd.read_csv(path, delimiter=delimiter, dtype=str)
        
    return mrn
