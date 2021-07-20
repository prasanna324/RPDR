
def load_RPDR_meds_multiple(dir_data, filename_meds, delimiter='|', make_lower=True):
    ''' load_RPDR_meds_multiple(dir_data, filename_meds, delimiter='|', make_lower=True):
        Sequentially loads all files from RPDR data dump when multiple files have the same name. 
        
        1. Starts in dir_data (should have trailing slash), grabs all sub-folders' names automatically, then sequentially loads: dir_data/[sub-folders]/path_labs (where path_labs is the name of the file)
        * Note for whatever reason, on a multiple-split file dump from RPDR the labs, demographics, etc files are all named the exact same, just in different zips
        2. Calls the traditional load meds function on each file
        3. Concatenates all results and returns 1 DF
        4. If make_lower==True, converts all medication names to lower case
        
        Warnings:
        1. Do not have any other subfolders besides the ones containing data in dir_data
        
        '''
    import os
    import pandas as pd
    
    # get list of subdirectories
    subdirectories = [x[0] for x in os.walk(dir_data)][1:]
    
    first=True
    # for each subdir, use the traditional load function to load data and concat
    for subdir in subdirectories:
        path_to_labs_full=subdir+'/'+filename_meds
        meds = load_RPDR_meds(path=path_to_labs_full, delimiter=delimiter)
        
        if first==True:
            concat_pd = meds
            first=False
        else:
            concat_pd=pd.concat([concat_pd, meds],ignore_index=True)
        
    if make_lower:
        concat_pd['Medication'] = concat_pd.Medication.apply(lambda x: x.lower())
        
    return concat_pd

def load_RPDR_meds(path, delimiter='|', datetime_col='Medication_Date', prune=True):
    import pandas as pd
    
    # load a medications record
    
    meds = pd.read_csv(path, delimiter=delimiter, dtype=str)
    # enforce the EMPI column is strings for later
    
    meds['datetime'] = pd.to_datetime(meds.loc[:,datetime_col])
    
    return meds[['EMPI', 'EPIC_PMRN', 'MRN_Type', 'MRN', 'datetime', 'Medication', 'Code_Type', 'Code', 'Quantity', 'Inpatient_Outpatient']]
