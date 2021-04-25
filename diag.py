
def load_RPDR_diag_multiple(dir_data, filename_diag, delimiter='|'):
    ''' load_RPDR_diag_multiple(dir_data, filename_diag, delimiter='|', make_lower=True):
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
        path_to_diag_full=subdir+'/'+filename_diag
        diag = load_RPDR_diag(path=path_to_diag_full, delimiter=delimiter)
        
        if first==True:
            concat_pd = diag
            first=False
        else:
            concat_pd=pd.concat([concat_pd, diag],ignore_index=True)
        
    return concat_pd

def load_RPDR_diag(path, delimiter='|', datetime_col='Date', prune=True):
    import pandas as pd
    
    diag = pd.read_csv(path, delimiter=delimiter, dtype=str)
    
    diag['Diagnosis_Name'] = diag['Diagnosis_Name'].str.lower()
    
    diag['datetime'] = pd.to_datetime(diag.loc[:,datetime_col])
    
    return diag

def tag_transplant(df, diag, remove=False):
    '''
    Accepts a dataframe with patient ID col name patient_ID, cross references the diagnosis database for that
    patient identifier for signature of transplant, and approximate date, adds two columns txplanted and txplanted_date
    
    WARNINGS:
    - assumes MRN, datetime are names of columns in both dataframes
    '''
    
    import pandas as pd
        
    # first identify transplant entries from diag
    
    fil1 = (diag.Diagnosis_Name.str.contains('transplant')) & (diag.Diagnosis_Name.str.contains('liver'))
    liver_dx1 = diag[fil1].copy()
    
    liver_dx2 = liver_dx1.sort_values(by='datetime', ascending=True).groupby('MRN').first()
    
    # reset index moves MRN to a column again, renumbers index
    liver_dx2 = liver_dx2.reset_index()
    
    # remove all columns except MRN and datetime
    liver_dx2 = liver_dx2[['MRN', 'datetime']]
    
    # rename datetime to txplant_earliest
    liver_dx2.rename(columns={'datetime': 'txplant_earliest'}, inplace=True)
    # make these all transplants
    liver_dx2['is_transplant'] = True
    
    # left merge the transplant database onto the passed dataframe
    out_df = df.merge(liver_dx2, how='left', on='MRN')
    
    # fill out NaN's in the is_transplant column with False if nan
    out_df['is_transplant'] = out_df.is_transplant.apply(lambda x: False if pd.isna(x) else x)
    
    return out_df, liver_dx2
