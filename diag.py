
def load_diag(path, delimiter='|', datetime_col='Date', prune=True):
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
