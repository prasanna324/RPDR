
def load_RPDR_endo(path,delimiter='|', datetime_col='Report_Date_Time'):
    ''' load_RPDR_endo(path, delimiter='|', datetime_col='Report_Date_Time'):
        Loads RPDR endoscopy notes file as pandas dataframe
        
        PARAMETERS:
        path: path to csv file or other text delimited file
        delimiter: delimiter for path file
        datetime_col: column name containing date/time information for each path report

        returns: pandas dataframe containing path information
        
        WARNINGS:
        1. Current function automatically searches for path + 'multiline_corrected', *if present* it assumes that is the correct 
            file. E.g., path='/data/path.txt', it searches for '/data/path_multiline_corrected.txt'.
        2. It will not overwrite this file if present
    
        '''
    import pandas as pd
    import os.path
    from os import path as os_path
    
    write_path = path.replace('.','_multiline_corrected.')
    if os_path.exists(write_path)==False:
        print('Reformatting path file to allow multi-line report text to be readable, saving as : {}'.format(write_path))
        
        with open(write_path,'w') as file_w:
            with open(path) as file_r:
                for i in range(1):
                    first_line = next(file_r)
                    file_w.write(first_line)
                for i,line in enumerate(file_r):
                    # Replace single quote with double quotes in all the lines
                    line = line.replace('"', '""')
                    # Find the right-most occurence of the character "|" in a line
                    index = line.rfind("|")
                    # Number of times "|" is present in a line
                    count = line.count("|")
                    
                    if index!=-1 and count==9:
                        # Replace the last occurence of '|' with '|"'
                        line = line[:index+1] + '"' + line[index+1:]
                    line = line.replace('[report_end]', '[report_end]"')
                    file_w.write(line)
        file_r.close()
        file_w.close()
        
    path = write_path
    
    # Read the processed .csv file from path location
    print('Reading from : ' + path)
    path_df = pd.read_csv(path, sep=delimiter, dtype=str)
    
    # Create unique_report_id by joining EMPI and Report_Number
    path_df['unique_report_id'] = path_df.apply(lambda x: str(x.EMPI) + '_' + str(x.Report_Number),axis=1)
    
    # Drop duplicates for 'unique_report_id' cases based on length of Result_Text instead of dropping the first observed case
    path_df['report_len'] = path_df['Report_Text'].str.len()
    
    path_df = (path_df
      .sort_values(['unique_report_id', 'report_len'])
      .drop_duplicates(subset=['unique_report_id'], keep='last', inplace=False, ignore_index=False)
      .drop(columns=['report_len'])
      .sort_index()
     )
    
    # now set index to unique_report_id
    path_df.set_index(keys='unique_report_id', inplace=True, verify_integrity=True)
    
    # Convert datetime column to pandas date time format
    path_df['datetime'] = pd.to_datetime(path_df.loc[:,datetime_col])
    
    return path_df




def truncate_dx_start(pathdf, update=True):
    ''' truncate_finaldx(pathdf, update=True)
    DESC: For MGH, BWH path reports, find the 'final diagnosis' line of the path report, remove everything preceding. Parse 
     to extract whether or not there is a final diagnosis line, what that line is, and the full report text (following this line)
    
    PARAMETERS:
    pathdf: pathology dataframe from load_RPDR_path
    update: return only MGH, BWH path as a new df or update pathdf with MGH, BWH results
    
    RETURNS: pathdf or new path dataframe with only MGH, BWH path reports (depending on update bool) with three new columns:
    ['has_final_diagnosis'] = did the function find a line of text that it thinks contains final diagnosis line?
    ['final_diagnosis_line'] = if has_final_diagnosis == True, then what is the final diagnosis line?
    ['Report_Text'] = report text after removing everything above final diagnosis
    '''
    import re
    
    # truncate to only final diagnosis (for MGH, BWH reports only)
    #print('Filtering only MGH, BWH path reports...')
    
    # first get only MGH values
    #fil_mgh = pathdf.MRN_Type == 'MGH'
    df_path = pathdf.copy()
    
    
    # FINAL DIAGNOSIS LINE FINDER
    num_reports = df_path.shape[0] # num rows of df_path
    has_final_diagnosis_col = []
    final_diagnosis_line_col = []
    trunc_path_col = []
    print('Truncating to only final diagnosis...')
    for i in range(0,num_reports):
        # extract path report for this entry
        report_text = df_path.iloc[i,:].Report_Text
        site = df_path.iloc[i,:].MRN_Type
        # split by newline character
        text_by_line = report_text.split('\n')

        has_final_diagnosis = False
        final_diagnosis_line = ''
        trunc_path_text = report_text

        # go line-by-line and perform some checks
        j=0
        
        if site in ['MGH','BWH','NWH','FH','NSM']:
            for line in text_by_line:
                lower_line = line.lower().strip()
                # capture situation where a line contains liver, biopsy; note will only grab first instance then short circuit
                
                if has_final_diagnosis==False and ('findings:' in lower_line or
                                                   'impression:' in lower_line):
                    has_final_diagnosis = True
                    final_diagnosis_line = line
                    trunc_path_text = '\n'.join(text_by_line[j:]) # should be a list of this line and all subsequent lines

                j=j+1
                
        has_final_diagnosis_col.append(has_final_diagnosis)
        final_diagnosis_line_col.append(final_diagnosis_line)
        # either returns the original report or the truncated form if it has a final diagnosis to truncate at
        trunc_path_col.append(trunc_path_text)

    df_path['has_dx_start'] = has_final_diagnosis_col
    df_path['dx_start_line'] = final_diagnosis_line_col
    df_path['Report_Text'] = trunc_path_col


    if update:
        # re-merge with original data
        #print('Updating input path dataframe with truncated MGH, BWH path reports')
        pathdf['has_dx_start'] = False
        pathdf['dx_start_line'] = ''
        pathdf.update(df_path)
        return_df = pathdf.copy()
    else:
        # return this mgh, bwh path only file
        #print('Returning only MGH, BWH entries with truncated path reports')
        return_df = df_path
        
        print('Done. | Status: ' + str(df_path[df_path.has_dx_start==True].shape[0]) + ' reports with a diagnosis beginning line, ' 
              + str(df_path[df_path.has_dx_start==False].shape[0]) + ' reports with no diagnosis beginning line')
    
    return return_df


def truncate_dx_end(pathdf, update=True, only_dx_start=True):
    
    #print('Filtering only MGH, BWH path reports...')
    fil_subset = pathdf.MRN_Type.isin(['MGH', 'BWH','NWH','FH','NSM'])
    df_path = pathdf[fil_subset].copy()
    
    if only_dx_start:
        # check the column exists first:
        if 'has_dx_start' in df_path.columns.tolist():
            fil_finaldx_trunc = df_path.has_dx_start == True
            df_path = df_path[fil_finaldx_trunc]
        else:
            print('The flag *only_dx_start=True* was passed, however truncate_finaldx() has not been called. Aborting...')
            return None

    num_reports = df_path.shape[0]
    has_lowersec_col = []
    lowersec_line_col = []
    lowersec_start_LAFD_col = []
    trunc_path_col = []
    
    for i in range(0,num_reports):
        # extract path report for this entry
        report_text = df_path.iloc[i,:].Report_Text
        # split by newline character
        text_by_line = report_text.split('\n')
        
        has_lowersec = False
        lowersec_line = ''
        lowersec_start_LAFD = -1
        trunc_path_text = report_text
        
        # go line-by-line and perform some checks
        j=0
        for line in text_by_line:
            lower_line = line.lower()
            
            if has_lowersec==False and ('recommendation:' in lower_line or
                                        'recommend:' in lower_line or
                                        'recommendations:' in lower_line or
                                        'procedure codes:' in lower_line or
                                        'procedure code(s):' in lower_line):
                
                has_lowersec = True
                lowersec_line = line
                lowersec_start_LAFD = j
                trunc_path_text = '\n'.join(text_by_line[:j])
            j=j+1

        has_lowersec_col.append(has_lowersec)
        lowersec_line_col.append(lowersec_line)
        lowersec_start_LAFD_col.append(lowersec_start_LAFD)
        trunc_path_col.append(trunc_path_text)
        
    df_path['has_dx_end'] = has_lowersec_col
    df_path['dx_end_line'] = lowersec_line_col
    df_path['dx_end_line_LAFD'] = lowersec_start_LAFD_col
    df_path['Report_Text'] = trunc_path_col
    
    
    if update:
        # re-merge with original data
        print('Updating input path dataframe with truncated MGH, BWH path reports')
        pathdf['has_dx_end'] = False
        pathdf['dx_end_line'] = ''
        pathdf['dx_end_line_LAFD'] = -1
        pathdf.update(df_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
        print('Returning MGH, BWH only entries with truncated path reports')
        return_df = df_path
    
    return return_df

    
