
def load_RPDR_hnp_multiple(dir_data, fname, delimiter='|', datetime_col='Report_Date_Time'):
    ''' load_RPDR_hnp_multiple(dir_data, fname, delimiter='|', datetime_col='Report_Date_Time'):
        Sequentially loads all files from RPDR data dump when output is split. 
        
        1. Starts in dir_data (should have trailing slash), grabs all sub-folders' names automatically, then sequentially loads: dir_data/[sub-folders]/fname (where fname is the name of the file)
        * Note for whatever reason, on a multiple-split file dump from RPDR the labs, demographics, etc files are all named the exact same, just in different zips
        2. Calls the traditional load path function on each file
        3. Concatenates all results and returns 1 DF
        
        See load_native_data for remainder of parameters which are passed to that function
        
        '''
    import os
    import pandas as pd
    
    # get list of subdirectories
    subdirectories = [x[0] for x in os.walk(dir_data)][1:]
    
    first=True
    # for each subdir, use the traditional load function to load data and concat
    for subdir in subdirectories:
        path_to_path=subdir+'/'+fname
        path = load_RPDR_dis(path=path_to_path,
                              delimiter=delimiter,
                              datetime_col=datetime_col)
        
        if first==True:
            concat_pd = path
            first=False
        else:
            concat_pd=pd.concat([concat_pd, path],ignore_index=True)
    
    return concat_pd

def load_RPDR_hnp(path,delimiter='|', datetime_col='Report_Date_Time'):
    ''' load_RPDR_hnp(path, delimiter='|', datetime_col='Report_Date_Time'):
        Loads RPDR history and physical file as pandas dataframe
        
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
