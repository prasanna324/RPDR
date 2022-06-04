
def load_RPDR_prc_multiple(dir_data, fname, delimiter='|', datetime_col='Date'):
    ''' load_RPDR_prc_multiple(dir_data, fname, delimiter='|', datetime_col='Date'):
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
        path = load_RPDR_prc(path=path_to_path,
                              delimiter=delimiter,
                              datetime_col=datetime_col)
        
        if first==True:
            concat_pd = path
            first=False
        else:
            concat_pd=pd.concat([concat_pd, path],ignore_index=True)
    
    return concat_pd


def load_RPDR_prc(path,delimiter='|', datetime_col='Date'):
    ''' load_RPDR_prc(path, delimiter='|', datetime_col='Date'):
        Loads RPDR procedures file as pandas dataframe
        
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
            
                    # Number of times "|" is present in a line
                    count = line.count("|")
                    
                    if count not in [9,14]:
                        # Replace the last occurence of '|' with '|"'
                        line = line.replace('\n', '')
                    file_w.write(line)
        file_r.close()
        file_w.close()
        
    path = write_path
    
    # Read the processed .csv file from path location
    print('Reading from : ' + path)
    path_df = pd.read_csv(path, sep=delimiter, dtype=str)
    
    # Convert datetime column to pandas date time format
    path_df['datetime'] = pd.to_datetime(path_df.loc[:,datetime_col])
    
    return path_df
