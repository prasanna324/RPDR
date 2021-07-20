
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
