
def load_RPDR_path_multiple(dir_data, fname, delimiter='|', datetime_col='Report_Date_Time'):
    ''' load_RPDR_path_multiple(dir_data, fname, delimiter='\t', datetime_col='Report_Date_Time'):
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
        path = load_RPDR_path(path=path_to_path,
                              delimiter=delimiter,
                              datetime_col=datetime_col)
        
        if first==True:
            concat_pd = path
            first=False
        else:
            concat_pd=pd.concat([concat_pd, path],ignore_index=True)
    
    return concat_pd

def load_RPDR_path(path,delimiter='|', datetime_col='Report_Date_Time'):
    ''' load_RPDR_path(path,string_format=False,delimiter='|', datetime_col='Report_Date_Time')
    DESC: loads an RPDR pathology file to pandas 
    1. removes deleted, canceled and in process path reports
    2. removes duplicate path reports, keeping the first entry
    3. resets index to 'Report_Number' column
    4. converts datetime_col to pd.DateTime format
    
    PARAMETERS:
    path: path to csv file or other text delimited file
    delimiter: delimiter for path file
    datetime_col: column name containing date/time information for each path report
    
    returns: pandas dataframe containing path information
    
    WARNINGS:
    1. Current function automatically searches for path + 'multiline_corrected', *if present* it assumes that is the correct 
        file. E.g., path='/data/path.txt', it searches for '/data/path_multiline_corrected.txt'.
    2. It will not overwrite this file if present
    
    TO-DO:
    1. Update path report selection when there are near-duplicates. These are almost always of the form 'Final' and 'Updated', 
        of which, when duplicated, we should take the 'Updated' report. I have looked at about 10 examples so far (Marc)
    '''
    import pandas as pd
    import os.path
    from os import path as os_path
    
    ## Mod 1 begins ##
    ## default Report_Text format is multi-line. If string_format==True, convert multi-line text by double-quoting all quotes
    ##  (which allows read_csv to read through by default), and enclosing full report in single-quotes
    write_path = path.replace('.','_multiline_corrected.')
    if os_path.exists(write_path)==False:
        print('Reformatting path file to allow multi-line report text to be readable, saving as : {}'.format(write_path))
        f = open(path,'r')
        filedata = f.read()
        f.close()
        newdata = filedata.replace('"', '""')
        newdata = newdata.replace('Accession Number:', '"Accession Number:')
        newdata = newdata.replace('[report_end]', '[report_end]"')
        f2 = open(write_path,'w')
        f2.write(newdata)
        f2.close()
        
        path=write_path
        
    ## Mod 1 ends ##
    
    print('Reading from : ' + path)
    path_df = pd.read_csv(path, sep=delimiter, dtype=str)
    
    # unique_identifier. Requires some explanation.. report_number SHOULD be unique (it is literally the accession
    #  number for finding a block of tissue), however it is NOT in some specific examples of NSMC and MGH. E.g., 
    #  in all_RPDR_path, report numbers  : 
    # 'S11-5510', 'S14-3070', 'S17-25856', 'S13-2400', 'S14-9841', 'S13-3071',
    #    'S12-6506', 'S14-218', 'S12-3414', 'S13-32', 'S13-6114', 'S12-2481',
    #    'S13-1077', 'S14-41', 'S13-8207', 'S12-10612', 'S12-5964', 'S10-9798',
    #    'S09-10295', 'S16-15842', 'S14-8374', 'S12-3350', 'S12-6495', 'S14-785',
    #    'S16-183'
    #  all give duplicates for NSMC & MGH. And they're *different* patients, totally diff reports
    #  solution: append EMPI (unique per patient) to report_number to screen out these cases, operate on this concatenated report id
    path_df['unique_report_id'] = path_df.apply(lambda x: str(x.EMPI) + '_' + str(x.Report_Number),axis=1)
    
    # remove deleted, cancelled, and in-process path reports
#     bad_reports = (path_df.Report_Status == 'Deleted') | (path_df.Report_Status == 'Cancelled') | (path_df.Report_Status == 'In Process') | (path_df.Report_Status == 'Preliminary') | (path_df.Report_Status == 'Hold')  | (path_df.Report_Status == 'Pending')
    bad_statuses=['Deleted', 'Cancelled', 'In Process', 'Preliminary', 'Hold', 'Pending','Unknown','In Revision']
    bad_reports = path_df.Report_Status.isin(bad_statuses)
    
    # drop rows with any of these features
    path_df2 = path_df[~bad_reports].copy()
    
    # re-index to 'Report_Number', dropping duplicates ONLY if duplicate on both report_number and Report_Text
    # (ie, there is a problem if they're not the same..)
    path_df2.drop_duplicates(subset=['unique_report_id','Report_Text'], keep='first', inplace=True, ignore_index=False)
    
    # print duplicated entries and delete 
    dup_entries=path_df2.duplicated(subset='unique_report_id', keep=False)
    print('Duplicated entries (will keep first only, please review):')
    print(path_df2[dup_entries].Report_Text)
    
    # drop these. Why twice? Because if the report text is identical, don't need to review. But these are non-identical duplicates,
    #  usually, these are path reports that were revised/updated
    path_df2.drop_duplicates(subset=['unique_report_id'], keep='first', inplace=True, ignore_index=False)
    
    # now set index to Report_Number
    path_df2.set_index(keys='unique_report_id', inplace=True, verify_integrity=True)
    
    # set date time as datetime
    path_df2['datetime'] = pd.to_datetime(path_df2.loc[:,datetime_col])

    return path_df2

def mgh_truncate_finaldx(pathdf, update=True):
    ''' mgh_truncate_finaldx(pathdf, update=True)
    DESC: For MGH path reports, find the 'final diagnosis' line of the path report, remove everything preceding. Parse 
     to extract whether or not there is a final diagnosis line, what that line is, and the full report text (following this line)
    
    PARAMETERS:
    pathdf: pathology dataframe from load_RPDR_path
    update: return only mgh path as a new df or update pathdf with MGH results
    
    RETURNS: pathdf or new path dataframe with only MGH path reports (depending on update bool) with three new columns:
    ['has_final_diagnosis'] = did the function find a line of text that it thinks contains final diagnosis line?
    ['final_diagnosis_line'] = if has_final_diagnosis == True, then what is the final diagnosis line?
    ['Report_Text'] = report text after removing everything above final diagnosis
    '''
    # truncate to only final diagnosis (for MGH reports only)
    print('Filtering only MGH path reports...')
    
    # first get only MGH values
    fil_mgh = pathdf.MRN_Type == 'MGH'
    mgh_path = pathdf[fil_mgh].copy()
    
    # FINAL DIAGNOSIS LINE FINDER
    num_reports = mgh_path.shape[0] # num rows of mgh_path
    has_final_diagnosis_col = []
    final_diagnosis_line_col = []
    trunc_path_col = []
    print('Truncating to only final diagnosis...')
    for i in range(0,num_reports):
        # extract path report for this entry
        report_text = mgh_path.iloc[i,:].Report_Text
        # split by newline character
        text_by_line = report_text.split('\n')

        has_final_diagnosis = False
        final_diagnosis_line = ''
        trunc_path_text = report_text

        # go line-by-line and perform some checks
        j=0
        for line in text_by_line:
            lower_line = line.lower()
            # capture situation where a line contains liver, biopsy; note will only grab first instance then short circuit
            if has_final_diagnosis==False and 'final' in lower_line and 'diagnosis' in lower_line and not 'amend' in lower_line:
                has_final_diagnosis = True
                final_diagnosis_line = line
                trunc_path_text = '\n'.join(text_by_line[j:]) # should be a list of this line and all subsequent lines
                
            j=j+1

        has_final_diagnosis_col.append(has_final_diagnosis)
        final_diagnosis_line_col.append(final_diagnosis_line)
        # either returns the original report or the truncated form if it has a final diagnosis to truncate at
        trunc_path_col.append(trunc_path_text)

    mgh_path['has_final_diagnosis'] = has_final_diagnosis_col
    mgh_path['final_diagnosis_line'] = final_diagnosis_line_col
    mgh_path['Report_Text'] = trunc_path_col


    if update:
        # re-merge with original data
        print('Updating input path dataframe with truncated MGH path reports')
        pathdf['has_final_diagnosis'] = False
        pathdf.update(mgh_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
        print('Returning MGH only entries with truncated path reports')
        return_df = mgh_path
        
        print('Done. | Status: ' + str(mgh_path[mgh_path.has_final_diagnosis==True].shape[0]) + ' reports with a final diagnosis, ' 
              + str(mgh_path[mgh_path.has_final_diagnosis==False].shape[0]) + ' reports with no final diagnosis')
    
    return return_df


def mgh_is_liver_biopsy(pathdf, update=True, only_finaldx=True): 
    ''' mgh_is_liver_biopsy(pathdf, update=True, only_finaldx=True)
    DESC: For MGH path reports, determine whether this pathology report is from a liver biopsy. 
    REQUIRES: a column named Report_Text containing the full text of the path report; if only_finaldx is true, must have run
     mgh_truncate_finaldx() first 
    
    PARAMETERS:
    pathdf: pathology dataframe from load_RPDR_path (or subsequent modification)
    update: return only mgh path as a new df or update pathdf with MGH results
    only_finaldx: bool, if True, only update rows where has_final_diagnosis==True
    
    RETURNS: pathdf or new path dataframe with MGH path reports (depending on update bool) with three new columns:
    ['is_liver_biopsy'] = bool, is this a liver biopsy or not?
    ['is_liver_biopsy_line'] = text of liver biopsy line if above true
    ['liver_biopsy_LAFD'] = LAFD=Lines After Final Diagnosis, ie if final diagnosis line on line 12, and liver biopsy line 16, =16-12
     This is to help screen for oddities, as there is a typical, rough number of lines between these entries.
    '''    
    print('Filtering only MGH path reports...')
    fil_mgh = pathdf.MRN_Type == 'MGH'
    mgh_path = pathdf[fil_mgh].copy()
    
    if only_finaldx:
        # check the column exists first:
        if 'has_final_diagnosis' in mgh_path.columns.tolist():
            fil_finaldx_trunc = mgh_path.has_final_diagnosis == True
            mgh_path = mgh_path[fil_finaldx_trunc]
        else:
            print('The flag *only_finaldx=True* was passed, however mgh_truncate_finaldx() has not been called. Aborting...')
            return None
    
    num_reports = mgh_path.shape[0]
    is_liver_biopsy_col = []
    is_liver_biopsy_line_col = []
    liver_biopsy_LAFD_col = []
    for i in range(0,num_reports):
        # extract path report for this entry
        report_text = mgh_path.iloc[i,:].Report_Text
        # split by newline character
        text_by_line = report_text.split('\n')

        is_liver_biopsy = False
        is_liver_biopsy_line = ''
        liver_biopsy_LAFD = -1

        # go line-by-line and perform some checks
        j=0
        for line in text_by_line:
            lower_line = line.lower()
            # capture situation where a line contains liver, biopsy; note will only grab first instance then short circuit
            if is_liver_biopsy==False and ('liver' in lower_line or 
                                           'hepatic' in lower_line) and ('biopsy' in lower_line or 
                                                                         'biopsies' in lower_line or 
                                                                         'lobe' in lower_line or 
                                                                         'non-focal' in lower_line or 
                                                                         'nonfocal' in lower_line or 
                                                                         'non focal' in lower_line or
                                                                         'tissue' in lower_line or 
                                                                         'random' in lower_line or 
                                                                         'core' in lower_line) and not ('colon' in lower_line or 
                                                                                                       'renal' in lower_line or
                                                                                                       'kidney' in lower_line):
                is_liver_biopsy = True
                is_liver_biopsy_line = line
                liver_biopsy_LAFD = j
                
            j=j+1

        is_liver_biopsy_col.append(is_liver_biopsy)
        is_liver_biopsy_line_col.append(is_liver_biopsy_line)
        liver_biopsy_LAFD_col.append(liver_biopsy_LAFD)

    mgh_path['is_liver_biopsy'] = is_liver_biopsy_col
    mgh_path['is_liver_biopsy_line'] = is_liver_biopsy_line_col
    mgh_path['liver_biopsy_LAFD'] = liver_biopsy_LAFD_col

    if update:
        # re-merge with original data
        print('Updating input path dataframe with truncated MGH path reports')
        pathdf['is_liver_biopsy'] = False
        pathdf['is_liver_biopsy_line'] = ''
        pathdf['liver_biopsy_LAFD'] = -1
        pathdf.update(mgh_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
        print('Returning MGH only entries with truncated path reports')
        return_df = mgh_path
        
    print('Done. | Status: ' + str(mgh_path[mgh_path.is_liver_biopsy==True].shape[0]) + ' reports with likely liver biopsy, ' 
          + str(mgh_path[mgh_path.has_final_diagnosis==False].shape[0]) + ' reports likely not a liver biopsy')
    
    return return_df

def mgh_find_note_start(pathdf, update=True, only_finaldx=True):
    
    if not 'is_liver_biopsy' in pathdf.columns.tolist():
        print('Function requires running mgh_is_biopsy() first (uses relative positioning of the biopsy line to call note start)')
        return None
    
    print('Filtering only MGH path reports...')
    fil_mgh = pathdf.MRN_Type == 'MGH'
    mgh_path = pathdf[fil_mgh].copy()
    
    if only_finaldx:
        # check the column exists first:
        if 'has_final_diagnosis' in mgh_path.columns.tolist():
            fil_finaldx_trunc = mgh_path.has_final_diagnosis == True
            mgh_path = mgh_path[fil_finaldx_trunc]
        else:
            print('The flag *only_finaldx=True* was passed, however mgh_truncate_finaldx() has not been called. Aborting...')
            return None
    
    num_reports = mgh_path.shape[0]
    has_note_start_col = []
    note_line_col = []
    note_start_LAFD_col = []
    for i in range(0,num_reports):
        # extract path report for this entry
        report_text = mgh_path.iloc[i,:].Report_Text
        # split by newline character
        text_by_line = report_text.split('\n')

        has_note = False
        has_note_after_biopsy = False
        note_line = ''
        note_start_LAFD = -1

        # go line-by-line and perform some checks
        j=0
        for line in text_by_line:
            lower_line = line.lower()
            # prep for finding the first non-space word in a line
            word_list = remove_extra_spaces(lower_line, return_as_list=True)
            # first conditional: make sure the word list isn't empty
            # second conditional: if we've already found note, don't look further
            # third and beyond conditionals: that the first or second word contains note is best signature of the start of the note section
            # the parentheses are to help it ignore catching note on a new line where it says '(see note)' or (see comment)
            if len(word_list) > 0 and has_note_after_biopsy==False and ((('note' in word_list[0] or 'comment' in word_list[0]) and not ')' in word_list[0]) or 
                                                           (len(word_list) > 1 and (('note' in word_list[1] or 'comment' in word_list[1]) and not ')' in word_list[1]))):
                
                
                # unlike other loop functions above, stipulate the 'best' note comes after the line declaring this to be a liver biopsy
                #  but if no such line materializes, still go with the FIRST line it was seen
                if j > mgh_path.iloc[i,:].liver_biopsy_LAFD:
                    # stop searching, this is a good candidate
                    has_note_after_biopsy = True
                    has_note = True
                    note_line = line
                    note_start_LAFD = j
                # this is the first encounter with a probable note line (has_note is False); keep it unless the if above is triggered with a better one
                elif has_note == False:
                    has_note = True
                    note_line = line
                    note_start_LAFD = j
                
            j=j+1

        has_note_start_col.append(has_note)
        note_line_col.append(note_line)
        note_start_LAFD_col.append(note_start_LAFD)

    mgh_path['has_note_start'] = has_note_start_col
    mgh_path['note_line'] = note_line_col
    mgh_path['note_start_LAFD'] = note_start_LAFD_col

    if update:
        # re-merge with original data
        print('Updating input path dataframe note starts and locations')
        pathdf['has_note_start'] = False
        pathdf['note_line'] = ''
        pathdf['note_start_LAFD'] = -1
        pathdf.update(mgh_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
        print('Returning MGH only entries annotated note starts')
        return_df = mgh_path
        
    print('Done. | Status: ' + str(mgh_path[mgh_path.has_note_start==True].shape[0]) + ' reports with an identifiable note entry, ' 
          + str(mgh_path[mgh_path.has_note_start==False].shape[0]) + ' reports without a note entry')
    
    return return_df

def mgh_extract_finaldx(pathdf, update=True):
    
    if not 'is_liver_biopsy' in pathdf.columns.tolist() and 'has_note_start' in pathdf.columns.tolist():
        print('Function requires running mgh_is_biopsy() and mgh_find_note_start first (uses relative positioning of the biopsy line to call note start)')
        return None
    
    print('Filtering only MGH path reports that are notated as liver biopsies...')
    fil_mgh = pathdf.MRN_Type == 'MGH'
    fil_liverbx = pathdf.is_liver_biopsy == True
    mgh_path = pathdf[fil_mgh & fil_liverbx].copy()
    
    num_reports = mgh_path.shape[0]
    has_finaldx_col = []
    finaldx_text_col = []
    finaldx_LAFD_col = []
    
    for i in range(0,num_reports):
        # extract path report for this entry
        report_text = mgh_path.iloc[i,:].Report_Text
        # split by newline character
        text_by_line = report_text.split('\n')

        has_finaldx = False
        finaldx_text = ''
        finaldx_LAFD = -1

        # grab markers from where biopsy line is noted and where the 'note' line is
        liverbx_LAFD = int(mgh_path.iloc[i,:].liver_biopsy_LAFD) # every entry should have a positive liver biopsy LAFD
        note_LAFD = int(mgh_path.iloc[i,:].note_start_LAFD) # NOT every entry will have a note; that's OK

        # if the order is liver biopsy line >> note, then the final diagnosis usually falls right in between
        if note_LAFD > liverbx_LAFD:
            finaldx_text = ' '.join(text_by_line[liverbx_LAFD+1:note_LAFD])
            # get rid of extra spaces
            finaldx_text = remove_extra_spaces(finaldx_text)
            
            has_finaldx = True
            finaldx_LAFD = liverbx_LAFD+1
        
        has_finaldx_col.append(has_finaldx)
        finaldx_text_col.append(finaldx_text)
        finaldx_LAFD_col.append(finaldx_LAFD)

    mgh_path['has_finaldx_text'] = has_finaldx_col
    mgh_path['finaldx_text'] = finaldx_text_col
    mgh_path['finaldx_LAFD'] = finaldx_LAFD_col

    if update:
        # re-merge with original data
        print('Updating input path dataframe with final diagnosis (short) text')
        pathdf['has_finaldx_text'] = False
        pathdf['finaldx_text'] = ''
        pathdf['finaldx_LAFD'] = -1
        pathdf.update(mgh_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
        print('Returning MGH only entries passing all filters above and annotated with final diagnosis (short) text')
        return_df = mgh_path
        
    print('Done. | Status: ' + str(mgh_path[mgh_path.has_finaldx_text==True].shape[0]) + ' reports with a probable final diagnosis, ' 
          + str(mgh_path[mgh_path.has_finaldx_text==False].shape[0]) + ' reports without a probable final diagnosis')
    
    return return_df

def remove_extra_spaces(input_as_str, return_as_list=False):
    word_list = input_as_str.split(' ')
    while '' in word_list: word_list.remove('')
    
    if return_as_list:
        out = word_list
    else:
        out = ' '.join(word_list)
    
    return out
