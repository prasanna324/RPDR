
def load_RPDR_path(path,delimiter='\t', datetime_col='Report_Date_Time'):
    ''' load_RPDR_path(path,delimiter='\t', datetime_col='Report_Date_Time')
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
    - 
    '''# FYIs:
    # - removes everything but MGH for now
    import pandas as pd
    
    path_df = pd.read_csv(path, sep=delimiter, dtype=str)
    
    # remove deleted, cancelled, and in-process path reports
    bad_reports = (path_df.Report_Status == 'Deleted') | (path_df.Report_Status == 'Cancelled') | (path_df.Report_Status == 'In Process')
    
    # drop rows with any of these features
    path_df2 = path_df[~bad_reports].copy()
    
    # re-index to 'Report_Number', dropping duplicates ONLY if duplicate on both report_number and Report_Text
    # (ie, there is a problem if they're not the same..)
    path_df2.drop_duplicates(subset=['Report_Number','Report_Text'], keep='first', inplace=True, ignore_index=False)
    
    # now set index to Report_Number
    path_df2.set_index(keys='Report_Number', inplace=True, verify_integrity=True)
    
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
