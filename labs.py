
def load_RPDR_labs_multiple(dir_data, path_labs, path_synonyms, datetime_col='Seq_Date_Time', result_col='Result', test_col='Group_Id', delim='|', clean_columns=True, return_discarded=False):
    ''' load_labs_multiple(dir_data, path_labs, path_synonyms, datetime_col='Seq_Date_Time', result_col='Result', test_col='Group_Id', delim='|', clean_columns=True):
        Sequentially loads all files from RPDR data dump when multiple files have the same name. 
        
        1. Starts in dir_data (should have trailing slash), grabs all sub-folders' names automatically, then sequentially loads: dir_data/[sub-folders]/path_labs (where path_labs is the name of the file)
        * Note for whatever reason, on a multiple-split file dump from RPDR the labs, demographics, etc files are all named the exact same, just in different zips
        2. Calls the traditional load LFTs function on each file
        3. Concatenates all results and returns 1 DF (may be VERY large. ~3-4 gb / 80,000 patients)
        
        See load_native_data for remainder of parameters which are passed to that function
        
        '''
    import os
    import pandas as pd
    
    # get list of subdirectories
    subdirectories = [x[0] for x in os.walk(dir_data)][1:]
    
    first=True
    # for each subdir, use the traditional load function to load data and concat
    for subdir in subdirectories:
        path_to_labs=subdir+'/'+path_labs
        lfts, discarded = load_RPDR_labs(path=path_to_labs,
                                           path_synonyms=path_synonyms,
                                           datetime_col=datetime_col,
                                           result_col=result_col,
                                           test_col=test_col, 
                                           delim=delim,
                                           clean_columns=clean_columns)
        
        if first==True:
            concat_pd = lfts
            if return_discarded:
                concat_discard = discarded
            first=False
        else:
            concat_pd=pd.concat([concat_pd, lfts],ignore_index=True)
            if return_discarded:
                concat_discard=pd.concat([concat_discard, discarded],ignore_index=True)
    
    return concat_pd, 

def load_RPDR_labs(path,path_synonyms,datetime_col='Seq_Date_Time', result_col='Result', test_col='Group_Id', delim='|', clean_columns=True):
    '''load_RPDR_labs(path,path_synonyms,datetime_col='Seq_Date_Time', result_col='Result', test_col='Group_Id', delim='|', clean_columns=True):
    
    DESC: This is the main labs loading function, for loading labs data from an RPDR data pull, and processing lab results. 
     Does MANY THINGS. 
     1. Homogenizes column names:
         rename result_col to 'result' if not already named result
         rename datetime_col to 'datetime' if not already named so; convert this to pd.DateTime format
         rename test_col to 'test_desc'
    path: path to labs file (in .csv or other text delimited file)
    path_synonyms: path to synonyms file (.csv). structure of synonyms file: first row is what you want to call each lab 
     defined in the test_col column; remainder of rows identify the exact text of labs that should be considered equivalent. 
     All of these will be homogenized to the row 1 name
    datetime_col: name of column with date/time of lab test
    result_col: name of column that contains the result of the lab test
    test_col: name of the column that contains the name of the lab test (this is what synonyms acts upon)
    delim: delimiter of the labs file (the file pointed to by path)
    clean_columns: bool. If true, removes all columns except those named the following:
     [['EMPI', 'test_desc', 'datetime', 'result_orig', 'result', 'MRN_Type', 'MRN', 'Test_Description', 'Result_Text', 'above_below_assay_limit']]
    
    RETURNS: lfts, discarded_rows
     lfts is the labs as panda database (should really be named labs, but initially developed for only liver function tests (LFTs))
     discarded_rows is a second pandas database containing anything removed by the processing
     
    Note: input result column will be returned as 'result_orig' and the processed result column will be returned as 'result'
    '''
    # read data
    # 
    # WARNINGS:
    # - default behavior is to remove < and >; could add a new column to mark these events if pertinent later
    import pandas as pd
    import re
    
    # read the data
    print('Loading data from ' + path)
    lfts = pd.read_csv(path,delimiter=delim, error_bad_lines=False)
    
#     if result_col != 'result':
#         lfts['result'] = lfts.loc[:,result_col]
#     else:
#         # result col is already named result. rename it
#         lfts['result_orig'] = lfts['result']

    ## Mod 1 ##
    if result_col in lfts.columns:
        lfts.rename(columns={result_col: 'result_orig'}, inplace=True)
        lfts['result'] = lfts.loc[:,'result_orig']
    else:
        raise ValueError('Incorrect name specified for result column')
    ## Mod 1 ##
    
    # convert datetime column to datetime format
    if datetime_col != 'datetime':
        lfts['datetime'] = pd.to_datetime(lfts.loc[:,datetime_col])
        lfts.drop(labels=datetime_col, axis=1, inplace=True)
    else:
        lfts['datetime'] = pd.to_datetime(lfts['datetime'])
    
    if test_col != 'test_desc':
        lfts['test_desc'] = lfts.loc[:,test_col]
    else:
        # test_desc col is already named test_desc. rename it
        lfts['test_desc_orig'] = lfts['test_desc']
    
    ######################
    # Identify synonyms of labs and homogenize based on a file called lab_synonyms.csv
    # structure: first row is what you want to call each lab; below the col are all the synonyms to change to the header name
    syn=pd.read_csv(path_synonyms)
    correct_lab_names=syn.columns.tolist()
    for correct_lab in correct_lab_names:
        alternate_names_to_replace=syn[correct_lab].tolist()
        for alt in alternate_names_to_replace:
            # irregular length columns, if this column has nan's at the end, don't cycle through them
            if isinstance(alt,str):
                lfts.test_desc.replace(alt,correct_lab,regex=False, inplace=True)
            
    # error checking
    if set(correct_lab_names) == set(lfts.test_desc.unique().tolist()):
        print('Successful homogenization of lab names')
        print(correct_lab_names)
    else:
        got_list=set(lfts.test_desc.unique().tolist())
        expected_list=set(correct_lab_names)
        print('FAILED TO HOMOGENIZE NAMES')
        print('Expected : ')
        print(expected_list)
        print('Of these, got only : ')
        print(got_list & expected_list)
        print('Got additionally : ')
        print(got_list-expected_list)
        print('...and we are missing : ')
        print(expected_list-got_list)
    
    ######################
    # CERTIFY INDIVIDUAL RESULTS THAT REQUIRE SPECIAL HANDLING
    # NUMERICS
    list_numeric_labs = ['ALT','AST','AKP','DBILI','TBILI', 'CR','INR','IGG']

    lfts['above_below_assay_limit'] = 0

    fil = lfts.test_desc.isin(list_numeric_labs)

    # def upper and lower bound finder
    def upper_lower_bound_finder(row):
        if not isinstance(row.result, float):
            if '<' in row.result:
                return -1
            elif '>' in row.result:
                return 1
            else:
                return 0

    # def upper and lower bound finder
    def remove_gg_ll_signs(row):
        if isinstance(row, str):
            if '<' in row:
                return row.replace('<','')
            elif '>' in row:
                return row.replace('>','')
            else:
                return row

    # first mark whether a result is at the upper or lower limit (otherwise NaN)
    lfts.loc[fil, 'above_below_assay_limit'] = lfts.loc[fil, ['result', 'above_below_assay_limit']].apply(upper_lower_bound_finder, axis=1)

    # remove greater than and less than signs
    lfts.loc[fil, 'result'] = lfts.loc[fil, 'result'].apply(remove_gg_ll_signs)
    
    # now that we've removed greater than and less than signs, data should be numeric. Anything non-numeric (in this category) is crap, 
    # e.g., "REFUSED", "HEMOLYZED" etc
    
    # filter for coercing data to numeric, catching NaNs (text)
    print('Removing non-numeric result values...')
    lfts.loc[fil, 'result'] = pd.to_numeric(lfts.loc[fil, 'result'], errors='coerce')
    # hold on to removed rows -- error checking
    removed_rows=lfts[(fil) & lfts.result.isnull()].copy()
    # remove these rows
    lfts.drop(lfts[(fil) & lfts.result.isnull()].index, inplace=True)
    
    
    ######################
    # CERTIFY INDIVIDUAL RESULTS THAT REQUIRE SPECIAL HANDLING
    # ANTIBODY TITERS
    
    list_titer_labs = ['LKM','ASMA', 'ANA', 'AMA', 'SLA']

    fil_ab = lfts.test_desc.isin(list_titer_labs)
    
    if not lfts[fil_ab].empty:

        def titer_interp(row):
            import math
            import re

            bool_switched_result_resultext=False

            if isinstance(row.result, float):
                if row.result < 20: # antibody titer cutoff for LKM and works for ANA, ASMA
                    return 0
                elif math.isnan(row.result):
                    if isinstance(row.Result_Text,str):
                        this_string=row.Result_Text
                        bool_switched_result_resultext=True
                    else:
                        return row.result
                else:
                    return row.result

            # function I'm going to use a bunch of times:
            def all_numbers_within_n_chars_of_word(text, word, n):
                idx_zero_pos_of_word=text.lower().find(word)
                if len(text[idx_zero_pos_of_word:]) < n:
                    sub_string_to_return=text[idx_zero_pos_of_word:]
                else:
                    sub_string_to_return=text[idx_zero_pos_of_word:idx_zero_pos_of_word+n]

                return find_all_num_in_string(sub_string_to_return)

            # CASES AT THIS POINT. Either this_string holds a swapped text or doesn't. Any string in result hasn't been accounted for yet
            # most of the logic will be shared, nevertheless keep them separate because, e.g., a single number in result is usually the value we want
            #  while a single number in result text might not be

            if isinstance(row.result, str):
                numbers = find_all_num_in_string(row.result)

                if '<' in row.result:
                    return 0
                elif 'ANA SCREEN POSITIVE, TITRE PERFORMED' in row.result or 'POSITIVE - SEE TITER' in row.result:
                    # this category will get removed; there's always a titer separately reported at same time point
                    return 'see titer'
                
                    
                ## FIRST IDENTIFY POSITIVELY OR NEGATIVELY IDENTIFIED EXAMPLES
                elif 'negative' in row.result.lower() and not 'positive' in row.result.lower():
                    return 0
                elif 'positive' in row.result.lower() and not 'negative' in row.result.lower():
                    
                    if len(numbers) == 1:
                        return numbers[0]
                    elif len(numbers) == 2 and numbers[0]==1:
                        return numbers[1]
                    
                    elif len(numbers) == 4 and numbers[2]==1:
                        return numbers[3]
                    
                    else: 
                        result_text=row.Result_Text
                        if isinstance(result_text,str):
                            
                            ## Some cases doenst contain the word 'positive' in Result_text which is generating a random incorrect result
                            if 'positive' in result_text.lower():
                                numbers=all_numbers_within_n_chars_of_word(result_text, 'positive', n=26)
                                if len(numbers) == 1:
                                    return numbers[0]
                                elif len(numbers) == 2 and numbers[0]==1:
                                    return numbers[1]
                                elif len(numbers) > 2:
                                    this_string=row.Result_Text.lower()
                                    if isinstance(this_string,str):
                                        idx_titer=this_string.find('titer')
                                        if len(row.Result_Text[idx_titer:]) <= 17:
                                            search_text=row.Result_Text[idx_titer:]
                                        else:
                                            search_text=row.Result_Text[idx_titer:idx_titer+17]
                                        numbers=find_all_num_in_string(search_text)
                                        if len(numbers) == 1:
                                            return numbers[0]
                                        elif len(numbers) == 2 and numbers[0]==1:
                                            return numbers[1]
                            else:
                                return 20

                elif 'positive' in row.result.lower() and 'negative' in row.result.lower():
                    # both pos and neg present; find the text following the 'positive' word (26 chars)
                    numbers=all_numbers_within_n_chars_of_word(row.result, 'positive', n=26)
                    if len(numbers) == 1:
                        return numbers[0]
                    elif len(numbers) == 2 and numbers[0]==1:
                        return numbers[1]
                    else:
                        # possible it's farther out?
                        text=row.result[index_positive:]
                        numbers=find_all_num_in_string(text)
                        print('pos & neg text present but no value within 26 chars; went out on a limb and returning : ' + str(numbers[1]))
                        return numbers[1]

                # okay, so the words 'positive' and 'negative' aren't present
                elif len(numbers) == 1:
                    return numbers[0]
                elif len(numbers) == 2 and numbers[0]==1:
                    # pair of numbers returned; make sure the first value is 1
                    return numbers[1]

            # CASES AT THIS POINT. Either this_string holds a swapped text or doesn't. Any interpretable string has been returned

            if bool_switched_result_resultext:
                numbers = find_all_num_in_string(this_string)

                ## FIRST IDENTIFY POSITIVELY OR NEGATIVELY IDENTIFIED EXAMPLES
                if 'negative' in this_string.lower() and not 'positive' in this_string.lower():
                    return 0
                elif 'positive' in this_string.lower() and not 'negative' in this_string.lower():
                    ## Positive at logic is not working when there are 4 numbers. Removed the if statement
                    ## Example: Positive at 1:40 and 1:160 (endpoint)
                    ## Finding all numbers within a close proximity of 'positive' seems to be working for all cases
                    numbers=all_numbers_within_n_chars_of_word(this_string, 'positive', n=26)
                    if len(numbers) == 1:
                        return numbers[0]
                    elif len(numbers) == 2 and numbers[0]==1:
                        return numbers[1]
                    elif len(numbers) == 4 and numbers[2]==1:
                        return numbers[3]
                    elif len(numbers)>0:
                        # first try to match the expression 'positive at 1:x'
                        m = re.search('positive at 1:(\d+)', this_string.lower())
                        if m:
                            return m.group(1)
                        else:
                            # no 'positive at' expressions.. get all the numbers and return the largest one, and print
                            #  result text because this is a weird case..
                            max_num = max(numbers)
                            print('positive only in result text but neither 1 nor 2 numbers, going on limb and taking : ' + str(max_num))
                            print(this_string)
                            return max_num
                    
                elif 'positive' in this_string.lower() and 'negative' in this_string.lower():
                    # both pos and neg present; find the text following the 'positive' word (15 chars)
                    index_positive=this_string.lower().find('positive')
                    if len(this_string[index_positive:]) <= 26:
                        text=this_string[index_positive:]
                    else:
                        text=this_string[index_positive:index_positive+26]
                    numbers=find_all_num_in_string(text)
        #             print(text)
        #             print(numbers)
                    if len(numbers) == 1:
                        return numbers[0]
                    elif len(numbers) == 2 and numbers[0]==1:
                        return numbers[1]
                    elif len(numbers) == 4:
                        return max(numbers)
                    elif 'negative at' in this_string.lower():
                        return 0
                # okay, so the words 'positive' and 'negative' aren't present
                elif len(numbers) == 1:
                    return numbers[0]
                elif len(numbers) == 2 and numbers[0]==1:
                    return numbers[1]


            return row.result

        # first mark whether a result is at the upper or lower limit (otherwise NaN)
        lfts.loc[fil_ab, 'result'] = lfts.loc[fil_ab, ['result', 'Result_Text']].apply(titer_interp, axis=1)

        print('Cleaning up antibody titer results...')
        lfts.loc[fil_ab, 'result'] = pd.to_numeric(lfts.loc[fil_ab, 'result'], errors='coerce')
        # hold on to removed rows -- error checking
        removed_rows_ab=lfts[(fil_ab) & lfts.result.isnull()].copy()
        # remove these reows
        lfts.drop(lfts[(fil_ab) & lfts.result.isnull()].index, inplace=True)
    
        ######################
        # COMBINE DISCARDED DATA
        removed_rows = pd.concat([removed_rows,removed_rows_ab],axis=0)
    
    ######################
    # REMOVE UNUSED/DISTRACTING/PRECURSOR COLUMNS
    if clean_columns:
        # easier to define as what we want to keep
        lfts = lfts[['EMPI', 'test_desc', 'datetime', 'result_orig', 'result', 'MRN_Type', 'MRN', 'Test_Description', 'Result_Text', 'above_below_assay_limit']].copy()
        
#     from IPython import embed; embed()
        
    # enforce the EMPI column is strings for later
    lfts['EMPI'] = lfts.EMPI.astype(str)
        
    print('...Done')
    
    return lfts, removed_rows

def find_all_num_in_string(sentence):
    # accepts a sentence; returns an array of all numbers in the sentence
    import re
    
    s = [float(s) for s in re.findall(r'-?\d+\.?\d*', sentence)]
    return s
