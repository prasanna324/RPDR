
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


def is_colitis(pathdf, corpus='en_core_sci_lg', term_set='en_clinical', update=True):

    '''
    corpus: en_ner_bc5cdr_md, en_core_sci_md, en_core_sci_lg
    termset: en, en_clinical, en_clinical_sensitive

    '''
    
    fil_subset = pathdf.MRN_Type.isin(['MGH', 'BWH','NWH','FH','NSM'])
    df_path = pathdf[fil_subset].copy()


    import spacy
    from negspacy.negation import Negex
    from negspacy.termsets import termset
    #from spacy.pipeline import EntityRuler

    ts = termset(term_set)

    config={
        "neg_termset":{
            "pseudo_negations": ts.terms['pseudo_negations'] + ['and stage', 'grade', 'active'],
            "preceding_negations": ts.terms['preceding_negations'] + ['negative'],
            "following_negations": ts.terms['following_negations'] + ['negative', 'is not', 'are not', 'does not', 'may not', 'have not', 'was not', 'were not', 'absent', 'not present'],
            "termination": ts.terms['termination'] + ['note:', 'moderate']
        }
    }


    nlp_2 = spacy.load(corpus) 

    # ruler = EntityRuler(nlp_2, overwrite_ents=True)
    # patterns = [
    #     {"label": "ENTITY", "pattern": [{"LOWER": "chronic inflammation"}]}
    #         ]
    # ruler.add_patterns(patterns)

    nlp_2.add_pipe(
        "negex",
        config = config
    )
    

    num_reports = df_path.shape[0]
    colitis_col = []
    chronic_colitis_col = []
    mild_colitis_col = []
    moderate_colitis_col = []
    severe_colitis_col = []
    inactive_colitis_col = []
    active_colitis_col= []
    acute_colitis_col = []
    disease_list_col = []

    for i in range(0,num_reports):

        # extract path report for this entry
        disease_list = []
        report_text = df_path.iloc[i,:].Report_Text
        result_text = entity_recognition_2(report_text, nlp=nlp_2)


        colitis = False
        chronic_colitis = False
        mild_colitis = False
        moderate_colitis = False
        severe_colitis = False
        inactive_colitis = False
        active_colitis = False
        acute_colitis = False


        for x in result_text.split('\n'):
            if 'colitis' in x and 'True' in x:
                colitis = True
            if 'chronic' in x and 'colitis' in x and 'True' in x:
                chronic_colitis = True
            if 'mild' in x and 'colitis' in x and 'True' in x:
                mild_colitis = True
            if 'moderate' in x and 'colitis' in x and 'True' in x:
                moderate_colitis = True
            if 'severe' in x and 'colitis' in x and 'True' in x:
                severe_colitis = True
            if 'active' in x and 'colitis' in x and 'True' in x:
                active_colitis = True
            if 'inactive' in x and 'colitis' in x and 'True' in x:
                inactive_colitis = True
            if 'acute' in x and 'colitis' in x and 'True' in x:
                acute_colitis = True

            if 'colitis' in x and 'True' in x:
                disease_list.append(x)
        
        colitis_col.append(colitis)
        chronic_colitis_col.append(chronic_colitis)
        mild_colitis_col.append(mild_colitis)
        moderate_colitis_col.append(moderate_colitis)
        severe_colitis_col.append(severe_colitis)
        active_colitis_col.append(active_colitis)
        inactive_colitis_col.append(inactive_colitis)
        acute_colitis_col.append(acute_colitis)
        disease_list_col.append(disease_list)
        
    df_path['colitis'] = colitis_col
    df_path['chronic_colitis'] = chronic_colitis_col
    df_path['mild_colitis'] = mild_colitis_col
    df_path['moderate_colitis'] = moderate_colitis_col
    df_path['severe_colitis'] = severe_colitis_col
    df_path['active_colitis'] = active_colitis_col
    df_path['inactive_colitis'] = inactive_colitis_col
    df_path['acute_colitis'] = acute_colitis_col
    df_path['disease_list'] = disease_list_col
   
    if update:
        # re-merge with original data
        print('Updating input path dataframe')
        pathdf['colitis'] = np.nan
        pathdf['chronic_colitis'] = np.nan
        pathdf['mild_colitis'] = np.nan
        pathdf['moderate_colitis'] = np.nan
        pathdf['severe_colitis'] = np.nan
        pathdf['active_colitis'] = np.nan
        pathdf['inactive_colitis'] = np.nan
        pathdf['acute_colitis'] = np.nan
        pathdf['disease_list'] = np.nan
        pathdf.update(df_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
        print('Returning MGH, BWH only entries with truncated path reports')
        return_df = df_path
        

    return return_df


def entity_recognition_2(text, nlp):
    
    text = text.lower()
    
    entity_result = ''
    
    for line in text.split('.'):
        
        #line = line.strip()
        line = " ".join(line.split())
        line = (line
#                 .replace('chronic active', 'chronic-active')
#                 .replace('active chronic', 'active-chronic')
#                 .replace('chronic inactive', 'chronic-inactive')
                .replace('severely', 'severe')
                .replace('moderately', 'moderate')
                .replace('mildly', 'mild').replace('mildl', 'mild')
                .replace('floridly', 'florid')
                .replace('severe pseudomembranous', 'severe-pseudomembranous')
                .replace('self limited', 'self-limited')
                .replace('moderate to severe', 'moderate&severe')
                .replace('moderate to focally severe', 'moderate&severe')
                .replace('mild to moderate', 'mild&moderate')
                .replace('mild to focally moderate', 'mild&moderate')
                .replace('mild to severe', 'mild&severe')
                .replace('mild to focally severe', 'mild&severe')
                .replace('chronic-inactive', 'chronic inactive')
                .replace('active chronic', 'active-chronic')
                .replace('acute and chronic', 'acute-chronic')
                .replace('acute on chronic', 'acute-chronic')
                .replace('severe active', 'severe-active')
                .replace('severe chronic', 'severe-chronic')
                .replace('active severe', 'active-severe')
                .replace('chronic severe', 'chronic-severe')
                .replace('pancolitis, moderate&severe', 'moderate&severe pancolitis')
                .replace('colitis, moderate&severe', 'moderate&severe colitis')
                .replace('colitis, severe', 'severe colitis')
                .replace('active-severe', 'severe-active')
                .replace('colitis, moderate', 'moderate colitis')
                .replace('active-moderate', 'moderate-active')
                .replace('severe ischemic', 'severe-ischemic')
                .replace('moderate active', 'moderate-active')
                .replace('moderate chronic', 'moderate-chronic')
                .replace('active moderate', 'active-moderate')
                .replace('chronic moderate', 'chronic-moderate')
                .replace('mild active', 'mild-active')
                .replace('mild chronic', 'mild-chronic')
                .replace(' active colitis', ' active-colitis')
                .replace('mild ', 'mild-')
                .replace('moderate ', 'moderate-')
                .replace('severe ', 'severe-')
                .replace('inactive ', 'inactive-')
                .replace('ulcerative colitis', 'ulcerative-colitis')
                .replace('healed colitis', 'healed-colitis')
                .replace('surveillance', 'not present')
#                 .replace('chronic active', 'chronic-active')
#                 .replace('inactive chronic', 'inactive-chronic')
# # #                 .replace('active colitis', 'active-colitis')
#                 .replace('inactive colitis', 'inactive-colitis')
               )
                
        
        global doc, e
            
        doc = nlp(line)
    
        for e in doc.ents:
            
            e_text = e.text
            e_text = re.sub(' +', ' ', e_text)
            e_bool = e._.negex
            
            # Replace negation words in the entity and adjust sentiment
            if e_text.startswith(('no ', 'non-', 'non ')):
                to_match = ['^no ', '^non-', '^non ']
                e_text = re.sub('|'.join(to_match), '', e_text)
                e_bool = not e_bool
            
#             chronic_col = bool(re.search(r'\b(?:chronic\W+(?:\w+\W+){0,1}?(colitis|pancolitis))\b', line))
#             active_col = bool(re.search(r'\b(?:active\W+(?:\w+\W+){0,1}?colitis)\b', line))
#             inactive_col = bool(re.search(r'\b(?:inactive\W+(?:\w+\W+){0,1}?colitis)\b', line))
            
#             if chronic_col and active_col and e_text=='colitis':
#                 entity_result = entity_result + 'chronic active ' + e_text + str(not e_bool) + '\n'
#             elif chronic_col and inactive_col and e_text=='colitis':
#                 entity_result = entity_result + 'chronic inactive ' + e_text + str(not e_bool) + '\n'
#             elif chronic_col and e_text=='colitis':
#                 entity_result = entity_result + 'chronic ' + e_text + str(not e_bool) + '\n'
#             elif active_col and e_text=='colitis':
#                 entity_result = entity_result + 'active ' + e_text + str(not e_bool) + '\n'
#             elif inactive_col and e_text=='colitis':
#                 entity_result = entity_result + 'inactive ' + e_text + str(not e_bool) + '\n'
#             else:

            e_text = " ".join(e_text.split())
            entity_result = entity_result + e_text + ' ' + str(not e_bool) + '\n'
                
        
    return entity_result
    
