
def is_colitis(pathdf, corpus='en_core_sci_lg', term_set='en_clinical', update=True, only_truncated=False):

    '''
    corpus: en_ner_bc5cdr_md, en_core_sci_md, en_core_sci_lg
    termset: en, en_clinical, en_clinical_sensitive

    '''
    
    fil_subset = pathdf.MRN_Type.isin(['MGH', 'BWH','NWH','FH','NSM'])
    df_path = pathdf[fil_subset].copy()
    
    if only_truncated:
        # check the column exists first:
        if 'has_final_diagnosis' in df_path.columns.tolist():
            fil_finaldx_trunc = df_path.has_final_diagnosis == True
            df_path = df_path[fil_finaldx_trunc]
        else:
            print('The flag *only_finaldx=True* was passed, however truncate_finaldx() has not been called. Aborting...')
            return None
    
    import spacy
    from negspacy.negation import Negex
    from negspacy.termsets import termset
    import numpy as np
    import pandas as pd
    import re
    #from spacy.pipeline import EntityRuler

    ts = termset(term_set)

    config={
        "neg_termset":{
            "pseudo_negations": ts.terms['pseudo_negations'] + ['and stage', 'grade', 'active'],
            "preceding_negations": ts.terms['preceding_negations'] + ['negative'],
            "following_negations": ts.terms['following_negations'] + ['negative', 'unremarkable', 'is not', 'are not', 'does not', 'may not', 'have not', 'was not', 'were not', 'absent', 'not present'],
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
    mayo_score_col = []
    disease_list_col = []

    for i in range(0,num_reports):

        # extract path report for this entry
        disease_list = []
        report_text = df_path.iloc[i,:].Report_Text
        result_text = entity_recognition_colon(report_text, nlp=nlp_2)


        colitis = False
        chronic_colitis = False
        mild_colitis = False
        moderate_colitis = False
        severe_colitis = False
        inactive_colitis = False
        active_colitis = False
        acute_colitis = False
        mayo_score = np.nan


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
            if 'colitis mayo-' in x:
                mayo_score = float(re.findall(r'.*?(\d+(?:,\d+)*(?:\.\d+)?)', x)[0])
                if mayo_score>3:
                    mayo_score = 3
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
        mayo_score_col.append(mayo_score)
        disease_list_col.append(disease_list)
        
    df_path['colitis'] = colitis_col
    df_path['chronic_colitis'] = chronic_colitis_col
    df_path['mild_colitis'] = mild_colitis_col
    df_path['moderate_colitis'] = moderate_colitis_col
    df_path['severe_colitis'] = severe_colitis_col
    df_path['active_colitis'] = active_colitis_col
    df_path['inactive_colitis'] = inactive_colitis_col
    df_path['acute_colitis'] = acute_colitis_col
    df_path['mayo_score'] = mayo_score_col
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
        pathdf['mayo_score'] = np.nan
        pathdf['disease_list'] = np.nan
        pathdf.update(df_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
        print('Returning MGH, BWH only entries with truncated path reports')
        return_df = df_path
        

    return return_df


def entity_recognition_colon(text, nlp):
    
    import re
    import numpy as np
    
    text = text.replace(' III ', '3').replace(' II ', '2').replace(' IV ', '4')
    text = text.lower()
    
    entity_result = ''
    mayo_score = -1
    mayo_bool = False
    
    for line in text.split('.'):
        
        #line = line.strip()
        line = " ".join(line.split())
        line = (line
                .replace('neither', 'no').replace('nor', 'no')
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
                .replace('mayo i ', 'mayo 1 ')
                .replace('grade i ', 'grade 1 ')
#                 .replace('mayo ', 'mayo-')
#                 .replace('chronic active', 'chronic-active')
#                 .replace('inactive chronic', 'inactive-chronic')
# # #                 .replace('active colitis', 'active-colitis')
#                 .replace('inactive colitis', 'inactive-colitis')
               )
                
        #global doc, e
            
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
        
        if 'mayo' in line:
            line = (line
                    .replace('0-1', '0.5').replace('1-2', '1.5').replace('2-3', '2.5').replace('3-4', '3.5')
                   )
        
            mayo_list = re.findall(r'mayo.*?(\d+(?:,\d+)*(?:\.\d+)?)', line)
            grade_list = re.findall(r'grade.*?(\d+(?:,\d+)*(?:\.\d+)?)', line)
            mayo_list = mayo_list + grade_list
            mayo_list = [z for z in mayo_list if 0<=float(z)<=4]
            
            if len(mayo_list)!=0 and float(max(mayo_list))>mayo_score:
                mayo_score = float(max(mayo_list))
                mayo_bool = True
    
    if mayo_bool==True:
        entity_result = entity_result + 'colitis mayo-' + str(mayo_score) + ' ' + str(True) + '\n'
        
    return entity_result


def is_liver_disease(pathdf, corpus='en_core_sci_lg', term_set='en_clinical', update=True, only_liv_biopsy=True):

    '''
    corpus: en_ner_bc5cdr_md, en_core_sci_md, en_core_sci_lg
    termset: en, en_clinical, en_clinical_sensitive

    '''
    
    fil_subset = pathdf.MRN_Type.isin(['MGH', 'BWH','NWH','FH','NSM'])
    df_path = pathdf[fil_subset].copy()
    
    
    if only_liv_biopsy:
        # check the column exists first:
        if 'is_liver_biopsy' in df_path.columns.tolist():
            fil_finaldx_trunc = df_path.is_liver_biopsy == True
            df_path = df_path[fil_finaldx_trunc]
        else:
            print('The flag *only_finaldx=True* was passed, however truncate_finaldx() has not been called. Aborting...')
            return None
        
    filter_keywords = df_path['Report_Text'].str.contains('steato|balloon|baloon|ballon|inflammat|hepatitis|hepatic', case=False, na=False)
    df_path = df_path[filter_keywords]

    import spacy
    from negspacy.negation import Negex
    from negspacy.termsets import termset
    import numpy as np
    import pandas as pd
    #from spacy.pipeline import EntityRuler

    ts = termset(term_set)

    config={
        "neg_termset":{
            "pseudo_negations": ts.terms['pseudo_negations'] + ['and stage', 'grade'],
            "preceding_negations": ts.terms['preceding_negations'] + ['negative'], #'grade 0'
            "following_negations": ts.terms['following_negations'] + ['negative', 'unremarkable', 'is not', 'are not', 'does not', 'may not', 'have not', 'was not', 'were not', 'absent', 'grade 0'],
            "termination": ts.terms['termination'] + ['note:', 'with', ';', ', negative', ',negative'] #'negative for'
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
    steatosis_col = []
    ballooning_col = []
    inflammation_col = []
    lobular_inflammation_col = []
#     portal_inflammation_col = []
    zone3_inflammation_col = []
    lobular_hepatitis_col= []
    zone3_hepatitis_col = []
    fibrosis_col = []
    
    fibrosis_stage_col = []
    bridging_fibrosis_col = []
    sinusoidal_fibrosis_col = []
    portal_fibrosis_col = []
    periportal_fibrosis_col = []
    pericellular_fibrosis_col = []
    perivenular_fibrosis_col = []
    septal_fibrosis_col = []
    central_fibrosis_col = []
    hepatitis_col = []
    autoimmune_hepatitis_col = []
    
    disease_list_col = []

    for i in range(0,num_reports):

        # extract path report for this entry
        disease_list = []
        report_text = df_path.iloc[i,:].Report_Text
        result_text = entity_recognition_liver(report_text, nlp=nlp_2)

        steatosis = False
        ballooning = False
        inflammation = False
        lobular_inflammation = False
#         portal_inflammation = False
        zone3_inflammation = False
        lobular_hepatitis = False
        zone3_hepatitis = False
        fibrosis = False
        
        fibrosis_stage = np.nan
        bridging_fibrosis = False
        sinusoidal_fibrosis = False
        portal_fibrosis = False
        periportal_fibrosis = False
        pericellular_fibrosis = False
        perivenular_fibrosis = False
        septal_fibrosis = False
        central_fibrosis = False
        hepatitis = False
        autoimmune_hepatitis = False
        
        steatosis_lt5 = sum([1 for ent_text in result_text.split('\n') if '<5%' in ent_text and 'steatosis' in ent_text])==0
        
        fib_stage = -1
        fib_ref = -1
        
        for x in result_text.split('\n'):
            if 'steatosis' in x and 'True' in x and steatosis_lt5:
                steatosis = True
            if ('balloon' in x or 'baloon' in x or 'ballon' in x) and 'True' in x:
                ballooning = True
            if 'inflammation' in x and 'True' in x and inflammation==False:
                inflammation = True
            if 'lobular' in x and 'inflammation' in x and 'True' in x:
                lobular_inflammation = True
#             if 'portal' in x and 'inflammation' in x and 'True' in x:
#                 portal_inflammation = True
            if 'zone-3' in x and 'inflammation' in x and 'True' in x:
                zone3_inflammation = True
            if 'lobular' in x and 'hepatitis' in x and 'True' in x:
                lobular_hepatitis = True
            if 'zone-3' in x and 'hepatitis' in x and 'True' in x:
                zone3_hepatitis = True
            
            if 'fibrosis' in x and 'True' in x:
                fibrosis = True
                
                if 'sinusoidal' in x:
                    sinusoidal_fibrosis = True
                if 'portal' in x and not 'peri-portal' in x:
                    portal_fibrosis = True
                if 'peri-portal' in x:
                    periportal_fibrosis = True
                if 'pericellular' in x:
                    pericellular_fibrosis = True
                if 'perivenular' in x:
                    perivenular_fibrosis = True
                if 'septal' in x:
                    septal_fibrosis = True
                if 'central' in x:
                    central_fibrosis = True
                    
                    
            if ('fibrosis' in x or 'bridging' in x or 'cirrhosis' in x) and 'True' in x and ' stage:' in x:
                fib_stage = float(x[-12:-9])
                fib_ref = float(x[-8:-5])
                
                if fib_ref<4.0:
                    fib_ref = 4.0
                
                fibrosis_stage = str(fib_stage) + '/' + str(fib_ref)
                if 'ishak' in x:
                    fibrosis_stage = 'ishak- ' + fibrosis_stage
            
            if 'bridging' in x and 'True' in x and not 'bridging-necrosis' in x:
                bridging_fibrosis = True
            
            if 'hepatitis' in x and 'True' in x:
                hepatitis = True
                
            if 'autoimmune hepatitis' in x and 'True' in x:
                autoimmune_hepatitis = True

            if (('steatosis' in x and not '<xf5%' in x and steatosis_lt5) or 'balloon' in x or 'baloon' in x or 'ballon' in x
                or 'inflammation' in x or 'hepatitis' in x or 'fibrosis' in x or 'bridging') and 'True' in x:
                disease_list.append(x)
        
        steatosis_col.append(steatosis)
        ballooning_col.append(ballooning)
        inflammation_col.append(inflammation)
        lobular_inflammation_col.append(lobular_inflammation)
#         portal_inflammation_col.append(portal_inflammation)
        zone3_inflammation_col.append(zone3_inflammation)
        lobular_hepatitis_col.append(lobular_hepatitis)
        zone3_hepatitis_col.append(zone3_hepatitis)
        fibrosis_col.append(fibrosis)
        
        fibrosis_stage_col.append(fibrosis_stage)
        bridging_fibrosis_col.append(bridging_fibrosis)
        sinusoidal_fibrosis_col.append(sinusoidal_fibrosis)
        portal_fibrosis_col.append(portal_fibrosis)
        periportal_fibrosis_col.append(periportal_fibrosis)
        pericellular_fibrosis_col.append(pericellular_fibrosis)
        perivenular_fibrosis_col.append(perivenular_fibrosis)
        septal_fibrosis_col.append(septal_fibrosis)
        central_fibrosis_col.append(central_fibrosis)
        hepatitis_col.append(hepatitis)
        autoimmune_hepatitis_col.append(autoimmune_hepatitis)
        
        disease_list_col.append(disease_list)
        
        
    df_path['steatosis'] = steatosis_col
    df_path['ballooning'] = ballooning_col
    df_path['inflammation'] = inflammation_col
    df_path['lobular_inflammation'] = lobular_inflammation_col
#     df_path['portal_inflammation'] = portal_inflammation_col
    df_path['zone3_inflammation'] = zone3_inflammation_col
    df_path['lobular_hepatitis'] = lobular_hepatitis_col
    df_path['zone3_hepatitis'] = zone3_hepatitis_col
    df_path['fibrosis'] = fibrosis_col
    
    df_path['fibrosis_stage'] = fibrosis_stage_col
    df_path['bridging_fibrosis'] = bridging_fibrosis_col
    df_path['sinusoidal_fibrosis'] = sinusoidal_fibrosis_col
    df_path['portal_fibrosis'] = portal_fibrosis_col
    df_path['periportal_fibrosis'] = periportal_fibrosis_col
    df_path['pericellular_fibrosis'] = pericellular_fibrosis_col
    df_path['perivenular_fibrosis'] = perivenular_fibrosis_col
    df_path['septal_fibrosis'] = septal_fibrosis_col
    df_path['central_fibrosis'] = central_fibrosis_col
    df_path['hepatitis'] = hepatitis_col
    df_path['autoimmune_hepatitis'] = autoimmune_hepatitis_col
    
    df_path['disease_list'] = disease_list_col
   
    if update:
        # re-merge with original data
        print('Updating input path dataframe')
        pathdf['steatosis'] = np.nan
        pathdf['ballooning'] = np.nan
        pathdf['inflammation'] = np.nan
        pathdf['lobular_inflammation'] = np.nan
#         pathdf['portal_inflammation'] = np.nan
        pathdf['zone3_inflammation'] = np.nan
        pathdf['lobular_hepatitis'] = np.nan
        pathdf['zone3_hepatitis'] = np.nan
        pathdf['fibrosis'] = np.nan
        
        pathdf['fibrosis_stage'] = np.nan
        pathdf['bridging_fibrosis'] = np.nan
        pathdf['sinusoidal_fibrosis'] = np.nan
        pathdf['portal_fibrosis'] = np.nan
        pathdf['periportal_fibrosis'] = np.nan
        pathdf['pericellular_fibrosis'] = np.nan
        pathdf['perivenular_fibrosis'] = np.nan
        pathdf['septal_fibrosis'] = np.nan
        pathdf['central_fibrosis'] = np.nan
        pathdf['hepatitis'] = np.nan
        pathdf['autoimmune_hepatitis'] = np.nan
        
        pathdf['disease_list'] = np.nan
        pathdf.update(df_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
        #print('Returning MGH, BWH only entries with truncated path reports')
        return_df = df_path
        

    return return_df


def entity_recognition_liver(text, nlp):
    
    import re
    
    text = text.lower().replace(' bridging.', ' active-bridging.') #.replace(';',' ')
    
    entity_result = ''
    
    for line in text.split('.'):
        
        line = " ".join(line.split())
        line = line.strip()
        line = (line
                .replace(' as well as', ', ')
                .replace('neither', 'no').replace('nor', 'no')
                .replace('very', '')
                .replace('inflammatory', 'inflammation')
                .replace('inflamed', 'inflammation')
                .replace('severely', 'severe')
                .replace('moderately', 'moderate')
                .replace('moderate ', 'moderate-')
                .replace('mild to moderate', 'mild&moderate')
                .replace('moderate to severe', 'moderate&severe')
                .replace('mild to severe', 'mild&severe')
                .replace('mild ', 'mild-').replace('mild', 'mild-')
                .replace('severe ', 'severe-')
                .replace('minimal ', 'minimal-')
                .replace('chronic ', 'chronic-')
                .replace('focal ', 'focal-')
                .replace(' areas', ' area')
                .replace(' area', '-area')
                .replace(' tracts', 'tract')
                .replace('-tracts', 'tract')
                .replace('portal tract', 'portaltract')
                .replace('centrilobular', 'centri-lobular')
                .replace('periportal', 'peri-portal')
                .replace('portal and lobular', 'portal&lobular')
                .replace('portal or lobular', 'portal&lobular')
                .replace('lobular and portal', 'lobular&portal')
                .replace('lobular or portal', 'lobular&portal')
                .replace('portal&lobular', 'lobular&portal')
                .replace('lobular inflammat', 'lobular-inflammat')
                .replace('mixed ', 'mixed-')
                .replace('kupffer cell', 'kupffer-cell')
                .replace('zone 3', 'zone-3')
                .replace('hepatic plate', 'hepatic-plate')
                .replace('hepatic ', 'hepatitis ')
                .replace('steatotic', 'steatosis')
                .replace('microvesicular ', 'microvesicular-')
                .replace('macrovesicular ', 'macrovesicular-')
                .replace('< 5%', '<5%')
                .replace('non classical', 'nonclassical')
                .replace('non-classical', 'nonclassical')
                .replace('ballooning', 'baloning')
                .replace('portal to portal', 'portal-portal')
                .replace('central to central', 'central-central')
                .replace('portal to central', 'portal-central')
                .replace('bridging necrosis', 'bridging-necrosis')
                .replace('fibrosis bridging', 'fibrosis-bridging')
                .replace('fibrous bridging', 'fibrosis-bridging')
                .replace('bridging fibrosis', 'bridging-fibrosis')
                .replace('sinusoidal fibrosis', 'sinusoidal-fibrosis')
                .replace('portal fibrosis', 'portal-fibrosis')
                .replace('pericellular fibrosis', 'pericellular-fibrosis')
                .replace('perivenular fibrosis', 'perivenular-fibrosis')
                .replace('septal fibrosis', 'septal-fibrosis')
                .replace('central fibrosis', 'central-fibrosis')
                .replace('ductal fibrosis', 'ductal-fibrosis')
                .replace('portal bridging', 'portal-bridging')
                .replace('central bridging', 'central-bridging')
                .replace(' bridging ', ' active-bridging ')
                .replace(' bridging;', ' active-bridging;')
                .replace(' bridging,', ' active-bridging,')
                .replace('(p-p)','')
                #.replace('chronic ', 'chronic-')
                #.replace('active ', 'active-')
               )
        
        if 'fibrosis' in line:
            line = (line
                    .replace('pericellular ', 'pericellular-fibrosis ').replace('pericellular,', 'pericellular-fibrosis ')
                    .replace('sinusoidal ', 'sinusoidal-fibrosis ').replace('sinusoidal,', 'sinusoidal-fibrosis ')
#                     .replace(' perisinusoidal ', ' perisinusoidal-fibrosis ').replace(' perisinusoidal,', ' perisinusoidal-fibrosis ')
                    .replace('portal ', 'portal-fibrosis ').replace('portal,', 'portal-fibrosis ')
                    .replace('central ', 'central-fibrosis ').replace('central,', 'central-fibrosis ')
                    .replace('septal ', 'septal-fibrosis ').replace('septal,', 'septal-fibrosis ')
                    .replace('ductal ', 'ductal-fibrosis ').replace('ductal,', 'ductal-fibrosis ')
                    .replace('perivenular ', 'perivenular-fibrosis ').replace('perivenular,', 'perivenular-fibrosis ')
                   )
        if 'bridging' in line:
            line = (line
                    .replace('portal-portal-fibrosis', 'portal-portal-fibrosis-bridging')
                    .replace('portal-portal ', 'portal-portal-bridging')
                    .replace('portal-central-fibrosis', 'portal-central-fibrosis-bridging')
                    .replace('portal-central ', 'portal-central-bridging')
                    .replace('central-central-fibrosis', 'central-central-fibrosis-bridging')
                    .replace('central-central ', 'central-central-bridging')
                   )
            
        
        #global doc, e
   
        doc = nlp(line)
    
#         print(line)
    
        for e in doc.ents:
            
            e_text = e.text
            e_text = re.sub(' +', ' ', e_text)
            e_bool = e._.negex
            
            # Replace negation words in the entity and adjust sentiment
            if e_text.startswith(('no ', 'non-', 'non ')):
                to_match = ['^no ', '^non-', '^non ']
                e_text = re.sub('|'.join(to_match), '', e_text)
                e_bool = not e_bool
            
#             if 'lobular' in e_text or 'portal' in e_text:
#                 e_text = re.sub('(mild|moderate|mixed)', '', e_text).strip()

            inf_count = line.count('inflammation')
            
            lobu_inf_1 = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,3}?inflammation)\b', line))
            lobu_inf_2 = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,3}?necroinflammatory)\b', line))
            lobu_inf = lobu_inf_1 or lobu_inf_2
            
            zone3_inf = bool(re.search(r'\b(?:zone-3\W+(?:\w+\W+){0,2}?inflammation|inflammation\W+(?:\w+\W+){0,1}?zone-3)\b', line))
            
            # port_inf = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,3}?inflammation|inflammation\W+(?:\w+\W+){0,1}?portal)\b', line))
            port_inf = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,3}?inflammation)\b', line))
            
            port_tract_inf = bool(re.search(r'\b(?:portaltract\W+(?:\w+\W+){0,5}?inflammation)\b', line))
            
            lob_port_inf = bool(re.search(r'\b(?:lobular&portal\W+(?:\w+\W+){0,2}?inflammation)\b', line))

            lobu_act = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,3}?activity)\b', line))
            port_act = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,3}?activity)\b', line))

            lobu_hep = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,2}?hepatitis)\b', line))
            
            zone3_hep = bool(re.search(r'\b(?:zone-3\W+(?:\w+\W+){0,3}?hepatitis)\b', line))
            
            lobu_dis = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,3}?disarray)\b', line)) 
            
#             stg_fib = bool(re.search(r'\b(?:fibrosis\W+(?:\w+\W+){0,8}?stage|stage\W+(?:\w+\W+){0,6}?fibrosis)\b', line))
#             stg_bri = bool(re.search(r'\b(?:bridging\W+(?:\w+\W+){0,8}?stage|stage\W+(?:\w+\W+){0,6}?bridging)\b', line))
#             stg_cir = bool(re.search(r'\b(?:cirrhosis\W+(?:\w+\W+){0,8}?stage|stage\W+(?:\w+\W+){0,6}?cirrhosis)\b', line))
            
            sin_fib = bool(re.search(r'\b(?:sinusoidal\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?sinusoidal)\b', line))
            perisin_fib = bool(re.search(r'\b(?:perisinusoidal\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?perisinusoidal)\b', line))
            periport_fib = bool(re.search(r'\b(?:peri-portal\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?peri-portal)\b', line))
            port_fib = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?portal)\b', line))
            bridg_fib = bool(re.search(r'\b(?:bridging\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?bridging)\b', line))
            cent_fib = bool(re.search(r'\b(?:central\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?central)\b', line))
            sept_fib = bool(re.search(r'\b(?:septal\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?septal)\b', line))
            periven_fib = bool(re.search(r'\b(?:perivenular\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?perivenular)\b', line))
            pericel_fib = bool(re.search(r'\b(?:pericellular\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?pericellular)\b', line))

            
            if 'steatosis' in line and ('<5%' in line or 'less than 5%' in line) and 'steatosis' in e_text: #('5-33%' not in line)
                e_text = '<5% ' + e_text
                
            if 'lobular&portal' in line and 'inflammation' in e_text and not 'lobular&portal' in e_text and lob_port_inf:
                e_text = 'lobular&portal ' + e_text
            
            if 'lobular' in line and 'inflammation' in e_text and not 'lobular' in e_text and lobu_inf and inf_count<=1:
                e_text = 'lobular ' + e_text
            
            if 'lobular' in line and 'inflammation' in e_text and not 'lobular' in e_text and not 'portal' in e_text and not 'lobular-inflammation' in line and lobu_inf and inf_count>1:
                e_text = 'lobular ' + e_text
            
            if 'portal' in line and 'inflammation' in e_text and not 'portal' in e_text and port_inf and inf_count<=1:
                e_text = 'portal ' + e_text
                
            if 'portal' in line and 'inflammation' in e_text and not 'portal' in e_text and not 'lobular' in e_text and not 'portal inflammation' in line and port_inf and inf_count>1:
                e_text = 'portal ' + e_text
            
            if ('portaltract' in line or 'portal-area' in line) and 'inflammation' in e_text and not 'portal' in e_text  and (port_inf or port_tract_inf):
                e_text = 'portal ' + e_text
                
            if 'lobular' in line and 'activity' in e_text and not 'lobular' in e_text and lobu_act:
                e_text = 'lobular ' + e_text
            if 'portal' in line and 'activity' in e_text and not 'portal' in e_text and port_act:
                e_text = 'portal ' + e_text
            
            if 'lobular' in line and 'disarray' in e_text and not 'lobular' in e_text and lobu_dis:
                e_text = 'lobular ' + e_text
            
            if 'lobular&portal' in line and 'hepatitis' in e_text and not 'lobular&portal' in e_text and lobu_hep:
                e_text = 'lobular&portal ' + e_text
            
            if 'lobular' in line and ('hepatitis' in e_text and not 'lobular' in e_text) and lobu_hep:
                e_text = 'lobular ' + e_text
                
            if 'zone-3' in line and not 'zone-3 injury' in line and 'hepatitis' in e_text and not 'zone-3' in e_text and zone3_hep:
                e_text = 'zone-3 ' + e_text
                
            if 'zone-3' in line and not 'zone-3 injury' in line and 'inflammation' in e_text and not 'zone-3' in e_text and zone3_inf:
                e_text = 'zone-3 ' + e_text
            
            
#             if 'fibrosis' in e_text:
            
#                 if 'sinusoidal' in line and 'sinusoidal' not in e_text and sin_fib:
#                     e_text = e_text + ' sinusoidal'
#                 if 'perisinusoidal' in line and 'perisinusoidal' not in e_text and perisin_fib:
#                     e_text = e_text + ' perisinusoidal'
#                 if 'peri-portal' in line and 'peri-portal' not in e_text and periport_fib:
#                     e_text = e_text + ' peri-portal'
#                 if 'portal' in line and 'portal' not in e_text and port_fib:
#                     e_text = e_text + ' portal'
#                 if 'bridging' in line and 'bridging' not in e_text and bridg_fib:
#                     e_text = e_text + ' bridging'
#                 if 'central' in line and 'central' not in e_text and cent_fib:
#                     e_text = e_text + ' central'
#                 if 'septal' in line and 'septal' not in e_text and sept_fib:
#                     e_text = e_text + ' septal'
#                 if 'perivenular' in line and 'perivenular' not in e_text and periven_fib:
#                     e_text = e_text + ' perivenular'
#                 if 'pericellular' in line and 'pericellular' not in e_text and pericel_fib:
#                     e_text = e_text + ' pericellular'
            
           
            if ('fibrosis' in e_text or 'bridging' in e_text or 'cirrhosis' in e_text) and 'stage' in line:
                
                line = (line
                        .replace('stage iii', 'stage 3').replace('stage ii', 'stage 2').replace('stage iv', 'stage 4')
                        .replace('stage vi', 'stage 6').replace('stage v', 'stage 5').replace('stage i', 'stage 1')
                        
                        .replace('0 to 1', '0-1').replace('1 to 2', '1-2').replace('2 to 3', '2-3')
                        .replace('3 to 4', '3-4').replace('4 to 5', '4-5').replace('5 to 6', '5-6')
                        .replace('0-1', '0.5').replace('1-2', '1.5').replace('2-3', '2.5')
                        .replace('3-4', '3.5').replace('4-5', '4.5').replace('5-6', '5.5')
                        .replace('1b-2', '1.5')
                        .replace(' 1a ', ' 1 ').replace(' 2a ', ' 2 ').replace(' 3a ', ' 3 ')
                        .replace(' 4a ', ' 4 ').replace(' 5a ', ' 5 ').replace(' 6a ', ' 6 ')
                   )
        
                stagelist_1 = re.findall(r'stage.*?(\d+(?:\.\d+)?).*?(\d+(?:\.\d+)?)', line)
                stagelist_2 = re.findall(r'stage.*?(\d+(?:\.\d+)?)', line)
                
                dist_fib = abs(line.find('fibrosis')-line.find('stage'))
                dist_pbc = abs(line.find('pbc')-line.find('stage'))
                
                
                if not (len(stagelist_2)==0 or ('pbc' in line and dist_pbc<dist_fib)):
                    
                    stage_val = float(stagelist_2[0])
                    ref_val = 4.0

                    if len(stagelist_1)==1:

                        bool_1 = '.' in stagelist_1[0][0] and '.' in stagelist_2[0]
                        bool_2 = '.' not in stagelist_2[0]

                        if bool_1 or bool_2:
                            stage_val = float(stagelist_1[0][0])
                            ref_val = float(stagelist_1[0][1])

                    if stage_val>=5:
                        ref_val = 6.0

                    if 'ishak' in line:
                        e_text = e_text + ' ishak'

                    e_text = e_text + ' stage: ' + str(stage_val) + '/' + str(ref_val)

                    if stage_val==0:
                        e_bool = True
                    else:
                        e_bool = False

            e_text = " ".join(e_text.split())
            
            
            entity_result = entity_result + e_text + ' ' + str(not e_bool) + '\n'
            
            entity_result = entity_result.replace('baloon', 'balloon').replace('ballon', 'balloon').replace('balon', 'balloon')
    
        
    return entity_result
