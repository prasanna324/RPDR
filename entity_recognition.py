
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
            "pseudo_negations": ts.terms['pseudo_negations'] + ['not limited to', 'not excluded', 'needs to be ruled out', 'although not apparent'],
            "preceding_negations": ts.terms['preceding_negations'] + ['negative', 'insufficient', 'without evidence of', 'rather than', 'history'],
            "following_negations": ts.terms['following_negations'] + ['negative', 'unremarkable', 'ruled out', 'less likely', 'is not', 'are not', 'does not', 'have not', 'was not', 'were not', 'absent', 'not present'],
            "termination": ts.terms['termination'] + ['note:', ';', ', negative', ',negative']
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
    
    erythema_col = []
    marked_erythema_col = []
    inflammation_col = []
    mild_inflammation_col = []
    moderate_inflammation_col = []
    severe_inflammation_col = []
    loss_vasculature_col = []
    dec_vasculature_col = []
    granularity_col = []
    ulceration_col = []
    friability_col = []
    mild_friability_col = []
    spont_bleeding_col = []
    adherent_blood_col = []
    erosion_col = []
    congestion_col = []
    edema_col = []
    pseudopolyp_col = []
    crohns_col = []
    
    superficial_ulcer_col = []
    shallow_ulcer_col = []
    aphthous_ulcer_col = []
    small_ulcer_col = []
    large_ulcer_col = []
    deep_ulcer_col = []
    mild_ulcer_col = []
    
    colitis_col = []
    chronic_colitis_col = []
    mild_colitis_col = []
    moderate_colitis_col = []
    severe_colitis_col = []
    inactive_colitis_col = []
    active_colitis_col= []
    acute_colitis_col = []
    ulcerative_colitis_col = []
    pan_colitis_col = []
    proctitis_col = []
    proctosigmoiditis_col = []
    left_sided_colitis_col = []
    
    active_ileitis_col = []
    chronic_ileitis_col = []
    arch_distortion_col = []
    basal_plasmacytosis_col = []
    active_enteritis_col = []
    chronic_enteritis_col = []
    crypt_abscess_col = []
    crypt_atrophy_col = []
    cryptitis_col = []
    lymphoid_agg_col = []
    lamina_propria_col = []
    granuloma_col = []
    noncaseating_gran_col = []
    nonnecrotizing_gran_col = []
    paneth_cell_col = []
    cdiff_col = []
    cmv_col = []
    
    mayo_score_col = []
    disease_list_col = []

    for i in range(0,num_reports):

        # extract path report for this entry
        disease_list = []
        report_text = df_path.iloc[i,:].Report_Text
        result_text = entity_recognition_colon(report_text, nlp=nlp_2)
        
        erythema = np.nan
        marked_erythema = np.nan
        inflammation = np.nan
        mild_inflammation = np.nan
        moderate_inflammation = np.nan
        severe_inflammation = np.nan
        loss_vasculature = np.nan
        dec_vasculature = np.nan
        granularity = np.nan
        ulceration = np.nan
        friability = np.nan
        mild_friability = np.nan
        spont_bleeding = np.nan
        adherent_blood = np.nan
        erosion = np.nan
        congestion = np.nan
        edema = np.nan
        pseudopolyp = np.nan
        crohns = np.nan
        
        superficial_ulcer = np.nan
        shallow_ulcer = np.nan
        aphthous_ulcer = np.nan
        small_ulcer = np.nan
        large_ulcer = np.nan
        deep_ulcer = np.nan
        mild_ulcer = np.nan
        
        colitis = np.nan
        chronic_colitis = np.nan
        mild_colitis = np.nan
        moderate_colitis = np.nan
        severe_colitis = np.nan
        inactive_colitis = np.nan
        active_colitis = np.nan
        acute_colitis = np.nan
        ulcerative_colitis = np.nan
        pan_colitis = np.nan
        proctitis = np.nan
        proctosigmoiditis = np.nan
        left_sided_colitis = np.nan
        
        active_ileitis = np.nan
        chronic_ileitis = np.nan
        arch_distortion = np.nan
        basal_plasmacytosis = np.nan
        active_enteritis = np.nan
        chronic_enteritis = np.nan
        crypt_abscess = np.nan
        crypt_atrophy = np.nan
        cryptitis = np.nan
        lymphoid_agg = np.nan
        lamina_propria = np.nan
        granuloma = np.nan
        noncaseating_gran = np.nan
        nonnecrotizing_gran = np.nan
        paneth_cell = np.nan
        cdiff = np.nan
        cmv = np.nan
        
        mayo_score = np.nan

        
        for x in result_text.split('\n'):
            
            is_disease = False
            
            if ('erythem' in x and not 'marked' in x) and (np.isnan(erythema) or 'True' in x):
                if 'True' in x: erythema, is_disease = True, True
                else: erythema = False
            if ('marked' in x and 'erythema' in x) and (np.isnan(marked_erythema) or 'True' in x):
                if 'True' in x: marked_erythema, is_disease = True, True
                else: marked_erythema = False
            if ('inflammat' in x or 'inflamed' in x) and (np.isnan(inflammation) or 'True' in x):
                if 'True' in x: inflammation, is_disease = True, True
                else: inflammation = False
            if ('mild' in x and 'inflammation' in x) and (np.isnan(mild_inflammation) or 'True' in x):
                if 'True' in x: mild_inflammation, is_disease = True, True
                else: mild_inflammation = False
            if ('moderate' in x and 'inflammation' in x) and (np.isnan(moderate_inflammation) or 'True' in x):
                if 'True' in x: moderate_inflammation, is_disease = True, True
                else: moderate_inflammation = False
            if ('severe' in x and 'inflammation' in x) and (np.isnan(severe_inflammation) or 'True' in x):
                if 'True' in x: severe_inflammation, is_disease = True, True
                else: severe_inflammation = False
            if ('loss-of-vasculature' in x or 'absent-vascular' in x) and (np.isnan(loss_vasculature) or 'True' in x):
                if 'True' in x: loss_vasculature, is_disease = True, True
                else: loss_vasculature = False
            if ('decreased-vasculature' in x or 'altered-vascular' in x) and (np.isnan(dec_vasculature) or 'True' in x):
                if 'True' in x: dec_vasculature, is_disease = True, True
                else: dec_vasculature = False
            if ('granular' in x or 'granulation' in x) and (np.isnan(granularity) or 'True' in x):
                if 'True' in x: granularity, is_disease = True, True
                else: granularity = False
            if ('ulcer' in x and not 'ulcerativecolitis' in x) and (np.isnan(ulceration) or 'True' in x):
                if 'True' in x: ulceration, is_disease = True, True
                else: ulceration = False
            if (('friable' in x or 'friability' in x) and not 'mild' in x) and (np.isnan(friability) or 'True' in x):
                if 'True' in x: friability, is_disease = True, True
                else: friability = False
            if (('friable' in x or 'friability' in x) and 'mild' in x) and (np.isnan(mild_friability) or 'True' in x):
                if 'True' in x: mild_friability, is_disease = True, True
                else: mild_friability = False
            if ('spontaneous' in x and ('bleed' in x or 'bled' in x)) and (np.isnan(spont_bleeding) or 'True' in x):
                if 'True' in x: spont_bleeding, is_disease = True, True
                else: spont_bleeding = False
            if 'adherent' in x and ('blood' in x or 'clot' in x) and (np.isnan(adherent_blood) or 'True' in x):
                if 'True' in x: adherent_blood, is_disease = True, True
                else: adherent_blood = False
            if 'erosion' in x and (np.isnan(erosion) or 'True' in x):
                if 'True' in x: erosion, is_disease = True, True
                else: erosion = False
            if 'congest' in x and (np.isnan(congestion) or 'True' in x):
                if 'True' in x: congestion, is_disease = True, True
                else: congestion = False
            if 'edema' in x and (np.isnan(edema) or 'True' in x):
                if 'True' in x: edema, is_disease = True, True
                else: edema = False
            if ('pseudopolyp' in x or 'pseudo-polyp' in x) and (np.isnan(pseudopolyp) or 'True' in x):
                if 'True' in x: pseudopolyp, is_disease = True, True
                else: pseudopolyp = False
            if 'crohn' in x and (np.isnan(crohns) or 'True' in x):
                if 'True' in x: crohns, is_disease = True, True
                else: crohns = False
                    
            if (('superficial' in x or 'superfacial' in x) and 'ulcer' in x) and (np.isnan(superficial_ulcer) or 'True' in x):
                if 'True' in x: superficial_ulcer, is_disease = True, True
                else: superficial_ulcer = False
            if ('shallow' in x and 'ulcer' in x) and (np.isnan(shallow_ulcer) or 'True' in x):
                if 'True' in x: shallow_ulcer, is_disease = True, True
                else: shallow_ulcer = False
            if ('aphthous' in x and 'ulcer' in x) and (np.isnan(aphthous_ulcer) or 'True' in x):
                if 'True' in x: aphthous_ulcer, is_disease = True, True
                else: aphthous_ulcer = False
            if ('small' in x and 'ulcer' in x) and (np.isnan(small_ulcer) or 'True' in x):
                if 'True' in x: small_ulcer, is_disease = True, True
                else: small_ulcer = False
            if ('large' in x and 'ulcer' in x) and (np.isnan(large_ulcer) or 'True' in x):
                if 'True' in x: large_ulcer, is_disease = True, True
                else: large_ulcer = False
            if ('deep' in x and 'ulcer' in x) and (np.isnan(deep_ulcer) or 'True' in x):
                if 'True' in x: deep_ulcer, is_disease = True, True
                else: deep_ulcer = False
            if ('mild' in x and (bool(re.search(r'\bulcer\b|\bulceration\b', x)))) and (np.isnan(mild_ulcer) or 'True' in x):
                if 'True' in x: mild_ulcer, is_disease = True, True
                else: mild_ulcer = False
                    
            
            if 'colitis' in x and (np.isnan(colitis) or 'True' in x):
                if 'True' in x: colitis, is_disease = True, True
                else: colitis = False
            if ('chronic' in x and 'colitis' in x) and (np.isnan(chronic_colitis) or 'True' in x):
                if 'True' in x: chronic_colitis, is_disease = True, True
                else: chronic_colitis = False
            if ('mild' in x and 'colitis' in x) and (np.isnan(mild_colitis) or 'True' in x):
                if 'True' in x: mild_colitis, is_disease = True, True
                else: mild_colitis = False
            if ('moderate' in x and 'colitis' in x) and (np.isnan(moderate_colitis) or 'True' in x):
                if 'True' in x: moderate_colitis, is_disease = True, True
                else: moderate_colitis = False
            if ('severe' in x and 'colitis' in x) and (np.isnan(severe_colitis) or 'True' in x):
                if 'True' in x: severe_colitis, is_disease = True, True
                else: severe_colitis = False
            if ('active' in x and 'colitis' in x and not 'inactive' in x) and (np.isnan(active_colitis) or 'True' in x):
                if 'True' in x: active_colitis, is_disease = True, True
                else: active_colitis = False
            if ('inactive' in x and 'colitis' in x) and (np.isnan(inactive_colitis) or 'True' in x):
                if 'True' in x: inactive_colitis, is_disease = True, True
                else: inactive_colitis = False
            if ('acute' in x and 'colitis' in x) and (np.isnan(acute_colitis) or 'True' in x):
                if 'True' in x: acute_colitis, is_disease = True, True
                else: acute_colitis = False
            if ('ulcerative' in x and 'colitis' in x) and (np.isnan(ulcerative_colitis) or 'True' in x):
                if 'True' in x: ulcerative_colitis, is_disease = True, True
                else: ulcerative_colitis = False
            if ('pancolitis' in x or 'pan colitis' in x or 'pan-colitis' in x) and (np.isnan(pan_colitis) or 'True' in x):
                if 'True' in x: pan_colitis, is_disease = True, True
                else: pan_colitis = False
            if 'proctitis' in x and (np.isnan(proctitis) or 'True' in x):
                if 'True' in x: proctitis, is_disease = True, True
                else: proctitis = False
            if 'proctosigmoiditis' in x and (np.isnan(proctosigmoiditis) or 'True' in x):
                if 'True' in x: proctosigmoiditis, is_disease = True, True
                else: proctosigmoiditis = False
            if 'left-sided colitis' in x and (np.isnan(left_sided_colitis) or 'True' in x):
                if 'True' in x: left_sided_colitis, is_disease = True, True
                else: left_sided_colitis = False
                    
            if ('active' in x and 'ileitis' in x and not 'inactive' in x) and (np.isnan(active_ileitis) or 'True' in x):
                if 'True' in x: active_ileitis, is_disease = True, True
                else: active_ileitis = False
            if ('chronic' in x and 'ileitis' in x) and (np.isnan(chronic_ileitis) or 'True' in x):
                if 'True' in x: chronic_ileitis, is_disease = True, True
                else: chronic_ileitis = False
            if ('architectural' in x and ('distortion' in x or 'disarray' in x or 'disorder' in x)) and (np.isnan(arch_distortion) or 'True' in x):
                if 'True' in x: arch_distortion, is_disease = True, True
                else: arch_distortion = False
            if ('basal' in x and 'plasmacytosis' in x) and (np.isnan(basal_plasmacytosis) or 'True' in x):
                if 'True' in x: basal_plasmacytosis, is_disease = True, True
                else: basal_plasmacytosis = False
            if ('active' in x and 'enteritis' in x and not 'inactive' in x) and (np.isnan(active_enteritis) or 'True' in x):
                if 'True' in x: active_enteritis, is_disease = True, True
                else: active_enteritis = False
            if ('chronic' in x and 'enteritis' in x) and (np.isnan(chronic_enteritis) or 'True' in x):
                if 'True' in x: chronic_enteritis, is_disease = True, True
                else: chronic_enteritis = False
            if ('crypt' in x and 'abscess' in x) and (np.isnan(crypt_abscess) or 'True' in x):
                if 'True' in x: crypt_abscess, is_disease = True, True
                else: crypt_abscess = False
            if ('crypt' in x and 'atrophy' in x) and (np.isnan(crypt_atrophy) or 'True' in x):
                if 'True' in x: crypt_atrophy, is_disease = True, True
                else: crypt_atrophy = False
            if ('cryptitis' in x) and (np.isnan(cryptitis) or 'True' in x):
                if 'True' in x: cryptitis, is_disease = True, True
                else: cryptitis = False
            if ('lymphoid' in x and 'aggregate' in x) and (np.isnan(lymphoid_agg) or 'True' in x):
                if 'True' in x: lymphoid_agg, is_disease = True, True
                else: lymphoid_agg = False
            if ('lamina-propria' in x and (('increase' in x and 'cellularity' in x) or ('inflammation' in x))) and (np.isnan(lamina_propria) or 'True' in x):
                if 'True' in x: lamina_propria, is_disease = True, True
                else: lamina_propria = False
            
            if ('granuloma' in x) and (np.isnan(granuloma) or 'True' in x):
                if 'True' in x: granuloma, is_disease = True, True
                else: granuloma = False
            if ('noncaseating' in x and 'granuloma' in x) and (np.isnan(noncaseating_gran) or 'True' in x):
                if 'True' in x: noncaseating_gran, is_disease = True, True
                else: noncaseating_gran = False
            if ('nonnecrotizing' in x and 'granuloma' in x) and (np.isnan(nonnecrotizing_gran) or 'True' in x):
                if 'True' in x: nonnecrotizing_gran, is_disease = True, True
                else: nonnecrotizing_gran = False
            if 'paneth-cell' in x and (np.isnan(paneth_cell) or 'True' in x):
                if 'True' in x: paneth_cell, is_disease = True, True
                else: paneth_cell = False
            if ('c-diff' in x or 'clostridium' in x) and (np.isnan(cdiff) or 'True' in x):
                if 'True' in x: cdiff, is_disease = True, True
                else: cdiff = False
            if ('cmv' in x) and (np.isnan(cmv) or 'True' in x):
                if 'True' in x: cmv, is_disease = True, True
                else: cmv = False     
            
            if 'colitis mayo-' in x:
                mayo_score = float(re.findall(r'.*?(\d+(?:,\d+)*(?:\.\d+)?)', x)[0])
                if mayo_score>3:
                    mayo_score = 3
            if is_disease:
                disease_list.append(x)
        
        erythema_col.append(erythema)
        marked_erythema_col.append(marked_erythema)
        inflammation_col.append(inflammation)
        mild_inflammation_col.append(mild_inflammation)
        moderate_inflammation_col.append(moderate_inflammation)
        severe_inflammation_col.append(severe_inflammation)4
        loss_vasculature_col.append(loss_vasculature)
        dec_vasculature_col.append(dec_vasculature)
        granularity_col.append(granularity)
        ulceration_col.append(ulceration)
        friability_col.append(friability)
        mild_friability_col.append(mild_friability)
        spont_bleeding_col.append(spont_bleeding)
        adherent_blood_col.append(adherent_blood)
        erosion_col.append(erosion)
        congestion_col.append(congestion)
        edema_col.append(edema)
        pseudopolyp_col.append(pseudopolyp)
        crohns_col.append(crohns)
        
        superficial_ulcer_col.append(superficial_ulcer)
        shallow_ulcer_col.append(shallow_ulcer)
        aphthous_ulcer_col.append(aphthous_ulcer)
        small_ulcer_col.append(small_ulcer)
        large_ulcer_col.append(large_ulcer)
        deep_ulcer_col.append(deep_ulcer)
        mild_ulcer_col.append(mild_ulcer)
    
        colitis_col.append(colitis)
        chronic_colitis_col.append(chronic_colitis)
        mild_colitis_col.append(mild_colitis)
        moderate_colitis_col.append(moderate_colitis)
        severe_colitis_col.append(severe_colitis)
        active_colitis_col.append(active_colitis)
        inactive_colitis_col.append(inactive_colitis)
        acute_colitis_col.append(acute_colitis)
        ulcerative_colitis_col.append(ulcerative_colitis)
        pan_colitis_col.append(pan_colitis)
        proctitis_col.append(proctitis)
        left_sided_colitis_col.append(left_sided_colitis)
        proctosigmoiditis_col.append(proctosigmoiditis)
        
        active_ileitis_col.append(active_ileitis)
        chronic_ileitis_col.append(chronic_ileitis)
        arch_distortion_col.append(arch_distortion)
        basal_plasmacytosis_col.append(basal_plasmacytosis)
        active_enteritis_col.append(active_enteritis)
        chronic_enteritis_col.append(chronic_enteritis)
        crypt_abscess_col.append(crypt_abscess)
        crypt_atrophy_col.append(crypt_atrophy)
        cryptitis_col.append(cryptitis)
        lymphoid_agg_col.append(lymphoid_agg)
        lamina_propria_col.append(lamina_propria)
        granuloma_col.append(granuloma)
        noncaseating_gran_col.append(noncaseating_gran)
        nonnecrotizing_gran_col.append(nonnecrotizing_gran)
        paneth_cell_col.append(paneth_cell)
        cdiff_col.append(cdiff)
        cmv_col.append(cmv)
        
        mayo_score_col.append(mayo_score)
        disease_list_col.append(disease_list)
        
    df_path['erythema'] = erythema_col
    df_path['marked_erythema'] = marked_erythema_col
    df_path['inflammation'] = inflammation_col
    df_path['mild_inflammation'] = mild_inflammation_col
    df_path['moderate_inflammation'] = moderate_inflammation_col
    df_path['severe_inflammation'] = severe_inflammation_col
    df_path['loss_vasculature'] = loss_vasculature_col
    df_path['dec_vasculature'] = dec_vasculature_col
    df_path['granularity'] = granularity_col
    df_path['ulceration'] = ulceration_col
    df_path['friability'] = friability_col
    df_path['mild_friability'] = mild_friability_col
    df_path['spont_bleeding'] = spont_bleeding_col
    df_path['adherent_blood'] = adherent_blood_col
    df_path['erosion'] = erosion_col
    df_path['congestion'] = congestion_col
    df_path['edema'] = edema_col
    df_path['pseudopolyp'] = pseudopolyp_col
    df_path['crohns'] = crohns_col
    
    df_path['superficial_ulcer'] = superficial_ulcer_col
    df_path['shallow_ulcer'] = shallow_ulcer_col
    df_path['aphthous_ulcer'] = aphthous_ulcer_col
    df_path['small_ulcer'] = small_ulcer_col
    df_path['large_ulcer'] = large_ulcer_col
    df_path['deep_ulcer'] = deep_ulcer_col
    df_path['mild_ulcer'] = mild_ulcer_col
    
    df_path['colitis'] = colitis_col
    df_path['chronic_colitis'] = chronic_colitis_col
    df_path['mild_colitis'] = mild_colitis_col
    df_path['moderate_colitis'] = moderate_colitis_col
    df_path['severe_colitis'] = severe_colitis_col
    df_path['active_colitis'] = active_colitis_col
    df_path['inactive_colitis'] = inactive_colitis_col
    df_path['acute_colitis'] = acute_colitis_col
    df_path['ulcerative_colitis'] = ulcerative_colitis_col
    df_path['pan_colitis'] = pan_colitis_col
    df_path['proctitis'] = proctitis_col
    df_path['proctosigmoiditis'] = proctosigmoiditis_col
    df_path['left_sided_colitis'] = left_sided_colitis_col
    
    df_path['active_ileitis'] = active_ileitis_col
    df_path['chronic_ileitis'] = chronic_ileitis_col
    df_path['arch_distortion'] = arch_distortion_col
    df_path['basal_plasmacytosis'] = basal_plasmacytosis_col
    df_path['active_enteritis'] = active_enteritis_col
    df_path['chronic_enteritis'] = chronic_enteritis_col
    df_path['crypt_abscess'] = crypt_abscess_col
    df_path['crypt_atrophy'] = crypt_atrophy_col
    df_path['cryptitis'] = cryptitis_col
    df_path['lymphoid_agg'] = lymphoid_agg_col
    df_path['lamina_propria'] = lamina_propria_col
    df_path['granuloma'] = granuloma_col
    df_path['noncaseating_gran'] = noncaseating_gran_col
    df_path['nonnecrotizing_gran'] = nonnecrotizing_gran_col
    df_path['paneth_cell'] = paneth_cell_col
    df_path['cdiff'] = cdiff_col
    df_path['cmv'] = cmv_col
    
    df_path['mayo_score'] = mayo_score_col
    df_path['disease_list'] = disease_list_col
   
    if update:
        # re-merge with original data
        print('Updating input path dataframe')
        
        pathdf['erythema'] = np.nan
        pathdf['marked_erythema'] = np.nan
        pathdf['inflammation'] = np.nan
        pathdf['mild_inflammation'] = np.nan
        pathdf['moderate_inflammation'] = np.nan
        pathdf['severe_inflammation'] = np.nan
        pathdf['loss_vasculature'] = np.nan
        pathdf['dec_vasculature'] = np.nan
        pathdf['granularity'] = np.nan
        pathdf['ulceration'] = np.nan
        pathdf['friability'] = np.nan
        pathdf['mild_friability'] = np.nan
        pathdf['spont_bleeding'] = np.nan
        pathdf['adherent_blood'] = np.nan
        pathdf['erosion'] = np.nan
        pathdf['congestion'] = np.nan
        pathdf['edema'] = np.nan
        pathdf['pseudopolyp'] = np.nan
        pathdf['crohns'] = np.nan
        
        pathdf['superficial_ulcer'] = np.nan
        pathdf['shallow_ulcer'] = np.nan
        pathdf['aphthous_ulcer'] = np.nan
        pathdf['small_ulcer'] = np.nan
        pathdf['large_ulcer'] = np.nan
        pathdf['deep_ulcer'] = np.nan
        pathdf['mild_ulcer'] = np.nan
        
        pathdf['colitis'] = np.nan
        pathdf['chronic_colitis'] = np.nan
        pathdf['mild_colitis'] = np.nan
        pathdf['moderate_colitis'] = np.nan
        pathdf['severe_colitis'] = np.nan
        pathdf['active_colitis'] = np.nan
        pathdf['inactive_colitis'] = np.nan
        pathdf['acute_colitis'] = np.nan
        pathdf['ulcerative_colitis'] = np.nan
        pathdf['pan_colitis'] = np.nan
        pathdf['proctitis'] = np.nan
        pathdf['proctosigmoiditis'] = np.nan
        pathdf['left_sided_colitis'] = np.nan
        
        pathdf['active_ileitis'] = np.nan
        pathdf['chronic_ileitis'] = np.nan
        pathdf['arch_distortion'] = np.nan
        pathdf['basal_plasmacytosis'] = np.nan
        pathdf['active_enteritis'] = np.nan
        pathdf['chronic_enteritis'] = np.nan
        pathdf['crypt_abscess'] = np.nan
        pathdf['crypt_atrophy'] = np.nan
        pathdf['cryptitis'] = np.nan
        pathdf['lymphoid_agg'] = np.nan
        pathdf['lamina_propria'] = np.nan
        pathdf['granuloma'] = np.nan
        pathdf['noncaseating_gran'] = np.nan
        pathdf['nonnecrotizing_gran'] = np.nan
        pathdf['paneth_cell'] = np.nan
        pathdf['cdiff'] = np.nan
        pathdf['cmv'] = np.nan
        
        pathdf['mayo_score'] = np.nan
        pathdf['disease_list'] = np.nan
        pathdf.update(df_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
#         print('Returning MGH, BWH only entries with truncated path reports')
        return_df = df_path
        
    return return_df


def entity_recognition_colon(text, nlp):
    
    import re
    import numpy as np
    
    text = (text
            .replace(' III ', '3').replace(' II ', '2').replace(' IV ', '4')
            .replace(' UC ', ' UlcerativeColitis ').replace(' UC.', ' UlcerativeColitis.')
            .replace(' UC,', ' UlcerativeColitis,').replace(' UC)', ' UlcerativeColitis)')
            .replace('(UC ', '(UlcerativeColitis ').replace('/UC ', '/UlcerativeColitis ')
            .replace(' UC-', ' UlcerativeColitis').replace(' UC\n', ' UlcerativeColitis\n')
            .replace(' CUC ', ' Chronic-UlcerativeColitis ').replace(' CUC.', ' Chronic-UlcerativeColitis.')
            .replace(' CUC,', ' Chronic-UlcerativeColitis,').replace(' CUC)', ' Chronic-UlcerativeColitis)')
            .replace('(CUC ', '(Chronic-UlcerativeColitis ').replace('/CUC ', '/Chronic-UlcerativeColitis ')
            .replace(' CUC', ' Chronic-UlcerativeColitis').replace(' CUC\n', ' Chronic-UlcerativeColitis\n')
            .replace(' U.C.', ' UlcerativeColitis')
            .replace('c. diff', 'c-diff').replace('C. diff', 'C-diff').replace('C. Diff', 'C-Diff').replace('c. Diff', 'c-Diff')
            .replace(' CD ', ' CD (Crohns) ').replace(' CD.', ' CD (Crohns).')
           )
    text = re.sub(' UC$', ' UlcerativeColitis', text)
    
    text = text.lower()
    
    entity_result = ''
    mayo_score = -1
    mayo_bool = False
    
    for line in text.split('.'):
                
        #line = line.strip()
        line = " ".join(line.split())
        line = (line
                .replace('+/-', ',')
                .replace(' no ', ' , no ')
                .replace('is not present', 'not present').replace('no present', 'not present').replace('not present', ' is not present')
                .replace(' minimal ', ' ,minimal ')
                .replace('ulcerations', 'ulceration')
                .replace('ulcers', 'ulcer')
                .replace('apthous', 'aphthous')
                .replace('noted in the', ' in ')
                .replace('is noted in', ' in ')
                .replace('neither', 'no').replace('nor', 'no')
                .replace('may not', 'will not')
#                 .replace('chronic active', 'chronic-active')
#                 .replace('active chronic', 'active-chronic')
#                 .replace('chronic inactive', 'chronic-inactive')
                .replace('severely', 'severe')
                .replace('moderately', 'moderate')
                .replace('mildly', 'mild').replace('mildl', 'mild').replace('midly', 'mild')
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
                .replace('mild ', 'mild_').replace('(mild) ', 'mild_').replace('(mild)', 'mild_')
                .replace('moderate ', 'moderate_').replace('(moderate) ', 'moderate_').replace('(moderate)', 'moderate_')
                .replace('severe ', 'severe_').replace('(severe) ', 'severe_').replace('(severe)', 'severe_')
                .replace('inactive ', 'inactive_')
                .replace('ulcerative colitis', 'ulcerativecolitis')
                .replace('healed colitis', 'healed-colitis')
#                 .replace('surveillance', 'not present')
                .replace('mayo i ', 'mayo 1 ')
                .replace('grade i ', 'grade 1 ')
                .replace('absent vascular', 'absent-vascular')
                .replace('decreased vasculature', 'decreased-vasculature')
                .replace('altered vascular', 'altered-vascular')
                .replace('spontaneous hemorrhage', 'spontaneous hemorrhage (spontaneous-bleeding)')
                .replace('spontaneous bleed ', 'spontaneous bleed (spontaneous-bleeding)')
                .replace('spontaneous bleed,', 'spontaneous bleed (spontaneous-bleeding),')
                .replace('spontaneously bleeding', 'spontaneously-bleeding (spontaneous-bleeding)')
                .replace('bleeding spontaneously', 'bleeding spontaneously (spontaneous-bleeding)')
                .replace('oozing and bleeding', 'oozing&bleeding')
                .replace('bleed actively', 'bleed actively (spontaneous-bleeding)')
                .replace('spontaneously bled ', 'spontaneously bled (spontaneous-bleeding)')
                .replace('bleed spontaneously', 'bleed spontaneously (spontaneous-bleeding)')
                .replace('non ulcer', 'no ulcer').replace('non-ulcer', 'no ulcer').replace('nonulcer', 'no ulcer')
                .replace('adherent blood', 'adherent-blood').replace('adherent clot', 'adherent-clot')
                .replace('pseudo- polyp', 'pseudopolyp')
                .replace('procto sigmoid', 'proctosigmoid').replace('procto-sigmoid', 'proctosigmoid')
                .replace('inflammation', 'inflammation')
                .replace('/budd', ' budd')
                .replace('architectural irregularity', 'architectural irregularity (architectural-distortion)')
                .replace('architectural irregularities', 'architectural irregularities (architectural-distortion)')
                .replace('architectural changes', 'architectural changes (architectural-distortion)')
                .replace('architectural alterations', 'architectural alterations (architectural-distortion)')
                .replace('architectural alteration', 'architectural alteration (architectural-distortion)')
                .replace('crypt distortion', 'crypt distortion (architectural-distortion)')
                .replace('crypt disarrary', 'crypt disarrary (architectural-distortion)')
                .replace('non-caseating', 'noncaseating').replace('non caseating', 'noncaseating')
                .replace('noncaseating granuloma','noncaseating-granuloma')
                .replace('non-necrotizing', 'nonnecrotizing').replace('non necrotizing', 'nonnecrotizing')
                .replace('nonnecrotizing granuloma','nonnecrotizing-granuloma')
                .replace('crypt abscess', 'crypt-abscess')
                .replace('crypt atrophy', 'crypt-atrophy')
                .replace('lymphoid aggregate', 'lymphoid-aggregate')
                .replace('lamina propria', 'lamina-propria')
                .replace('inflammatory', 'inflammation')
                .replace('paneth cell', 'paneth-cell')
                .replace('c diff', 'c-diff')
                .replace('marked erythema', 'marked-erythema')
                .replace('superficial erosion', 'superficial_erosion')
                .replace('superficial ulcer', 'superficial_ulcer')
                .replace('shallow ulcer', 'shallow_ulcer')
                .replace('aphthous ulcer', 'aphthous_ulcer')
                .replace('aphthous lesion', 'aphthous_lesion')
                .replace('small bowel', 'small_bowel').replace('small intestine', 'small_intestine')
                .replace('small ulcer', 'small_ulcer')
                .replace('large ulcer', 'large_ulcer')
                .replace('deep ulcer', 'deep_ulcer')
                .replace('mild ulcer', 'mild_ulcer')
                
                
#                 .replace('mayo ', 'mayo-')
#                 .replace('chronic active', 'chronic-active')
#                 .replace('inactive chronic', 'inactive-chronic')
# # #                 .replace('active colitis', 'active-colitis')
#                 .replace('inactive colitis', 'inactive-colitis')
               )
                
        #global doc, e
            
        doc = nlp(line)
        
        loss_vasc_1 = bool(re.search(r'\b(?:loss\W+(?:\w+\W+){0,2}?vasculature)\b', line))
        loss_vasc_2 = bool(re.search(r'\b(?:loss\W+(?:\w+\W+){0,2}?vascular)\b', line))
        loss_vasc = (loss_vasc_1 or loss_vasc_2)
        
        dec_vasc_1 = bool(re.search(r'\b(?:decrease\W+(?:\w+\W+){0,2}?vasculature)\b', line))
        dec_vasc_2 = bool(re.search(r'\b(?:decrease\W+(?:\w+\W+){0,2}?vascular)\b', line))
        dec_vasc_3 = bool(re.search(r'\b(?:decreased\W+(?:\w+\W+){0,2}?vasculature)\b', line))
        dec_vasc_4 = bool(re.search(r'\b(?:decreased\W+(?:\w+\W+){0,2}?vascular)\b', line))
        dec_vasc = (dec_vasc_1 or dec_vasc_2 or dec_vasc_3 or dec_vasc_4)
        
        mild_infl = bool(re.search(r'\b(?:.*mild\W+(?:\w+\W+){0,2}?inflammation|inflammation\W+(?:\w+\W+){0,2}?mild)\b', line))
        mod_infl = bool(re.search(r'\b(?:.*moderate\W+(?:\w+\W+){0,2}?inflammation|inflammation\W+(?:\w+\W+){0,2}?moderate)\b', line))
        sev_infl = bool(re.search(r'\b(?:.*severe\W+(?:\w+\W+){0,2}?inflammation|inflammation\W+(?:\w+\W+){0,2}?severe)\b', line))
        
        adhr_blood_1 = bool(re.search(r'\b(?:adherent\W+(?:\w+\W+){0,1}?blood)\b', line))
        adhr_clot_1 = bool(re.search(r'\b(?:adherent\W+(?:\w+\W+){0,1}?clot)\b', line))
        adhr_clot_2 = bool(re.search(r'\b(?:adherent\W+(?:\w+\W+){0,1}?clots)\b', line))
        adhr_blood = (adhr_blood_1 or adhr_clot_1 or adhr_clot_2)
        
        left_colitis_1 = bool(re.search(r'\b(?:left\W+(?:\w+\W+){0,4}?.*colitis|.*colitis\W+(?:\w+\W+){0,3}?left)\b', line))
        left_colitis_2 = bool(re.search(r'\b(?:left-sided\W+(?:\w+\W+){0,4}?.*colitis|.*colitis\W+(?:\w+\W+){0,3}?left-sided)\b', line))
        left_colitis_3 = bool(re.search(r'\b(?:left-side\W+(?:\w+\W+){0,4}?.*colitis|.*colitis\W+(?:\w+\W+){0,3}?left-side)\b', line))
        left_colitis = (left_colitis_1|left_colitis_2|left_colitis_3)
        
        mild_col = bool(re.search(r'\b(?:.*mild\W+(?:\w+\W+){0,2}?.*colitis|.*colitis\W+(?:\w+\W+){0,2}?mild)\b', line))
        mod_col = bool(re.search(r'\b(?:.*moderate\W+(?:\w+\W+){0,2}?.*colitis|.*colitis\W+(?:\w+\W+){0,2}?moderate)\b', line))
        sev_col = bool(re.search(r'\b(?:.*severe\W+(?:\w+\W+){0,2}?.*colitis|.*colitis\W+(?:\w+\W+){0,2}?severe)\b', line))
        
        pan_col_1 = bool(re.search(r'\b(?:.*pancolonic\W+(?:\w+\W+){0,2}?.*colitis)\b', line))
        pan_col_2 = bool(re.search(r'\b(?:pan\W+(?:\w+\W+){0,2}?.*colitis)\b', line))
        pan_col = (pan_col_1 or pan_col_2)
        
        act_col = bool(re.search(r'\b(?:.*active\W+(?:\w+\W+){0,2}?.*colitis|.*colitis\W+(?:\w+\W+){0,2}?active)\b', line))
        inact_col = bool(re.search(r'\b(?:.*inactive\W+(?:\w+\W+){0,2}?.*colitis|.*colitis\W+(?:\w+\W+){0,2}?inactive)\b', line))
                
        prcsg_col = bool(re.search(r'\b(?:proctosigmoid\W+(?:\w+\W+){0,2}?.*colitis|.*colitis\W+(?:\w+\W+){0,2}?proctosigmoid)\b', line))

        act_ile = bool(re.search(r'\b(?:.*active\W+(?:\w+\W+){0,2}?.*ileitis|.*ileitis\W+(?:\w+\W+){0,2}?active)\b', line))
        act_ent = bool(re.search(r'\b(?:.*active\W+(?:\w+\W+){0,2}?.*enteritis|.*enteritis\W+(?:\w+\W+){0,2}?active)\b', line))

        noncas_gran_1 = bool(re.search(r'\b(?:noncaseating\W+(?:\w+\W+){0,2}?granuloma|granuloma\W+(?:\w+\W+){0,2}?noncaseating)\b', line))
        noncas_gran_2 = bool(re.search(r'\b(?:noncaseating\W+(?:\w+\W+){0,2}?granulomas|granulomas\W+(?:\w+\W+){0,2}?noncaseating)\b', line)) 
        noncas_gran = (noncas_gran_1 or noncas_gran_2)
        
        nonnec_gran_1 = bool(re.search(r'\b(?:nonnecrotizing\W+(?:\w+\W+){0,2}?granuloma|granuloma\W+(?:\w+\W+){0,2}?nonnecrotizing)\b', line))
        nonnec_gran_2 = bool(re.search(r'\b(?:nonnecrotizing\W+(?:\w+\W+){0,2}?granulomas|granulomas\W+(?:\w+\W+){0,2}?nonnecrotizing)\b', line)) 
        nonnec_gran = (nonnec_gran_1 or nonnec_gran_2)
        
        lamprop_inf = bool(re.search(r'\b(?:lamina-propria\W+(?:\w+\W+){0,2}?inflammation|inflammation\W+(?:\w+\W+){0,3}?lamina-propria)\b', line))
        
        sup_ulcer_1 = bool(re.search(r'\b(?:superficial\W+(?:\w+\W+){0,4}?ulcer|ulcer\W+(?:\w+\W+){0,2}?superficial)\b', line))
        sup_ulcer_2 = bool(re.search(r'\b(?:superficial\W+(?:\w+\W+){0,4}?ulceration|ulceration\W+(?:\w+\W+){0,2}?superficial)\b', line)) 
        sup_ulcer = (sup_ulcer_1 or sup_ulcer_2)
        
        shal_ulcer_1 = bool(re.search(r'\b(?:shallow\W+(?:\w+\W+){0,4}?ulcer|ulcer\W+(?:\w+\W+){0,2}?shallow)\b', line))
        shal_ulcer_2 = bool(re.search(r'\b(?:shallow\W+(?:\w+\W+){0,4}?ulceration|ulceration\W+(?:\w+\W+){0,2}?shallow)\b', line)) 
        shal_ulcer = (shal_ulcer_1 or shal_ulcer_2)
        
        aph_ulcer_1 = bool(re.search(r'\b(?:aphthous\W+(?:\w+\W+){0,3}?ulcer|ulcer\W+(?:\w+\W+){0,1}?aphthous)\b', line))
        aph_ulcer_2 = bool(re.search(r'\b(?:aphthous\W+(?:\w+\W+){0,3}?ulceration|ulceration\W+(?:\w+\W+){0,1}?aphthous)\b', line)) 
        aph_ulcer = (aph_ulcer_1 or aph_ulcer_2)
        
        sml_ulcer_1 = bool(re.search(r'\b(?:small\W+(?:\w+\W+){0,3}?ulcer)\b', line))
        sml_ulcer_2 = bool(re.search(r'\b(?:small\W+(?:\w+\W+){0,3}?ulceration)\b', line))
        sml_ulcer = (sml_ulcer_1 or sml_ulcer_2)
        
        lrg_ulcer_1 = bool(re.search(r'\b(?:large\W+(?:\w+\W+){0,4}?ulcer|ulcer\W+(?:\w+\W+){0,1}?large)\b', line))
        lrg_ulcer_2 = bool(re.search(r'\b(?:large\W+(?:\w+\W+){0,4}?ulceration|ulceration\W+(?:\w+\W+){0,1}?large)\b', line)) 
        lrg_ulcer = (lrg_ulcer_1 or lrg_ulcer_2)
        
        dp_ulcer_1 = bool(re.search(r'\b(?:deep\W+(?:\w+\W+){0,5}?ulcer|ulcer\W+(?:\w+\W+){0,2}?deep)\b', line))
        dp_ulcer_2 = bool(re.search(r'\b(?:deep\W+(?:\w+\W+){0,5}?ulceration|ulceration\W+(?:\w+\W+){0,2}?deep)\b', line)) 
        dp_ulcer = (dp_ulcer_1 or dp_ulcer_2)
        
        
        lamprop_inccel = False
        if 'lamina-propria' in line and 'increase' in line and 'cellularity' in line:
            lamprop_inccel = True
            
#         print(line, end='\n')
    
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
            
            if ('vascular' in e_text or 'vasculature' in e_text) and loss_vasc:
                e_text = e_text + ' (loss-of-vasculature)'
            if ('vascular' in e_text or 'vasculature' in e_text) and dec_vasc:
                e_text = e_text + ' (decreased-vasculature)'
                
            if 'inflammation' in e_text and mild_infl and not 'mild' in e_text:
                e_text = e_text + ' (mild-inflammation)'
            if 'inflammation' in e_text and mod_infl and not 'moderate' in e_text:
                e_text = e_text + ' (moderate-inflammation)'
            if 'inflammation' in e_text and sev_infl and not 'severe' in e_text:
                e_text = e_text + ' (severe-inflammation)'
                
            if 'adherent' in e_text and adhr_blood and (not ('clot' in e_text or 'blood' in e_text)):
                e_text = e_text + ' (adherent-blood)'
                
            if 'colitis' in e_text and pan_col and not 'pancolitis' in e_text:
                 e_text = e_text + ' (pancolitis)'
                    
#             if 'colitis' in e_text and left_colitis:
#                 print('e_text')
                    
            if 'colitis' in e_text and left_colitis and ('surveillance' in line or not (('biops' in line or ' bx ' in line) and 'taken' in line)):
                e_text = e_text + ' (left-sided colitis)'
                
            if 'ileitis' in e_text and act_ile and not 'active' in e_text:
                 e_text = e_text + ' (active)'
                    
            if 'enteritis' in e_text and act_ent and not 'active' in e_text:
                 e_text = e_text + ' (active)'
                    
            if 'granuloma' in e_text and noncas_gran and not 'noncaseating' in e_text:
                e_text = e_text + ' (noncaseating-granuloma)'
                
            if 'granuloma' in e_text and nonnec_gran and not 'nonnecrotizing' in e_text:
                e_text = e_text + ' (nonnecrotizing-granuloma)'
                
            if 'lamina-propria' in e_text and lamprop_inf and not 'inflammation' in e_text:
                e_text = e_text + ' (inflammation)'
                
            if 'lamina-propria' in e_text and lamprop_inccel:
                e_text = e_text + ' (increased-cellularity)'
                
            if 'superficial' in e_text and not 'ulcer' in e_text and sup_ulcer:
                e_text = e_text + ' (superficial-ulcer)'
                
            if 'shallow' in e_text and not 'ulcer' in e_text and shal_ulcer:
                e_text = e_text + ' (shallow-ulcer)'
                
            if 'aphthous' in e_text and not 'ulcer' in e_text and aph_ulcer:
                e_text = e_text + ' (aphthous-ulcer)'
                
            if 'small' in e_text and not 'ulcer' in e_text and sml_ulcer:
                e_text = e_text + ' (small-ulcer)'
                
            if 'large' in e_text and not 'ulcer' in e_text and lrg_ulcer:
                e_text = e_text + ' (large-ulcer)'
                
            if 'deep' in e_text and not 'ulcer' in e_text and dp_ulcer:
                e_text = e_text + ' (deep-ulcer)'
                
    
                
            if 'colitis' in e_text:
                if mild_col and 'mild' not in e_text:
                    e_text = e_text + ' (mild)'
                if mod_col and 'moderate' not in e_text:
                    e_text = e_text + ' (moderate)'
                if sev_col and 'severe' not in e_text:
                    e_text = e_text + ' (severe)'
                if act_col and 'active' not in e_text:
                    e_text = e_text + ' (active)'
                if inact_col and 'inactive' not in e_text:
                    e_text = e_text + ' (inactive)'
                if prcsg_col and 'proctosigmoid' not in e_text:
                    e_text = e_text + ' (proctosigmoiditis)'
                
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
    
    fil_subset = pathdf.MRN_Type.isin(['MGH','BWH','NWH','FH','NSM'])
    df_path = pathdf[fil_subset].copy()
    
    
    if only_liv_biopsy:
        # check the column exists first:
        if 'is_liver_biopsy' in df_path.columns.tolist():
            fil_finaldx_trunc = df_path.is_liver_biopsy == True
            df_path = df_path[fil_finaldx_trunc]
        else:
            print('The flag *only_finaldx=True* was passed, however truncate_finaldx() has not been called. Aborting...')
            return None
        
#     filter_keywords = df_path['Report_Text'].str.contains('steato|balloon|baloon|ballon|inflam|hepatitis|hepatic|fibrosis|bridging|cirrhosis|aih', case=False, na=False)
#     df_path = df_path[filter_keywords]

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
            "pseudo_negations": ts.terms['pseudo_negations'] + ['not limited to', 'not excluded', 'needs to be ruled out', 'although not apparent'],
            "preceding_negations": ts.terms['preceding_negations'] + ['negative', 'insufficient', 'without evidence of', 'rather than'], #'grade 0'
            "following_negations": ts.terms['following_negations'] + ['negative', 'unremarkable', 'ruled out', 'less likely', 'is not', 'are not', 'does not', 'have not', 'was not', 'were not', 'absent', 'grade 0', 'not present'],
            "termination": ts.terms['termination'] + ['note:', ';', ', negative', ',negative'] #'negative for', 'with'
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
    bridging_fibrosis_col = []
    sinusoidal_fibrosis_col = []
    portal_fibrosis_col = []
    periportal_fibrosis_col = []
    pericellular_fibrosis_col = []
    perivenular_fibrosis_col = []
    septal_fibrosis_col = []
    central_fibrosis_col = []
    zone3_fibrosis_col = []
    zone1_fibrosis_col = []
    centrilob_fibrosis_col = []
    hepatitis_col = []
    autoimmune_hepatitis_col = []
    mallory_col = []
    
    pbc_col = []
    cirrhosis_col = []
    steatohepatitis_col = []
    hepatitisa_col = []
    hepatitisb_col = []
    hepatitisc_col = []
    drug_hepatitis_col = []
    interface_hepatitis_col = []
    viral_hepatitis_col = []
    granulomatous_hepatitis_col = []
    hepatic_parenchyma_col = []
    
    hemochromatosis_col = []
    antitrypsin_col = []
    cholangitis_col = []
    wilsons_col = []
    drug_ind_liv_inj_col = []
    budd_chiari_col = []
    alcoholic_col = []
    carcinoma_col = []
    
    nafld_col = []
    nash_col = []
    
    fibrosis_stage_4_col = []
    fibrosis_stage_6_col = []
    
    disease_list_col = []

    for i in range(0,num_reports):

        # extract path report for this entry
        disease_list = []
        report_text = df_path.iloc[i,:].Report_Text
        result_text = entity_recognition_liver(report_text, nlp=nlp_2)

        steatosis = np.nan
        ballooning = np.nan
        inflammation = np.nan
        lobular_inflammation = np.nan
#         portal_inflammation = False
        zone3_inflammation = np.nan
        lobular_hepatitis = np.nan
        zone3_hepatitis = np.nan
        
        fibrosis = np.nan
        bridging_fibrosis = np.nan
        sinusoidal_fibrosis = np.nan
        portal_fibrosis = np.nan
        periportal_fibrosis = np.nan
        pericellular_fibrosis = np.nan
        perivenular_fibrosis = np.nan
        septal_fibrosis = np.nan
        central_fibrosis = np.nan
        zone3_fibrosis = np.nan
        zone1_fibrosis = np.nan
        centrilob_fibrosis = np.nan
        
        hepatitis = np.nan
        autoimmune_hepatitis = np.nan
        mallory = np.nan
        
        pbc = np.nan
        cirrhosis = np.nan
        steatohepatitis = np.nan
        hepatitisa = np.nan
        hepatitisb = np.nan
        hepatitisc = np.nan
        drug_hepatitis = np.nan
        interface_hepatitis = np.nan
        viral_hepatitis = np.nan
        granulomatous_hepatitis = np.nan

        hepatic_parenchyma = np.nan
        hemochromatosis = np.nan
        antitrypsin = np.nan
        cholangitis = np.nan
        wilsons = np.nan
        drug_ind_liv_inj = np.nan
        budd_chiari = np.nan
        alcoholic = np.nan
        carcinoma = np.nan
        
        nafld = np.nan
        nash = np.nan
        
        fibrosis_stage_4 = np.nan
        fibrosis_stage_6 = np.nan
        
        other_liv_diseases = np.nan
        
        steatosis_lt5 = sum([1 for ent_text in result_text.split('\n') if '<5%' in ent_text and 'steatosis' in ent_text])==0
        
        fib_stage = -1
        fib_ref = -1
        
        for x in result_text.split('\n'):
            
            is_disease = False
            
            if 'steatosis' in x and (np.isnan(steatosis) or 'True' in x): #and steatosis_lt5
                if 'True' in x: steatosis, is_disease = True, True
                else: steatosis = False
            if ('balloon' in x or 'baloon' in x or 'ballon' in x) and (np.isnan(ballooning) or 'True' in x):
                if 'True' in x: ballooning, is_disease = True, True
                else: ballooning = False
            if 'inflammation' in x and (np.isnan(inflammation) or 'True' in x):
                if 'True' in x: inflammation, is_disease = True, True
                else: inflammation = False
            if 'lobular' in x and ('inflammation' in x or 'activity' in x or 'infiltrate' in x) and (np.isnan(lobular_inflammation) or 'True' in x):
                if 'True' in x: lobular_inflammation, is_disease = True, True
                else: lobular_inflammation = False
            if 'zone-3' in x and 'inflammation' in x and (np.isnan(zone3_inflammation) or 'True' in x):
                if 'True' in x: zone3_inflammation, is_disease = True, True
                else: zone3_inflammation = False
            if ('lobular' in x and 'hepatitis' in x and not 'steato' in x) and (np.isnan(lobular_hepatitis) or 'True' in x):
                if 'True' in x: lobular_hepatitis, is_disease = True, True
                else: lobular_hepatitis = False
            if ('zone-3' in x and 'hepatitis' in x and not 'steato' in x) and (np.isnan(zone3_hepatitis) or 'True' in x):
                if 'True' in x: zone3_hepatitis, is_disease = True, True
                else: zone3_hepatitis = False
            
            if 'fibrosis' in x and (np.isnan(fibrosis) or 'True' in x):
                if 'True' in x: fibrosis, is_disease = True, True
                else: fibrosis = False
            if ('bridging' in x and not 'bridging-necrosis' in x) and (np.isnan(bridging_fibrosis) or 'True' in x):
                if 'True' in x: bridging_fibrosis, is_disease = True, True
                else: bridging_fibrosis = False
            if ('fibrosis' in x and 'sinusoidal' in x) and (np.isnan(sinusoidal_fibrosis) or 'True' in x):
                if 'True' in x: sinusoidal_fibrosis, is_disease = True, True
                else: sinusoidal_fibrosis = False
            if ('fibrosis' in x and 'portal' in x and not 'peri-portal' in x) and (np.isnan(portal_fibrosis) or 'True' in x):
                if 'True' in x: portal_fibrosis, is_disease = True, True
                else: portal_fibrosis = False
            if ('fibrosis' in x and 'peri-portal' in x) and (np.isnan(periportal_fibrosis) or 'True' in x):
                if 'True' in x: periportal_fibrosis, is_disease = True, True
                else: periportal_fibrosis = False
            if ('fibrosis' in x and 'pericellular' in x) and (np.isnan(pericellular_fibrosis) or 'True' in x):
                if 'True' in x: pericellular_fibrosis, is_disease = True, True
                else: pericellular_fibrosis = False
            if ('fibrosis' in x and 'perivenular' in x) and (np.isnan(perivenular_fibrosis) or 'True' in x):
                if 'True' in x: perivenular_fibrosis, is_disease = True, True
                else: perivenular_fibrosis = False
            if ('fibrosis' in x and 'septal' in x) and (np.isnan(septal_fibrosis) or 'True' in x):
                if 'True' in x: septal_fibrosis, is_disease = True, True
                else: septal_fibrosis = False
            if ('fibrosis' in x and 'central' in x) and (np.isnan(central_fibrosis) or 'True' in x):
                if 'True' in x: central_fibrosis, is_disease = True, True
                else: central_fibrosis = False
            if ('fibrosis' in x and 'zone-3' in x) and (np.isnan(zone3_fibrosis) or 'True' in x):
                if 'True' in x: zone3_fibrosis, is_disease = True, True
                else: zone3_fibrosis = False
            if ('fibrosis' in x and 'zone-1' in x) and (np.isnan(zone1_fibrosis) or 'True' in x):
                if 'True' in x: zone1_fibrosis, is_disease = True, True
                else: zone1_fibrosis = False
            if ('fibrosis' in x and 'centrilobular' in x) and (np.isnan(centrilob_fibrosis) or 'True' in x):
                if 'True' in x: centrilob_fibrosis, is_disease = True, True
                else: centrilob_fibrosis = False
                

            if ('hepatitis' in x and not 'steato' in x) and (np.isnan(hepatitis) or 'True' in x):
                if 'True' in x: hepatitis, is_disease = True, True
                else: hepatitis = False
            if ('autoimmune hepatitis' in x or bool(re.search(r'\baih\b', x)) or 'auto-immune hepatitis' in x) and (np.isnan(autoimmune_hepatitis) or 'True' in x):
                if 'True' in x: autoimmune_hepatitis, is_disease = True, True
                else: autoimmune_hepatitis = False
            if 'mallory' in x and (np.isnan(mallory) or 'True' in x):
                if 'True' in x: mallory, is_disease = True, True
                else: mallory = False
            
            if ('biliary' in x and 'cirrhosis' in x) and (np.isnan(pbc) or 'True' in x):
                if 'True' in x: pbc, is_disease = True, True
                else: pbc = False
            if ('cirrhosis' in x and 'biliary' not in x) and (np.isnan(cirrhosis) or 'True' in x):
                if 'True' in x: cirrhosis, is_disease = True, True
                else: cirrhosis = False
            if 'steatohepatitis' in x and (np.isnan(steatohepatitis) or 'True' in x):
                if 'True' in x: steatohepatitis, is_disease = True, True
                else: steatohepatitis = False
                    
#             if 'portal' in x and 'inflammation' in x and 'True' in x:
#                 portal_inflammation = True
#             if ('chronic' in x and 'hepatitis' in x) and (np.isnan(chronic_hepatitis) or 'True' in x):
#                 if 'True' in x: chronic_hepatitis, is_disease = True, True
#                 else: chronic_hepatitis = False

            if (bool(re.search(r'\bhepatitis a\b', x)) or bool(re.search(r'\bhepatitis-a\b', x)) or bool(re.search(r'\bhep a\b', x))) and (np.isnan(hepatitisa) or 'True' in x):
                if 'True' in x: hepatitisa, is_disease = True, True
                else: hepatitisa = False
            if (bool(re.search(r'\bhepatitis b\b', x)) or bool(re.search(r'\bhepatitis-b\b', x)) 
                or bool(re.search(r'\bhep b\b', x)) or bool(re.search(r'\bhbv\b', x))) and (np.isnan(hepatitisb) or 'True' in x):
                if 'True' in x: hepatitisb, is_disease = True, True
                else: hepatitisb = False
            # Need to add bool(re.search(r'\bchronichepatitis c\b', x)) bool(re.search(r'\bhepatitis c\b', x))
            # bool(re.search(r'\bhepatitis c\b', x))
            if (bool(re.search(r'\bhepatitis c\b', x)) or bool(re.search(r'\bhepatitis-c\b', x)) 
                    or bool(re.search(r'\bhcv\b', x)) or bool(re.search(r'\bhep c\b', x))
                    or 'ishak' in x) and (np.isnan(hepatitisc) or 'True' in x):
                if 'True' in x: hepatitisc, is_disease = True, True
                else: hepatitisc = False
            if ('drug' in x and 'hepatitis' in x) and (np.isnan(drug_hepatitis) or 'True' in x):
                if 'True' in x: drug_hepatitis, is_disease = True, True
                else: drug_hepatitis = False
            if ('interface' in x and 'hepatitis' in x) and (np.isnan(interface_hepatitis) or 'True' in x):
                if 'True' in x: interface_hepatitis, is_disease = True, True
                else: inflammation = False
            if ('viral' in x and 'hepatitis' in x) and (np.isnan(viral_hepatitis) or 'True' in x):
                if 'True' in x: viral_hepatitis, is_disease = True, True
                else: viral_hepatitis = False
            if ('granulomatous' in x and 'hepatitis' in x) and (np.isnan(granulomatous_hepatitis) or 'True' in x):
                if 'True' in x: granulomatous_hepatitis, is_disease = True, True
                else: granulomatous_hepatitis = False
            
            if ('hepatic' in x and 'parenchyma' in x) and (np.isnan(hepatic_parenchyma) or 'True' in x):
                if 'True' in x: hepatic_parenchyma, is_disease = True, True
                else: hepatic_parenchyma = False
            if 'hemochromatosis' in x and (np.isnan(hemochromatosis) or 'True' in x):
                if 'True' in x: hemochromatosis, is_disease = True, True
                else: hemochromatosis = False
            if 'antitrypsin' in x and (np.isnan(antitrypsin) or 'True' in x):
                if 'True' in x: antitrypsin, is_disease = True, True
                else: antitrypsin = False
            if 'cholangitis' in x and (np.isnan(cholangitis) or 'True' in x):
                if 'True' in x: cholangitis, is_disease = True, True
                else: cholangitis = False
            if "wilson's" in x and (np.isnan(wilsons) or 'True' in x):
                if 'True' in x: wilsons, is_disease = True, True
                else: wilsons = False
            if ('drug-induced' in x or bool(re.search(r'\bdili\b', x)) or 'drug-related' in x) and (np.isnan(drug_ind_liv_inj) or 'True' in x):
                if 'True' in x: drug_ind_liv_inj, is_disease = True, True
                else: drug_ind_liv_inj = False
            if 'budd-chiari' in x and (np.isnan(budd_chiari) or 'True' in x):
                if 'True' in x: budd_chiari, is_disease = True, True
                else: budd_chiari = False
            if bool(re.search(r'\balcoholic\b', x)) and (np.isnan(alcoholic) or 'True' in x):
                if 'True' in x: alcoholic, is_disease = True, True
                else: alcoholic = False
            if ('metastatic' in x or 'metastases' in x or 'metastasis' in x or 'carcinoma' in x or 'lymphoma' in x
                       or bool(re.search(r'\bhcc\b', x)) or 'malign' in x or 'cancer' in x or 'carcinoid' in x
                       or 'angiosarcoma' in x) and (np.isnan(carcinoma) or 'True' in x):
                if 'True' in x: carcinoma, is_disease = True, True
                else: carcinoma = False
            
            if (bool(re.search(r'\bnafld\b', x)) or 'nonalcoholic fatty liver disease' in x) and (np.isnan(nafld) or 'True' in x):
                if 'True' in x: nafld, is_disease = True, True
                else: nafld = False
            if (bool(re.search(r'\bnash\b', x)) or 'nonalcoholic steatohepatitis' in x) and (np.isnan(nash) or 'True' in x):
                if 'True' in x: nash, is_disease = True, True
                else: nash = False
                    
            if ('fibrosis' in x or 'bridging' in x or 'cirrhosis' in x) and ' stage:' in x:
                try:
                    if 'True' in x:
                        fib_stage = float(x[-12:-9])
                        fib_ref = float(x[-8:-5])
                    elif 'False' in x:
                        fib_stage = float(x[-13:-10])
                        fib_ref = float(x[-9:-6])

                    if fib_ref<4.0:
                        fib_ref = 4.0
                        
                    if fib_ref==4 and np.isnan(fibrosis_stage_4):
                        fibrosis_stage_4 = fib_stage
                    elif fib_ref==6 and np.isnan(fibrosis_stage_6):
                        fibrosis_stage_6 = fib_stage
                    
                    is_disease = True
                    
                except:
                    pass
                
            if is_disease:
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
        
        
        bridging_fibrosis_col.append(bridging_fibrosis)
        sinusoidal_fibrosis_col.append(sinusoidal_fibrosis)
        portal_fibrosis_col.append(portal_fibrosis)
        periportal_fibrosis_col.append(periportal_fibrosis)
        pericellular_fibrosis_col.append(pericellular_fibrosis)
        perivenular_fibrosis_col.append(perivenular_fibrosis)
        septal_fibrosis_col.append(septal_fibrosis)
        central_fibrosis_col.append(central_fibrosis)
        zone3_fibrosis_col.append(zone3_fibrosis)
        zone1_fibrosis_col.append(zone1_fibrosis)
        centrilob_fibrosis_col.append(centrilob_fibrosis)
        hepatitis_col.append(hepatitis)
        autoimmune_hepatitis_col.append(autoimmune_hepatitis)
        mallory_col.append(mallory)
        
        pbc_col.append(pbc)
        cirrhosis_col.append(cirrhosis)
        steatohepatitis_col.append(steatohepatitis)
        hepatitisa_col.append(hepatitisa)
        hepatitisb_col.append(hepatitisb)
        hepatitisc_col.append(hepatitisc)
        drug_hepatitis_col.append(drug_hepatitis)
        interface_hepatitis_col.append(interface_hepatitis)
        viral_hepatitis_col.append(viral_hepatitis)
        granulomatous_hepatitis_col.append(granulomatous_hepatitis)
        hepatic_parenchyma_col.append(hepatic_parenchyma)
        
        hemochromatosis_col.append(hemochromatosis)
        antitrypsin_col.append(antitrypsin)
        cholangitis_col.append(cholangitis)
        wilsons_col.append(wilsons)
        drug_ind_liv_inj_col.append(drug_ind_liv_inj)
        budd_chiari_col.append(budd_chiari)
        alcoholic_col.append(alcoholic)
        carcinoma_col.append(carcinoma)
        
        nafld_col.append(nafld)
        nash_col.append(nash)
        
        fibrosis_stage_4_col.append(fibrosis_stage_4)
        fibrosis_stage_6_col.append(fibrosis_stage_6)
        
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
    df_path['bridging_fibrosis'] = bridging_fibrosis_col
    df_path['sinusoidal_fibrosis'] = sinusoidal_fibrosis_col
    df_path['portal_fibrosis'] = portal_fibrosis_col
    df_path['periportal_fibrosis'] = periportal_fibrosis_col
    df_path['pericellular_fibrosis'] = pericellular_fibrosis_col
    df_path['perivenular_fibrosis'] = perivenular_fibrosis_col
    df_path['septal_fibrosis'] = septal_fibrosis_col
    df_path['central_fibrosis'] = central_fibrosis_col
    df_path['zone3_fibrosis'] = zone3_fibrosis_col
    df_path['zone1_fibrosis'] = zone1_fibrosis_col
    df_path['centrilob_fibrosis'] = centrilob_fibrosis_col
    df_path['hepatitis'] = hepatitis_col
    df_path['autoimmune_hepatitis'] = autoimmune_hepatitis_col
    df_path['mallory'] = mallory_col
    
    df_path['pbc'] = pbc_col
    df_path['cirrhosis'] = cirrhosis_col
    df_path['steatohepatitis'] = steatohepatitis_col
    df_path['hepatitisa'] = hepatitisa_col
    df_path['hepatitisb'] = hepatitisb_col
    df_path['hepatitisc'] = hepatitisc_col
    df_path['drug_hepatitis'] = drug_hepatitis_col
    df_path['interface_hepatitis'] = interface_hepatitis_col
    df_path['viral_hepatitis'] = viral_hepatitis_col
    df_path['granulomatous_hepatitis'] = granulomatous_hepatitis_col
    df_path['hepatic_parenchyma'] = hepatic_parenchyma_col
    
    df_path['hemochromatosis'] = hemochromatosis_col
    df_path['antitrypsin'] = antitrypsin_col
    df_path['cholangitis'] = cholangitis_col
    df_path['wilsons'] = wilsons_col
    df_path['drug_ind_liv_inj'] = drug_ind_liv_inj_col
    df_path['budd_chiari'] = budd_chiari_col
    df_path['alcoholic'] = alcoholic_col
    df_path['carcinoma'] = carcinoma_col
    
    df_path['nafld'] = nafld_col
    df_path['nash'] = nash_col
    
    df_path['fibrosis_stage_4'] = fibrosis_stage_4_col
    df_path['fibrosis_stage_6'] = fibrosis_stage_6_col
        
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
        pathdf['bridging_fibrosis'] = np.nan
        pathdf['sinusoidal_fibrosis'] = np.nan
        pathdf['portal_fibrosis'] = np.nan
        pathdf['periportal_fibrosis'] = np.nan
        pathdf['pericellular_fibrosis'] = np.nan
        pathdf['perivenular_fibrosis'] = np.nan
        pathdf['septal_fibrosis'] = np.nan
        pathdf['central_fibrosis'] = np.nan
        pathdf['zone3_fibrosis'] = np.nan
        pathdf['zone1_fibrosis'] = np.nan
        pathdf['centrilob_fibrosis'] = np.nan
        pathdf['hepatitis'] = np.nan
        pathdf['autoimmune_hepatitis'] = np.nan
        pathdf['mallory'] = np.nan
        
        pathdf['pbc'] = np.nan
        pathdf['cirrhosis'] = np.nan
        pathdf['steatohepatitis'] = np.nan
        pathdf['hepatitisa'] = np.nan
        pathdf['hepatitisb'] = np.nan
        pathdf['hepatitisc'] = np.nan
        pathdf['drug_hepatitis'] = np.nan
        pathdf['interface_hepatitis'] = np.nan
        pathdf['viral_hepatitis'] = np.nan
        pathdf['granulomatous_hepatitis'] = np.nan
        pathdf['hepatic_parenchyma'] = np.nan
        
        pathdf['hemochromatosis'] = np.nan
        pathdf['antitrypsin'] = np.nan
        pathdf['cholangitis'] = np.nan
        pathdf['wilsons'] = np.nan
        pathdf['drug_ind_liv_inj'] = np.nan
        pathdf['budd_chiari'] = np.nan
        pathdf['alcoholic'] = np.nan
        pathdf['carcinoma'] = np.nan
        
        pathdf['nafld'] = np.nan
        pathdf['nash'] = np.nan
        
        pathdf['fibrosis_stage_4'] = np.nan
        pathdf['fibrosis_stage_6'] = np.nan
            
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
    
    text = (text
            .replace('Hep. A', 'Hep A').replace('Hep. B', 'Hep B').replace('Hep. C', 'Hep C')
            .replace('/IV', '/4').replace('/VI', '/6')
            .replace('PBC', 'PBC (primary biliary cirrhosis)')
            .replace('PSC', 'PSC (primary sclerosing cholangitis)')
           )
    
    text = text.lower().replace(' bridging.', ' active-bridging.') #.replace(';',' ')
    
    entity_result = ''
    
    for line in text.split('.'):
        
        line = " ".join(line.split())
        line = line.strip()
        line = (line
                .replace('+/-', ',')
                .replace(' no ', ' , no ')
                .replace('no present', 'not present')
                .replace('is not present', 'not present').replace('not present', ' is not present')
                .replace(' minimal ', ' ,minimal ')
                .replace('noted in the', ' in ')
                .replace('is noted in', ' in ')
                .replace(' as well as', ', ')
                .replace('may not', 'will not')
                .replace('neither', 'no').replace('nor', 'no')
                .replace('very', '')
                .replace('mildly', 'mild').replace('mildl', 'mild')
                .replace('non alcoholic', 'nonalcoholic')
                .replace('non-alcoholic', 'nonalcoholic')
                .replace('steato-hepatitis', 'steatohepatitis')
                .replace('steato hepatitis', 'steatohepatitis')
                .replace('nonalcoholic steatohepatitis', 'nonalcoholic-steatohepatitis')
                .replace('inflammatory infiltrate', 'inflammatory-infiltrate')
                .replace('necroinflammatory', 'necro-inflammation')
                .replace('centric inflammation', 'centric-inflammation')
                .replace('inflammatory', 'inflammation')
                .replace('inflamed', 'inflammation')
                .replace('severely', 'severe')
                .replace('moderately', 'moderate')
                .replace('moderate ', 'moderate_').replace('(moderate) ', 'moderate_').replace('(moderate)', 'moderate_')
                .replace('mild to moderate', 'mild&moderate')
                .replace('moderate to severe', 'moderate&severe')
                .replace('mild to severe', 'mild&severe')
                .replace('mild active', 'mild-active')
                .replace('mild chronic', 'mild-chronic')
                .replace('mild ', 'mild_').replace('(mild) ', 'mild_').replace('(mild)', 'mild_')
                .replace('severe ', 'severe_').replace('(severe) ', 'severe_').replace('(severe)', 'severe_')
                .replace('minimal ', 'minimal_')
                .replace('chronic ', 'chronic_')
                .replace('focal ', 'focal_')
                .replace(' areas', ' area')
                .replace('-area ', ' area ')
                .replace('portal area', 'portalarea')
                .replace(' tracts', 'tract')
                .replace('-tracts', 'tract')
                .replace('portal tract', 'portaltract')
#                 .replace('centrilobular', 'centri-lobular')
                .replace('periportal', 'peri-portal')
                .replace('portal and lobular', 'portal&lobular')
                .replace('portal or lobular', 'portal&lobular')
                .replace('lobular and portal', 'lobular&portal')
                .replace('lobular or portal', 'lobular&portal')
                .replace('portal&lobular', 'lobular&portal')
                .replace('lobular inflammat', 'lobular-inflammat')
                .replace('mixed ', 'mixed-')
                .replace('kupffer cell', 'kupffer-cell')
                .replace('zone 3', 'zone-3').replace('zone3', 'zone-3')
                .replace('zone 1', 'zone-1').replace('zone1', 'zone-1')
                .replace('hepatic plate', 'hepatic-plate')
                .replace('hepatic parenchyma', 'hepatic-parenchyma')
#                 .replace('hepatic ', 'hepatitis ')
                .replace('steatotic', 'steatosis')
#                 .replace('microvesicular ', 'microvesicular-')
#                 .replace('macrovesicular ', 'macrovesicular-')
                .replace('< 5%', '<5%')
                .replace('non classical', 'nonclassical')
                .replace('non-classical', 'nonclassical')
                .replace('ballooning degeneration', 'ballooning_degeneration')
                .replace('hepatocyte ballooning', 'hepatocyte-ballooning')
                .replace('hepatocytic ballooning', 'hepatocytic-ballooning')
                .replace('heptocellular ballooning', 'heptocellular-ballooning')
                .replace('ballooned', 'ballooning')
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
                .replace('centrilobular fibrosis', 'centrilobular-fibrosis')
                .replace('septal fibrosis', 'septal-fibrosis')
                .replace('fibrous septa', 'septal-fibrosis')
                .replace('central fibrosis', 'central-fibrosis')
                .replace('ductal fibrosis', 'ductal-fibrosis')
                .replace('portal bridging', 'portal-bridging')
                .replace('central bridging', 'central-bridging')
                .replace(' bridging ', ' active-bridging ')
                .replace(' bridging;', ' active-bridging;')
                .replace(' bridging,', ' active-bridging,')
                .replace('(p-p)','')
                .replace(' hep ', ' hepatitis ')
                .replace('a1at', 'alpha-1-antitrypsin')
#                 .replace('alpha-1 antitrypsin', 'alpha-1-antitrypsin')
#                 .replace('alpha 1 antitrypsin', 'alpha-1-antitrypsin')
#                 .replace('alpha - 1 antitrypsin', 'alpha-1-antitrypsin')
#                 .replace('alpha 1-antitrypsin', 'alpha-1-antitrypsin')
                .replace('wilsons disease', "wilson's disease")
                .replace('wilson disease', "wilson's disease")
#                 .replace('drug induced liver injury', 'drug-induced-liver-injury')
#                 .replace('drug-induced liver injury', 'drug-induced-liver-injury')
#                 .replace('drug induced cholestatic liver injury', 'drug-induced-liver-injury')
#                 .replace('drug induced injury', 'drug-induced-liver-injury')
#                 .replace('drug induced hepatocellular injury', 'drug-induced-liver-injury')
#                 .replace('drug induced lesion', 'drug-induced-liver-injury')
                .replace('drug induced', 'drug-induced').replace('drug related', 'drug-related')
                .replace('budd chiari', 'budd-chiari')
                .replace('hepatocellular carcinoma', 'hepatocellular-carcinoma')
                .replace('cirrhoses', 'cirrhosis')
                .replace('hepatitis-a', 'hepatitis a')
                .replace('hepatitis-b', 'hepatitis b')
                .replace('hepatitis-c', 'hepatitis c')
                .replace('steato-hepatitis', 'steatohepatitis').replace('steato hepatitis', 'steatohepatitis')
                .replace('centri-lobular', 'centrilobular')
                .replace('centrilobular necrosis', 'centrilobular-necrosis')
                .replace('centrilobular hepatic necrosis', 'centrilobular-hepatic-necrosis')
                .replace('droplet steatosis', 'droplet-steatosis')
#                 .replace('_inflammation', ' iniflammation')
                .replace('_lobular', ' lobular')
                .replace('_portal', ' portal')
                .replace('early cirrhotic', 'early cirrhotic (cirrhosis)')
                .replace('primary biliary cirrhosis', 'primary_biliary_cirrhosis')
                .replace('biliary cirrhosis', 'biliary_cirrhosis')
                .replace('billiary cirrhosis', 'biliary_cirrhosis')
                .replace('AMA negative', 'AMA-negative')
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
            
            lobu_infil = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,4}?infiltrate)\b', line))
            
            
            zone3_inf = bool(re.search(r'\b(?:zone-3\W+(?:\w+\W+){0,2}?inflammation|inflammation\W+(?:\w+\W+){0,1}?zone-3)\b', line))
            
            # port_inf = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,3}?inflammation|inflammation\W+(?:\w+\W+){0,1}?portal)\b', line))
            port_inf = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,3}?inflammation)\b', line))
            
            port_tract_inf = bool(re.search(r'\b(?:portaltract\W+(?:\w+\W+){0,5}?inflammation)\b', line))
            
            lob_port_inf = bool(re.search(r'\b(?:lobular&portal\W+(?:\w+\W+){0,3}?inflammation)\b', line))

            lobu_act = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,3}?activity)\b', line))
            port_act = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,3}?activity)\b', line))

            lobu_hep = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,2}?hepatitis)\b', line))
            
            zone3_hep = bool(re.search(r'\b(?:zone-3\W+(?:\w+\W+){0,3}?hepatitis)\b', line))
            
            lobu_dis = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,3}?disarray)\b', line)) 
            
            zone3_fib = bool(re.search(r'\b(?:zone-3\W+(?:\w+\W+){0,3}?fibrosis|fibrosis\W+(?:\w+\W+){0,3}?zone-3)\b', line))
            
            zone1_fib = bool(re.search(r'\b(?:zone-1\W+(?:\w+\W+){0,3}?fibrosis|fibrosis\W+(?:\w+\W+){0,3}?zone-1)\b', line))
            
            cent_scar_1 = bool(re.search(r'\b(?:central\W+(?:\w+\W+){0,2}?scarring|scarring\W+(?:\w+\W+){0,4}?central)\b', line))
            cent_scar_2 = bool(re.search(r'\b(?:centrilobular\W+(?:\w+\W+){0,3}?scarring|scarring\W+(?:\w+\W+){0,4}?centrilobular)\b', line))
            cent_scar_3 = bool(re.search(r'\b(?:pericentral\W+(?:\w+\W+){0,3}?scarring|scarring\W+(?:\w+\W+){0,4}?pericentral)\b', line))
            cent_scar = cent_scar_1 or cent_scar_3 or cent_scar_3
            
            port_scar_1 = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,3}?scarring|scarring\W+(?:\w+\W+){0,1}?portal)\b', line))
            port_scar_2 = bool(re.search(r'\b(?:portaltract\W+(?:\w+\W+){0,3}?scarring|scarring\W+(?:\w+\W+){0,1}?portaltract)\b', line))
            port_scar_3 = bool(re.search(r'\b(?:portalarea\W+(?:\w+\W+){0,3}?scarring|scarring\W+(?:\w+\W+){0,1}?portalarea)\b', line))
            port_scar_4 = bool(re.search(r'\b(?:periportal\W+(?:\w+\W+){0,3}?scarring|scarring\W+(?:\w+\W+){0,1}?periportal)\b', line))
            port_scar = port_scar_1 or port_scar_2 or port_scar_3 or port_scar_4
            
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
            centrilob_fib = bool(re.search(r'\b(?:centrilobular\W+(?:\w+\W+){0,2}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?centrilobular)\b', line))

            
            if 'steatosis' in line and ('<5%' in line or 'less than 5%' in line) and 'steatosis' in e_text: #('5-33%' not in line)
                e_text = '<5% ' + e_text
                
            if 'lobular&portal' in line and 'inflammation' in e_text and not 'lobular&portal' in e_text and lob_port_inf:
                e_text = 'lobular&portal ' + e_text
            
            if 'inflammation' in e_text and not 'lobular' in e_text and lobu_inf and inf_count<=1:
                e_text = 'lobular ' + e_text
            
            if 'inflammation' in e_text and not 'lobular' in e_text and not 'portal' in e_text and not 'lobular-inflammation' in line and lobu_inf and inf_count>1:
                e_text = 'lobular ' + e_text
            
            if 'inflammation' in e_text and not 'portal' in e_text and port_inf and inf_count<=1:
                e_text = 'portal ' + e_text
                
            if 'inflammation' in e_text and not 'portal' in e_text and not 'lobular' in e_text and not 'portal inflammation' in line and port_inf and inf_count>1:
                e_text = 'portal ' + e_text
            
            if ('portaltract' in line or 'portalarea' in line) and 'inflammation' in e_text and not 'portal' in e_text  and (port_inf or port_tract_inf):
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
                
            if 'fibrosis' in e_text and zone3_fib and not 'zone-3' in e_text:
                e_text = 'zone-3 ' + e_text
            
            if 'fibrosis' in e_text and zone1_fib and not 'zone-1' in e_text:
                e_text = 'zone-1 ' + e_text
                
            if 'fibrosis' in e_text and centrilob_fib and not 'centrilobular' in e_text:
                e_text = 'centrilobular ' + e_text
                
            if 'scarring' in e_text and cent_scar and not 'central' in e_text:
                e_text = 'central ' + e_text
                
            if 'scarring' in e_text and port_scar and not 'portal' in e_text:
                e_text = 'portal ' + e_text
                
            if 'infiltrate' in e_text and lobu_infil and not 'lobular' in e_text:
                e_text = 'lobular ' + e_text
            
            
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
                
                stagelist_ishak = re.findall(r'ishak.*?(\d+(?:\.\d+)?).*?(\d+(?:\.\d+)?)', line)
                
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
                    
                    if len(stagelist_ishak)==1:
                        stage_val = float(stagelist_ishak[0][0])
                        ref_val = float(stagelist_ishak[0][1])

                    if stage_val>=5:
                        ref_val = 6.0
                    
                    if ref_val<=6 and stage_val<=ref_val and (ref_val in [4,6]):
                        
                        if 'ishak' in line:
                            e_text = e_text + ' ishak'
                            ref_val = 6.0

                        e_text = e_text + ' stage: ' + str(stage_val) + '/' + str(ref_val)

                    if stage_val==0:
                        e_bool = True
                    else:
                        e_bool = False

            e_text = " ".join(e_text.split())
            
            
            entity_result = entity_result + e_text + ' ' + str(not e_bool) + '\n'
            
            entity_result = entity_result.replace('baloon', 'balloon').replace('ballon', 'balloon').replace('balon', 'balloon')
    
        
    return entity_result
