
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
            "preceding_negations": ts.terms['preceding_negations'] + ['negative', 'insufficient', 'without evidence of', 'rather than', 'history', 'precludes'],
            "following_negations": ts.terms['following_negations'] + ['negative', 'ruled out', 'less likely', 'is not', 'are not', 'does not', 'have not', 'was not', 'were not', 'absent', 'not present'], #unremarkable
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
    microscopic_colitis_col = []
    collagenous_colitis_col = []
    lymphocytic_colitis_col = []
    pan_colitis_col = []
    proctitis_col = []
    ulc_proctitis_col = []
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
        microscopic_colitis = np.nan
        collagenous_colitis = np.nan
        lymphocytic_colitis = np.nan
        pan_colitis = np.nan
        proctitis = np.nan
        ulc_proctitis = np.nan
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
            if (('mild' in x or 'minimal' in x) and ('inflammat' in x or 'inflamed' in x)) and (np.isnan(mild_inflammation) or 'True' in x):
                if 'True' in x: mild_inflammation, is_disease = True, True
                else: mild_inflammation = False
            if ('moderate' in x and ('inflammat' in x or 'inflamed' in x)) and (np.isnan(moderate_inflammation) or 'True' in x):
                if 'True' in x: moderate_inflammation, is_disease = True, True
                else: moderate_inflammation = False
            if (('severe' in x or 'very' in x) and ('inflammat' in x or 'inflamed' in x)) and (np.isnan(severe_inflammation) or 'True' in x):
                if 'True' in x: severe_inflammation, is_disease = True, True
                else: severe_inflammation = False
            if (('loss' in x or 'absent' in x) and ('vasculature' in x or 'vascular' in x)) and (np.isnan(loss_vasculature) or 'True' in x):
                if 'True' in x: loss_vasculature, is_disease = True, True
                else: loss_vasculature = False
            if (('decreased' in x or 'altered' in x) and ('vasculature' in x or 'vascular' in x)) and (np.isnan(dec_vasculature) or 'True' in x):
                if 'True' in x: dec_vasculature, is_disease = True, True
                else: dec_vasculature = False
            if ('granular' in x or 'granulation' in x) and (np.isnan(granularity) or 'True' in x):
                if 'True' in x: granularity, is_disease = True, True
                else: granularity = False
            if ('ulcer' in x and not 'ulcerativecolitis' in x) and (np.isnan(ulceration) or 'True' in x):
                if 'True' in x: ulceration, is_disease = True, True
                else: ulceration = False
            if (('friable' in x or 'friability' in x) and not ('mild' in x or 'minimal' in x)) and (np.isnan(friability) or 'True' in x):
                if 'True' in x: friability, is_disease = True, True
                else: friability = False
            if (('friable' in x or 'friability' in x) and ('mild' in x or 'minimal' in x)) and (np.isnan(mild_friability) or 'True' in x):
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
            if (('mild' in x or 'minimal' in x) and (bool(re.search(r'\bulcer\b|\bulceration\b', x)))) and (np.isnan(mild_ulcer) or 'True' in x):
                if 'True' in x: mild_ulcer, is_disease = True, True
                else: mild_ulcer = False
                    
            
            if 'colitis' in x and (np.isnan(colitis) or 'True' in x):
                if 'True' in x: colitis, is_disease = True, True
                else: colitis = False
            if ('chronic' in x and 'colitis' in x) and (np.isnan(chronic_colitis) or 'True' in x):
                if 'True' in x: chronic_colitis, is_disease = True, True
                else: chronic_colitis = False
            if (('mild' in x or 'minimal' in x) and 'colitis' in x) and (np.isnan(mild_colitis) or 'True' in x):
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
            if ('microscopic' in x and 'colitis' in x) and (np.isnan(microscopic_colitis) or 'True' in x):
                if 'True' in x: microscopic_colitis, is_disease = True, True
                else: microscopic_colitis = False
            if ('collagenous' in x and 'colitis' in x) and (np.isnan(collagenous_colitis) or 'True' in x):
                if 'True' in x: collagenous_colitis, is_disease = True, True
                else: collagenous_colitis = False
            if ('lymphocytic' in x and 'colitis' in x) and (np.isnan(lymphocytic_colitis) or 'True' in x):
                if 'True' in x: lymphocytic_colitis, is_disease = True, True
                else: lymphocytic_colitis = False
            if ('pancolitis' in x) and (np.isnan(pan_colitis) or 'True' in x):
                if 'True' in x: pan_colitis, is_disease = True, True
                else: pan_colitis = False
            if 'proctitis' in x and (np.isnan(proctitis) or 'True' in x):
                if 'True' in x: proctitis, is_disease = True, True
                else: proctitis = False
            if ('proctitis' in x and 'ulcerative' in x) and (np.isnan(ulc_proctitis) or 'True' in x):
                if 'True' in x: ulc_proctitis, is_disease = True, True
                else: ulc_proctitis = False
            if (('proctosigmoiditis' in x) or ('pro' in x and 'sigmoiditis' in x)) and (np.isnan(proctosigmoiditis) or 'True' in x):
                if 'True' in x: proctosigmoiditis, is_disease = True, True
                else: proctosigmoiditis = False
            if ('left-side' in x and 'colitis' in x) and (np.isnan(left_sided_colitis) or 'True' in x):
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
        severe_inflammation_col.append(severe_inflammation)
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
        microscopic_colitis_col.append(microscopic_colitis)
        collagenous_colitis_col.append(collagenous_colitis)
        lymphocytic_colitis_col.append(lymphocytic_colitis)
        pan_colitis_col.append(pan_colitis)
        proctitis_col.append(proctitis)
        ulc_proctitis_col.append(ulc_proctitis)
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
    df_path['microscopic_colitis'] = microscopic_colitis_col
    df_path['collagenous_colitis'] = collagenous_colitis_col
    df_path['lymphocytic_colitis'] = lymphocytic_colitis_col
    df_path['pan_colitis'] = pan_colitis_col
    df_path['proctitis'] = proctitis_col
    df_path['ulc_proctitis'] = ulc_proctitis_col
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
        pathdf['microscopic_colitis'] = np.nan
        pathdf['collagenous_colitis'] = np.nan
        pathdf['lymphocytic_colitis'] = np.nan
        pathdf['pan_colitis'] = np.nan
        pathdf['proctitis'] = np.nan
        pathdf['ulc_proctitis'] = np.nan
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
            .replace('Impression:', '.Impression:')
           )
    text = re.sub(' UC$', ' UlcerativeColitis', text)
    text = re.sub('  - ', '.', text)
    
    text = text.lower()
        
    entity_result = ''
    mayo_score = -1
    mayo_bool = False
    pointer=0
    inflam_pointer=0
    
    for line in text.split('.'):
                
        #line = line.strip()
        line = " ".join(line.split())
#         line = re.sub(r'^-', r'', line)
#         line = re.sub(r'^impression: -', r'impression: ', line)
        line = (line
                .replace('+/-', ',')
                .replace(' no ', ' , no ')
                .replace('is not present', 'not present').replace('no present', 'not present').replace('not present', ' is not present')
                .replace(' minimal ', ' ,minimal ')
                .replace('ulcerations', 'ulceration')
                .replace('ulcers', 'ulcer')
                .replace('apthous', 'aphthous')
                .replace('colitis-', 'colitis ')
                .replace('ulcerative colitis', 'ulcerativecolitis')
                .replace('microscopic colitis', 'microscopic_colitis')
                .replace('collagenous colitis', 'collagenous_colitis')
                .replace('lymphocytic colitis', 'lymphocytic_colitis')
#                 .replace('active ulcerativecolitis', 'ulcerativecolitis (active)')
                .replace('healing ulcerativecolitis', 'healing-ulcerativecolitis')
                .replace('pan colitis', 'pancolitis').replace('pan-colitis', 'pancolitis')
                .replace('noted in the', ' in ')
                .replace('is noted in', ' in ')
                .replace(' neither ', ' no ').replace(' nor ', ' no ')
                .replace('may not', 'will not')
#                 .replace('chronic active', 'chronic-active')
#                 .replace('active chronic', 'active-chronic')
#                 .replace('chronic inactive', 'chronic-inactive')
                .replace('severely', 'severe')
                .replace('moderately severe', 'moderate')
                .replace('moderately', 'moderate')
                .replace('mildly', 'mild').replace('mildl', 'mild').replace('midly', 'mild')
                .replace('minimally', 'minimal')
                .replace('floridly', 'florid')
                .replace('severe pseudomembranous', 'severe-pseudomembranous')
                .replace('self limited', 'self-limited')
#                 .replace('moderate to severe', 'moderate&severe')
#                 .replace('moderate to focally severe', 'moderate&severe')
#                 .replace('mild to moderate', 'mild&moderate')
#                 .replace('mild to focally moderate', 'mild&moderate')
#                 .replace('mild to severe', 'mild&severe')
#                 .replace('mild to focally severe', 'mild&severe')
                .replace('chronic-inactive', 'chronic inactive')
#                 .replace('active chronic', 'active-chronic')
#                 .replace('acute and chronic', 'acute-chronic')
#                 .replace('acute on chronic', 'acute-chronic')
#                 .replace('severe active', 'severe-active')
#                 .replace('severe chronic', 'severe-chronic')
#                 .replace('active severe', 'active-severe')
#                 .replace('chronic severe', 'chronic-severe')
#                 .replace('pancolitis, moderate&severe', 'moderate&severe pancolitis')
#                 .replace('colitis, moderate&severe', 'moderate&severe colitis')
#                 .replace('colitis, severe', 'severe colitis')
#                 .replace('active-severe', 'severe-active')
#                 .replace('colitis, moderate', 'moderate colitis')
#                 .replace('active-moderate', 'moderate-active')
#                 .replace('severe ischemic', 'severe-ischemic')
#                 .replace('moderate active', 'moderate-active')
#                 .replace('moderate chronic', 'moderate-chronic')
#                 .replace('active moderate', 'active-moderate')
#                 .replace('chronic moderate', 'chronic-moderate')
#                 .replace('mild active', 'mild-active')
#                 .replace('mild chronic', 'mild-chronic')
#                 .replace(' active colitis', ' active-colitis')
                .replace('(mild)', 'mild')
                .replace('mild and ', 'mild ').replace('mild ', 'mild and ')
                .replace('(moderate)', 'moderate')
                .replace('(severe)', 'severe')
#                 .replace('mild ', 'mild_').replace('(mild) ', 'mild_').replace('(mild)', 'mild_')
#                 .replace('moderate ', 'moderate_').replace('(moderate) ', 'moderate_').replace('(moderate)', 'moderate_')
#                 .replace('severe ', 'severe_').replace('(severe) ', 'severe_').replace('(severe)', 'severe_')
#                 .replace('inactive ', 'inactive_')
#                 .replace('healed colitis', 'healed-colitis')
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
                .replace('inflamed', 'inflammation')
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
                .replace('ulcerated mucosa', 'ulcerated-mucosa')
                .replace('very inflam', 'very_inflam')
                .replace('friable', 'friability')
                .replace('left side', 'left-side')
                .replace('leftside', 'left-side')
                .replace('chonic', 'chronic')
                .replace('active', 'actv')
#                 .replace('procto-sigmoiditis', 'proctosigmoiditis')
                .replace('extensive colitis', 'extensive colitis (pancolitis)')
#                 .replace('mild friability', 'mild_friability')
#                 .replace('mayo ', 'mayo-')
#                 .replace('chronic active', 'chronic-active')
#                 .replace('inactive chronic', 'inactive-chronic')
# # #                 .replace('active colitis', 'active-colitis')
#                 .replace('inactive colitis', 'inactive-colitis')
               )
    
        if 'inflam' in line:
            inflam_pointer = pointer
        
        # term 'inflam' should be atmost 2 sentences prior to be carried forward to current sentence
        if (pointer-inflam_pointer<=2) and line.startswith(('this was')):
            line = (line
                    .replace('this was mild ', 'inflammation was mild ')
                    .replace('this was moderate ', 'inflammation was moderate ')
                    .replace('this was severe ', 'inflammation was severe ')
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
        
        mild_infl = bool(re.search(r'\b(?:mild\W+(?:\w+\W+){0,4}?inflammation|inflammation\W+(?:\w+\W+){0,3}?mild)\b', line))
        min_infl = bool(re.search(r'\b(?:minimal|minor\W+(?:\w+\W+){0,2}?inflammation|inflammation\W+(?:\w+\W+){0,2}?minimal|minor)\b', line))
        mod_infl = bool(re.search(r'\b(?:moderate\W+(?:\w+\W+){0,3}?inflammation|inflammation\W+(?:\w+\W+){0,3}?moderate)\b', line))
        sev_infl = bool(re.search(r'\b(?:severe\W+(?:\w+\W+){0,3}?inflammation|inflammation\W+(?:\w+\W+){0,3}?severe)\b', line))
        
        adhr_blood_1 = bool(re.search(r'\b(?:adherent\W+(?:\w+\W+){0,1}?blood)\b', line))
        adhr_clot_1 = bool(re.search(r'\b(?:adherent\W+(?:\w+\W+){0,1}?(clot\w*))\b', line))
        adhr_blood = (adhr_blood_1 or adhr_clot_1)
        
        left_colitis_1 = bool(re.search(r'\b(?:left\W+(?:\w+\W+){0,4}?(\w*colitis)|(\w*colitis)\W+(?:\w+\W+){0,3}?left)\b', line))
        left_colitis_2 = bool(re.search(r'\b(?:left-sided\W+(?:\w+\W+){0,4}?(\w*colitis)|(\w*colitis)\W+(?:\w+\W+){0,3}?left-sided)\b', line))
        left_colitis_3 = bool(re.search(r'\b(?:left-side\W+(?:\w+\W+){0,4}?(\w*colitis)|(\w*colitis)\W+(?:\w+\W+){0,3}?left-side)\b', line))
        left_colitis = (left_colitis_1|left_colitis_2|left_colitis_3)
        
        mild_col = bool(re.search(r'\b(?:mild\W+(?:\w+\W+){0,4}?(\w*colitis)|(\w*colitis)\W+(?:\w+\W+){0,3}?mild)\b', line))
        mod_col = bool(re.search(r'\b(?:moderate\W+(?:\w+\W+){0,3}?(\w*colitis)|(\w*colitis)\W+(?:\w+\W+){0,3}?moderate)\b', line))
        sev_col = bool(re.search(r'\b(?:severe\W+(?:\w+\W+){0,3}?(\w*colitis)|(\w*colitis)\W+(?:\w+\W+){0,3}?severe)\b', line))
        
        pan_col_1 = bool(re.search(r'\b(?:(\w*pancolonic)\W+(?:\w+\W+){0,2}?(\w*colitis))\b', line))
        pan_col_2 = bool(re.search(r'\b(?:pan\W+(?:\w+\W+){0,2}?(\w*colitis))\b', line))
        pan_col = (pan_col_1 or pan_col_2)
        
        act_col = bool(re.search(r'\b(?:actv\W+(?:\w+\W+){0,2}?(\w*colitis)|(\w*colitis)\W+(?:\w+\W+){0,2}?actv)\b', line))
        inact_col = bool(re.search(r'\b(?:inactv\W+(?:\w+\W+){0,2}?(\w*colitis)|(\w*colitis)\W+(?:\w+\W+){0,2}?inactv)\b', line))
                
        prcsg_col = bool(re.search(r'\b(?:proctosigmoid\W+(?:\w+\W+){0,2}?(\w*colitis)|(\w*colitis)\W+(?:\w+\W+){0,2}?proctosigmoid)\b', line))
        
        ulcproc_col = bool(re.search(r'\b(?:(ulcerative\w*)\W+(?:\w+\W+){0,3}?(\w*proctitis)|(\w*proctitis)\W+(?:\w+\W+){0,2}?(ulcerative\w*))\b', line))

        act_ile = bool(re.search(r'\b(?:actv\W+(?:\w+\W+){0,2}?(\w*ileitis)|(\w*ileitis)\W+(?:\w+\W+){0,2}?actv)\b', line))
        act_ent = bool(re.search(r'\b(?:actv\W+(?:\w+\W+){0,2}?(\w*enteritis)|(\w*enteritis)\W+(?:\w+\W+){0,2}?actv)\b', line))

        noncas_gran_1 = bool(re.search(r'\b(?:noncaseating\W+(?:\w+\W+){0,2}?granuloma|granuloma\W+(?:\w+\W+){0,2}?noncaseating)\b', line))
        noncas_gran_2 = bool(re.search(r'\b(?:noncaseating\W+(?:\w+\W+){0,2}?granulomas|granulomas\W+(?:\w+\W+){0,2}?noncaseating)\b', line)) 
        noncas_gran = (noncas_gran_1 or noncas_gran_2)
        
        nonnec_gran_1 = bool(re.search(r'\b(?:nonnecrotizing\W+(?:\w+\W+){0,2}?granuloma|granuloma\W+(?:\w+\W+){0,2}?nonnecrotizing)\b', line))
        nonnec_gran_2 = bool(re.search(r'\b(?:nonnecrotizing\W+(?:\w+\W+){0,2}?granulomas|granulomas\W+(?:\w+\W+){0,2}?nonnecrotizing)\b', line)) 
        nonnec_gran = (nonnec_gran_1 or nonnec_gran_2)
        
        lamprop_inf = bool(re.search(r'\b(?:lamina-propria\W+(?:\w+\W+){0,2}?inflammation|inflammation\W+(?:\w+\W+){0,3}?lamina-propria)\b', line))
        
        sup_ulcer = bool(re.search(r'\b(?:superficial\W+(?:\w+\W+){0,4}?(ulcer\w*)|(ulcer\w*)\W+(?:\w+\W+){0,2}?superficial)\b', line))
        
        shal_ulcer = bool(re.search(r'\b(?:shallow\W+(?:\w+\W+){0,4}?(ulcer\w*)|(ulcer\w*)\W+(?:\w+\W+){0,2}?shallow)\b', line))
        
        aph_ulcer = bool(re.search(r'\b(?:aphthous\W+(?:\w+\W+){0,3}?(ulcer\w*)|(ulcer\w*)\W+(?:\w+\W+){0,1}?aphthous)\b', line))
        
        sml_ulcer = bool(re.search(r'\b(?:small\W+(?:\w+\W+){0,3}?(ulcer\w*))\b', line))
        
        lrg_ulcer = bool(re.search(r'\b(?:large\W+(?:\w+\W+){0,4}?(ulcer\w*)|(ulcer\w*)\W+(?:\w+\W+){0,1}?large)\b', line))
        
        dp_ulcer = bool(re.search(r'\b(?:deep\W+(?:\w+\W+){0,5}?(ulcer\w*)|(ulcer\w*)\W+(?:\w+\W+){0,2}?deep)\b', line))
        
        mild_friab = bool(re.search(r'\b(?:mild\W+(?:\w+\W+){0,1}?friability|friability\W+(?:\w+\W+){0,1}?mild)\b', line))
        
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
            if 'inflammation' in e_text and min_infl and not ('minimal' in e_text or 'minor' in e_text):
                e_text = e_text + ' (minimal-inflammation)'
            if 'inflammation' in e_text and mod_infl and not 'moderate' in e_text:
                e_text = e_text + ' (moderate-inflammation)'
            if 'inflammation' in e_text and sev_infl and not 'severe' in e_text:
                e_text = e_text + ' (severe-inflammation)'
             
            if 'friability' in e_text and mild_friab and not 'mild' in e_text:
                e_text = e_text + ' (mild)'
                
            if 'adherent' in e_text and adhr_blood and (not ('clot' in e_text or 'blood' in e_text)):
                e_text = e_text + ' (adherent-blood)'
                
            if 'colitis' in e_text and pan_col and not 'pancolitis' in e_text:
                 e_text = e_text + ' (pancolitis)'
                    
#             if 'colitis' in e_text and left_colitis:
#                 print('e_text')
                    
            if 'colitis' in e_text and left_colitis and ('surveillance' in line or not (('biops' in line or ' bx ' in line) and 'taken' in line)):
                e_text = e_text + ' (left-sided colitis)'
                
            if 'ileitis' in e_text and act_ile and not 'actv' in e_text:
                 e_text = e_text + ' (actv)'
                    
            if 'enteritis' in e_text and act_ent and not 'actv' in e_text:
                 e_text = e_text + ' (actv)'
                    
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
                
            if 'ulcerative' in e_text and not 'proctitis' in e_text and ulcproc_col:
                e_text = e_text + ' (proctitis)'
                
    
                
            if 'colitis' in e_text:
                if mild_col and 'mild' not in e_text:
                    e_text = e_text + ' (mild)'
                if mod_col and 'moderate' not in e_text:
                    e_text = e_text + ' (moderate)'
                if sev_col and 'severe' not in e_text:
                    e_text = e_text + ' (severe)'
                if act_col and 'actv' not in e_text:
                    e_text = e_text + ' (actv)'
                if inact_col and 'inactv' not in e_text:
                    e_text = e_text + ' (inactv)'
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
                    
#         print(line, end='\n')
        
        pointer+=1
    
    if mayo_bool==True:
        entity_result = entity_result + 'colitis mayo-' + str(mayo_score) + ' ' + str(True) + '\n'
    
    entity_result = entity_result.replace('actv', 'active')
        
    return entity_result
