from ScrumPy.Bioinf import Biomass
from ScrumPy.Util import Set
from importlib import reload
import ExtraMolWts_fosB
reload(ExtraMolWts_fosB)
Biomass.MolWts.MolWts.update(ExtraMolWts_fosB.MolWts) # update Mol weights for 'undefined' comps in model (identified in ExtraMolWts)


DefaultMu = 0.432 # Observed growthrate on complete MHHW media
GAM = 60 #growth assotiated maintanance cost
NGAM=8 #non-growth associated maintanance cost

'''
keep all these values as they are for Staph. Epidermidis
'''
Prot = Biomass.Composition(
    0.43,{
        "L-ALPHA-ALANINE" : 0.044020536,
        "ARG" : 0.044936584,
        "ASN" : 0.053067709,
        "L-ASPARTATE" : 0.080208356,
        "CYS" : 0.005394561,
        "GLN" : 0.068050796,
        "GLT" : 0.067775302,
        "GLY" : 0.033141282,
        "HIS" : 0.026470545,
        "ILE" : 0.082115006,
        "LEU" : 0.086000908,
        "LYS" : 0.080585505,
        "MET" : 0.029034192,
        "PHE" : 0.052790901,
        "PRO" : 0.031106174,
        "SER" : 0.048301381,
        "THR" : 0.049528989,
        "TRP" : 0.010911318,
        "TYR" : 0.050087580,
        "VAL" : 0.056472366
    }
)


DNAnuc = Biomass.Composition(
    0.03,{
        "DATP" : 0.329314752,
        "DCTP" : 0.15283191,
        "DGTP" : 0.16604016,
        "TTP" : 0.323221964,
   }
)

RNAnuc = Biomass.Composition(
    0.12,{
        "ATP" : 0.254594912,
        "CTP" : 0.237646992,
        "GTP" : 0.257498896,
        "UTP" : 0.242936672
    }
)

cellmembrane = Biomass.Composition(
    0.07,{
        "DIACYLGLYCEROL" : 0.159977540, #DAG_bm_tx
        "L-1-PHOSPHATIDYL-GLYCEROL" : 0.481503385, #PG_bm_tx
        "CARDIOLIPIN" : 0.113131574, #CL_bm_tx
        "Type-I-LTA" : 0.06831944,  #one of the substrate UDP-GLC (CPD-12575),LTA_bm_tx
        "diacyl-3-O-glucl-1-6-gluc-sn-glycerol" : 0.088583010, #one of the substrate UDP-GLC (CPD-12575), GLC2-DAG_bm_tx
        "D-Glucosyl-12-diacyl-glycerols" : 0.012214416, #one of the substrate UDP-GLC, GLC-DAG_bm_tx
        "MENAQUINONE" : 0.076270631 #Menaquinone_bm_tx
    }
)
#Precursor UDP-N-ACETYL-D-GLUCOSAMINE
cellwall =  Biomass.Composition(
    0.24,{
        "CPD-12230" : 0.234311640,#PepGly_bm_tx
        "Rbo-P-Teichoic-aureus-peptidoglycan" : 0.765688359#WTA-PG_bm_tx
    }
)

sol = Biomass.Composition(
    0.11, {
        "glycogen" : 0.864808552,
        "ACETYL-COA" : 0.002672369, 
        "SUC-COA" : 0.000173428,
        "CO-A" : 0.003070192,
        "FAD" : 0.005270756,
        "NAD" : 0.095428126,
        "NADH" : 0.002200921,
        "NADP" : 0.00647552,
        "NADPH" : 0.019900124
    }
)

#Precursor UDP-N-ACETYL-D-GLUCOSAMINE
Biofilm = Biomass.Composition( # hexoses plus hexosamines
    1.0, {
        "PIA1" : 0.740000000, 
        "PIA2" : 0.200000000,
        "PIA3" : 0.060000000  
    }
)

WholeCell = Biomass.Composition(
    1.0, 
    {
        "Prot" : Prot,
        "DNAnuc" : DNAnuc,
        "RNAnuc" : RNAnuc,
        "cellwall" : cellwall,
        "cellmembrane" : cellmembrane,
        "sol" : sol
    }
)

All = Biomass.Composition(
    1.0,{
        "WholeCell": WholeCell,
        "Biofilm" :Biofilm
    }
)


StripComp = lambda s: s.rsplit("_",1)[0] #keep only the name of the metabolite

def WhatExports(m, compound): #find biomass producing reactions associated to a given compound
    Externs = m.md.GetExtMetNames() #external x_ metabolites from md = model description
    for reac in m.sm.InvolvedWith(compound).keys(): 
        if "bm_tx" in reac and len(Set.Intersect(Externs, m.smx.InvolvedWith(reac).keys())) != 0: 
            return reac #if the bm reaction involved with the compound is also involved with the external metabolite


def FluxDic(m,mu=DefaultMu, comp=All): # All = WholeCell + Biofilm (above), mu can be modified according to the growth conditions
    rv = {}
    for met in comp.GetLeaves(): # leaves are the end nodes in BM(final subcomponents)
        xp = WhatExports(m,met)
        if xp == None:
            print ("no tx for ", met)
        else:
            rv[xp] = -All.MolsOf(met,StripComp)*mu*1000  # set flux of exporters (- cause is exporting) so the corresponding number of 'x' moles of each comp in BM description for prod of 1 arbitrary unit of BM are exported (*1000 to work on mmols)                                                   
    return rv


