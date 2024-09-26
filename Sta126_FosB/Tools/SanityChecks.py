#Corrections made for Python3.6

from ScrumPy.Bioinf import PyoCyc
from ScrumPy.Util import Set
import ModelCuration
#from ScrumPy.Structural import StoCons

path = '.'


def Updatedb(db,path=path,ExtraCompounds='ExtraCompounds.dat'):
    '''pre: True
    Post: Update db with mol formula of compounds in ExtraCompounds.dat'''

    if ExtraCompounds !=None:
        updatedb = PyoCyc.Compound.DB(path,ExtraCompounds)
        for met in updatedb.keys():
            db.dbs['Compound'][met] = updatedb[met]
    return db
   

def CheckImbals(m,db,path='.', ExtraCompounds='ExtraCompounds.dat',reacs=[]):
    '''Pre:True
    Post: return dictionary of imbalanced reactions'''

    rv = {}
    db = Updatedb(db,path,ExtraCompounds)

    if reacs == []:
        reacs = filter(lambda x: not "_tx" in x, m.smx.cnames)
    for reac in reacs:
        stod = m.smx.InvolvedWith(reac)
        imbal = db.dbs['Compound'].AtomImbal(stod)

        if len(imbal) >0:
                rv[reac] = imbal

    return rv   


def FindUnknownCompounds(imbaldic):
    """ pre: imbaldic = CheckImbals(..)
    post: list of Compound names in imbaldic that were prefixed as "Unknown compound"
    """

    rvd = {}

    for ib in imbaldic.values():

        names = ib.keys()
        for name in names:
            if "Unknown" in name:
                name = name.rsplit(" ",1)[1]
                rvd[name] =1
    return rvd.keys()



def  SeperateUnknowns(imbaldic):
    '''Pre:imbaldic = CheckImbals(..)
    Post: remove reactions involedwith "Unknown compound" from imbalanced dictionary'''
    rv = {}

    for reac in list(imbaldic.keys()):
        imbal = imbaldic[reac]
        names = "".join(imbal.keys())
        if "Unknown" in names:
            rv[reac] = imbal
            del imbaldic[reac]
    return rv
        

                

def FindAtoms(imbaldic, atom="C"):
    """ pre: imbaldic = CheckImbals(..)
    post: list of reactions with an atomic imbalance in atom
   	"""

    rv = {}
    for reac,ib in imbaldic.items():
            if ib.has_key(atom):
                    rv[reac] = ib
    return rv
 

def IdentifyEmpForm(imbaldic):
    """ pre: imbaldic = CheckImbals(..)
        post: return empirical formula for unknown compounds in imbalanced dictionary
    """
    UBDict = imbaldic.copy()
    rv = {}
    uks = list(filter(lambda s: "Unknown compound " in s, UBDict))
    if len(uks)==1:
        uk = uks[0]
        del UBDict[uk]
        uk = uk.replace("Unknown compound ","")
        for k in UBDict:
                rv[k] = int(abs(UBDict[k]))
        return uk, rv
    return None  
   

def EmpFormAsStr(name, empform):
    """ pre: Name of metabolite and empirical formula
        post: return empirical formula as in BioCyc entry
    """

    uid = "UNIQUE-ID - " + name
    rvl = [uid]
    form = "CHEMICAL-FORMULA - (atom val)"

    for atom in empform:
        atomval = form.replace("atom", atom).replace("val", str(empform[atom]))
        rvl.append(atomval)
    rvl.append("//")

    return "\n".join(rvl)

        
def AssignUnknown(imbal):
    """ pre: imbaldic = CheckImbals(..)
        post: prints empirical formula as in BioCyc entry
    """
    mets = []
    Extras = ''
    for r in imbal.keys():
        emp=IdentifyEmpForm(imbal[r])
        if emp != None and emp[0] not in mets:
            tt = EmpFormAsStr(emp[0], emp[1])
            mets.append(emp[0])
            Extras += tt + '\n'#why print one at a time when we can print in one big block? why not add it directly to the ExtraCoumpounds.dat file?
    return Extras

'''
new function to append all compounds to an ExtraCompounds file
'''
def AddUnk(imbal): #open the file and add the new compunds
    comp = AssignUnknown(imbal)
    manual_revision = input('want to see what has been added? [Y/N]')
    if manual_revision == 'Y' or manual_revision == 'y':
        print(comp)
    if len(comp)!=0:
        comp = '\n'+comp
        file = open('ExtraCompounds.dat','a')
        file.write(comp)
        file.close()
        print('new compounds added')
    else:
        print('no further addition of unknown compounds')
        
'''
new function to remove unknwon orphan metabolites and tRNA
'''
def del_orph_unk(model, imbalances):
    orph = model.OrphanMets()
    unk = FindUnknownCompounds(imbalances)
    rna = ModelCuration.FindtRNA(model)
    del_reactions = []
    for met in Set.Intersect(orph,unk)+rna:
        for r in model.smx.InvolvedWith(met).keys():
            if r not in del_reactions:
                del_reactions.append(r)
    print(f'{len(del_reactions)} reactions involving unknown orphan metabolites and tRNA were removed')
    print(f'{del_reactions} were removed')
    model.DelReactions(del_reactions)
    
    

#def StoCon(m):
#    '''Pre:
#    Post: Unconserved metabolite'''
#    smt = m.smx.Copy(tx=1)
#    nc = StoCons.UnconservedMets(smt)      
#    return nc   




def GetCharge(db,met):
    rv = {}
    charge = 0

    if "ATOM-CHARGES" in db.dbs['Compound'][met].keys():
        for c in range(0,len(db.dbs['Compound'][met]["ATOM-CHARGES"])):
            charge += int(db.dbs['Compound'][met]['ATOM-CHARGES'][c].split(' ')[-1].replace(')',''))
        ncharge = charge
        #print s,ncharge
    elif "SMILES" in db.dbs['Compound'][met].keys() and "ATOM-CHARGES" not in db.dbs['Compound'][met].keys():
        ncharge = 0
    else: 
        ncharge = 'Unknown charge ' + met

    return ncharge

def ChargeDic(db,StoDic):
    rv = {}
    for s in StoDic:
        rv[s] = (StoDic[s],GetCharge(db,s))
    return rv


def CheckChargeImbals(m,db,path=path,ExtraCompounds='ExtraCompounds.dat',reacs=[]):
    '''Pre:True
    Post: return dictionary of charge imbalanced reactions'''

    rv = {}
    db = Updatedb(db,path,ExtraCompounds)

    if reacs == []:
        reacs = filter(lambda x: not "_tx" in x, m.smx.cnames)
    for reac in reacs:
        dic = {}
        charge = 0.0
        stod = m.smx.InvolvedWith(reac)
        chargedic = ChargeDic(db,stod)
        for e in filter(lambda s:type(s[1])!=str,chargedic.values()):
            charge += int(e[0])*e[1] #coeff * charge
        dic['ImbalCharge'] = charge

        unknowns = filter(lambda s:type(s[1])==str,chargedic.values())

        for unk in unknowns:
            dic[unk[1]] =  unk[0]	

        if dic['ImbalCharge'] != 0.0 or unknowns != []:
            rv[reac] = dic

    return rv





