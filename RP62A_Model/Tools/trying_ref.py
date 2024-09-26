import BuildLP
import Analysis

def solLP(m):
    lp=BuildLP.BoundedUptakeLP(m)
    lp.Solve()
    sol=lp.GetPrimSol()
    '''
    for r in sol:
        print(f'{r}\t{sol[r]}')
    '''
    print(lp.GetObjVal())
    return sol

def matching(m):
    sol = solLP(m)
    pres_reac = ['GLU6PDEHYDROG-RXN','RXN-12878','VALINE-PYRUVATE-AMINOTRANSFER-RXN','6PGLUCONOLACT-RXN',
     'MALIC-NADP-RXN',  'RXN-14207','FORMATE_bp_tx', 'RXN-14192']
    abs_reac = ['DGDPKIN-RXN','BRANCHED-CHAINAMINOTRANSFERVAL-RXN','DADPKIN-RXN','2-OXOGLUTARATE-SYNTHASE-RXN',
     'GARTRANSFORMYL2-RXN', 'GLUCOSE-1-DEHYDROGENASE-RXN-(NADP)']
    p = []
    for r in pres_reac:
        if r in sol:
            p.append(r)
    a = []
    for r in abs_reac:
        if r not in sol:
            a.append(r)
    return p+a

def find_sol(m):
    rec = matching(m)
    if rec == ['RXN-14207', 'RXN-14192','DGDPKIN-RXN','DADPKIN-RXN']:
        print('I')
    elif rec == ['RXN-12878', 'MALIC-NADP-RXN', '2-OXOGLUTARATE-SYNTHASE-RXN', 'GLUCOSE-1-DEHYDROGENASE-RXN-(NADP)']:
        print('II')
    elif rec == ['GLU6PDEHYDROG-RXN','6PGLUCONOLACT-RXN','RXN-14207', 'RXN-14192', 'DGDPKIN-RXN', 'DADPKIN-RXN', 'GLUCOSE-1-DEHYDROGENASE-RXN-(NADP)']:
        print('III')
    elif rec == []:
        print('IV')
    elif rec == ['GLU6PDEHYDROG-RXN','VALINE-PYRUVATE-AMINOTRANSFER-RXN','6PGLUCONOLACT-RXN','BRANCHED-CHAINAMINOTRANSFERVAL-RXN', 'GLUCOSE-1-DEHYDROGENASE-RXN-(NADP)']:
        print('V')
    elif rec == ['RXN-12878', 'MALIC-NADP-RXN','FORMATE_bp_tx', '2-OXOGLUTARATE-SYNTHASE-RXN', 'GARTRANSFORMYL2-RXN', 'GLUCOSE-1-DEHYDROGENASE-RXN-(NADP)']:
        print('VI')
    elif rec == ['FORMATE_bp_tx', 'GARTRANSFORMYL2-RXN', ]:
        print('VII')
    

    
    
        
