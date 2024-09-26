from ScrumPy.Util import Set
from importlib import reload
import BuildLP_fosB
reload(BuildLP_fosB)
import re
import MinBiomass_fosB
reload(MinBiomass_fosB)
import pandas as pd

def CheckUnitFlux(m,reac,neg=False,lp=None):
    """ pre: reac in m.sm.cnames
       post: the lp solution for a unit flux in reac or {} if none exists"""

    if lp ==None:
        lp = BuildLP_fosB.BuildLP(m)

    j = -1 if neg else 1
    lp.SetFixedFlux({reac:j})
    lp.Solve()  #automatically print optimal for each biomass-producing reaction
    lp.ClearFluxConstraint(reac)
    return lp.GetPrimSol()


def EnergyCheck(m,lp=None,check='ATPASE-RXN'):
    """pre: reac in m.sm.cnames
       post: the lp solution for a unit flux in reac in absence of external energy source"""
    if lp ==None:
        lp = BuildLP_fosB.BlockAAUptakeLP(m)
    if 'HYDROGEN-MOLECULE_tx' in m.sm.cnames:
        lp.SetFixedFlux({'HYDROGEN-MOLECULE_tx':0})
    sol = CheckUnitFlux(m,check,lp=lp)
    return sol

			
def CheckBMProduction(m,neg=True,lp=None,fd=None):
	"""pre: True
	post: the lp solution for production of each biomass components"""
	rv = {}
	if fd == None:
		bm = list(filter(lambda s:'_bm_tx' in s, m.sm.cnames))
	else:
		fd = FluxDic.fd
		bm = list(fd.keys())
	if lp == None:
        	lp = BuildLP_fosB.BuildLP(m)
	for reac in bm:
		sol = CheckUnitFlux(m,reac,neg,lp)
		if sol == {}:
			print(reac)
			rv[reac] = 'NF'
		else:
			rv[reac] = sol

	return rv
 