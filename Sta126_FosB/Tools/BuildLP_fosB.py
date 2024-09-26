from importlib import reload
import MinBiomass_fosB
reload(MinBiomass_fosB)
import Media_fosB
reload(Media_fosB)

media = Media_fosB.media

def MaintCost(mu=1,GAM=MinBiomass_fosB.GAM, NGAM=MinBiomass_fosB.NGAM ):
    """ Growth associated and non-growth associated maintenance cost """
    return mu*GAM+NGAM

  
def BuildLP(m): #basic LP problem
    lp = m.GetLP() 
    lp.SetObjective(m.sm.cnames)   # minimise total flux
    return lp


def BiomassLP(m, mu=MinBiomass_fosB.DefaultMu, GAM=MinBiomass_fosB.GAM, NGAM=MinBiomass_fosB.NGAM): #lp with constrain biomass and cell maintenance
    fd = MinBiomass_fosB.FluxDic(m, mu)
    fd["ATPASE-RXN"] = MaintCost(mu,GAM, NGAM)
    lp = BuildLP(m)
    lp.SetFixedFlux(fd)
    return lp


def ClosedLP(m): #lp with constraints influx to 0
    lp = BuildLP(m)
    for tx in filter(lambda s: "_tx" in s, m.sm.cnames): # all tx in model
        lp.SetFixedFlux({tx:0.0}) # set flux to 0.0
    return lp

def BlockAAUptakeLP(m,Biomass=False, mu=MinBiomass_fosB.DefaultMu): #lp with constaints influx of amino acids to 0
	if Biomass == False:
		lp = BuildLP(m)
	else:
		lp = BiomassLP(m,mu)
	for tx in filter(lambda s:"_AA_mm_tx" in s, m.sm.cnames):
		lp.SetFixedFlux({tx:0.0})
	return lp

def BoundedUptakeLP(m,mu=MinBiomass_fosB.DefaultMu,lp=None): #lp with upper bound on substrate uptake
	if lp == None:
		lp = BiomassLP(m,mu)
	for Compound in media.keys():
		txs = [tx for tx in m.smx.InvolvedWith('x_'+Compound) if not "_bm_" in tx and not "_bp_" in tx]
		if len(txs)==1.0:
			lp.SetFluxBounds({txs[0]:(0.0,media[Compound])})
		elif len(txs)>1.0:
			lp.SetSumFluxConstraint(txs,media[Compound],"Total"+Compound)
			lp.SetRowBounds("Total"+Compound,0,media[Compound])
		else:
			print("no tx found for", Compound) 
	return lp
