from importlib import reload
import SepiBiomass
reload(SepiBiomass)
import Media
reload(Media)

media = Media.media

def MaintCost(mu=1,GAM=SepiBiomass.GAM, NGAM=SepiBiomass.NGAM ):
    """ Growth associated and non-growth associated maintenance cost """
    return mu*GAM+NGAM

  
def BuildLP(m):
    """ build lp with objective as minimisation of total flux """
    
    lp = m.GetLP() 
    lp.SetObjective(m.sm.cnames)   # minimise total flux

    return lp


def BiomassLP(m, mu=SepiBiomass.DefaultMu, GAM=SepiBiomass.GAM, NGAM=SepiBiomass.NGAM):
    """build lp with constraint on biomass and cell maintenance and objective as minimisation of total flux"""

    fd = SepiBiomass.FluxDic(m, mu)
    fd["ATPASE-RXN"] = MaintCost(mu,GAM, NGAM)

    lp = BuildLP(m)
    lp.SetFixedFlux(fd)

    return lp


def ClosedLP(m):

    """ lp object contraining flux through tx reacs to zero """

    lp = BuildLP(m)
    
    for tx in filter(lambda s: "_tx" in s, m.sm.cnames): # all tx in model
        lp.SetFixedFlux({tx:0.0}) # set flux to 0.0

    return lp

def BlockAAUptakeLP(m,Biomass=False, mu=SepiBiomass.DefaultMu):
	"""Pre:True
	Post: lp object, constrained  AA media transporter fluxes to 0  """
	if Biomass == False:
		lp = BuildLP(m)
	else:
		lp = BiomassLP(m,mu)
	for tx in filter(lambda s:"_AA_mm_tx" in s, m.sm.cnames):
		lp.SetFixedFlux({tx:0.0})
	return lp

def BoundedUptakeLP(m,mu=SepiBiomass.DefaultMu,lp=None):
	'''Pre: True
	Post: return lp with upper bound on substrate uptake''' 
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
