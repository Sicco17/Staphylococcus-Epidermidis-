from ScrumPy.Util import Set
import re

def ETCReacs(db,Type='Electron-Transfer'):
	'''pre:True
	post:list of electron transfer reactions in database'''
	rv = []
	for p in db.dbs['Pathway'].keys():
		if 'TYPES' in list(db[p].keys()) and Type in db[p]['TYPES']:
			for r in db[p]['REACTION-LIST']:
				if r not in rv:
					rv.append(r)
	return rv


def FindtRNA(m,pattern='RNA'):
	'''pre:True
	post:list of met associated with RNA in database'''
	rv = []
	for met in m.sm.rnames:
		z = re.search(pattern, met)
		if z != None:
			rv.append(met)
	return rv


def ReacToPW(m,db):
	rv = {}
	for r in m.sm.cnames:
		if r in db.dbs['REACTION'].keys() and 'IN-PATHWAY' in db.dbs['REACTION'][r].keys():
			rv[r] = db.dbs['REACTION'][r]['IN-PATHWAY']
	return rv

def FindSubReac(m,db):
	'''pre:True
	post:dictionary of subreaction with main reaction as value'''
	rv = {}
	r2p = ReacToPW(m,db)
	reactions = m.sm.cnames
	for r in r2p.keys():
		common = Set.Intersect(r2p[r],reactions)
		if common != []:  
			rv[r] = common
	return rv 

def FindMainReac(db):
	'''pre:True
	post:dictionary of main reactions with sub reaction as value'''
	rv = {}
	for r in db.dbs['REACTION'].keys():
		if 'REACTION-LIST' in db.dbs['REACTION'][r].keys():
			rv[r] = db.dbs['REACTION'][r]['REACTION-LIST']
	return rv

def DelSubReac(m,db):
	'''pre:True
	post:list of subreaction that causes enzyme inconsistent set'''
	main = FindMainReac(db)
	ess = m.EnzSubsets()
	esss = map(lambda s: ess[s].keys(), ess.Inconsistants)
	subreacs = map(lambda s: main[s], main)
	return Set.Intersect(list(set().union(*esss)),list(set().union(*subreacs)))


def GeneToTx(db):
	rv = {}
	for p in db.dbs['PROTEIN'].keys():
		if 'COMMON-NAME' in db.dbs['PROTEIN'][p].keys():
			if re.search('port',db.dbs['PROTEIN'][p]['COMMON-NAME'][0]):
				#print(p)
				if 'GENE' in db.dbs['PROTEIN'][p].keys():
					gene = db.dbs['PROTEIN'][p]['GENE'][0]
					#print (gene)
					if 'COMMON-NAME' in db.dbs['GENE'][gene].keys():
						rv[gene] = db.dbs['GENE'][gene]['COMMON-NAME'][0]
					else:
						rv[gene] = gene
	return rv


def Paths(r2p):
	reacs = r2p.keys()
	path = []
	for r in reacs:
		for p in r2p[r]:
			if p not in path:
				path.append(p)
	return path

def PathToReac(r2p):
	rv = {}
	path = Paths(r2p)
	for p in path:
		reac = []
		for r in r2p.keys():
			if p in r2p[r]:
				reac.append(r)
		rv[p] = reac
	return rv
	

def WhyDead(m,db,dead=[]):
	rv = {}
	if dead == []:
		dead = m.DeadReactions()
	orp = m.OrphanMets()
	for r in dead:
		if r in db.dbs['REACTION'].keys() and 'IN-PATHWAY' in db.dbs['REACTION'][r].keys():
			p = Set.Complement(db.dbs['REACTION'][r]['IN-PATHWAY'],db.dbs['REACTION'].keys())[0]
			print (p,r)
			reacs = []
			orph = []
			not_model = []
			for r in db.dbs['Pathway'][p]['REACTION-LIST']:
				if r in m.sm.cnames:
					for met in m.sm.InvolvedWith(r).keys():
						if met in orp:
							reacs.append(r)
							print ("Orphan",met)
							orph.append(met)
						
				else:
					print ('not in model',r)
					not_model.append(r)
			
			print ('*******************************')
			rv[p] = (reacs,orph,not_model)
	return rv


