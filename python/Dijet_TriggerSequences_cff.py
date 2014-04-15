import FWCore.ParameterSet.Config as cms

###
### trigger paths of interest
### Trigger evolution in 2011 cf. http://fwyzard.web.cern.ch/fwyzard/hlt/summary
###
def getTriggerPaths(version=2012) :

    #################
    # 2012 triggers #
    #################
    if(version==2012) :
	mcTrigs = ['HLT_Mu17_v',
		   'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v']
        
	Muon = ['HLT_Mu17_v']

	Electron = ['HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v']

    
    return Muon, Electron, mcTrigs

