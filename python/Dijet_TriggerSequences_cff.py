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
        
	Muon = ['HLT_Mu8_v',
		'HLT_Mu12_v',
		'HLT_Mu17_v',
		'HLT_Mu24_',
		'HLT_IsoMu24_']

	Electron = ['HLT_Ele8_',
		    'HLT_Ele8_CaloIdL_CaloIsoVL_v',
		    'HLT_Ele17_CaloIdL_CaloIsoVL_v']

    
    return Muon, Electron

