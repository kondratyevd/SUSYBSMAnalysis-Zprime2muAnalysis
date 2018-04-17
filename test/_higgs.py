#!/usr/bin/env python


miniAOD = True
Electrons = False

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD

process.source.fileNames =[#'file:./pat.root'
'/store/mc/RunIISummer16MiniAODv2/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/04BC77C2-4BE0-E611-88B7-0090FAA57CC4.root'
#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/00312D7A-FEBD-E611-A713-002590DB923E.r\oot'                                                                                                                                                                                                              
#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/00355459-F4BD-E611-B5E8-D4AE526A11F3.r\oot'                                                                                                                                                                                                              
#'/store/mc/RunIISummer16MiniAODv2/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/36967CD0-3CC1-E611-A615-D8D385FF1996.root'
]			   

process.maxEvents.input = 10000
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2_v1' # for MC
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # default 1000

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match, trigger_paths, prescaled_trigger_paths, overall_prescale, offline_pt_threshold, prescaled_offline_pt_threshold
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
electrons_miniAOD(process)

import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionNew_cff as OurSelectionNew
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionDec2012_cff as OurSelectionDec2012
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2016_cff as OurSelection2016

dils = [('MuonsPlusMuonsMinus',          '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
	('MuonsPlusMuonsPlus',           '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+',         'daughter(0).pdgId() + daughter(1).pdgId() == -26'),
	('MuonsMinusMuonsMinus',         '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 26'),
	('MuonsSameSign',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
	('MuonsAllSigns',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
	]


cut_name = 'Simple'
Selection = OurSelectionDec2012


	# Keep track of modules to put in the path for this set of cuts.
path_list = []
path_list.append(process.egmGsfElectronIDSequence)
	    
leptons_name = cut_name + 'Leptons'
#    if cut_name == 'Simple':
muon_cuts = ''
leptons = process.leptonsMini.clone(muon_cuts = muon_cuts)
    
setattr(process, leptons_name, leptons)
path_list.append(leptons)

    # Make all the combinations of dileptons we defined above.
for dil_name, dil_decay, dil_cut in dils:

    name = cut_name + dil_name
    allname = 'all' + name

    alldil = Selection.allDimuons.clone(decay = dil_decay % locals(), cut = dil_cut)
    if 'AllSigns' in dil_name:
        alldil.checkCharge = cms.bool(False)

    dil = Selection.dimuons.clone(src = cms.InputTag(allname))

        #if cut_name == 'Simple':
    alldil.electron_cut_mask = cms.uint32(0)
        #alldil.loose_cut = 'isGlobalMuon && pt > 20.'#to be changed for first runs
    alldil.loose_cut = 'isGlobalMuon && pt > 20.'
    alldil.tight_cut = ''
    dil.max_candidates = 100
    dil.sort_by_pt = True
    dil.do_remove_overlap = False
    if hasattr(dil, 'back_to_back_cos_angle_min'):
        delattr(dil, 'back_to_back_cos_angle_min')
    if hasattr(dil, 'vertex_chi2_max'):
        delattr(dil, 'vertex_chi2_max')
    if hasattr(dil, 'dpt_over_pt_max'):
        delattr(dil, 'dpt_over_pt_max')
       
        # Add all these modules to the process and the path list.
    setattr(process, allname, alldil)
    setattr(process, name, dil)
    path_list.append(alldil * dil)


       #define the list of MC samples to be read here. be careful that if WWinclusive or tautau sample are not commented it will apply the filters when running locally.

    samples = [
            ('dyInclusive50','/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')	
        ]

    # Finally, make the path for this set of cuts.
    pathname = 'path' + cut_name
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DileptonPreselector_cfi')
    pobj = process.dileptonPreseletor * reduce(lambda x,y: x*y, path_list)
 
    path = cms.Path(pobj)
    setattr(process, pathname, path)


def ntuplify(process, fill_gen_info=True):
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
    obj = process.prunedMCLeptons
    obj.src = cms.InputTag('prunedGenParticles')

    process.HiggsToMuMu = cms.EDAnalyzer('HiggsToMuMu',
                                           dimu_src = cms.InputTag('SimpleMuonsAllSigns'),
						met_src = cms.InputTag("slimmedMETs"),
						jet_src = cms.InputTag("slimmedJets"),
                                           beamspot_src = cms.InputTag('offlineBeamSpot'),
                                           vertices_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
								TriggerResults_src = cms.InputTag('TriggerResults', '', 'PAT'),	#mc
# 								TriggerResults_src = cms.InputTag('TriggerResults', '', 'RECO'),	#data
                                           genEventInfo = cms.untracked.InputTag('generator'),
                                           metFilter = cms.VInputTag( cms.InputTag("Flag_HBHENoiseFilter"), cms.InputTag("Flag_HBHENoiseIsoFilter"), 
                                                                      cms.InputTag("Flag_EcalDeadCellTriggerPrimitiveFilter"), cms.InputTag("Flag_eeBadScFilter"), 
                                                                      cms.InputTag("Flag_globalTightHalo2016Filter"))
                                           )        
    if hasattr(process, 'pathSimple'):
	if miniAOD and fill_gen_info:
        	process.pathSimple *=obj * process.HiggsToMuMu 
#
ntuplify(process) #to have ntuples also running in interactive way

def printify(process):
    process.MessageLogger.categories.append('PrintEvent')

    process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')
    process.triggerSummaryAnalyzerAOD.inputTag = cms.InputTag('hltTriggerSummaryAOD','','HLT')
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.triggerSummaryAnalyzerAOD

    process.PrintOriginalMuons = cms.EDAnalyzer('PrintEvent', muon_src = cms.InputTag('cleanPatMuonsTriggerMatch'), trigger_results_src = cms.InputTag('TriggerResults','','HLT'))
    process.pathSimple *= process.PrintOriginalMuons

    pe = process.PrintEventSimple = cms.EDAnalyzer('PrintEvent', dilepton_src = cms.InputTag('SimpleMuonsPlusMuonsMinus'))
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.PrintEventSimple

    #- 2011-2012 selection (Nlayers > 8)
    #process.PrintEventOurNew = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsPlusMuonsMinus'))
    #process.PrintEventOurNewSS = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsSameSign'))
    #process.PrintEventOurNewEmu = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsElectronsOppSign'))
    #process.pathOurNew *= process.PrintEventOurNew * process.PrintEventOurNewSS * process.PrintEventOurNewEmu

    #- December 2012 selection (Nlayers > 5, re-tuned TuneP, dpT/pT < 0.3)
    if hasattr(process, 'pathOur2012'):
        process.PrintEventOur2012    = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsPlusMuonsMinus'))
        process.PrintEventOur2012SS  = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsSameSign'))
        process.PrintEventOur2012Emu = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsElectronsOppSign'))
        process.pathOur2012 *= process.PrintEventOur2012 * process.PrintEventOur2012SS * process.PrintEventOur2012Emu

def check_prescale(process, trigger_paths, hlt_process_name='HLT'):
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.CheckPrescale_cfi')
    process.CheckPrescale.trigger_paths = cms.vstring(*trigger_paths)
    process.pCheckPrescale = cms.Path(process.CheckPrescale)

def for_data(process):
    process.GlobalTag.globaltag ='80X_dataRun2_2016SeptRepro_v6' #RunB to RunG   #change line 52
#     process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v14' #RunH              #change line 52
    ntuplify(process)
    #check_prescale(process, trigger_paths) ####### Now it seams that there are no prescaled path ########

def for_mc(process, hlt_process_name, fill_gen_info):
    ntuplify(process, fill_gen_info)
    switch_hlt_process_name(process, hlt_process_name) # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)

def get_dataset(run):
    #JMTBAD common with dataset_details in submit below, make a DataSamples.py?
    run = int(run)
    if 190450 <= run <= 191284:
        return '/SingleMu/tucker-datamc_SingleMuRun2012A_Prompt_190450_191284_20120418134612-57b19813ab8f2ab142c4566dc6738156/USER'
    else:
        raise ValueError('dunno how to do run %i' % run)

if 'int_data' in sys.argv:
    for_data(process)
    printify(process)
    
if 'int_mc' in sys.argv:
    for_mc(process, 'HLT', False)
    printify(process)
    
if 'gogo' in sys.argv:
    for_data(process)
    printify(process)
    
    n = sys.argv.index('gogo')
    run, lumi, event = sys.argv[n+1], sys.argv[n+2], sys.argv[n+3]
    print run, lumi, event
    run = int(run)
    lumi = int(lumi)
    event = int(event)
    filename = [x for x in sys.argv if x.endswith('.root')]
    if filename:
        filename = filename[0]
    else:
        dataset = get_dataset(run)
        print dataset
        output = os.popen('dbs search --url https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet --query="find file where dataset=%s and run=%s and lumi=%s"' % (dataset, run, lumi)).read()
        print repr(output)
        filename = [x for x in output.split('\n') if x.endswith('.root')][0]
    print filename
    process.source.fileNames = [filename]
    from SUSYBSMAnalysis.Zprime2muAnalysis.cmsswtools import set_events_to_process
    set_events_to_process(process, [(run, event)])

#f = file('outfile_histos1', 'w')
#f.write(process.dumpPython())
#f.close()

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'ana_datamc_%(name)s'
config.General.workArea = 'crab'
#config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_crab.py'   
#config.JobType.priority = 1
config.Data.inputDataset =  '%(ana_dataset)s'
config.Data.inputDBS = 'global'
job_control
config.Data.publication = False
config.Data.outputDatasetTag = 'ana_datamc_%(name)s'
config.Data.outLFNDirBase = '/store/user/dkondrat'
#config.Data.ignoreLocality = True 
#config.Site.whitelist = ["T2_IT_Bari"]
config.Site.storageSite = 'T2_RU_JINR'
'''
    
    just_testing = 'testing' in sys.argv
        
    # Run on data.
    if 'no_data' not in sys.argv:
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import *

        dataset_details = [

						('SingleMuonRun2016B-ReReco-v3', '/SingleMuon/Run2016B-23Sep2016-v3/MINIAOD'),
						('SingleMuonRun2016C-ReReco-v1', '/SingleMuon/Run2016C-23Sep2016-v1/MINIAOD'),
						('SingleMuonRun2016D-ReReco-v1', '/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD'),
						('SingleMuonRun2016E-ReReco-v1', '/SingleMuon/Run2016E-23Sep2016-v1/MINIAOD'),
						('SingleMuonRun2016F-ReReco-v1', '/SingleMuon/Run2016F-23Sep2016-v1/MINIAOD'),
						('SingleMuonRun2016G-ReReco-v1', '/SingleMuon/Run2016G-23Sep2016-v1/MINIAOD'),
# 							('SingleMuonRun2016H-PromptReco-v3', '/SingleMuon/Run2016H-PromptReco-v3/MINIAOD'), #change global tag: 26 and 408
# 							('SingleMuonRun2016H-PromptReco-v2', '/SingleMuon/Run2016H-PromptReco-v2/MINIAOD'),  ##change global tag: 26 and 408

				# 							('SingleMuonRun2016G-PromptReco-v1',  '/SingleMuon/Run2016G-PromptReco-v1/MINIAOD')

            ]

        lumi_lists = [
		# 'NoLumiMask'
		#           'DCSOnly',
		#            'Run2012PlusDCSOnlyMuonsOnly',
		'Run2016MuonsOnly',
		# 'Run2015',
		]

        jobs = []
        for lumi_name in lumi_lists:
            ll = eval(lumi_name + '_ll') if lumi_name != 'NoLumiMask' else None
            for dd in dataset_details:
                jobs.append(dd + (lumi_name, ll))
                
        for dataset_name, ana_dataset, lumi_name, lumi_list in jobs:
            json_fn = 'tmp.json'
            lumi_list.writeJSON(json_fn)
            lumi_mask = json_fn

            name = '%s_%s' % (lumi_name, dataset_name)
            print name

            new_py = open('histos.py').read()
            new_py += "\nfor_data(process)\n"
            open('histos_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()

            job_control = '''
config.Data.splitting = 'LumiBased'
#config.Data.runRange = '256843-257490'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 200
config.Data.lumiMask = '%(lumi_mask)s' #######
''' % locals()

            new_crab_cfg = new_crab_cfg.replace('job_control', job_control)
            open('crabConfig.py', 'wt').write(new_crab_cfg)

            if not just_testing:
                os.system('crab submit -c crabConfig.py')
            else:
                cmd = 'diff histos.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crab.py'
                print cmd
                os.system(cmd)

        if not just_testing:
            #os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')
            os.system('rm crabConfig.py histos_crab.py histos_crab.pyc tmp.json')

    if 'no_mc' not in sys.argv:
        # Set crab_cfg for MC.
        crab_cfg = crab_cfg.replace('job_control','''
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 1000000
    ''')

       
        for name, ana_dataset in samples:
            print name

            new_py = open('histos.py').read()
            
            open('histos_crab.py', 'wt').write(new_py)

            open('crabConfig.py', 'wt').write(crab_cfg % locals())
            if not just_testing:
                os.system('crab submit -c crabConfig.py')
            else:
                cmd = 'diff histos.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crabConfig.py'
                print cmd
                os.system(cmd)

        if not just_testing:
		os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')

