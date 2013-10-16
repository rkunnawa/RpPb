// macro to calculate RpPb and Rcp

#include <iostream>
#include <stdio.h>

#include <TRandom.h>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <cstdlib>
#include <cmath>

#include "RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBayes.h"

#include "RooUnfold-1.1.1/src/RooUnfoldSvd.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"

#include "headers/utilities.h"
#include "headers/bayesianUnfold.h"
#include "headers/prior.h"




void Unfold_RpPb(){
	
#ifdef __CINT__
	gSystem->Load("RooUnfold-1.1.1/src/libRooUnfold");
#endif
	
	// start the timer
	TStopwatch timer;
	timer.Start();
	
	gStyle->SetErrorX(0.5);
	gStyle->SetPaintTextFormat("3.2f");
	gStyle->SetOptLogz(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadRightMargin(0.13);
	
	cout<<" ------------------- RpPb - measured and unfolded, Raghav 09/16/2013 --------------------"<<endl;
	cout<<" Program start"<<endl;
	/*
	int nBayesianITer = 4;
	char chmet1[100];
	
	if(method==1){
		sprintf(chmet1,"Bayes_unfo");
	}
	
	printf("Method : %s \n",chmet1);
	printf("AlgopPb: %s \n",algonamePP[algo]);
	*/
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	
	// load the required pp signal sample 5.02 GeV files as a cross check for pPb.
	/*
	const int nbins_pthat_pptrack = 8;
	Double_t boundaries_pthat_pptrack[nbins_pthat_pptrack+1];
	char *fileName_pthat_pptrack[nbins_pthat_pptrack+1];
	Double_t xsection_pptrack[nbins_pthat_pptrack+1]
	
	boundaries_pthat_pptrack[0] = 15;
	fileName_pthat_pptrack[0] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt15/HiForest_v77_v2_merged01/pt15_HP04_hiforest77_hiSignal.root";
	xsection_pptrack[0] = 5.335e-01;
	
	boundaries_pthat_pptrack[1] = 30;
	fileName_pthat_pptrack[1] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt30/HiForest_v77_v2_merged01/pt30_HP04_hiforest77_hiSignal.root";
	xsection_pptrack[1] = 3.378e-02;
	
	boundaries_pthat_pptrack[2] = 50;
	fileName_pthat_pptrack[2] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt50/HiForest_v77_v2_merged01/pt50_HP04_hiforest77_hiSignal.root";
	xsection_pptrack[2] = 3.778e-03;
	
	boundaries_pthat_pptrack[3] = 80;
	fileName_pthat_pptrack[3] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt80/HiForest_v77_v2_merged01/pt80_HP04_hiforest77_hiSignal.root";
	xsection_pptrack[3] = 4.412e-04;
	
	boundaries_pthat_pptrack[4] = 120;
	fileName_pthat_pptrack[4] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt120/HiForest_v77_v2_merged01/pt120_HP04_hiforest77_hiSignal.root";
	xsection_pptrack[4] = 6.147e-05;
	
	boundaries_pthat_pptrack[5] = 170;
	fileName_pthat_pptrack[5] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt170/HiForest_v77_v2_merged01/pt170_HP04_hiforest77_hiSignal.root";
	xsection_pptrack[5] = 1.018e-05;
	
	boundaries_pthat_pptrack[6] = 220;
	fileName_pthat_pptrack[6] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt220/HiForest_v77_v2_merged02/pt220_HP04_hiforest77_hiSignal.root";
	xsection_pptrack[6] = 2.477e-06;
	
	boundaries_pthat_pptrack[7] = 280;
	fileName_pthat_pptrack[7] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt280/HiForest_v77_v2_merged01/pt280_HP04_hiforest77_hiSignal.root";
	xsection_pptrack[7] = 6.160e-07;
	
	boundaries_pthat_pptrack[8] = 1000;
	xsection_pptrack[8] = 0;
	
	//these are the files required for the unfolding. 
	const int nbins_pthat_ppb = 9;
	Double_t boundaries_pthat_ppb[nbins_pthat_ppb+1];
	char *fileName_pthat_ppb[nbins_pthat_ppb+1];
	Double_t xsection_ppb[nbins_pthat_ppb+1]
	
	boundaries_pthat_ppb[0] = 15;
	fileName_pthat_ppb[0] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt15/HiForest_v77_merged01/pt15_HP04_prod16_v77_merged_forest_0.root";
	xsection_ppb[0] = 5.335e-01;// make sure that you have the correct one here.
	
	boundaries_pthat_ppb[1] = 30;
	fileName_pthat_ppb[1] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt30/HiForest_v77_merged01/pt30_HP04_prod16_v77_merged_forest_0.root";
	xsection_ppb[1] = 3.378e-02;
	
	boundaries_pthat_ppb[2] = 50;
	fileName_pthat_ppb[2] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt50/HiForest_v77_merged01/pt50_HP04_prod16_v77_merged_forest_0.root";
	xsection_ppb[2] = 3.778e-03;
	
	boundaries_pthat_ppb[3] = 80;
	fileName_pthat_ppb[3] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt80/HiForest_v77_merged01/pt80_HP04_prod16_v77_merged_forest_0.root";
	xsection_ppb[3] = 4.412e-04;
	
	boundaries_pthat_ppb[4] = 120;
	fileName_pthat_ppb[4] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt120/HiForest_v77_merged01/pt120_HP04_prod16_v77_merged_forest_0.root";
	xsection_ppb[4] = 6.147e-05;
	
	boundaries_pthat_ppb[5] = 170;
	fileName_pthat_ppb[5] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt170/HiForest_v77_merged01/pt170_HP04_prod16_v77_merged_forest_0.root";
	xsection_ppb[5] = 1.018e-05;
	
	boundaries_pthat_ppb[6] = 220;
	fileName_pthat_ppb[6] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt220/HiForest_v77_merged01/pt220_HP04_prod16_v77_merged_forest_0.root";
	xsection_ppb[6] = 2.477e-06;
	
	boundaries_pthat_ppb[7] = 280;
	fileName_pthat_ppb[7] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt280/HiForest_v77_merged01/pt280_HP04_prod16_v77_merged_forest_0.root";
	xsection_ppb[7] = 6.160e-07;
	
	boundaries_pthat_ppb[8] = 370;
	fileName_pthat_ppb[8] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt370/HiForest_v77_merged01/pt370_HP04_prod16_v77_merged_forest_0.root";
	xsection_ppb[8] = 1.088e-07;
	
	boundaries_pthat_ppb[8] = 1000;
	xsection_ppb[8] = 0;
	*/
	
	// lumi number for the sample 
	float ppblumi=30.9;
	//float pplumi=5300; //dont need this right now since comparison is with the pythia. 
	
	
	//output file
	TFile *output = new TFile("root://eoscms//eos/cms/store/group/phys_heavyions/rkunnawa/rpPb_calculation_data_w_mc_pp.root","RECREATE");
	
	// declare the histograms that are required for data and mc collection. do the unfolding later. also having only one centrality (0-100) for the preliminary analysis. 
	TH1F *meas_ppb_data_ak3pf = new TH1F("meas_ppb_data_ak3pf","measured ppb data ak3pf jt pt spectrum",nbins_rec,boundaries_rec);
	TH1F *meas_ppb_data_akpu3pf = new TH1F("meas_ppb_data_akpu3pf","measured ppb data akpu3pf jt pt spectrum",nbins_rec,boundaries_rec);

	//TH1F *meas_ppb_mc = new TH1F("meas_ppb_mc","measured ppb mc pythia + hijing jtpt spectrum",nbins_truth,boundaries_truth);
	//TH1F *meas_pp_mc = new TH1F("meas_pp_mc","measured pp signal sample 5.02 GeV jtpt spectrum",nbins_truth,boundaries_truth);
	
	//ppb data file:
	TFile *infData = new TFile();
	infData = TFile::Open("root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/pPb2013/promptReco/PA2013_HiForest_PromptReco_JSonPPb_forestv77.root");
	TTree *tDataEvt = (TTree*)infData->Get("hiEvtAnalyzer/HiTree");
	TTree *tDataSkim = (TTree*)infData->Get("skimanalysis/HltTree");
	TTree *tDataHlt = (TTree*)infData->Get("hltanalysis/HltTree");
	TTree *tDataJet1 = (TTree*)infData->Get("ak3PFJetAnalyzer/t");
	TTree *tDataJet2 = (TTree*)infData->Get("akPu3PFJetAnalyzer/t");
	
	tDataJet1->AddFriend(tDataEvt);
	tDataJet1->AddFriend(tDataHlt);
	tDataJet1->AddFriend(tDataSkim);
	tDataJet2->AddFriend(tDataEvt);
	tDataJet2->AddFriend(tDataHlt);
	tDataJet2->AddFriend(tDataSkim);
	
	TCut dataSelectionpPb = "abs(vz)<15&&pPAcollisionEventselectionPA&&pHBHENoiseFilter&&abs(jteta)<2";
	TCut TriggerSelectionpPb = "HLT_PAJet80__NoJetID_v1&&HLT_PAJet100_NoJetID_v1";
	
	tDataJet1->Project("meas_ppb_data_ak3pf","jtpt",dataSelectionpPb&&TriggerSelectionpPb);
	tDataJet2->Project("meas_ppb_data_akpu3pf","jtpt",dataSelectionpPb&&TriggerSelectionpPb);
	
	TCanvas *c1 = new TCanvas("c1","Pu difference in ak3 for ppb",600,600);
	meas_ppb_data_ak3pf->SetLineColor(kBlue);
	meas_ppb_data_akpu3pf->SetLineColor(kRed);
	meas_ppb_data_ak3pf->Draw();
	meas_ppb_data_akpu3pf->Draw("same");
	meas_ppb_data_akpu3pf->Write();
	meas_ppb_data_ak3pf->Write();
	output->Write();
	output->Close();
	timer.Stop();
	cout<<"cpu time = "<<timer.CpuTime()<<endl;
	cout<<"real time = "<<timer.RealTime()<<endl;

}



















































