// macro to merge the different hlt data files.
// here we are merging HLT_100 with HLT_80 and HLT_60.
// first lets normalize by the d\eta
// scale by the fraction of events which passed your eta cut (say from -1 to +1) by GetEntries()/GetEntries(abs(eta)<1)
// Procedure for the merging is nerged spectrum  = n3*prescl3 + n2 + n1 where
// n1 = HLT_100 from the HLT_80 file
// n2 = HLT_80 !HLT_100 from the HLT_80 file
// n3 = HLT_60 !HLT_80 !HLT_100 from the HLT_60 file
// prescl3 = GetEntries()/GetEntries(HLT_60) from the HLT_80 file
// Normalization for the data spectrum is 1/N_mb d^2N/(dp_T * d\eta)
// once the merged spectra is available then apply divide by bin width to get dp_T
// then we can normalize by the N_mb  = sigma_inelastic * Integrated Lumi.  



#include <TChain.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom.h>
#include <TCut.h>
#include <cstdlib>
#include <TCanvas.h>
#include <cmath>
#include <TLegend.h>
#include <TLatex.h>
#include <iostream>
using namespace std;

// divide by bin width
void divideBinWidth(TH1 *h)
{
	h->Sumw2();
	for (int i=0;i<=h->GetNbinsX();i++)
	{
		Float_t val = h->GetBinContent(i);
		Float_t valErr = h->GetBinError(i);
		val/=h->GetBinWidth(i);
		valErr/=h->GetBinWidth(i);
		h->SetBinContent(i,val);
		h->SetBinError(i,valErr);
	}
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
}

TLegend *myLegend(double x1,double y1,double x2, double y2)
{
	TLegend *leg = new TLegend(x1,y1,x2,y2);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	return leg;
	
}
void putCMSPrel(double x, double y, double size){
	TLatex *tex=0;
	tex = new TLatex(x,y,"CMS Preliminary");
	tex->SetTextSize(size);
	tex->SetLineWidth(2);
	tex->SetNDC();
	tex->Draw();
}
void drawText(const char *text, float xp, float yp, int size){
	TLatex *tex = new TLatex(xp,yp,text);
	tex->SetTextFont(63);
	tex->SetTextSize(size);
	tex->SetTextColor(kBlack);
	tex->SetLineWidth(1);
	//tex->SetTextFont(42);
	tex->SetNDC();
	tex->Draw();
}



void ppb_merge(){
	
	TStopwatch timer;
	timer.Start();
	
	
	//Float_t N_mb = Lumi_ppb*sigma_inelastic*1000;
	//Float_t Lumi_ppb = 30.9;// inverse micro barns
	//Float_t sigma_inelastic = 70.0;//milli barns

	TH1::SetDefaultSumw2();
	
	
	
	TString inname1 = "root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/pPb2013/promptReco/PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
	TString inname2 = "root://eoscms//eos/cms/store/group/phys_heavyions/krajczar/inbound/mnt/hadoop/cms/store/user/krajczar/pPb_Jet40Jet60_Full_v1/mergedJet40Jet60_KK.root";
	
	TFile *f1 = TFile::Open(inname1);
	TFile *f2 = TFile::Open(inname2);
	
	cout<<" File for HLT_100 HLT_80 = "<<inname1<<endl;
	cout<<" File for HLT_60 HLT_40 = "<<inname2<<endl;
	
	TFile *outfile = new TFile("pPbmerged_output.root","RECREATE");
	
	TTree* jet_80 = (TTree*)f1->Get("akPu3PFJetAnalyzer/t");
	TTree* jet_80_hlt = (TTree*)f1->Get("hltanalysis/HltTree");
	TTree* jet_80_skim = (TTree*)f1->Get("skimanalysis/HltTree");
	TTree* jet_80_evt = (TTree*)f1->Get("hiEvtAnalyzer/HiTree");
	jet_80->AddFriend(jet_80_hlt);
	jet_80->AddFriend(jet_80_skim);
	jet_80->AddFriend(jet_80_evt);
	
	TTree* jet_60 = (TTree*)f2->Get("akPu3PFJetAnalyzer/t");
	TTree* jet_60_hlt = (TTree*)f2->Get("hltanalysis/HltTree");
	TTree* jet_60_skim = (TTree*)f2->Get("skimanalysis/HltTree");
	TTree* jet_60_evt = (TTree*)f2->Get("hiEvtAnalyzer/HiTree");
	jet_60->AddFriend(jet_60_hlt);
	jet_60->AddFriend(jet_60_skim);
	jet_60->AddFriend(jet_60_evt);
	
	TCut Sel = "abs(vz)<15&&pHBHENoiseFilter&&pPAcollisionEventSelectionPA";
	TCut Trig_100 = "HLT_PAJet100_NoJetID_v1";
	TCut Trig_80 = "HLT_PAJet80_NoJetID_v1";
	TCut Trig_60 = "HLT_PAJet60_NoJetID_v1";
	TCut eta = "abs(jteta)<1";

	Float_t N_mb = 7.71e13;
	
	Float_t prescl3 = (Float_t)jet_80->GetEntries()/jet_80->GetEntries(Trig_60);
	Float_t prescl3_test = (Float_t)jet_80->GetEntries(Sel&&eta&&Trig_100)/jet_80->GetEntries(Sel&&eta&&Trig_60);

	
	const Int_t nbins = 16;
	const Double_t bound[nbins+1] = {30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 140., 150., 160., 200., 220., 260., 300.};
	
	TH1F *hMeas_100 = new TH1F("hMeas_100","PPb HLT_100 Measured histo",nbins,bound);
	TH1F *hMeas_80 = new TH1F("hMeas_80","PPb HLT_80 Measured histo",nbins,bound);
	TH1F *hMeas_60 = new TH1F("hMeas_60","PPb HLT_60 Measured histo",nbins,bound);
	TH1F *hCombined = new TH1F("hCombined","PPb Combined spectra",nbins,bound);
	
	hMeas_100->Sumw2();
	hMeas_60->Sumw2();
	hMeas_80->Sumw2();
	
	jet_80->Draw("jtpt>>hMeas_100",Sel&&Trig_100&&eta);
	jet_80->Draw("jtpt>>hMeas_80",Sel&&!Trig_100&&Trig_80&&eta);
	jet_60->Draw("jtpt>>hMeas_60","1.58606"*Sel&&eta&&Trig_60&&!Trig_80&&!Trig_100);
	
	
	divideBinWidth(hMeas_100);
	divideBinWidth(hMeas_80);
	divideBinWidth(hMeas_60);
	
	hMeas_100->Scale(1./2); //scaling by d eta  - from -1 to +1 is 2. 
	hMeas_80->Scale(1./2);
	hMeas_60->Scale(1./2);
	
	hCombined->Add(hMeas_100,1);
	hCombined->Add(hMeas_80,1);
	hCombined->Add(hMeas_60,1);
	
	hCombined->Scale(1./N_mb);
	hMeas_100->Scale(1./N_mb);
	hMeas_80->Scale(1./N_mb);
	hMeas_60->Scale(1./N_mb);
	
	TCanvas *cMerged = new TCanvas("cMerged","Merged PPb spectra",800,600);
	cMerged->SetLogy();

	hCombined->SetXTitle("Jet p_{T} [GeV/c]");
	hCombined->SetYTitle("1/N_mb d^2N/dp_t d eta");
	hCombined->SetMarkerColor(kBlack);
	hCombined->SetMarkerStyle(20);
	hCombined->SetAxisRange(30,300,"X");
	hCombined->Draw();
	hMeas_100->SetMarkerColor(kRed);
	hMeas_100->SetMarkerStyle(21);
	hMeas_100->Draw("same");
	hMeas_80->SetMarkerColor(kBlue);
	hMeas_80->SetMarkerStyle(22);
	hMeas_80->Draw("same");
	hMeas_60->SetMarkerColor(kGreen);
	hMeas_60->SetMarkerStyle(23);
	hMeas_60->Draw("same");
	
	TLegend *leg_PPb = myLegend(0.6,0.65,0.95,0.9);
	leg_PPb->SetTextSize(0.05);
	leg_PPb->AddEntry(hCombined,"Merged PPb spectra ","pl");
	leg_PPb->AddEntry(hMeas_100,"HLT_100 spectra","pl");
	leg_PPb->AddEntry(hMeas_80,"HLT_80 spectra","pl");
	leg_PPb->AddEntry(hMeas_60,"HLT_60 spectra","pl");
	leg_PPb->Draw();
	
	putCMSPrel(0.2,0.83,0.06);
	drawText("PPb AKPu3PF |eta|<1 |vz|<15",0.2,0.23,20);
	
	cMerged->SaveAs("pPb_merged.pdf","RECREATE");
	
	hMeas_100->Write();
	hCombined->Write();
	hMeas_80->Write();
	hMeas_60->Write();
	

	TFile *fYaxian = TFile::Open("AkPu3PFJetRpA.root");
	TH1F *Yaxian = (TH1F*)fYaxian->Get("DataJetWideBin;3");
	TH1F *test = (TH1F*)Yaxian->Clone("test");

	cout<<"hi"<<endl;
	//outfile->cd();
	Yaxian->Print("base");
	test->Print("base");
	Yaxian->Divide(hCombined);
	TCanvas *yaxian = new TCanvas("yaxian","",800,600);
	yaxian->Divide(2,1);
	yaxian->cd(1);
	Yaxian->SetTitle("ratio of Yaxian's measured pPb spectra to Mine");
	Yaxian->SetXTitle("Jet p_{T} [GeV/c]");
	Yaxian->SetMarkerColor(kBlack);
	Yaxian->SetMarkerStyle(23);
	Yaxian->Draw();
	
	yaxian->cd(2);
	yaxian->cd(2)->SetLogy();
	test->SetMarkerStyle(22);
	test->SetMarkerColor(kBlack);
	test->SetXTitle("Jet p_{T} [GeV/c]");
	test->SetYTitle("1/N_mb d^2N/dp_t d eta");
	test->Draw();
	hCombined->Draw("same");
	
	TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
	leg->SetTextSize(0.05);
	leg->AddEntry(hCombined,"Merged PPb spectra ","pl");
	leg->AddEntry(test,"Yaxian's spectra","pl");
	leg->Draw();
	
	putCMSPrel(0.2,0.83,0.06);
	drawText("PPb AKPu3PF |eta|<1 |vz|<15",0.2,0.23,20);
	
	yaxian->SaveAs("Yaxian_Comparison_pPb_pt_spectra.root","RECREATE");
	
	outfile->cd();
	Yaxian->Write();
	test->Write();
	
	outfile->Write();
	outfile->Close();
	
	timer.Stop();
    float rtime  = timer.RealTime();
    float ctime  = timer.CpuTime();
    
    std::cout<<"\t"<<std::endl;
    std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
    std::cout<<"\t"<<std::endl;
    std::cout<<"Good bye : " <<"\t"<<std::endl;
    
	//define the required histograms
	//static const Int_t nbins = 22;
	//static const Double_t bound[nbins+1] = {30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,180.,200.,220.,260.,300.,350.,400.,450.,500.};
	
}
