#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TAxis.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>

using namespace std;

using std::cout;
using std::endl;


void BRFinalAna(){

	TString infile = "infiles/Final.root";


	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();


	float Value;
	float ValueError;

	const int NHFQA = 9;
	int QAVtxPDGID[NHFQA] = {411,421,431,4122,511,521,531,443,553};
	TString HFName[NHFQA] = {"D^{+}","D^{0}","D^{+}_{s}","#Lambda_{c}^{+}","B^{0}","B^{+}","B^{0}_{s}","J/#psi","#Upsilon(1S)"};
	TString HFNameSave[NHFQA] = {"D+","D0","Ds","LambdaC","B0","B+","Bs","Jpsi","Upsilon"};

	float PDGcTau[NHFQA] = {1.040,0.4101,0.500,0.200,1.519,1.638,1.515,7.2/1000000000,1.21/100000000};
	float PDGcTauErr[NHFQA] = {0.007,0.0015,0.007,0.006,0.004,0.004,0.004,0.1/1000000000.1/100000000};

	const int NChannels = 10;

	TString DecayChannelLabels[NHFQA][NChannels] = {{"#pi^{0}#pi^{+}","#pi^{+}#pi^{-}#pi^{+}","#pi^{+}#pi^{0}K^{0}_{S}","#pi^{+}#pi^{+}K^{-}","#pi^{+}#pi^{+}#pi^{+}#pi^{-}#pi^{-}","#pi^{+}#phi^{0}","#mu^{+}#nu_{#mu}","K^{+}K^{0}_{S}","e^{+}X","#phi^{0}X"}, //D+
		{"#pi^{+}K^{-}","#pi^{+}#pi^{0}K^{-}","#pi^{-}K^{+}","K^{+}K^{-}K^{0}_{S}","3K^{0}_{S}","2#pi^{0}","#pi^{+}#pi^{-}","K^{0}_{S}X","#pi^{0}X","e^{+}X"}, //D0
		{"#mu^{+}#nu_{#mu}#phi^{0}","#mu^{+}#nu_{#mu}","#pi^{+}K^{+}K^{-}","K^{+}K^{0}_{S}","2pi^{0}#pi^{+}","p #bar p","#pi^{+}#phi^{0}","K^{0}_{S}X","#phi^{0}X","e^{+}X"}, //Ds
		{"pK^{-}#pi^{+}","pK^{0}_{S}","p#phi","p#pi^{0}#pi^{+}K^{-}","p#pi^{+}#pi^{-}#bar K^{0}","p#pi^{+}#pi^{-}#pi^{+}#pi^{-}","p#pi^{+}#pi^{-}","K^{+}K^{0}_{S}","e^{+}X","#phi^{0}X"}, //Lambda_C
		{"#pi^{0}#pi^{+}","#pi^{+}#pi^{-}#pi^{+}","#pi^{+}#pi^{0}K^{0}_{S}","#pi^{+}#pi^{+}K^{-}","#pi^{+}#pi^{+}#pi^{+}#pi^{-}#pi^{-}","#pi^{+}#phi^{0}","#mu^{+}#nu_{#mu}","#bar K^{0} X","pX","e^{+}X"},  //B+
		{"#pi^{-}K^{+}","#pi^{-}K^{+}#pi^{0}","#pi^{+}#pi^{-}","2#pi^{0}","D^{-}#pi^{+}","D^{-}D_{s}^{+}","J/#psiK^{+}#pi^{-}","J/#psi#pi^{+}#pi^{-}","DX","J/#psiX"}, //B0
		{"#pi^{0}#pi^{+}","#pi^{+}#pi^{-}#pi^{+}","#pi^{+}#pi^{0}K^{0}_{S}","#pi^{+}#pi^{+}K^{-}","#pi^{+}#pi^{+}#pi^{+}#pi^{-}#pi^{-}","#pi^{+}#phi^{0}","#mu^{+}#nu_{#mu}","K^{+}K^{0}_{S}","e^{+}X","#phi^{0}X"},  //Bs
		{"#pi^{0}#pi^{+}","#pi^{+}#pi^{-}#pi^{+}","#pi^{+}#pi^{0}K^{0}_{S}","#pi^{+}#pi^{+}K^{-}","#pi^{+}#pi^{+}#pi^{+}#pi^{-}#pi^{-}","#pi^{+}#phi^{0}","#mu^{+}#nu_{#mu}","K^{+}K^{0}_{S}","e^{+}X","#phi^{0}X"},  //Jpsi
		{"#pi^{0}#pi^{+}","#pi^{+}#pi^{-}#pi^{+}","#pi^{+}#pi^{0}K^{0}_{S}","#pi^{+}#pi^{+}K^{-}","#pi^{+}#pi^{+}#pi^{+}#pi^{-}#pi^{-}","#pi^{+}#phi^{0}","#mu^{+}#nu_{#mu}","K^{+}K^{0}_{S}","e^{+}X","#phi^{0}X"}};  //Upsilon(1S)



	float Light = 2.9979 * 100000000;
	float unit = 1.0/1000000;
	//Units
	for(int i = 0; i < NHFQA; i++){
		PDGcTau[i] = PDGcTau[i] * Light * unit;
		PDGcTauErr[i] = PDGcTauErr[i] * Light * unit;

	}

	const int NChannel = 10; 

	TFile * fin = new TFile(infile.Data());
	TH1D * HFHadronStat = (TH1D *) fin->Get("HFHadronStat");
	//HFHadronStat->SetCanExtend(TH1::kAllAxes);	


	TH1D * BR1DHis[NHFQA];

	TH1D * EvtGenBR[NHFQA];
	TH1D * PDGBR[NHFQA];

	TH1D * EvtGenEncoded[NHFQA];


	TH1D * Diff[NHFQA];
	TH1D * Diff2[NHFQA];
	TH1D * Diff3[NHFQA];
	int TotalEvents[NHFQA];

	TH1D * LifeTime[NHFQA];
	TF1 * func[NHFQA];


	TLegend * leg[NHFQA]; 
	TLegend * leg2[NHFQA]; 
	TLegend * leg3[NHFQA]; 
	TLegend * leg4[NHFQA]; 


	TH1D * QAE[NHFQA];
	TH1D * QAPx[NHFQA];
	TH1D * QAPy[NHFQA];
	TH1D * QAPz[NHFQA];
	TH1D * QACosTheta[NHFQA];


	for(int i = 0; i < NHFQA; i++){

		cout << "i = " << i << endl;


		LifeTime[i] = (TH1D *) fin->Get(Form("ProperLifeTime_%d",i));
		LifeTime[i]->GetXaxis()->SetTitle(Form("%s Decay Lifetime c #tau (cm)",HFName[i].Data()));
		LifeTime[i]->GetYaxis()->SetTitle("Number of Events");
		LifeTime[i]->GetXaxis()->CenterTitle();
		LifeTime[i]->GetYaxis()->CenterTitle();
		LifeTime[i]->GetXaxis()->SetTitleOffset(1.2);
		LifeTime[i]->GetYaxis()->SetTitleOffset(1.4);
		//	LifeTime[i]->GetYaxis()->SetTitleSize(0.14);

		LifeTime[i]->SetMarkerStyle(20);
		LifeTime[i]->SetMarkerSize(1);		
		LifeTime[i]->SetLineColor(kBlack);
		LifeTime[i]->SetMarkerColor(kBlack);

		LifeTime[i]->SetMaximum(200000);
		LifeTime[i]->SetMinimum(0);


		//QA E
		QAE[i] = (TH1D *) fin->Get(Form("QAE_%d",i));
		QAE[i]->GetXaxis()->SetTitle(Form("%s Energy: (Mother - All Daughters)/Mother",HFName[i].Data()));
		QAE[i]->GetYaxis()->SetTitle("Number of Events");
		QAE[i]->GetXaxis()->CenterTitle();
		QAE[i]->GetYaxis()->CenterTitle();
		QAE[i]->GetXaxis()->SetTitleOffset(1.2);
		QAE[i]->GetYaxis()->SetTitleOffset(1.4);

		QAE[i]->SetMarkerStyle(20);
		QAE[i]->SetMarkerSize(1);		
		QAE[i]->SetLineColor(kBlack);
		QAE[i]->SetMarkerColor(kBlack);
		QAE[i]->SetMinimum(0);


		//QA Px
		QAPx[i] = (TH1D *) fin->Get(Form("QAPx_%d",i));
		QAPx[i]->GetXaxis()->SetTitle(Form("%s p_{x}: (Mother - All Daughters)/Mother",HFName[i].Data()));
		QAPx[i]->GetYaxis()->SetTitle("Number of Events");
		QAPx[i]->GetXaxis()->CenterTitle();
		QAPx[i]->GetYaxis()->CenterTitle();
		QAPx[i]->GetXaxis()->SetTitleOffset(1.2);
		QAPx[i]->GetYaxis()->SetTitleOffset(1.4);

		QAPx[i]->SetMarkerStyle(20);
		QAPx[i]->SetMarkerSize(1);		
		QAPx[i]->SetLineColor(kBlack);
		QAPx[i]->SetMarkerColor(kBlack);
		QAPx[i]->SetMinimum(0);

		//QA Py
		QAPy[i] = (TH1D *) fin->Get(Form("QAPy_%d",i));
		QAPy[i]->GetXaxis()->SetTitle(Form("%s p_{y}: (Mother - All Daughters)/Mother",HFName[i].Data()));
		QAPy[i]->GetYaxis()->SetTitle("Number of Events");
		QAPy[i]->GetXaxis()->CenterTitle();
		QAPy[i]->GetYaxis()->CenterTitle();
		QAPy[i]->GetXaxis()->SetTitleOffset(1.2);
		QAPy[i]->GetYaxis()->SetTitleOffset(1.4);

		QAPy[i]->SetMarkerStyle(20);
		QAPy[i]->SetMarkerSize(1);		
		QAPy[i]->SetLineColor(kBlack);
		QAPy[i]->SetMarkerColor(kBlack);
		QAPy[i]->SetMinimum(0);



		//QA Pz
		QACosTheta[i] = (TH1D *) fin->Get(Form("QACosTheta_%d",i));
		QACosTheta[i]->GetXaxis()->SetTitle(Form("%s cos(#theta)",HFName[i].Data()));
		QACosTheta[i]->GetYaxis()->SetTitle("Number of Events");
		QACosTheta[i]->GetXaxis()->CenterTitle();
		QACosTheta[i]->GetYaxis()->CenterTitle();
		QACosTheta[i]->GetXaxis()->SetTitleOffset(1.2);
		QACosTheta[i]->GetYaxis()->SetTitleOffset(1.4);

		QACosTheta[i]->SetMarkerStyle(20);
		QACosTheta[i]->SetMarkerSize(1);		
		QACosTheta[i]->SetLineColor(kBlack);
		QACosTheta[i]->SetMarkerColor(kBlack);
		QACosTheta[i]->SetMinimum(0);



		//QA CosTheta
		QAPz[i] = (TH1D *) fin->Get(Form("QAPy_%d",i));
		QAPz[i]->GetXaxis()->SetTitle(Form("%s p_{z}: (Mother - All Daughters)/Mother",HFName[i].Data()));
		QAPz[i]->GetYaxis()->SetTitle("Number of Events");
		QAPz[i]->GetXaxis()->CenterTitle();
		QAPz[i]->GetYaxis()->CenterTitle();
		QAPz[i]->GetXaxis()->SetTitleOffset(1.2);
		QAPz[i]->GetYaxis()->SetTitleOffset(1.4);

		QAPz[i]->SetMarkerStyle(20);
		QAPz[i]->SetMarkerSize(1);		
		QAPz[i]->SetLineColor(kBlack);
		QAPz[i]->SetMarkerColor(kBlack);
		QAPz[i]->SetMinimum(0);


		EvtGenBR[i] = new TH1D(Form("EvtGenBR_%d",i),"",NChannel,-0.5,NChannel-0.5);

		EvtGenBR[i]->GetXaxis()->SetTitle(Form("EvtGen Simulation: %s Channel ID",HFName[i].Data()));
		EvtGenBR[i]->GetYaxis()->SetTitle("Branching Ratio");
		EvtGenBR[i]->GetXaxis()->CenterTitle();
		EvtGenBR[i]->GetYaxis()->CenterTitle();
		EvtGenBR[i]->GetYaxis()->SetTitleOffset(0.7);
		EvtGenBR[i]->GetYaxis()->SetTitleSize(0.06);

		EvtGenBR[i]->SetMarkerStyle(20);
		EvtGenBR[i]->SetMarkerSize(1);		
		EvtGenBR[i]->SetLineColor(kRed);
		EvtGenBR[i]->SetMarkerColor(kRed);

		EvtGenBR[i]->SetMaximum(1.5);
		EvtGenBR[i]->SetMinimum(5e-04);







		EvtGenEncoded[i] = new TH1D(Form("EvtGenEncoded_%d",i),"",NChannel,-0.5,NChannel-0.5);

		EvtGenEncoded[i]->GetXaxis()->SetTitle(Form("EvtGen Encoded: %s Channel ID",HFName[i].Data()));
		EvtGenEncoded[i]->GetYaxis()->SetTitle("Branching Ratio");
		EvtGenEncoded[i]->GetXaxis()->CenterTitle();
		EvtGenEncoded[i]->GetYaxis()->CenterTitle();
		EvtGenEncoded[i]->GetYaxis()->SetTitleOffset(0.7);
		EvtGenEncoded[i]->GetYaxis()->SetTitleSize(0.06);

		EvtGenEncoded[i]->SetMarkerStyle(20);
		EvtGenEncoded[i]->SetMarkerSize(1);		
		EvtGenEncoded[i]->SetLineColor(kGreen);
		EvtGenEncoded[i]->SetMarkerColor(kGreen);
		EvtGenEncoded[i]->SetMaximum(1.6);
		EvtGenEncoded[i]->SetMinimum(5e-04);




		




		PDGBR[i] = new TH1D(Form("PDGBR_%d",i),"",NChannel,-0.5,NChannel-0.5);
		PDGBR[i]->GetXaxis()->SetTitle(Form("PDG 2022: %s Channel ID",HFName[i].Data()));
		PDGBR[i]->GetYaxis()->SetTitle("Branching Ratio");
		PDGBR[i]->GetXaxis()->CenterTitle();
		PDGBR[i]->GetYaxis()->CenterTitle();
		PDGBR[i]->GetYaxis()->SetTitleOffset(1.2);
		PDGBR[i]->GetYaxis()->SetTitleSize(0.13);

		PDGBR[i]->SetMarkerStyle(20);
		PDGBR[i]->SetMarkerSize(1);		
		PDGBR[i]->SetLineColor(kBlue);
		PDGBR[i]->SetMarkerColor(kBlue);
		PDGBR[i]->SetMaximum(1.6);
		PDGBR[i]->SetMinimum(5e-04);






		Diff[i] = new TH1D(Form("Diff_%d",i),"",NChannel,-0.5,NChannel-0.5);
		//Diff[i]->GetXaxis()->SetTitle(Form("%s Channel ID",HFName[i].Data()));
		Diff[i]->GetYaxis()->SetTitle("Dev From PDG");
		Diff[i]->GetXaxis()->CenterTitle();
		Diff[i]->GetYaxis()->CenterTitle();
		Diff[i]->GetYaxis()->SetTitleOffset(0.4);
		Diff[i]->GetYaxis()->SetTitleSize(0.10);

		Diff[i]->GetXaxis()->SetTitleOffset(1.0);
		Diff[i]->GetXaxis()->SetTitleSize(0.13);


		Diff[i]->SetMarkerStyle(20);
		Diff[i]->SetMarkerSize(1);		
		Diff[i]->SetLineColor(kBlack);
		Diff[i]->SetMarkerColor(kBlack);
		Diff[i]->GetXaxis()->SetLabelSize(0.15);
		Diff[i]->GetYaxis()->SetLabelSize(0.07);

		Diff[i]->SetMaximum(1.6);
		Diff[i]->SetMinimum(-1.6);






		Diff2[i] = new TH1D(Form("Diff2_%d",i),"",NChannel,-0.5,NChannel-0.5);
		//Diff2[i]->GetXaxis()->SetTitle(Form("%s Channel ID",HFName[i].Data()));
		Diff2[i]->GetYaxis()->SetTitle("Dev From Decay Table");
		Diff2[i]->GetXaxis()->CenterTitle();
		Diff2[i]->GetYaxis()->CenterTitle();
		Diff2[i]->GetYaxis()->SetTitleOffset(0.4);
		Diff2[i]->GetYaxis()->SetTitleSize(0.10);

		Diff2[i]->GetXaxis()->SetTitleOffset(1.0);
		Diff2[i]->GetXaxis()->SetTitleSize(0.13);


		Diff2[i]->SetMarkerStyle(20);
		Diff2[i]->SetMarkerSize(1);		
		Diff2[i]->SetLineColor(kBlack);
		Diff2[i]->SetMarkerColor(kBlack);
		Diff2[i]->GetXaxis()->SetLabelSize(0.15);
		Diff2[i]->GetYaxis()->SetLabelSize(0.07);

		Diff2[i]->SetMaximum(1.6);
		Diff2[i]->SetMinimum(-1.6);






		Diff3[i] = new TH1D(Form("Diff3_%d",i),"",NChannel,-0.5,NChannel-0.5);
		//Diff3[i]->GetXaxis()->SetTitle(Form("%s Channel ID",HFName[i].Data()));
		Diff3[i]->GetYaxis()->SetTitle("Dev From PDG");
		Diff3[i]->GetXaxis()->CenterTitle();
		Diff3[i]->GetYaxis()->CenterTitle();
		Diff3[i]->GetYaxis()->SetTitleOffset(0.4);
		Diff3[i]->GetYaxis()->SetTitleSize(0.10);

		Diff3[i]->GetXaxis()->SetTitleOffset(1.0);
		Diff3[i]->GetXaxis()->SetTitleSize(0.13);


		Diff3[i]->SetMarkerStyle(20);
		Diff3[i]->SetMarkerSize(1);		
		Diff3[i]->SetLineColor(kBlack);
		Diff3[i]->SetMarkerColor(kBlack);
		Diff3[i]->GetXaxis()->SetLabelSize(0.15);
		Diff3[i]->GetYaxis()->SetLabelSize(0.07);

		Diff3[i]->SetMaximum(1.6);
		Diff3[i]->SetMinimum(-1.6);


		//Redefintion of Axis//
		for(int r = 0; r < NChannels; r++){

			Value = EvtGenBR[i]->GetBinContent(r+1);
			ValueError = EvtGenBR[i]->GetBinError(r+1);
	
			EvtGenBR[i]->Fill(DecayChannelLabels[i][r],Value);
			EvtGenBR[i]->SetBinContent(r+1,Value);
			EvtGenBR[i]->SetBinError(r+1,ValueError);

			Value = EvtGenEncoded[i]->GetBinContent(r+1);
			ValueError = EvtGenEncoded[i]->GetBinError(r+1);
		
			EvtGenEncoded[i]->Fill(DecayChannelLabels[i][r],Value);
			EvtGenEncoded[i]->SetBinContent(r+1,Value);
			EvtGenEncoded[i]->SetBinError(r+1,ValueError);


			Value = PDGBR[i]->GetBinContent(r+1);
			ValueError = PDGBR[i]->GetBinError(r+1);
		
			PDGBR[i]->Fill(DecayChannelLabels[i][r],Value);
			PDGBR[i]->SetBinContent(r+1,Value);
			PDGBR[i]->SetBinError(r+1,ValueError);

	

			Value = Diff[i]->GetBinContent(r+1);
			ValueError = Diff[i]->GetBinError(r+1);
		
			Diff[i]->Fill(DecayChannelLabels[i][r],Value);
			Diff[i]->SetBinContent(r+1,Value);
			Diff[i]->SetBinError(r+1,ValueError);

			Value = Diff2[i]->GetBinContent(r+1);
			ValueError = Diff2[i]->GetBinError(r+1);
		
			Diff2[i]->Fill(DecayChannelLabels[i][r],Value);
			Diff2[i]->SetBinContent(r+1,Value);
			Diff2[i]->SetBinError(r+1,ValueError);


			Value = Diff3[i]->GetBinContent(r+1);
			ValueError = Diff3[i]->GetBinError(r+1);
		
			Diff3[i]->Fill(DecayChannelLabels[i][r],Value);
			Diff3[i]->SetBinContent(r+1,Value);
			Diff3[i]->SetBinError(r+1,ValueError);


			EvtGenBR[i]->GetYaxis()->SetTitle("Branching Ratio");
			EvtGenBR[i]->GetYaxis()->CenterTitle();
			EvtGenBR[i]->GetXaxis()->SetLabelSize(0.08);
			EvtGenBR[i]->GetXaxis()->SetLabelOffset(HFHadronStat->GetXaxis()->GetLabelOffset() * 3.2);

			EvtGenEncoded[i]->GetYaxis()->SetTitle("Branching Ratio");
			EvtGenEncoded[i]->GetYaxis()->CenterTitle();
			EvtGenEncoded[i]->GetXaxis()->SetLabelSize(0.08);
			EvtGenEncoded[i]->GetXaxis()->SetLabelOffset(HFHadronStat->GetXaxis()->GetLabelOffset() * 3.2);


			PDGBR[i]->GetYaxis()->SetTitle("Branching Ratio");
			PDGBR[i]->GetYaxis()->CenterTitle();
			PDGBR[i]->GetXaxis()->SetLabelSize(0.08);
			PDGBR[i]->GetXaxis()->SetLabelOffset(HFHadronStat->GetXaxis()->GetLabelOffset() * 3.2);



			Diff[i]->GetYaxis()->SetTitle("Dev From PDG");
			Diff[i]->GetYaxis()->CenterTitle();
			Diff[i]->GetXaxis()->SetLabelSize(0.20);
			Diff[i]->GetXaxis()->SetLabelOffset(HFHadronStat->GetXaxis()->GetLabelOffset() * 6.2);

			Diff2[i]->GetYaxis()->SetTitle("Dev From Decay Table");
			Diff2[i]->GetYaxis()->CenterTitle();
			Diff2[i]->GetXaxis()->SetLabelSize(0.20);
			Diff2[i]->GetXaxis()->SetLabelOffset(HFHadronStat->GetXaxis()->GetLabelOffset() * 6.2);


			Diff3[i]->GetYaxis()->SetTitle("Dev From PDG");
			Diff3[i]->GetYaxis()->CenterTitle();
			Diff3[i]->GetXaxis()->SetLabelSize(0.20);
			Diff3[i]->GetXaxis()->SetLabelOffset(HFHadronStat->GetXaxis()->GetLabelOffset() * 6.2);



		}



		/*
		   Decay Channel Explanations


		   D+

		   0. Decay D+ -> pi0 pi+ with all intermediate resonances 
		   1. Decay D+ -> pi+ pi- pi+ with all intermediate resonances 
		   2. Decay D+ -> pi+ pi0 K0s 
		   3. Decay D+ -> pi+ pi+ K- 
		   4. Decay D+ -> pi+ pi+ pi+ pi- pi-
		   5. Decay D+ -> pi+ phi(1020)
		   6. Decay D+ -> mu+ + v_mu
		   7. Decay D+ -> K+ K0s 
		   8. Decay D+ -> e+ + X
		   9. Decay D+ -> phi + X

		   D0

		   0. Decay D0 -> pi+ K- with all intermediate resonances 
		   1. Decay D0 -> pi+ pi0 K-  
		   2. Decay D0 -> pi- K+ with all intermediate resonances 
		   3. Decay D0 -> K+ K- K0s 
		   4. Decay D0 -> K0s K0s K0s 
		   5. Decay D0 -> pi0 pi0
		   6. Decay D0 -> pi+ pi- 
		   7. Decay D0 -> K0s + X
		   8. Decay D0 -> pi0 + X
		   9. Decay D0 -> e+ + X


		   Ds+

		   0. Decay Ds+ -> mu+ v_mu phi 
		   1. Decay Ds+ -> mu+ v_mu 
		   2. Decay Ds+ -> pi+ K+ K-  
		   3. Decay Ds+ -> K+ K0s 
		   4. Decay Ds+ -> pi+ eta 
		   5. Decay Ds+ -> p n bar
		   6. Decay Ds+ -> pi+ phi
		   7. Decay Ds+ -> Ks + X
		   8. Decay Ds+ -> phi + X
		   9. Decay Ds+ -> e+ + X



		   LambdaC+

		   0. Decay LambdaC+ -> p K- pi+ 
		   1. Decay LambdaC+ -> p K0 bar
		   2. Decay LambdaC+ -> p phi  
		   3. Decay LambdaC+ -> p K- pi+ pi0 
		   4. Decay LambdaC+ -> p pi+ pi- K0bar
		   5. Decay LambdaC+ -> p pi+ pi- pi+ pi-
		   6. Decay LambdaC+ -> p pi+ pi-
		   7. Decay LambdaC+ -> K0bar + X
		   8. Decay LambdaC+ -> p + X
		   9. Decay LambdaC+ -> e+ + X

		   B0

		   0. Decay B0 -> pi- K+ with all intermediate resonances 
		   1. Decay B0 -> pi- K+ pi0 with all intermediate resonances 
		   2. Decay B0 -> pi+ pi- 
		   3. Decay B0 -> pi0 pi0  
		   4. Decay B0 -> D- pi+
		   5. Decay B0 -> D- Ds+ 
		   6. Decay B0 -> Jpsi K+ pi-
		   7. Decay B0 -> Jpsi pi+ pi-
		   8. Decay B0 -> D + X
		   9. Decay B0 -> Jpsi + X 

		   B+

		   0. Decay B+ -> K+ K- pi+ with all intermediate resonances 
		   1. Decay B+ -> K+ K- K+ with all intermediate resonances 
		   2. Decay B+ -> K+ phi(1020) 
		   3. Decay B+ -> K+ J/psi  
		   4. Decay B+ -> K+ Dbar
		   5. Decay B+ -> Ds+ Dbar
		   6. Decay B+ -> K+ K*0Bar(1430) (K*0Bar(1430) -> K- pi+)
		   7. Decay B+ -> e ve X
		   8. Decay B+ -> D + X
		   9. Decay B+ -> Jpsi + X



			Bs0

			0. Decay Bs -> pi+ K-
			1. Decay Bs -> K+ K- 
			2. Decay Bs -> p bar p K+ K-
			3. Decay Bs -> phi gamma
			4. Decay Bs -> D- pi+
			5. Decay Bs -> Ds+ Ds-
			6. Decay Bs -> Jpsi K+ K-
			7. Decay Bs -> Jpsi phi
			8. Decay Bs -> Ds- + X
			9. Decay Bs -> Jpsi + X


			Jpsi

			0. Decay Jpsi -> e+ e- 
			1. Decay Jpsi -> mu+ mu-
			2. Decay Jpsi -> pi+ pi-
			3. Decay Jpsi -> pi+ pi- pi0
			4. Decay Jpsi -> K+ K-
			5. Decay Jpsi -> K+ K- pi+ pi-	
			6. Decay Jpsi -> p+ p-
			7. Decay Jpsi -> n bar n
			8. Decay Jpsi -> gamma gamma gamma
			9. Decay Jpsi -> gamma pi0


			Upsilon(1S)

			0. Decay Upsilon -> e+ e- 
			1. Decay Upsilon -> mu+ mu-
			2. Decay Upsilon -> gamma pi+ pi-
			3. Decay Upsilon -> pi+ pi- pi0
			4. Decay Upsilon -> gamma K+ K-
			5. Decay Upsilon -> gamma pi0 pi0	
			6. Decay Upsilon -> phi K+ K-
			7. Decay Upsilon -> pi+ pi- pi0	pi0
			8. Decay Upsilon -> gamma pi+ pi- p+ p-
			9. Decay Upsilon -> Jpsi + X


			*/
			//  {0.001247,0.00327,0.01562,0.0983,0.00570,0.00169,0,0.1607,0.0112,0.063}
			//  {0.00033,0.00018,0.00031,0.0016,0.00014,0.00011,0,0.003,0.0004,0.007}

			//float PDGValue[NHFQA][NChannel] = {{0.001247,0.00327,0.0736,0.0938,0.00570,0.00166,3.74e-04,0.00304,0.1607,0.0112},{0.03947,0.0164,0.000150,0.00442,7.5e-04,8.26e-04,0.00515,0.3,0.2,0.0649},{0.019,0.00543,0.0538,0.01453,0.0168,0.00122,0.045,0.190,0.157,0.0633},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{5.2e-06,4.88e-05,8.8e-06,0.00102,3.69e-04,0.009,3.8e-07,0.108,0.95,0.01},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}};

			//	float PDGValue[NHFQA][NChannel] = {{0.001247,0.00327,0.0736,0.0938,0.00570,0.00166,3.74e-04,0.00304,0.1607,0.0112},{0.03947,0.144,0.000150,0.00442,7.5e-04,8.26e-04,0.001454,0.3,0.2,0.0649},{0.019,0.00543,0.0538,0.01453,0.0168,0.00122,0.045,0.190,0.157,0.0633},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{5.2e-06,4.88e-05,8.8e-06,0.00102,3.69e-04,0.009,3.8e-07,0.108,0.95,0.01},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}};


			//	float PDGValueErr[NHFQA][NChannel] = {{0.00033,0.00018,0.0021,0.0016,0.00014,0.00016,0.17e-04,0.00009,0.003,0.0004},{0.00030,0.0050,0.000007,0.00032,0.7e-04,0.25e-04,0.000024,0.01,0.005,0.0011},{0.005,0.00015,0.0010,0.0035,0.0009,0.00011,0.004,0.011,0.010,0.0015},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0.4e-06,3.24293e-05,0.65e-06,0.000019,0.16e-04,0.0009,1.3e-07,0.004,0.002,0.002},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}};


			/*
			   Upsilon 

			   0. Decay Upsilon -> e+ e- 
			   1. Decay Upsilon -> mu+ mu-
			   2. Decay Upsilon -> gamma pi+ pi-
			   3. Decay Upsilon -> pi+ pi- pi0
			   4. Decay Upsilon -> gamma K+ K-
			   5. Decay Upsilon -> gamma pi0 pi0	
			   6. Decay Upsilon -> phi K+ K-
			   7. Decay Upsilon -> pi+ pi- pi0	pi0
			   8. Decay Upsilon -> gamma pi+ pi- p+ p-
			   9. Decay Upsilon -> Jpsi + X


*/
			cout << "Pass 1" << endl;


		float PDGValue[NHFQA][NChannel] = {{0.001247,0.00327,0.0736,0.0938,0.00570,0.00166,3.74e-04,0.00304,0.1607,0.0112},    //D+
			{0.03947,0.0164,0.000150,0.00442,7.5e-04,8.26e-04,0.00515,0.3,0.2,0.0649},  //D0 
			{0.019,0.00543,0.0538,0.01453,0.0168,0.00122,0.045,0.190,0.157,0.0633},  //Ds
			{0.0628,0.0159,0.00106,0.000111,0.00010,0.0182,0.00461,0.099,0.50,0.0395},   //LambdaC+
			{1.96e-05,3.78e-05,5.12e-06,1.59e-06,2.51e-03,7.2e-03,1.15e-03,4.00e-05,1.027,0.0192},   //B0
			{5.2e-06,4.88e-05,8.8e-06,0.00102,3.69e-04,0.009,3.8e-07,0.108,0.95,0.01}, //B+ 
			{2.2e-06,2.66e-05,1.85e-05,3.4e-05,2.98e-03,4.4e-03,7.9e-04,1.04e-03,0.62,0.0192},   //Bs
			{0.05971,0.05961,1.47e-04,0.0210,2.88e-03,6.86e-03,2.120e-03,2.09e-03,1.16e-05,3.56e-05},   //Jpsi
			{0.0238,0.0248,6.3e-05,2.1e-06,1.14e-05,1.7e-05,2.4e-06,1.28e-05,1.5e-04,5.4e-04}    //Upsilon(1S)
		};  


		float PDGValueErr[NHFQA][NChannel] = {
			{0.00033,0.00018,0.0021,0.0016,0.00014,0.00016,0.17e-04,0.00009,0.003,0.0004},   //D+ 
			{0.00030,0.0050,0.000007,0.00032,0.7e-04,0.25e-04,0.000024,0.01,0.005,0.0011},  //D0
			{0.005,0.00015,0.0010,0.0035,0.0009,0.00011,0.004,0.011,0.010,0.0015},   //Ds
			{0.0032,0.0008,0.00014,0.000018,0.00004,0.0025,0.00028,0.0070,0.16,0.0035},   //LambdaC+
			{0.05e-05,0.32e-05,0.19e-06,0.26e-06,0.08e-03,0.8e-03,0.05e-03,0.15e-05,0.010,0.002},   //B0
			{0.4e-06,3.24293e-05,0.65e-06,0.000019,0.16e-04,0.0009,1.3e-07,0.004,0.002,.002},   //B+
			{0.7e-06,0.22e-05,0.14e-05,0.4e-05,0.14e-03,0.5e-03,0.7e-04,0.04e-03,0.06,0.002},   //Bs
			{0.032,0.033,0.14e-04,0.0008,0.12e-03,0.28e-03,0.029e-03,0.16e-03,0.22e-05,0.17e-05},   //Jpsi
			{0.0011,0.0005,1.8e-05,0.8e-06,0.13e-05,0.7e-05,0.5e-06,0.30e-05,0.6e-04,0.4e-04}     //Upsilon(1S)
		}; 

		/*

		   Upsilon(1S)


		   0. Decay Upsilon -> e+ e- 
		   1. Decay Upsilon -> mu+ mu-
		   2. Decay Upsilon -> gamma pi+ pi-
		   3. Decay Upsilon -> pi+ pi- pi0
		   4. Decay Upsilon -> gamma K+ K-
		   5. Decay Upsilon -> gamma pi0 pi0	
		   6. Decay Upsilon -> phi K+ K-
		   7. Decay Upsilon -> pi+ pi- pi0	pi0
		   8. Decay Upsilon -> gamma pi+ pi- p+ p-
		   9. Decay Upsilon -> Jpsi + X
		   */


		float EncodedValue[NHFQA][NChannel] = {{0.00126,0.00244,0.0699,0.09400,0.001660,0.004660214,0.000382,0.00286,0.16313194,0.027660214},   //D+
			{0.0389,0.1390,0.00013745,0.00465,0.00095,0.0008,0.001397,0.1529059,0.37957735,0.0651},   //D0
			{0.018309605,0.005800000,0.003811325,0.01490,0.0065,0.00130,0.04500,0.02944,0.19589603,0.06432038},   //Ds+
			{0.0254,0.0230,0.00082,0.02300,0.02600,0.00180,0.0007,0.10370,0.20462,0.0480},   //LambdaC 
			{0.0000194,0.0000275,0.00000513,0.00000162,0.002680,0.00720,0.00133,0.0000270,0.108,0.0071881},   //B0
			{0.000005000,0.000025400,0.000008300,0.001014000,0.000368000,0.010000000,0.000001,0.108,1.09,0.0192},   //B+
			{0.0000049,0.000033,0.000014,0.000057,0.00320,0.0104,0.0007,0.00130,0.1144,0.0265},   //Bs
			{0.05940,0.05930,0.000147,0.02070,0.000237000,0.0003400,0.002170,0.00220,0.000012,0.0000349},   //Jpsi
			{0.024800,0.024800,0.0000630,2.1e-06,0.0000114,0.0000170,0.014959973,1.28e-05,0.000250000,0.044879919}    //Upsilon(1S)
		};





		for(int j = 0; j < NChannel; j++){

			PDGBR[i]->SetBinContent(j+1,PDGValue[i][j]);
			PDGBR[i]->SetBinError(j+1,PDGValueErr[i][j]);

			EvtGenEncoded[i]->SetBinContent(j+1,EncodedValue[i][j]);
			EvtGenEncoded[i]->SetBinError(j+1,EncodedValue[i][j] * 0.01);


		}




		float ChannelStat;
		float TotalStat;
		float BR;


		float ChannelStatErr;
		float BRErr;



		float Dev;
		float DevErr;


		float Dev2;
		float Dev2Err;


		float Dev3;
		float Dev3Err;


		TotalEvents[i] = HFHadronStat->GetBinContent(i+1);		

		cout << "Pass 2" << endl;


		BR1DHis[i] = (TH1D *) fin->Get(Form("BR1DHis_%d",i));
		TotalStat = HFHadronStat->GetBinContent(i+1);

		//		cout << "TotalStat = " << TotalStat << endl;
		for(int j = 0; j < NChannel; j++){

			ChannelStat = BR1DHis[i]->GetBinContent(j+1);
			ChannelStatErr = BR1DHis[i]->GetBinError(j+1);
			BR = ChannelStat/TotalStat;
			BRErr = ChannelStatErr/TotalStat;

			EvtGenBR[i]->SetBinContent(j+1,BR);
			EvtGenBR[i]->SetBinError(j+1,BRErr);

			Dev =  (BR - PDGValue[i][j])/PDGValue[i][j];
			DevErr = (BR/PDGValue[i][j]) * sqrt(BRErr/BR * BRErr/BR + PDGValueErr[i][j]/PDGValue[i][j] * PDGValueErr[i][j]/PDGValue[i][j]);

			Diff[i]->SetBinContent(j+1,Dev);
			Diff[i]->SetBinError(j+1,DevErr);


			Dev2 =  (BR - EncodedValue[i][j])/EncodedValue[i][j];
			Dev2Err = BRErr/EncodedValue[i][j];


			Diff2[i]->SetBinContent(j+1,Dev2);
			Diff2[i]->SetBinError(j+1,Dev2Err);


			Dev3 =  (EncodedValue[i][j] - PDGValue[i][j])/PDGValue[i][j];
			Dev3Err = PDGValueErr[i][j]/EncodedValue[i][j];



			Diff3[i]->SetBinContent(j+1,Dev3);
			Diff3[i]->SetBinError(j+1,Dev3Err);


			//	if(i == PartID) cout << "Dev = " << Dev << endl;
			if(i == 2) cout << "Dev2 = " << Dev2 << endl;

		}









		TLatex *lat = new TLatex();
		lat->SetNDC();
		lat->SetTextSize(0.1);


		TLine* lPull = new TLine(-0.5, 0, 9.5, 0);
		lPull->SetLineWidth(1);
		lPull->SetLineStyle(2);
		lPull->SetLineColor(2);

		c->cd();

		leg[i] = new TLegend(0.12,0.60,0.50,0.75,NULL,"brNDC");
		leg[i]->SetBorderSize(0);
		leg[i]->SetTextSize(0.045);
		leg[i]->SetTextFont(42);
		leg[i]->SetFillStyle(0);
		leg[i]->SetLineWidth(2);

		TPad* p1 = new TPad("p1","p1",0,0.3,1,1);
		p1->SetBottomMargin(0);

		p1->Draw();
		p1->cd();

		p1->SetLogy();


		EvtGenBR[i]->Draw("ep");
		PDGBR[i]->Draw("epSAME");

		lat->DrawLatex(0.15,0.78,HFName[i].Data());
		lat->SetTextSize(0.05);	
		lat->DrawLatex(0.40,0.78,Form("Total Statistics: %d",TotalEvents[i]));
		lat->SetTextSize(0.1);	

		leg[i]->AddEntry(EvtGenBR[i],"EvtGen","PL");
		leg[i]->AddEntry(PDGBR[i],"PDG 2022","PL");

		leg[i]->Draw("SAME");

		c->cd();

		TPad* p2 = new TPad("p2","p2",0,0,1,0.3);
		p2->SetTopMargin(0);
		p2->SetBottomMargin(0.3);

		p2->Draw();
		p2->cd();


		Diff[i]->Draw("ep");

		lat->DrawLatex(0.15,0.8,"(EvtGen Sim - PDG)/PDG");

		lPull->Draw("SAME");

		c->cd();

		c->SaveAs(Form("Plots/%s/SimPDG.png",HFNameSave[i].Data()));



		c->cd();

		leg2[i] = new TLegend(0.12,0.60,0.50,0.75,NULL,"brNDC");
		leg2[i]->SetBorderSize(0);
		leg2[i]->SetTextSize(0.045);
		leg2[i]->SetTextFont(42);
		leg2[i]->SetFillStyle(0);
		leg2[i]->SetLineWidth(2);

		TPad* p3 = new TPad("p1","p1",0,0.3,1,1);
		p3->SetBottomMargin(0);

		p3->Draw();
		p3->cd();

		p3->SetLogy();


		EvtGenBR[i]->Draw("ep");
		EvtGenEncoded[i]->Draw("epSAME");

		lat->DrawLatex(0.15,0.78,HFName[i].Data());

		lat->SetTextSize(0.05);			
		lat->DrawLatex(0.40,0.78,Form("Total Statistics: %d",TotalEvents[i]));
		lat->SetTextSize(0.1);	

		leg2[i]->AddEntry(EvtGenBR[i],"EvtGen Sim","PL");
		leg2[i]->AddEntry(EvtGenEncoded[i],"EvtGen Decay Table","PL");

		leg2[i]->Draw("SAME");

		c->cd();

		TPad* p4 = new TPad("p2","p2",0,0,1,0.3);
		p4->SetTopMargin(0);
		p4->SetBottomMargin(0.3);

		p4->Draw();PDGcTau[i] * Light * unit;

		p4->cd();

		Diff2[i]->Draw("ep");

		lat->DrawLatex(0.15,0.8,"EvtGen: (Sim - Decay Table)/Decay Table");

		lPull->Draw("SAME");

		c->cd();

		c->SaveAs(Form("Plots/%s/SimEvtGen.png",HFNameSave[i].Data()));





		c->cd();

		leg3[i] = new TLegend(0.12,0.60,0.50,0.75,NULL,"brNDC");
		leg3[i]->SetBorderSize(0);
		leg3[i]->SetTextSize(0.045);
		leg3[i]->SetTextFont(42);
		leg3[i]->SetFillStyle(0);
		leg3[i]->SetLineWidth(2);

		TPad* p5 = new TPad("p1","p1",0,0.3,1,1);
		p5->SetBottomMargin(0);

		p5->Draw();
		p5->cd();

		p5->SetLogy();


		EvtGenEncoded[i]->Draw("ep");
		PDGBR[i]->Draw("epSAME");

		lat->DrawLatex(0.15,0.78,HFName[i].Data());

		lat->SetTextSize(0.05);			
		lat->DrawLatex(0.40,0.78,Form("Total Statistics: %d",TotalEvents[i]));
		lat->SetTextSize(0.1);	


		leg3[i]->AddEntry(EvtGenEncoded[i],"EvtGen Decay Table","PL");		
		leg3[i]->AddEntry(PDGBR[i],"PDG 2022","PL");

		leg3[i]->Draw("SAME");

		c->cd();

		TPad* p6 = new TPad("p2","p2",0,0,1,0.3);
		p6->SetTopMargin(0);
		p6->SetBottomMargin(0.3);

		p6->Draw();
		p6->cd();

		Diff3[i]->Draw("ep");

		lat->DrawLatex(0.15,0.8,"(EvtGen Decay Table - PDG)/PDG");

		lPull->Draw("SAME");

		c->cd();
		c->SaveAs(Form("Plots/%s/EvtGenPDG.png",HFNameSave[i].Data()));


		c->cd();

		func[i] = new TF1(Form("Func_%d",i),"[1] * TMath::Exp(-[0] * x)",0.01,0.04);

		LifeTime[i]->Fit(func[i],"R");
		LifeTime[i]->Draw("ep");

		leg4[i] = new TLegend(0.12,0.60,0.50,0.75,NULL,"brNDC");
		leg4[i]->SetBorderSize(0);
		leg4[i]->SetTextSize(0.045);
		leg4[i]->SetTextFont(42);
		leg4[i]->SetFillStyle(0);
		leg4[i]->SetLineWidth(2);

		lat->DrawLatex(0.15,0.78,HFName[i].Data());
		lat->SetTextSize(0.038);			
		lat->DrawLatex(0.32,0.78,Form("Total Statistics: %d",TotalEvents[i]));

		leg4[i]->AddEntry(LifeTime[i],"EvtGen Simulation","PL");		
		leg4[i]->AddEntry(func[i],"Fit","PL");

		float cTau =  1.0/func[i]->GetParameter(0) * 10000;
		float cTauErr =  func[i]->GetParError(0) * func[i]->GetParError(0)/func[i]->GetParameter(0) * 10000;

		lat->DrawLatex(0.32,0.70,Form("EvtGen: c #tau = %.1f #pm %.1f (#mu m)",cTau,cTauErr));
		lat->DrawLatex(0.32,0.62,Form("PDG: c #tau = %.1f #pm %.1f (#mu m)",PDGcTau[i],PDGcTauErr[i]));

		c->SaveAs(Form("Plots/%s/LifeTime.png",HFNameSave[i].Data()));

		QAE[i]->Draw("ep");
		float EnergyRMS = QAE[i]->GetRMS();

		lat->SetTextSize(0.1);			
		lat->DrawLatex(0.15,0.78,HFName[i].Data());
		lat->SetTextSize(0.038);			
		lat->DrawLatex(0.32,0.78,Form("Total Statistics: %d",TotalEvents[i]));
		lat->DrawLatex(0.32,0.70,Form("RMS: %.3f",EnergyRMS));

		c->SaveAs(Form("Plots/%s/Energy.png",HFNameSave[i].Data()));
			
		QAPx[i]->Draw("ep");

		lat->SetTextSize(0.1);			
		lat->DrawLatex(0.15,0.78,HFName[i].Data());
		lat->SetTextSize(0.038);			
		lat->DrawLatex(0.32,0.78,Form("Total Statistics: %d",TotalEvents[i]));

		c->SaveAs(Form("Plots/%s/Px.png",HFNameSave[i].Data()));

		QAPy[i]->Draw("ep");

		lat->SetTextSize(0.1);			
		lat->DrawLatex(0.15,0.78,HFName[i].Data());
		lat->SetTextSize(0.038);			
		lat->DrawLatex(0.32,0.78,Form("Total Statistics: %d",TotalEvents[i]));

		c->SaveAs(Form("Plots/%s/Py.png",HFNameSave[i].Data()));



		QAPz[i]->Draw("ep");

		lat->SetTextSize(0.1);			
		lat->DrawLatex(0.15,0.78,HFName[i].Data());
		lat->SetTextSize(0.038);			
		lat->DrawLatex(0.32,0.78,Form("Total Statistics: %d",TotalEvents[i]));

		c->SaveAs(Form("Plots/%s/Pz.png",HFNameSave[i].Data()));

		QACosTheta[i]->Draw("ep");
		lat->SetTextSize(0.1);			
		lat->DrawLatex(0.15,0.78,HFName[i].Data());
		lat->SetTextSize(0.038);			
		lat->DrawLatex(0.32,0.78,Form("Total Statistics: %d",TotalEvents[i]));

		c->SaveAs(Form("Plots/%s/CosTheta.png",HFNameSave[i].Data()));



	}



	//Redefintion of Axis//
	for(int r = 0; r < NHFQA; r++){

		Value = HFHadronStat->GetBinContent(r+1);
		ValueError = HFHadronStat->GetBinError(r+1);
		
		HFHadronStat->Fill(HFName[r],Value);
		HFHadronStat->SetBinContent(r+1,Value);
		HFHadronStat->SetBinError(r+1,ValueError);


	}

	c->cd();
	c->SetLeftMargin(0.15);

	HFHadronStat->GetYaxis()->SetTitle("Counts");
	HFHadronStat->GetYaxis()->CenterTitle();
	HFHadronStat->GetXaxis()->SetLabelSize(0.08);
	HFHadronStat->GetXaxis()->SetLabelOffset(HFHadronStat->GetXaxis()->GetLabelOffset() * 3.2);

	HFHadronStat->Draw("ep");
	//HFHadronStat->Draw("text");
	c->SaveAs("Plots/HFHadronStat.png");

}
