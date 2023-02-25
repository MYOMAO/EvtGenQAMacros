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


void PHOTOSExtra(){


	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.1);





	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	
	c->SetLeftMargin(0.13);	

	TFile * fin = new TFile("Extras/ZZ_qa.root");
	fin->cd();


	TH1D * D0KPiExtraPhotons = (TH1D *) fin->Get("D0KPiExtraPhotons");
	D0KPiExtraPhotons->GetXaxis()->SetTitle("Number of Photons Per Event");
	D0KPiExtraPhotons->GetYaxis()->SetTitle("Number of Events");
	D0KPiExtraPhotons->GetXaxis()->CenterTitle();
	D0KPiExtraPhotons->GetYaxis()->CenterTitle();
	D0KPiExtraPhotons->GetXaxis()->SetTitleOffset(1.2);
	D0KPiExtraPhotons->GetYaxis()->SetTitleOffset(1.7);

	D0KPiExtraPhotons->SetMarkerStyle(20);
	D0KPiExtraPhotons->SetMarkerSize(1);		
	D0KPiExtraPhotons->SetLineColor(kBlack);
	D0KPiExtraPhotons->SetMarkerColor(kBlack);
	D0KPiExtraPhotons->SetMinimum(0);

//	c->SetLogy();

	D0KPiExtraPhotons->Draw("ep");

	c->SaveAs("ExtraPlots/D0KPiExtraPhotons.png");




	TH1D * PhotonE = (TH1D *) fin->Get("PhotonE");
	PhotonE->GetXaxis()->SetTitle("Photon Energy (GeV)");
	PhotonE->GetYaxis()->SetTitle("Number of Photons");
	PhotonE->GetXaxis()->CenterTitle();
	PhotonE->GetYaxis()->CenterTitle();
	PhotonE->GetXaxis()->SetTitleOffset(1.2);
	PhotonE->GetYaxis()->SetTitleOffset(1.7);

	PhotonE->SetMarkerStyle(20);
	PhotonE->SetMarkerSize(1);		
	PhotonE->SetLineColor(kBlack);
	PhotonE->SetMarkerColor(kBlack);
	PhotonE->SetMinimum(0);

	PhotonE->Draw("hist");
//	c->SetLogy();

	c->SaveAs("ExtraPlots/PhotonE.png");


	TH1D * PiKMass = (TH1D *) fin->Get("InvMass_1");
	PiKMass->GetXaxis()->SetTitle("m_{#pi K} (GeV/c^{2})");
	PiKMass->GetYaxis()->SetTitle("Number of Events");
	PiKMass->GetXaxis()->CenterTitle();
	PiKMass->GetYaxis()->CenterTitle();
	PiKMass->GetXaxis()->SetTitleOffset(1.2);
	PiKMass->GetYaxis()->SetTitleOffset(1.7);

	PiKMass->SetMarkerStyle(20);
	PiKMass->SetMarkerSize(1);		
	PiKMass->SetLineColor(kBlack);
	PiKMass->SetMarkerColor(kBlack);
	PiKMass->SetMinimum(0);

	PiKMass->Draw("hist");
//	c->SetLogy();

	c->SaveAs("ExtraPlots/PiKMass.png");

}
