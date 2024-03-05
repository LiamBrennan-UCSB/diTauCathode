#include <iostream>
#include <cstring>
#include <vector>
#include <cmath>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <cstdlib>
//#include "ExRootTreeReader.h"
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "../../classes/DelphesClasses.h"
#include "../../external/ExRootAnalysis/ExRootTreeReader.h"
#include "../../external/ExRootAnalysis/ExRootResult.h"
#endif

/*
 *  example running:
 *   root -l create_dataset.C
 *    */

// class TFile;

using namespace std;

void Delphes_Preprocessing(string file_n, int label, TString newfilename) {


  int debug = 0;
  gSystem->Load("libDelphes");
  int isSig = label;

  char infile[200], outfile[200];
  // string csv_path = "/uscms/home/lbrennan/work/CATHODE/delphes/";

  string in_path = "root://cmseos.fnal.gov//store/user/tvami/diTauCathode/";
  string file_name = file_n.c_str();
  sprintf(infile,"%s%s.root",in_path.c_str(),file_name.c_str());
  TFile * fin = TFile::Open(infile);
  // FILE *fout;
  // sprintf(outfile,"%s/diTauCathode/src/%s.csv",csv_path.c_str(),file_name.c_str());
  // fout = fopen(outfile, "w");
  
  cout << "Sample used is " << file_name.c_str() << endl;
  

  TChain chain("Delphes");
  chain.Add(infile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  // TClonesArray *branchsubJet = treeReader->UseBranch("Jet.SoftDroppedP4");

  float tau1_pt, tau1_eta, tau1_phi, tau2_pt, tau2_eta, tau2_phi, tau1_m, tau2_m, m_tau1tau2;
  
  
  cout << "Running on " << numberOfEntries << " events" << endl;


  auto new_f = std::make_unique<TFile>(newfilename,"RECREATE");
  new_f->SetCompressionAlgorithm(ROOT::kLZMA);
  auto newtree = std::make_unique<TTree>("Events", "ntuple from nanoAOD"); //mySnapshot or Events


  // std::vector<Int_t> Jet_pt;
  std::vector<Float_t> tau1_m_tree;
  // newtree->Branch("Jet_pt",  &Jet_pt);
  newtree->Branch("tau1_m_tree",  &tau1_m_tree);


  std::vector<Float_t> SubJet_pt;
  std::vector<Float_t> SubJet_phi;
  std::vector<Float_t> SubJet_eta;
  std::vector<Float_t> SubJet_mass;
  std::vector<Float_t> SubJet_px;
  std::vector<Float_t> SubJet_py;
  std::vector<Float_t> SubJet_pz;

  newtree->Branch("SubJet_pt",  &SubJet_pt);
  newtree->Branch("SubJet_phi",  &SubJet_phi);
  newtree->Branch("SubJet_eta",  &SubJet_eta);
  newtree->Branch("SubJet_mass",  &SubJet_mass);

  newtree->Branch("SubJet_px",  &SubJet_px);
  newtree->Branch("SubJet_py",  &SubJet_py);
  newtree->Branch("SubJet_pz",  &SubJet_pz);


  newtree->Branch("Number_of_Subjets",  &Number_of_Subjets);


  for (Long64_t entry = 0; entry < numberOfEntries; ++entry) {
    // Jet_pt.push_back(1);
    if (entry % 2000 == 0) {
      cout << "Processing event " << entry << endl;
      // Jet_pt.push_back(entry);
    }
    treeReader->ReadEntry(entry);
    bool filled = false;
    bool filledTau = false;
    int numTauJets = 0;
    int numJets = 0;


    TLorentzVector jetMom[5]; // Example TLorentzVector array

    for (int i = 0; i < branchJet->GetEntries(); ++i) {


      Jet *jet = (Jet*) branchJet->At(i);
      if (!jet or jet->TauTag != 1) continue;


      // Iterate over the TLorentzVector array
      int NSubJets = 0;
      for (Int_t subjet = 0; subjet < 5; ++subjet) {

        TLorentzVector& subjetMom = jet->SoftDroppedP4[subjet];
        
        // Access the transverse momentum (Pt) directly from TLorentzVector
        Float_t subjetPt = subjetMom.Pt();
        if (subjetPt != 0) {NSubJets++;}

        Float_t subjetPhi = subjetMom.Phi();
        Float_t subjetEta = subjetMom.Eta();
        Float_t subjetMass = subjetMom.M();

        Float_t subjetPx = subjetMom.Px();
        Float_t subjetPy = subjetMom.Py();
        Float_t subjetPz = subjetMom.Pz();


        // Print the transverse momentum of the subjet
        std::cout << "Subjet " << subjet << " PT: " << jetPt << " Phi:" << jetPhi " Eta:" << subjetEta << " Mass: " << subjetMass << " Px:" << subjetPx << " Py:" << subjetPy << " Pz:" << subjetPz << std::endl;

        SubJet_pt.push_back(subjetPt);
        SubJet_phi.push_back(subjetPhi);
        SubJet_eta.push_back(subjetEta);
        SubJet_mass.push_back(subjetMass);
        SubJet_px.push_back(subjetPx);
        SubJet_py.push_back(subjetPy);
        SubJet_pz.push_back(subjetPz);



        // Perform more analysis based on jet properties
        // Add your analysis code here
      }

      cout << "Number of Subjets: " << NSubJets << endl;
      Number_of_Subjets.push_back(NSubJets)

      if (jet->TauTag == 1) {
        ++numTauJets;
        tau1_pt = jet->PT;
        tau1_eta = jet->Eta;
        tau1_phi = jet->Phi;
        tau1_m = (jet->P4()).M();
        // cout << "tau1_m " << tau1_m << endl;
        tau1_m_tree.push_back(tau1_m);
        // cout << "tau1_m_tree: " << len(tau1_m_tree) << endl;


      }

      for (int j = 0; j < branchJet->GetEntries(); ++j) {
        auto* jet2 = static_cast<Jet*>(branchJet->At(j));
        if (!jet2 || jet == jet2) continue;


        if (jet2->TauTag == 1) {
          if (!filledTau) {
            m_tau1tau2 = (jet->P4() + jet2->P4()).M();
            tau2_pt = jet2->PT;
            tau2_eta = jet2->Eta;
            tau2_phi = jet2->Phi;
            tau2_m = (jet2->P4()).M();
              }
          filledTau = true;
        }
      }
    }
    // if(filledTau) fprintf(fout,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%i\n", tau1_pt, tau1_eta, tau1_phi, tau2_pt, tau2_eta, tau2_phi, tau1_m, tau2_m, m_tau1tau2, isSig);
    newtree->Fill();
    // Jet_pt.clear();
    tau1_m_tree.clear();

    SubJet_pt.clear();
    SubJet_phi.clear();
    SubJet_eta.clear();
    SubJet_mass.clear();
    SubJet_px.clear();
    SubJet_py.clear();
    SubJet_pz.clear();

  }
  newtree->Write("", TObject::kOverwrite);

// return 1;
}

