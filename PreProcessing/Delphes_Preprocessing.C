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


  std::vector<Float_t> Tau_Jet_pt;
  std::vector<Float_t> Tau_Jet_eta;
  std::vector<Float_t> Tau_Jet_phi;
  std::vector<Float_t> Tau_Jet_mass;

  std::vector<Float_t> Tau_SubJet_pt;
  std::vector<Float_t> Tau_SubJet_phi;
  std::vector<Float_t> Tau_SubJet_eta;
  std::vector<Float_t> Tau_SubJet_mass;
  std::vector<Float_t> Tau_SubJet_px;
  std::vector<Float_t> Tau_SubJet_py;
  std::vector<Float_t> Tau_SubJet_pz;

  std::vector<Int_t> Tau_Number_of_Jets;
  std::vector<Int_t> Tau_Number_of_Subjets;




  std::vector<Float_t> Jet_pt;
  std::vector<Float_t> Jet_eta;
  std::vector<Float_t> Jet_phi;
  std::vector<Float_t> Jet_mass;




  std::vector<Float_t> B_Jet_pt;
  std::vector<Float_t> B_Jet_eta;
  std::vector<Float_t> B_Jet_phi;
  std::vector<Float_t> B_Jet_mass;

  std::vector<Float_t> B_SubJet_pt;
  std::vector<Float_t> B_SubJet_phi;
  std::vector<Float_t> B_SubJet_eta;
  std::vector<Float_t> B_SubJet_mass;
  std::vector<Float_t> B_SubJet_px;
  std::vector<Float_t> B_SubJet_py;
  std::vector<Float_t> B_SubJet_pz;


  std::vector<Int_t> B_Number_of_Jets;
  std::vector<Int_t> B_Number_of_Subjets;




  std::vector<Int_t> n_Jets;
  std::vector<Int_t> Number_of_Jets;
  std::vector<Int_t> Number_of_Subjets;


  newtree->Branch("Tau_Jet_pt",  &Tau_Jet_pt);
  newtree->Branch("Tau_Jet_phi",  &Tau_Jet_phi);
  newtree->Branch("Tau_Jet_eta",  &Tau_Jet_eta);
  newtree->Branch("Tau_Jet_mass",  &Tau_Jet_mass);

  newtree->Branch("Tau_SubJet_pt",  &Tau_SubJet_pt);
  newtree->Branch("Tau_SubJet_phi",  &Tau_SubJet_phi);
  newtree->Branch("Tau_SubJet_eta",  &Tau_SubJet_eta);
  newtree->Branch("Tau_SubJet_mass",  &Tau_SubJet_mass);

  newtree->Branch("Tau_SubJet_px",  &Tau_SubJet_px);
  newtree->Branch("Tau_SubJet_py",  &Tau_SubJet_py);
  newtree->Branch("Tau_SubJet_pz",  &Tau_SubJet_pz);

  newtree->Branch("Tau_Number_of_Jets",  &Tau_Number_of_Jets);
  newtree->Branch("Tau_Number_of_Subjets",  &Tau_Number_of_Subjets);





  newtree->Branch("B_Jet_pt",  &B_Jet_pt);
  newtree->Branch("B_Jet_phi",  &B_Jet_phi);
  newtree->Branch("B_Jet_eta",  &B_Jet_eta);
  newtree->Branch("B_Jet_mass",  &B_Jet_mass);

  newtree->Branch("B_SubJet_pt",  &B_SubJet_pt);
  newtree->Branch("B_SubJet_phi",  &B_SubJet_phi);
  newtree->Branch("B_SubJet_eta",  &B_SubJet_eta);
  newtree->Branch("B_SubJet_mass",  &B_SubJet_mass);

  newtree->Branch("B_SubJet_px",  &B_SubJet_px);
  newtree->Branch("B_SubJet_py",  &B_SubJet_py);
  newtree->Branch("B_SubJet_pz",  &B_SubJet_pz);

  newtree->Branch("B_Number_of_Jets",  &B_Number_of_Jets);
  newtree->Branch("B_Number_of_Subjets",  &B_Number_of_Subjets);



  newtree->Branch("Jet_pt",  &Jet_pt);
  newtree->Branch("Jet_phi",  &Jet_phi);
  newtree->Branch("Jet_eta",  &Jet_eta);
  newtree->Branch("Jet_mass",  &Jet_mass);



  newtree->Branch("Number_of_Jets",  &Number_of_Jets);
  newtree->Branch("n_Jets",  &n_Jets);
  newtree->Branch("Number_of_Subjets",  &Number_of_Subjets);



  for (Long64_t entry = 0; entry < numberOfEntries; ++entry) {
    if (entry % 2000 == 0) {
      cout << "Processing event " << entry << endl;
    }
    treeReader->ReadEntry(entry);
    bool filled = false;
    bool filledTau = false;
    int numTauJets = 0;
    int numJets = 0;
    Int_t Number_of_Tau_Jets = 0;
    Int_t Number_of_B_Jets = 0;

    TLorentzVector jetMom[5]; // Example TLorentzVector array

    Int_t nJets = 0;
    // cout << "branchJet: "<< branchJet->GetEntries() << endl;
    Int_t numberJets = branchJet->GetEntries();
    Number_of_Jets.push_back(numberJets);
    for (int i = 0; i < branchJet->GetEntries(); ++i) {
      nJets +=1;

      Jet *jet = (Jet*) branchJet->At(i);
      // if (!jet or jet->TauTag != 1) continue;
      if (!jet) continue;


      // Iterate over the TLorentzVector array
      Int_t Tau_NSubJets = 0;
      if (jet->TauTag == 1 && jet->BTag == 0) {

        Number_of_Tau_Jets++;

        float_t Tau_Jetpt = jet->PT;
        float_t Tau_Jetphi = jet->Phi;
        float_t Tau_Jeteta = jet->Eta;
        float_t Tau_Jetmass = jet->Mass;



        Tau_Jet_pt.push_back(Tau_Jetpt);
        Tau_Jet_phi.push_back(Tau_Jetphi);
        Tau_Jet_eta.push_back(Tau_Jeteta);
        Tau_Jet_mass.push_back(Tau_Jetmass);

        for (Int_t subjet = 0; subjet < 5; ++subjet) {

          TLorentzVector& Tau_subjetMom = jet->SoftDroppedP4[subjet];
          
          // Access the transverse momentum (Pt) directly from TLorentzVector
          Float_t Tau_subjetPt = Tau_subjetMom.Pt();
          if (Tau_subjetPt != 0) {Tau_NSubJets++;}



          Float_t Tau_subjetPhi = Tau_subjetMom.Phi();
          Float_t Tau_subjetEta = Tau_subjetMom.Eta();
          Float_t Tau_subjetMass = Tau_subjetMom.M();

          Float_t Tau_subjetPx = Tau_subjetMom.Px();
          Float_t Tau_subjetPy = Tau_subjetMom.Py();
          Float_t Tau_subjetPz = Tau_subjetMom.Pz();


          // Print the transverse momentum of the subjet
          // cout << "Subjet " << subjet << " PT: " << subjetPt << " Phi:" << subjetPhi << " Eta: " << subjetEta << endl;
          // cout << " Mass: " << subjetMass << " Px: " << subjetPx << " Py:" << subjetPy << " Pz:" << subjetPz << endl;

          if (Tau_subjetPt != 0){
            Tau_SubJet_pt.push_back(Tau_subjetPt);
            Tau_SubJet_phi.push_back(Tau_subjetPhi);
            Tau_SubJet_eta.push_back(Tau_subjetEta);
            Tau_SubJet_mass.push_back(Tau_subjetMass);
            Tau_SubJet_px.push_back(Tau_subjetPx);
            Tau_SubJet_py.push_back(Tau_subjetPy);
            Tau_SubJet_pz.push_back(Tau_subjetPz);
          }


        // Perform more analysis based on jet properties
        // Add your analysis code here
        }
        Tau_Number_of_Subjets.push_back(Tau_NSubJets);
      }
      Int_t B_NSubJets = 0;
      if (jet->BTag == 1 && jet->TauTag == 0) {

        Number_of_B_Jets++;

        float_t B_Jetpt = jet->PT;
        float_t B_Jetphi = jet->Phi;
        float_t B_Jeteta = jet->Eta;
        float_t B_Jetmass = jet->Mass;



        B_Jet_pt.push_back(B_Jetpt);
        B_Jet_phi.push_back(B_Jetphi);
        B_Jet_eta.push_back(B_Jeteta);
        B_Jet_mass.push_back(B_Jetmass);

        for (Int_t subjet = 0; subjet < 5; ++subjet) {

          TLorentzVector& B_subjetMom = jet->SoftDroppedP4[subjet];
          
          // Access the transverse momentum (Pt) directly from TLorentzVector
          Float_t B_subjetPt = B_subjetMom.Pt();
          if (B_subjetPt != 0) {B_NSubJets++;}

          Float_t B_subjetPhi = B_subjetMom.Phi();
          Float_t B_subjetEta = B_subjetMom.Eta();
          Float_t B_subjetMass = B_subjetMom.M();

          Float_t B_subjetPx = B_subjetMom.Px();
          Float_t B_subjetPy = B_subjetMom.Py();
          Float_t B_subjetPz = B_subjetMom.Pz();


          // Print the transverse momentum of the subjet
          // cout << "Subjet " << subjet << " PT: " << subjetPt << " Phi:" << subjetPhi << " Eta: " << subjetEta << endl;
          // cout << " Mass: " << subjetMass << " Px: " << subjetPx << " Py:" << subjetPy << " Pz:" << subjetPz << endl;

          if (B_subjetPt != 0){
            B_SubJet_pt.push_back(B_subjetPt);
            B_SubJet_phi.push_back(B_subjetPhi);
            B_SubJet_eta.push_back(B_subjetEta);
            B_SubJet_mass.push_back(B_subjetMass);
            B_SubJet_px.push_back(B_subjetPx);
            B_SubJet_py.push_back(B_subjetPy);
            B_SubJet_pz.push_back(B_subjetPz);
          }



        // Perform more analysis based on jet properties
        // Add your analysis code here
        }
        B_Number_of_Subjets.push_back(B_NSubJets);

      }


      // cout << "Number of Subjets: " << NSubJets << endl;
      // Tau_Number_of_Subjets.push_back(Tau_NSubJets);
      // B_Number_of_Subjets.push_back(B_NSubJets);

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

      float_t Jetpt = jet->PT;
      float_t Jetphi = jet->Phi;
      float_t Jeteta = jet->Eta;
      float_t Jetmass = jet->Mass;

      Jet_pt.push_back(Jetpt);
      Jet_phi.push_back(Jetphi);
      Jet_eta.push_back(Jeteta);
      Jet_mass.push_back(Jetmass);


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
      // cout << "nJets: " << nJets << endl;
      n_Jets.push_back(nJets);
      Tau_Number_of_Jets.push_back(Number_of_Tau_Jets);
      B_Number_of_Jets.push_back(Number_of_B_Jets);
    }
    // if(filledTau) fprintf(fout,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%i\n", tau1_pt, tau1_eta, tau1_phi, tau2_pt, tau2_eta, tau2_phi, tau1_m, tau2_m, m_tau1tau2, isSig);
    newtree->Fill();
    // Jet_pt.clear();
    tau1_m_tree.clear();

    Jet_pt.clear();
    Jet_phi.clear();
    Jet_eta.clear();
    Jet_mass.clear();


    Tau_Jet_pt.clear();
    Tau_Jet_phi.clear();
    Tau_Jet_eta.clear();
    Tau_Jet_mass.clear();


    B_Jet_pt.clear();
    B_Jet_phi.clear();
    B_Jet_eta.clear();
    B_Jet_mass.clear();



    Tau_SubJet_pt.clear();
    Tau_SubJet_phi.clear();
    Tau_SubJet_eta.clear();
    Tau_SubJet_mass.clear();
    Tau_SubJet_px.clear();
    Tau_SubJet_py.clear();
    Tau_SubJet_pz.clear();

    B_SubJet_pt.clear();
    B_SubJet_phi.clear();
    B_SubJet_eta.clear();
    B_SubJet_mass.clear();
    B_SubJet_px.clear();
    B_SubJet_py.clear();
    B_SubJet_pz.clear();

    n_Jets.clear();
    Number_of_Jets.clear();
    B_Number_of_Subjets.clear();
    Tau_Number_of_Subjets.clear();

    Tau_Number_of_Jets.clear();
    B_Number_of_Jets.clear();
  }
  newtree->Write("", TObject::kOverwrite);

// return 1;
}

