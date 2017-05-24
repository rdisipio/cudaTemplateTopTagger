#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
//#include "external/ExRootAnalysis/ExRootResult.h"
#endif

struct MyEvent {
  Int_t           eventNumber;

  vector< float > mc_pt;
  vector< float > mc_eta;
  vector< float > mc_phi;
  vector< float > mc_E;
  vector< float > mc_m;
  vector< int >   mc_pdgId;
  vector< int >   mc_status;

  vector< float > jet_pt;
  vector< float > jet_eta;
  vector< float > jet_phi;
  vector< float > jet_E;
  vector< float > jet_m;
  vector< float > jet_tau32;

  vector< int >   cl_matched_jet;
  vector< float > cl_pt;
  vector< float > cl_eta;
  vector< float > cl_phi;
  vector< float > cl_E;
};

void read_delphes( const char * inputFile )
{
   gSystem->Load("libDelphes");
   TChain chain("Delphes");
   chain.Add(inputFile);

   ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
   Long64_t numberOfEntries = treeReader->GetEntries();

//   TClonesArray *branchEvent = treeReader->UseBranch("Event");
   TClonesArray *branchJet = treeReader->UseBranch("Jet");
   TClonesArray *branchParticle = treeReader->UseBranch("Particle");
   TClonesArray *branchTower = treeReader->UseBranch("Tower");

   // create output file
   MyEvent event;
   TFile * outfile = TFile::Open( "output.root", "RECREATE" );
   TTree * outtree = new TTree( "physics", "Physics" );

   outtree->Branch( "eventNumber", &event.eventNumber, "eventNumber/I" );

   outtree->Branch( "jet_pt",  &event.jet_pt );
   outtree->Branch( "jet_eta", &event.jet_eta );
   outtree->Branch( "jet_phi", &event.jet_phi );
   outtree->Branch( "jet_E",   &event.jet_E );
   outtree->Branch( "jet_m",   &event.jet_m );
   outtree->Branch( "jet_tau32", &event.jet_tau32 );

   outtree->Branch( "cl_matched_jet", &event.cl_matched_jet );
   outtree->Branch( "cl_pt",  &event.cl_pt );
   outtree->Branch( "cl_eta", &event.cl_eta );
   outtree->Branch( "cl_phi", &event.cl_phi );
   outtree->Branch( "cl_E",   &event.cl_E );

   outtree->Branch( "mc_pt",  &event.mc_pt );
   outtree->Branch( "mc_eta", &event.mc_eta );
   outtree->Branch( "mc_phi", &event.mc_phi );
   outtree->Branch( "mc_E",   &event.mc_E );
   outtree->Branch( "mc_m",   &event.mc_m );
   outtree->Branch( "mc_pdgId", &event.mc_pdgId );
   outtree->Branch( "mc_status", &event.mc_status );

   for(Int_t entry = 0; entry < numberOfEntries; ++entry ) {
     treeReader->ReadEntry(entry);

 //    event.eventNumber = branchEvent->Number;
     event.eventNumber = entry;

     event.mc_pt.clear();
     event.mc_phi.clear();
     event.mc_eta.clear();
     event.mc_E.clear();
     event.mc_m.clear();
     event.mc_pdgId.clear();
     event.mc_status.clear();

     // JETS
     const int jets_n = branchJet->GetEntries();

     event.jet_pt.resize(jets_n);
     event.jet_eta.resize(jets_n);
     event.jet_phi.resize(jets_n);
     event.jet_E.resize(jets_n);
     event.jet_m.resize(jets_n);
     event.jet_tau32.resize(jets_n);

     event.cl_matched_jet.clear();
     event.cl_pt.clear();
     event.cl_eta.clear();
     event.cl_phi.clear();
     event.cl_E.clear();

     Int_t n_good_jets = 0;
     for( Int_t i = 0 ; i < jets_n ; ++i ) {
        Jet *jet = (Jet*) branchJet->At( i );
        TLorentzVector p = jet->P4();

        if( p.Pt() < 200. ) continue;
        if( fabs(p.Eta()) > 2.0 ) continue;
        if( fabs( p.M() - 172.5 ) > 50. ) continue;

        n_good_jets++;

        event.jet_pt[i]  = p.Pt();
        event.jet_phi[i] = p.Phi();
        event.jet_eta[i] = p.Eta();
        event.jet_E[i]   = p.E();
        event.jet_m[i]   = p.M();

        event.jet_tau32[i] = ( jet->Constituents.GetEntriesFast() > 3 ) ?
                               jet->Tau[2] / jet->Tau[1] : 0.;

        for( Int_t j = 0; j < jet->Constituents.GetEntriesFast(); ++j) {
           TObject *object = jet->Constituents.At(j);
           if( object == 0 ) continue;

           if( object->IsA() != Tower::Class() ) continue;

           Tower * tower = (Tower*)object;

           TLorentzVector p = tower->P4();
           event.cl_matched_jet.push_back( i );
           event.cl_pt.push_back( p.Pt() );
           event.cl_eta.push_back( p.Eta() );
           event.cl_phi.push_back( p.Phi() );
           event.cl_E.push_back( p.E() );
        } // end loop clusters

     } // end loop jets

     if( n_good_jets == 0 ) continue;

     // Look for top quarks
     const int mc_n = branchParticle->GetEntriesFast();
     for( Int_t i = 0 ; i < mc_n ; ++i ) {
        if( i > 200 ) continue;

        GenParticle * particle = (GenParticle*) branchParticle->At(i);

        const int pid = particle->PID;
        const int apid = abs(pid);
        if( apid != 6 ) continue; 
        const int status = particle->Status;

        TLorentzVector p = particle->P4();
        event.mc_pt.push_back( p.Pt() );
        event.mc_eta.push_back( p.Eta() );
        event.mc_phi.push_back( p.Phi() );
        event.mc_E.push_back( p.E() );
        event.mc_m.push_back( p.M() );
        event.mc_pdgId.push_back( pid );
        event.mc_status.push_back( status );
     }

     // Finally, save event
     outtree->Fill();
   }

   outtree->Write();
   outfile->Close();
}
