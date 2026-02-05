#include "AnalysisManager.h" 
#include "HepMCElectro.h"        
#include "KinematicsProcElectro.h"        
#include "ReactionKinematics.h" 
#include "CommonDefines.h"
#include "RDFUtils.h"
#include <TBenchmark.h>

void CombiJpsi(){
   using namespace rad;
   ROOT::EnableImplicitMT(8);
 
   using Reaction = HepMCElectro;
   using Processor = KinematicsProcElectro;
   // =================================================================================
  // 1. SETUP & INJECTION
  // =================================================================================
   AnalysisManager<Reaction,Processor>  mgr{
    "Jpsi_Analysis", //name
    "hepmc3_tree",  //treename
    "jpac_jpsi_106_pomeron.root" //filename
  };
  mgr.SetOutputDir("output");
  auto& reaction = mgr.Reaction();
  
  reaction.SetupMC();
 
  // =================================================================================
  // 2. PARTICLE DEFINITIONS
  // =================================================================================
  reaction.SetBeamIonIndex(1);
  reaction.SetBeamElectronIndex(0);
  reaction.SetScatElectronIndex(2); 
  reaction.SetParticleIndex("ele", 4); 
  reaction.SetParticleIndex("pos", 5); 
  reaction.SetParticleIndex("p", 3);   

  reaction.MakeCombinations();
  rad::PrintDefinedColumnNames(reaction.CurrFrame());
  mgr.SetTypes(MC());

  // [A] SHARED KINEMATICS (Topology)
  // Applied to both Rec and Truth streams. Defines the decay chain.
  auto topology_recipe = [](Processor& p) {
    // 1. Reconstruct J/psi -> e+ e-
    p.Creator().Sum("Jpsi", {{"ele", "pos"}});         
    p.SetMesonParticles({"Jpsi"});
    p.SetBaryonParticles({"p"});
    
    // 3. Calculate Invariant Masses
    p.Mass("MassJ",     {"Jpsi"});             
    
    4. Calculate Mandelstam t (requires beam definition)
    p.RegisterCalc("tb", rad::physics::TBot);
    p.RegisterCalc("DeltaPhiJxP", rad::DeltaPhi,{{"Jpsi","p"}});

    p.ParticleTheta({{"scat_ele","Jpsi","p"}});
    p.ParticlePhi({{"scat_ele","Jpsi","p"}});
    p.ParticleP({{"scat_ele","Jpsi","p"}});
    
  };
  
  //HISTOGRAMS
  // Shared definitions for Rec and Truth.
  auto histogram_recipe = [](histo::Histogrammer& h) {
    h.Create("MassJ",     "Invariant Mass J/psi; Mass [GeV]", 100, 2.0, 4.0, "MassJ");
  };
  
  mgr.ConfigureKinematics(topology_recipe);
  mgr.ConfigureHistograms(histogram_recipe);
  mgr.Snapshot();

 
  gBenchmark->Start("analysis");
  mgr.Run();
  gBenchmark->Stop("analysis");
  gBenchmark->Print("analysis");

}
