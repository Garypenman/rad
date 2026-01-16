#include "CommonDefines.h"
#include "HepMCElectro.h"
#include "KinematicsProcElectro.h"
#include "KineCalculation.h"
#include "Indicing.h"
#include "ElectronScatterKinematics.h"
#include <TBenchmark.h>

/**
 * @brief Example Script: Z_c(3900) Analysis with Combinatorics and Missing Mass.
 * Updated to use SetMesonParticles shortcut and CloneLinked.
 */
void ProcessHepMCZCombi(){
  
  using namespace rad;
  using namespace rad::consts::data_type; 

  gBenchmark->Start("df");

  // =================================================================================
  // 1. SETUP & INJECTION
  // =================================================================================
  HepMCElectro hepmc{
      "hepmc3_tree", 
      "/home/dglazier/Dropbox/EIC/EventGenerators/elSpectro/examples/out/jpac_z3900_10x100.hepmc.root"
  };
  hepmc.SetupMC();
    
  // =================================================================================
  // 2. PARTICLE DEFINITIONS
  // =================================================================================
  hepmc.SetBeamIonIndex(1);
  hepmc.SetBeamElectronIndex(0);
  hepmc.SetScatElectronCandidates({6, 2}); 
  hepmc.SetParticleCandidates("ele", {2, 6}); 
  hepmc.SetParticleIndex("pos", 3); 
  hepmc.SetParticleIndex("pip", 4); 
  hepmc.SetParticleIndex("n", 5);   

  hepmc.MakeCombinations();

  // =================================================================================
  // 3. KINEMATICS PROCESSOR (Standard Topology)
  // =================================================================================
  KinematicsProcElectro kine{&hepmc, MC()}; 

  // A. Define Particles FIRST (Topology Construction)
  kine.Creator().Sum("Jpsi", {{"ele", "pos"}});       
  kine.Creator().Sum("Z",    {{"Jpsi", "pip"}});      
  
  // Missing Mass: n_calc = (Beam + Target) - (ScatEle + Jpsi + pi+)
  kine.Creator().Diff("n_calc", 
      {{consts::BeamIon(), consts::BeamEle()},  
       {"Jpsi", "pip", consts::ScatEle()}}      
  );

  // B. Define Groups NEXT (Lazy Definition)
  // This uses the Processor's SetGroup wrapper, ensuring "Z" exists before the group is built.
  kine.SetMesonParticles({"Z"});
  kine.SetBaryonParticles({"n"});

  // =================================================================================
  // 4. CALCULATIONS (Registered on Master)
  // =================================================================================
  // These will be calculated for Master AND automatically copied to the Linked clone
  
  kine.Q2();          
  kine.CosThetaCM(); 
  kine.PhiCM();       
  kine.Mass("ZMass", {"Z"});
  
  kine.Mass2("MissMass2", 
      {consts::BeamIon(), consts::BeamEle()}, 
      {"Jpsi", "pip", consts::ScatEle(), "n"}
  );

  // Custom Calc (t_prime)
  kine.RegisterCalc("tp", rad::physics::TPrimeTop);

  // =================================================================================
  // 5. LINKED PROCESSOR (Cloned Hypothesis)
  // =================================================================================
  
  // 1. Clone: Copies all the calculations registered above (Q2, Phi, Mass, tp...)
  //    Creates a new processor for the "mc_" stream but with suffix "_ncalc"
  // auto kine_miss = kine.CloneLinked("_ncalc");

  // // 2. Customize: Override "Baryons" to be "n_calc" (missing neutron) for this hypothesis
  // kine_miss->SetBaryonParticles({"n_calc"});

  // =================================================================================
  // 6. INITIALIZATION & EXECUTION
  // =================================================================================

  // Init Master: Executes definitions and calculations for Standard Topology
  kine.Init();

  // Init Linked: Binds to 'kine' topology and executes copied calculations for Missing Topology
  //kine_miss->InitLinked(kine);

  gBenchmark->Start("snapshot");
  hepmc.Snapshot("HepMCZCombi.root");
  gBenchmark->Stop("snapshot");
  gBenchmark->Print("snapshot");
}
