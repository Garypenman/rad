#include "CommonDefines.h"
#include "HepMCElectro.h"
#include "KinematicsProcElectro.h"
#include "KineCalculation.h"
#include "Indicing.h"
#include "ElectronScatterKinematics.h"
#include "gammaN_2_Spin0Spin0SpinHalf.h"
#include <TBenchmark.h>

/**
 * @brief Example Script: TCS Analysis with Combinatorics and Missing Mass.
 * Updated to use SetMesonParticles shortcut and CloneLinked.
 */
void ProcessHepMCTCSCombi(){
  
  using namespace rad;
  using namespace rad::consts::data_type; 

  gBenchmark->Start("df");

  // =================================================================================
  // 1. SETUP & INJECTION
  // =================================================================================
  HepMCElectro hepmc{
      "hepmc3_tree", 
      "/w/work5/home/garyp/eic/Farm/data/EpIC_DDVCS_ee_18x275/rootfiles/18x275_ddvcs_edecay_hplus.root"
  };
  hepmc.SetupMC();
    
  // =================================================================================
  // 2. PARTICLE DEFINITIONS
  // =================================================================================
  hepmc.SetBeamIonIndex(3);
  hepmc.SetBeamElectronIndex(0);
  //hepmc.SetScatElectronCandidates({1, 6}); 
  //hepmc.SetParticleCandidates("ele", {1, 6}); 
  hepmc.SetScatElectronIndex(1);
  hepmc.SetParticleIndex("ele",6);
  hepmc.SetParticleIndex("pos", 7); 
  hepmc.SetParticleIndex("pprime", 5);   

  hepmc.MakeCombinations();

  // =================================================================================
  // 3. KINEMATICS PROCESSOR (Standard Topology)
  // =================================================================================
  KinematicsProcElectro kine{&hepmc, MC()}; 

  // A. Define Particles FIRST (Topology Construction)
  kine.Creator().Sum("gprime", {{"ele", "pos"}});       
  
  // Missing Mass: n_calc = (Beam + Target) - (ScatEle + Jpsi + pi+)
  kine.Creator().Diff("pprime_calc",{
        {consts::BeamIon(), consts::BeamEle()},  
	{"gprime", consts::ScatEle()}
    });

  // B. Define Groups NEXT (Lazy Definition)
  // This uses the Processor's SetGroup wrapper, ensuring "Z" exists before the group is built.
  kine.SetMesonParticles({"ele","pos"});
  kine.SetBaryonParticles({"pprime"});

  // =================================================================================
  // 4. CALCULATIONS (Registered on Master)
  // =================================================================================
  // These will be calculated for Master AND automatically copied to the Linked clone
  
  kine.Q2();         
  kine.xbj();
  kine.y();         
  kine.nu();
  kine.tau();
  kine.tauprime();

  kine.RegisterCalc("GammaPol",rad::physics::ElS_PolVirtPhot);
  kine.RegisterCalc("GammaPolCirc",rad::physics::ElS_CircPolVirtPhot);

  kine.RegisterCalc("CosThetaHel",rad::gn2s0s0s12::CosThetaHel);
  kine.RegisterCalc("ThetaHel",rad::gn2s0s0s12::ThetaHel);
  kine.RegisterCalc("PhiHel",rad::gn2s0s0s12::PhiHel);
  
  kine.CosThetaCM(); 
  kine.PhiCM();       
  
  kine.Mass("GMass", {"gprime"});
  kine.Mass2("MissMass2", 
	     {consts::BeamIon(), consts::BeamEle()}, 
	     {"gprime", consts::ScatEle(), "pprime"});
  
  
  // Custom Calc (t_prime)
  kine.RegisterCalc("t_top", rad::physics::TTop);
  kine.RegisterCalc("t_bot", rad::physics::TBot);
  kine.RegisterCalc("tp_top", rad::physics::TPrimeTop);
  kine.RegisterCalc("tp_bot", rad::physics::TPrimeBot);

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
  hepmc.Snapshot("/w/work5/home/garyp/combirad_trees/HepMC_ddvcs_ee_18x275_hplus.root");
  gBenchmark->Stop("snapshot");
  gBenchmark->Print("snapshot");
}
