#pragma once

#include "ElectroIonReaction.h"
#include "ParticleInjector.h" 

namespace rad{

  using rad::consts::data_type::MC;

  /**
   * @class HepMCElectro
   * @brief HepMC3-specific configuration for Electro-Ion reactions.
   * * This class handles the mapping of raw HepMC tree branches to standard framework components.
   * It assumes a specific naming convention (e.g. `particles.momentum.m_v1`) common in HepMC3 ROOT files.
   */
  class HepMCElectro : public ElectroIonReaction {

  public:
    // --- Constructors ---
    HepMCElectro(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={}) 
      : ElectroIonReaction{treeName,fileNameGlob,columns} {}

    HepMCElectro(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} ) 
      : ElectroIonReaction{treeName,filenames,columns} {}

    /**
     * @brief Maps raw HepMC momentum branches to standard framework components.
     * * Pipeline:
     * 1. Registers the "MC" data type.
     * 2. Configures `ParticleInjector` with standard suffixes.
     * 3. Creates unified vectors (e.g., `mc_px`, `mc_phi`).
     */
    void SetupMC();

    /**
     * @brief Sets up aliases with additional filtering for stable particles.
     * * In addition to `SetupMC`, defines `mc_final_pid` which masks out unstable particles
     * (GenStatus != 1).
     */
    void DefineStableMomentumComponents();
  };


  // =================================================================================
  // IMPLEMENTATION
  // =================================================================================

  inline void HepMCElectro::SetupMC(){
      AddType(MC());
      
      ParticleInjector injector(this);
      
      // Define the target standardized suffixes
      injector.DefineParticleInfo({"px", "py", "pz", "m", "pid"});
      
      // Map the specific HepMC branches to these suffixes
      injector.AddSource(MC(), 
			 {"particles.momentum.m_v1","particles.momentum.m_v2","particles.momentum.m_v3","particles.mass","particles.pid"}); 
      
      // Execute the aliasing/creation logic
      injector.CreateUnifiedVectors();
 
      // Define standard kinematic variables for all types
      DefineForAllTypes("phi", Form("rad::ThreeVectorPhi(components_p3)"));
      DefineForAllTypes("theta", Form("rad::ThreeVectorTheta(components_p3)"));
      DefineForAllTypes("eta", Form("rad::ThreeVectorEta(components_p3)"));
      DefineForAllTypes("pmag", Form("rad::ThreeVectorMag(components_p3)"));
  }

  inline void HepMCElectro::DefineStableMomentumComponents(){
      SetupMC();
       
      // Map the status branch
      setBranchAlias("particles.status",MC()+"genStat");

      // Define 'final_pid' which retains PID only for stable particles (status == 1)
      Define(MC()+"final_pid",[](const ROOT::RVecI& pid,const ROOT::RVecI& stat){
	auto n = pid.size();
	auto final_pid = pid;
	// Vectorized masking: pid becomes 0 if stat != 1
	final_pid*=(stat==1);
	return final_pid;

      },{MC()+"pid",MC()+"genStat"});
  }

}
