#pragma once

#include "KinematicsProcessor.h"
#include "ElectronScatterKinematics.h"

namespace rad {

  /**
   * @brief Physics Processor specialized for Electro-production reactions (e.g. ep -> e' X).
   * * This class extends the generic engine to define the standard Electron Scattering topology:
   * 1. Defines a "Scattered Electron" group (containing the 'scat_ele' particle).
   * 2. Automatically creates the Virtual Photon (Beam - Scattered).
   * 3. Provides Q2, CosThetaCM, and PhiCM shortcuts.
   */
  class KinematicsProcElectro : public KinematicsProcessor {
  public:
    
    /**
     * @brief Standard Constructor.
     * @param cr The configuration interface.
     * @param prefix The data type prefix (e.g., "rec_", "tru_").
     * @param suffix Optional suffix for output columns (e.g. "_miss").
     */
    KinematicsProcElectro(ConfigReaction* cr, const std::string& prefix, const std::string& suffix = "");

    /**
     * @brief Copy Constructor (The Fork).
     * Creates a new processor inheriting the configuration of the other, but with a new suffix.
     * Use this to create a "Missing Neutron" hypothesis from a standard processor.
     */
    KinematicsProcElectro(const KinematicsProcElectro& other, const std::string& new_suffix);

    // =================================================================================
    // Electro-Specific Shortcuts
    // =================================================================================

    /** @brief Calculates Q2 (Negative Squared 4-Momentum Transfer). */
    void Q2();
    
    /** @brief Calculates Cos(Theta) in the Center of Mass frame. */
    void CosThetaCM();
    
    /** @brief Calculates Phi (Azimuthal Angle) in the Center of Mass frame. */
    void PhiCM();
  };

  // =================================================================================
  // IMPLEMENTATION
  // =================================================================================

  inline KinematicsProcElectro::KinematicsProcElectro(ConfigReaction* cr, const std::string& prefix, const std::string& suffix) 
      : KinematicsProcessor(cr, prefix, suffix) 
  {
      // 1. Define Topology: Create the "Scattered Electron Group" in RDF
      // We must use the specific prefix (e.g. "rec_") 
      Reaction()->setGroupParticles(consts::ScatGroup(), GetPrefix(), {consts::ScatEle()});

      // 2. Register this group with the Creator
      // Creator() now knows that "scat_ele" logic relies on this RDF group.
      Creator().DefineGroup(consts::OrderScatEle(), consts::ScatGroup());

      // 3. Define the Virtual Photon: Beam Electron - Scattered Electron
      // This is a "Created" particle that exists in the combinatorial arrays.
      Creator().Diff(consts::VirtGamma(), {{consts::BeamEle()}, {consts::ScatEle()}});
  }

  inline KinematicsProcElectro::KinematicsProcElectro(const KinematicsProcElectro& other, const std::string& new_suffix)
      : KinematicsProcessor(other, new_suffix) 
  {
      // The Base Class copy constructor handles copying the ParticleCreator state 
      // (including the VirtGamma definition and Group mappings).
  }

  inline void KinematicsProcElectro::Q2() {
       KineCalculation Q2(*this, "Q2", physics::ElS_Q2);
  }
    
  inline void KinematicsProcElectro::CosThetaCM() {
       KineCalculation CosTheta(*this, "CosThetaCM", physics::ElS_CosThetaCM);
  }
    
  inline void KinematicsProcElectro::PhiCM() {
       KineCalculation Phi(*this, "PhiCM", physics::ElS_PhiCM);
  }

} // namespace rad
