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
     * @param suffix Optional suffix for output columns (e.g. "_miss").
     */
    KinematicsProcElectro(config::ConfigReaction* cr, const std::string& suffix = "") 
      : KinematicsProcessor(cr, suffix) 
    {
      // 1. Define Topology: Create the "Scattered Electron Group" in RDF
      // We take the single particle name "scat_ele" and put it in a group "scat_ele_group".
      // This allows ParticleCreator to treat it uniformly as a group lookup.
      Reaction()->setGroupParticles(names::ScatGroup(), {names::ScatEle()});

      // 2. Register this group with the Creator
      // Tells the engine: "The ScatEle slot (index 4) is filled by the group 'scat_ele_group'"
      Creator().DefineGroup(names::OrderScatEle(), names::ScatGroup());

      // 3. Define Topology: Create the Virtual Photon
      // Rule: VirtGamma = BeamEle - ScatEle
      Creator().Diff(names::VirtGamma(), {{names::BeamEle()}, {names::ScatEle()}});
    }

    /**
     * @brief Copy Constructor (The Fork).
     * Creates a new processor inheriting the configuration of the other, but with a new suffix.
     * Use this to create a "Missing Neutron" hypothesis from a standard processor.
     */
    KinematicsProcElectro(const KinematicsProcElectro& other, const std::string& new_suffix)
      : KinematicsProcessor(other, new_suffix) 
    {
        // The Base Class copy constructor handles copying the ParticleCreator state 
        // (including the VirtGamma definition and Group mappings), so no extra setup is needed here.
    }

    // =================================================================================
    // Electro-Specific Shortcuts
    // =================================================================================

    /**
     * @brief Calculates Q2 (Negative Squared 4-Momentum Transfer).
     * Uses the ElS_Q2 kernel.
     */
    void Q2() {
       KineCalculation Q2(*this, "Q2", rad::ElS_Q2);
    }
    
    /**
     * @brief Calculates Cos(Theta) in the Center of Mass frame.
     */
    void CosThetaCM() {
       KineCalculation CosTheta(*this, "CosThetaCM", rad::ElS_CosThetaCM);
    }
    
    /**
     * @brief Calculates Phi (Azimuthal Angle) in the Center of Mass frame.
     */
    void PhiCM() {
       KineCalculation Phi(*this, "PhiCM", rad::ElS_PhiCM);
    }

  };

} // namespace rad


/* #pragma once */

/* #include "KinematicsProcessor.h" */
/* #include "ElectronScatterKinematics.h" // Physics Kernels */

/* namespace rad { */

/*   class KinematicsProcElectro : public KinematicsProcessor { */
/*   public: */
    
/*     KinematicsProcElectro(config::ConfigReaction* cr) : KinematicsProcessor(cr) { */
      
/*       // 1. Define Topology: Register the Scattered Electron Group */
/*       // We take the single particle name "scat_ele" and put it in a group "scat_ele_group" */
/*       // This allows ParticleCreator to treat it uniformly as a group. */
/*       Reaction()->setGroupParticles(names::ScatGroup(), {names::ScatEle()}); */

/*       // 2. Register this group with the Creator */
/*       Creator().DefineGroup(names::OrderScatEle(), names::ScatGroup()); */

/*       // 2. Define Topology: Create the Virtual Photon */
/*       // Diff: BeamEle - ScatEle */
/*       Creator().Diff(names::VirtGamma(), {{names::BeamEle()}, {names::ScatEle()}}); */
/*     } */

/*     // --- Electro-Specific Shortcuts --- */

/*     void Q2() { */
/*        KineCalculation Q2(*this, "Q2", rad::ElS_Q2); */
/*     } */
    
/*     void CosThetaCM() { */
/*        KineCalculation CosTheta(*this, "CosThetaCM", rad::ElS_CosThetaCM); */
/*     } */
    
/*     void PhiCM() { */
/*        KineCalculation Phi(*this, "PhiCM", rad::ElS_PhiCM); */
/*     } */

/*     // Note: Generic shortcuts (Mass, Pt) are inherited from KinematicsProcessor */
/*   }; */

/* } // namespace rad */
