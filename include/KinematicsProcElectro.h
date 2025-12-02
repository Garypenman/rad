#pragma once

#include "KinematicsProcessor.h"
#include "ElectronScatterKinematics.h" // Physics Kernels

namespace rad {

  class KinematicsProcElectro : public KinematicsProcessor {
  public:
    
    KinematicsProcElectro(config::ConfigReaction* cr) : KinematicsProcessor(cr) {
      
      // 1. Define Topology: Register the Scattered Electron Group
      // We take the single particle name "scat_ele" and put it in a group "scat_ele_group"
      // This allows ParticleCreator to treat it uniformly as a group.
      Reaction()->setGroupParticles(names::ScatGroup(), {names::ScatEle()});

      // 2. Register this group with the Creator
      Creator().DefineGroup(names::OrderScatEle(), names::ScatGroup());

      // 2. Define Topology: Create the Virtual Photon
      // Diff: BeamEle - ScatEle
      Creator().Diff(names::VirtGamma(), {{names::BeamEle()}, {names::ScatEle()}});
    }

    // --- Electro-Specific Shortcuts ---

    void Q2() {
       KineCalculation Q2(*this, "Q2", rad::ElS_Q2);
    }
    
    void CosThetaCM() {
       KineCalculation CosTheta(*this, "CosThetaCM", rad::ElS_CosThetaCM);
    }
    
    void PhiCM() {
       KineCalculation Phi(*this, "PhiCM", rad::ElS_PhiCM);
    }

    // Note: Generic shortcuts (Mass, Pt) are inherited from KinematicsProcessor
  };

} // namespace rad
