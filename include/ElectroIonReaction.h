#pragma once

#include "ConfigReaction.h"
#include "RVecHelpers.h"
#include "BasicKinematics.h"
#include "Beams.h"
#include "CommonDefines.h"
#include "ParticleInjector.h" 

namespace rad {
  namespace electroion {
    /// @brief Returns a comma-separated string of beam indices for configuration.
    inline const std::string BeamIndices() { 
      return Form("%s,%s", consts::BeamIon().data(), consts::BeamEle().data()); 
    }
  }
}
using rad::electroion::BeamIndices;

namespace rad {
    
    /**
     * @class ElectroIonReaction
     * @brief Configuration base class for Electron-Ion scattering experiments.
     * * This class extends `ConfigReaction` to provide specific utilities for:
     * 1. **Beam Definitions:** Setting up Electron and Ion beam 4-vectors.
     * 2. **Scattered Electron:** Helper methods to define the scattered electron candidate.
     */
    class ElectroIonReaction : public ConfigReaction {

    public:
      // --- Constructors ---
      ElectroIonReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns ={} );
      ElectroIonReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t& columns ={} );
      ElectroIonReaction(ROOT::RDataFrame rdf);

      // --- Beam Component Definition ---
      /**
       * @brief Standardizes the Beam 4-vectors in the output tree.
       * Creates columns like `BeamEle_px`, `BeamIon_pz`, etc.
       */
      void DefineBeamComponents();

      // --- Scattered Electron Setters ---
      
      /** @brief Set a fixed single index for the scattered electron (e.g. index 0). */
      void setScatElectronIndex(const int idx, const std::string& type = "");
      
      /** @brief Set a list of potential candidates for the scattered electron. */
      void setScatElectronCandidates(const Indices_t& idx, const std::string& type = "");
      
      /** @brief Define scattered electron candidates using a lambda function. */
      template<typename Lambda>
      void setScatElectronCandidates(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns);
      
      // Convenience Overloads (Default Type)
      template<typename Lambda>
      void setScatElectronCandidates(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns);

      template<typename Lambda>
      void setScatElectronIndex(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns);
      
      template<typename Lambda>
      void setScatElectronIndex(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns);

      // --- Beam Setters ---
      
      void setBeamElectronIndex(const int idx, const std::string& type = "");
      
      template<typename Lambda>
      void setBeamElectron(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns);
      
      template<typename Lambda>
      void setBeamElectron(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns);

      void setBeamIonIndex(const int idx, const std::string& type = "");
      
      template<typename Lambda>
      void setBeamIon(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns);
      
      template<typename Lambda>
      void setBeamIon(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns);


      // --- Beam Kinematics Storage ---

      /** @brief Set fixed beam electron momentum (GeV). */
      void setBeamElectron(double x, double y, double z);
     
      /** @brief Set fixed beam ion momentum (GeV). Default mass is Proton. */
      void setBeamIon(double x, double y, double z, double m=0.938272);
      
      void DefineBeamElectron();
      void DefineBeamIon();

      /** @brief Force using fixed beam values, ignoring MC truth info. */
      void FixBeamElectronMomentum(double x,double y,double z);
      
      /** @brief Force using fixed beam values, ignoring MC truth info. */
      void FixBeamIonMomentum(double x,double y,double z,double m=0.938272);
   
      PxPyPzMVector P4BeamIon()const {return _p4ion_beam;}
      PxPyPzMVector P4BeamEle()const {return _p4el_beam;}

    protected:
      PxPyPzMVector _p4el_beam;
      PxPyPzMVector _p4ion_beam;
      bool _useBeamsFromMC = false;

    }; // ElectroIonReaction

  // =================================================================================
  // IMPLEMENTATION
  // =================================================================================

  // --- Constructors ---
  inline ElectroIonReaction::ElectroIonReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns) 
    : ConfigReaction{treeName,fileNameGlob,columns} {}

  inline ElectroIonReaction::ElectroIonReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t& columns) 
    : ConfigReaction{treeName,filenames,columns} {}

  inline ElectroIonReaction::ElectroIonReaction(ROOT::RDataFrame rdf) 
    : ConfigReaction{rdf} {}


  // --- Beam Logic ---
  inline void ElectroIonReaction::DefineBeamComponents() {
      if(ColumnExists("BeamEle_px")) return;

      auto def_vec = [&](std::string name, double val) { Define(name, [val](){ return ROOT::RVecD{val}; }, {}); };
      auto def_pid = [&](std::string name, int val) { Define(name, [val](){ return ROOT::RVecI{val}; }, {}); };

      // Electron Beam
      if(_useBeamsFromMC && ColumnExists("MCParticles.momentum.z")) {
         Define("BeamEle_px", "ROOT::RVecD{MCParticles.momentum.x[0]}"); 
         Define("BeamEle_py", "ROOT::RVecD{MCParticles.momentum.y[0]}");
         Define("BeamEle_pz", "ROOT::RVecD{MCParticles.momentum.z[0]}");
         Define("BeamEle_m",  "ROOT::RVecD{MCParticles.mass[0]}");
         def_pid("BeamEle_pid", 11);
      } else {
         def_vec("BeamEle_px", _p4el_beam.Px());
         def_vec("BeamEle_py", _p4el_beam.Py());
         def_vec("BeamEle_pz", _p4el_beam.Pz());
         def_vec("BeamEle_m",  _p4el_beam.M());
         def_pid("BeamEle_pid", 11);
      }

      // Ion Beam
      if(_useBeamsFromMC && ColumnExists("MCParticles.momentum.z")) {
         Define("BeamIon_px", "ROOT::RVecD{MCParticles.momentum.x[1]}"); 
         Define("BeamIon_py", "ROOT::RVecD{MCParticles.momentum.y[1]}");
         Define("BeamIon_pz", "ROOT::RVecD{MCParticles.momentum.z[1]}");
         Define("BeamIon_m",  "ROOT::RVecD{MCParticles.mass[1]}");
         def_pid("BeamIon_pid", 2212); 
      } else {
         def_vec("BeamIon_px", _p4ion_beam.Px());
         def_vec("BeamIon_py", _p4ion_beam.Py());
         def_vec("BeamIon_pz", _p4ion_beam.Pz());
         def_vec("BeamIon_m",  _p4ion_beam.M());
         def_pid("BeamIon_pid", 2212);
      }
  }

  // --- Scattered Electron ---
  inline void ElectroIonReaction::setScatElectronIndex(const int idx, const std::string& type){
    setParticleIndex(consts::ScatEle().data(), type.empty() ? GetDefaultType() : type, idx);
  }

  inline void ElectroIonReaction::setScatElectronCandidates(const Indices_t& idx, const std::string& type){
    setParticleCandidates(consts::ScatEle().data(), type.empty() ? GetDefaultType() : type, idx);
  }

  template<typename Lambda>
  inline void ElectroIonReaction::setScatElectronCandidates(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns){
     setParticleCandidates(consts::ScatEle().data(), type.empty() ? GetDefaultType() : type, func, columns);
  }
  
  template<typename Lambda>
  inline void ElectroIonReaction::setScatElectronCandidates(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns){
     setParticleCandidates(consts::ScatEle().data(), GetDefaultType(), func, columns);
  }

  template<typename Lambda>
  inline void ElectroIonReaction::setScatElectronIndex(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns){
    setParticleIndex(consts::ScatEle().data(), type.empty() ? GetDefaultType() : type, func, columns);
  }
  
  template<typename Lambda>
  inline void ElectroIonReaction::setScatElectronIndex(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns){
    setParticleIndex(consts::ScatEle().data(), GetDefaultType(), func, columns);
  }

  // --- Beam Electron ---
  inline void ElectroIonReaction::setBeamElectronIndex(const int idx, const std::string& type){
    setParticleIndex(consts::BeamEle().data(), type.empty() ? GetDefaultType() : type, idx);
  }
  template<typename Lambda>
  inline void ElectroIonReaction::setBeamElectron(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns){
    setParticleIndex(consts::BeamEle().data(), type.empty() ? GetDefaultType() : type, func, columns);
  }
  template<typename Lambda>
  inline void ElectroIonReaction::setBeamElectron(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns){
    setParticleIndex(consts::BeamEle().data(), GetDefaultType(), func, columns);
  }

  // --- Beam Ion ---
  inline void ElectroIonReaction::setBeamIonIndex(const int idx, const std::string& type){
    setParticleIndex(consts::BeamIon().data(), type.empty() ? GetDefaultType() : type, idx);
  }
  template<typename Lambda>
  inline void ElectroIonReaction::setBeamIon(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns){
    setParticleIndex(consts::BeamIon().data(), type.empty() ? GetDefaultType() : type, func, columns);
  }
  template<typename Lambda>
  inline void ElectroIonReaction::setBeamIon(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns){
    setParticleIndex(consts::BeamIon().data(), GetDefaultType(), func, columns);
  }

  // --- Fixed Beam Values ---
  inline void ElectroIonReaction::setBeamElectron(double x, double y, double z){
    _p4el_beam = PxPyPzMVector{x,y,z,0.000510999};
  }
 
  inline void ElectroIonReaction::setBeamIon(double x, double y, double z, double m){
    _p4ion_beam = PxPyPzMVector{x,y,z,m};
  }
  
  inline void ElectroIonReaction::DefineBeamElectron(){
    if( !_useBeamsFromMC ) return;
    auto p4=_p4el_beam;
    Define(rad::consts::P4BeamEle(),[p4](){return p4;},{});
  }
  inline void ElectroIonReaction::DefineBeamIon(){
    if( !_useBeamsFromMC ) return;
    auto p4=_p4ion_beam;
    Define(rad::consts::P4BeamIon(),[p4](){return p4;},{});
  }

  inline void ElectroIonReaction::FixBeamElectronMomentum(double x,double y,double z){
    setBeamElectron(x,y,z);
    _useBeamsFromMC=false; 
    DefineBeamElectron();
  }
  inline void ElectroIonReaction::FixBeamIonMomentum(double x,double y,double z,double m){
    setBeamIon(x,y,z,m);
    _useBeamsFromMC=false; 
    DefineBeamIon();
  }


  namespace electroion{
    /**
     * @brief Calculates the Virtual Photon 4-vector (q = k - k').
     */
    template<typename Tp, typename Tm>
    inline PxPyPzMVector PhotoFourVector(const RVecIndexMap& react, const Tp &px, const Tp &py, const Tp &pz, const Tm &m){
      // OrderVirtGamma() is usually 0 within the Created Particles group.
      return FourVector(react[consts::OrderCreated()][consts::OrderVirtGamma()],px,py,pz,m);
    }
  }
}
using rad::electroion::PhotoFourVector;
