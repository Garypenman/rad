/**
 * @file ElectroIonReaction.h
 * @brief Configuration class for Electron-Ion scattering experiments.
 * @details
 * This class acts as the bridge between raw data and physics analysis for EIC-like topologies.
 * It manages:
 * 1. **Beam Definitions:** Event-by-Event, or MC-derived beam 4-vectors.
 * 2. **Scattered Electron:** Helper methods to identify the scattered lepton.
 * 3. **Component Definition:** Automatic generation of standard beam columns.
 */

#pragma once

#include "ConfigReaction.h"
#include "RVecHelpers.h"
#include "BasicKinematics.h"
#include "CommonDefines.h"     // Defines ResultType_t, RVecResultType
#include "DefineNames.h"       // Defines NamePx(), NamePy()...
#include "Constants.h"         // Defines PDG_ele(), M_ele()...
#include "ParticleInjector.h" 

#include <TString.h>           // For Form()


//set BeamIndices for electroion reaction
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
     * @details
     * Extends the generic `ConfigReaction` to handle the specific needs of DIS/Photoproduction events.
     * It provides a standardized way to define the initial state (Beams) and the primary
     * scattered lepton, which are prerequisites for calculating $Q^2$, $W$, etc.
     */
    class ElectroIonReaction : public ConfigReaction {

    public:
        // =================================================================================
        // Constructors
        // =================================================================================

        /** * @brief Constructor for globbed filenames.
         * @param treeName Name of the input TTree (e.g. "events").
         * @param fileNameGlob File pattern (e.g. "data/XXX.root").
         * @param columns Optional list of columns to read (optimization).
         */
        ElectroIonReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns = {});

        /** * @brief Constructor for a vector of filenames.
         * @param treeName Name of the input TTree.
         * @param filenames Vector of explicit file paths.
         * @param columns Optional list of columns to read.
         */
        ElectroIonReaction(const std::string_view treeName, const ROOT::RVec<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns = {});

        /** * @brief Constructor for wrapping an existing RDataFrame.
         * @param rdf Existing RDataFrame node.
         */
        ElectroIonReaction(ROOT::RDataFrame rdf);

        /** * @brief Constructor for wrapping an existing RDataFrame.
         * @param rdf Existing RDataFrame node.
         */
        ElectroIonReaction(ROOT::RDF::RNode rdf);

        // =================================================================================
        // Scattered Electron Interface
        // =================================================================================

        /** * @brief Sets a fixed single index for the scattered electron.
         * @details Useful when the scattered electron is always at a known position (e.g., index 0).
         * @param idx The index in the unified vector.
         * @param type The data prefix (e.g., "rec_"). Defaults to the primary type.
         */
        void SetScatElectronIndex(const int idx, const std::string& type = "");

        /** * @brief Sets a list of potential candidates for the scattered electron.
         * @param idx List of indices to consider.
         * @param type The data prefix.
         */
        void SetScatElectronCandidates(const Indices_t& idx, const std::string& type = "");

        /** * @brief Defines scattered electron candidates using a generic lambda.
         * @tparam Lambda Functor returning `Indices_t` (RVecI).
         * @param func The filter function.
         * @param columns Input columns required by the function.
         */
        template<typename Lambda>
        void SetScatElectronCandidates(Lambda&& func, const ROOT::RDF::ColumnNames_t& columns);

        /** * @brief Defines scattered electron candidates using a lambda (Explicit Type).
         * @tparam Lambda Functor returning `Indices_t` (RVecI).
         * @param func The filter function.
         * @param type Explicit type prefix (e.g., "mc_").
         * @param columns Input columns required by the function.
         */
        template<typename Lambda>
        void SetScatElectronCandidates(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t& columns);

        // =================================================================================
        // Beam Configuration
        // =================================================================================

        /**
         * @brief Sets fixed beam momentum (Collider/Fixed Target mode).
         * @details 
         * Sets the internal 4-vector and ensures the unified Beam Index is set to 0.
         * **Note:** calling this disables the automatic MC beam lookup.
         * @param x Px (GeV)
         * @param y Py (GeV)
         * @param z Pz (GeV)
         */
        void SetBeamElectron(double x, double y, double z);

        /**
         * @brief Sets fixed beam ion momentum.
         * @details 
         * Sets the internal 4-vector and ensures the unified Beam Index is set to 0.
         * **Note:** calling this disables the automatic MC beam lookup.
         * @param x Px (GeV)
         * @param y Py (GeV)
         * @param z Pz (GeV)
         * @param m Mass (GeV), defaults to Proton mass.
         */
        void SetBeamIon(double x, double y, double z, double m = consts::M_pro());

        /**
         * @brief Sets the indices of the beams within the MCParticles array.
         * @details 
         * By default, RAD assumes Electron=[0] and Ion=[1] in the `MCParticles` branch.
         * Use this method if your generator stores the beams at different indices.
         * @param eleIdx The index of the beam electron in MC.
         * @param ionIdx The index of the beam ion in MC.
         */
        void SetMCBeamIndices(int eleIdx, int ionIdx);

        /**
         * @brief Maps input data columns to the standardized Beam 4-vectors (Data Mode).
         * @details 
         * Use this when beam properties vary event-by-event in the input file.
         * This creates aliases mapping your input columns to the standard `BeamEle_px`, etc.
         * @param px Input Column Px
         * @param py Input Column Py
         * @param pz Input Column Pz
         * @param m  Input Column Mass
         * @param type Prefix (e.g. "rec_"). Defaults to primary type.
         */
        void SetBeamElectronColumns(const std::string& px, const std::string& py, 
                                    const std::string& pz, const std::string& m, 
                                    const std::string& type = "");

        /**
         * @brief Maps input data columns to the standardized Beam Ion 4-vectors.
         * @see SetBeamElectronColumns
         */
        void SetBeamIonColumns(const std::string& px, const std::string& py, 
                               const std::string& pz, const std::string& m, 
                               const std::string& type = "");

        // =================================================================================
        // Internal Machinery
        // =================================================================================

        /**
         * @brief Creates the RDataFrame columns for beam components.
         * @details 
         * This method standardizes the beam interface regardless of source:
         * 1. **MC Mode:** Extracts single particles from `MCParticles` using configured indices.
         * 2. **Fixed Mode:** Defines constant columns from stored P4 vectors.
         * 3. **Data Mode:** Uses aliases defined by `SetBeam...Columns`.
         * @param prefix  data event type (Rec() or Truth()...)
         */
        void DefineBeamComponents(const std::string& prefix);

        /**
         * @brief Overrides ConfigReaction::MakeCombinations to ensure beams are defined.
         * @details Triggers `DefineBeamComponents()` if not already called, ensuring
         * beam columns exist before combinatorial logic begins.
         */
        void MakeCombinations();

        /** @return The configured Ion Beam 4-vector (valid only in Fixed/MC mode). */
        PxPyPzMVector P4BeamIon() const;

        /** @return The configured Electron Beam 4-vector (valid only in Fixed/MC mode). */
        PxPyPzMVector P4BeamEle() const;

          /** @brief Helper to set the internal index for the Beam Electron (Always 0). */
        void SetBeamElectronIndex(const int idx, const std::string& type = "");

        /** @brief Helper to set the internal index for the Beam Ion (Always 0). */
        void SetBeamIonIndex(const int idx, const std::string& type = "");

    protected:
      
      PxPyPzMVector _p4el_beam;   ///< Internal storage for Electron Beam P4
      PxPyPzMVector _p4ion_beam;  ///< Internal storage for Ion Beam P4
        
 
      bool _useBeamsFromMC = true; ///< Flag to determine if beams should be read from MC
        
     private:

      // Indices in the source MCParticles array (Used for Extraction)
        int _mcBeamEleIdx = 0; ///< Index of Beam Electron in MC array
        int _mcBeamIonIdx = 1; ///< Index of Beam Ion in MC array

    }; // class ElectroIonReaction


    namespace electroion {
        /**
         * @brief Helper to generate a Virtual Photon 4-vector (q = k - k').
         * @details Requires the reaction map to locate the "Virtual Gamma" in the created particles list.
         */
        template<typename Tp, typename Tm>
        inline PxPyPzMVector PhotoFourVector(const RVecIndexMap& react, const Tp &px, const Tp &py, const Tp &pz, const Tm &m);
    }


    // =================================================================================
    // IMPLEMENTATION
    // =================================================================================

    // --- Constructors ---

    inline ElectroIonReaction::ElectroIonReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns) 
        : ConfigReaction(treeName, fileNameGlob, columns) {}

    inline ElectroIonReaction::ElectroIonReaction(const std::string_view treeName, const ROOT::RVec<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns) 
        : ConfigReaction(treeName, filenames, columns) {}

    inline ElectroIonReaction::ElectroIonReaction(ROOT::RDataFrame rdf) 
        : ConfigReaction(rdf) {}

    inline ElectroIonReaction::ElectroIonReaction(ROOT::RDF::RNode rdf) 
        : ConfigReaction(rdf) {}

    // --- Scattered Electron Interface ---

    inline void ElectroIonReaction::SetScatElectronIndex(const int idx, const std::string& type) {
        SetParticleIndex(consts::ScatEle().data(), type.empty() ? GetDefaultType() : type, idx);
    }

    inline void ElectroIonReaction::SetScatElectronCandidates(const Indices_t& idx, const std::string& type) {
        SetParticleCandidates(consts::ScatEle().data(), type.empty() ? GetDefaultType() : type, idx);
    }

    template<typename Lambda>
    inline void ElectroIonReaction::SetScatElectronCandidates(Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
        SetParticleCandidates(consts::ScatEle().data(), GetDefaultType(), std::forward<Lambda>(func), columns);
    }

    template<typename Lambda>
    inline void ElectroIonReaction::SetScatElectronCandidates(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t& columns) {
        SetParticleCandidates(consts::ScatEle().data(), type.empty() ? GetDefaultType() : type, std::forward<Lambda>(func), columns);
    }

    // --- Beam Configuration ---

    inline void ElectroIonReaction::SetBeamElectron(double x, double y, double z) {
        _p4el_beam = PxPyPzMVector(x, y, z, consts::M_ele());
        _useBeamsFromMC = false;
        // The beam is now a fixed scalar column of size 1, so index is always 0
        SetBeamElectronIndex(0); 
    }

    inline void ElectroIonReaction::SetBeamIon(double x, double y, double z, double m) {
        _p4ion_beam = PxPyPzMVector(x, y, z, m);
        _useBeamsFromMC = false;
        SetBeamIonIndex(0); 
    }

    inline void ElectroIonReaction::SetMCBeamIndices(int eleIdx, int ionIdx) {
        _mcBeamEleIdx = eleIdx;
        _mcBeamIonIdx = ionIdx;
    }

    inline void ElectroIonReaction::SetBeamElectronColumns(const std::string& px, const std::string& py, 
                                                           const std::string& pz, const std::string& m, 
                                                           const std::string& type) 
    {
        std::string p = type.empty() ? GetDefaultType() : type;
        std::string name = p + consts::BeamEle() + "_src_";
        
        // Define aliases mapping standard names to input columns
        if(!ColumnExists(name + consts::NamePx())) Define(name + consts::NamePx(), px);
        if(!ColumnExists(name + consts::NamePy())) Define(name + consts::NamePy(), py);
        if(!ColumnExists(name + consts::NamePz())) Define(name + consts::NamePz(), pz);
        if(!ColumnExists(name + consts::NameM()))  Define(name + consts::NameM(),  m);
        
        // Aliased column is treated as the primary source, index is 0
        SetBeamElectronIndex(0, p); 
        _useBeamsFromMC = false; 
    }

    inline void ElectroIonReaction::SetBeamIonColumns(const std::string& px, const std::string& py, 
                                                      const std::string& pz, const std::string& m, 
                                                      const std::string& type) 
    {
        std::string p = type.empty() ? GetDefaultType() : type;
        std::string name = p + consts::BeamIon() + "_src_";

        if(!ColumnExists(name + consts::NamePx())) Define(name + consts::NamePx(), px);
        if(!ColumnExists(name + consts::NamePy())) Define(name + consts::NamePy(), py);
        if(!ColumnExists(name + consts::NamePz())) Define(name + consts::NamePz(), pz);
        if(!ColumnExists(name + consts::NameM()))  Define(name + consts::NameM(),  m);
        
        SetBeamIonIndex(0, p); 
        _useBeamsFromMC = false; 
    }

    // --- Internal Machinery ---

    inline void ElectroIonReaction::MakeCombinations() {
         ConfigReaction::MakeCombinations();
    }

    inline void ElectroIonReaction::DefineBeamComponents(const std::string& prefix) {
      // Helpers to define constant columns
        auto def_vec = [&](std::string name, double val) { 
            if(!ColumnExists(name)) Define(name, [val](){ return rad::RVecResultType{val}; }, {}); 
        };
        auto def_pid = [&](std::string name, int val) { 
            if(!ColumnExists(name)) Define(name, [val](){ return rad::Indices_t{val}; }, {}); 
        };

        // 1. Electron Beam
        std::string eleName = prefix + consts::BeamEle() + "_src_"; 
        
        if(_useBeamsFromMC && ColumnExists("MCParticles.momentum.z")) {
            // MC Source: Use Configurable Indices via Form()
            // Creates a vector of Size 1 containing the value extracted from _mcBeamEleIdx
            Define(eleName + consts::NamePx(), Form("rad::RVecResultType{static_cast<rad::ResultType_t>(MCParticles.momentum.x[%d])}", _mcBeamEleIdx)); 
            Define(eleName + consts::NamePy(), Form("rad::RVecResultType{static_cast<rad::ResultType_t>(MCParticles.momentum.y[%d])}", _mcBeamEleIdx));
            Define(eleName + consts::NamePz(), Form("rad::RVecResultType{static_cast<rad::ResultType_t>(MCParticles.momentum.z[%d])}", _mcBeamEleIdx));
            Define(eleName + consts::NameM(),  Form("rad::RVecResultType{static_cast<rad::ResultType_t>(MCParticles.mass[%d])}",       _mcBeamEleIdx));
            def_pid(eleName + consts::NamePid(), consts::PDG_ele_beam());
        } else {
            // Fixed/Column Source: Check if columns were already aliased by SetBeamElectronColumns
            // If not, use the internal Fixed values (_p4el_beam)
            if(!ColumnExists(eleName + consts::NamePx())) {
                def_vec(eleName + consts::NamePx(), _p4el_beam.Px());
                def_vec(eleName + consts::NamePy(), _p4el_beam.Py());
                def_vec(eleName + consts::NamePz(), _p4el_beam.Pz());
                def_vec(eleName + consts::NameM(),  _p4el_beam.M());
                def_pid(eleName + consts::NamePid(), consts::PDG_ele_beam());
            }
        }

        // 2. Ion Beam
        std::string ionName = prefix + consts::BeamIon() + "_src_"; 
        
        if(_useBeamsFromMC && ColumnExists("MCParticles.momentum.z")) {
            Define(ionName + consts::NamePx(), Form("rad::RVecResultType{static_cast<rad::ResultType_t>(MCParticles.momentum.x[%d])}", _mcBeamIonIdx)); 
            Define(ionName + consts::NamePy(), Form("rad::RVecResultType{static_cast<rad::ResultType_t>(MCParticles.momentum.y[%d])}", _mcBeamIonIdx));
            Define(ionName + consts::NamePz(), Form("rad::RVecResultType{static_cast<rad::ResultType_t>(MCParticles.momentum.z[%d])}", _mcBeamIonIdx));
            Define(ionName + consts::NameM(),  Form("rad::RVecResultType{static_cast<rad::ResultType_t>(MCParticles.mass[%d])}",       _mcBeamIonIdx));
            def_pid(ionName + consts::NamePid(), consts::PDG_pro_beam()); 
        } else {
            if(!ColumnExists(ionName + consts::NamePx())) {
                def_vec(ionName + consts::NamePx(), _p4ion_beam.Px());
                def_vec(ionName + consts::NamePy(), _p4ion_beam.Py());
                def_vec(ionName + consts::NamePz(), _p4ion_beam.Pz());
                def_vec(ionName + consts::NameM(),  _p4ion_beam.M());
                def_pid(ionName + consts::NamePid(), consts::PDG_pro_beam());
            }
        }
    }

    inline PxPyPzMVector ElectroIonReaction::P4BeamIon() const { return _p4ion_beam; }
    inline PxPyPzMVector ElectroIonReaction::P4BeamEle() const { return _p4el_beam; }

    inline void ElectroIonReaction::SetBeamElectronIndex(const int idx, const std::string& type) {
        SetParticleIndex(consts::BeamEle().data(), type.empty() ? GetDefaultType() : type, idx);
    }

    inline void ElectroIonReaction::SetBeamIonIndex(const int idx, const std::string& type) {
        SetParticleIndex(consts::BeamIon().data(), type.empty() ? GetDefaultType() : type, idx);
    }

    // --- Physics Helpers ---

    namespace electroion {
        template<typename Tp, typename Tm>
        inline PxPyPzMVector PhotoFourVector(const RVecIndexMap& react, const Tp &px, const Tp &py, const Tp &pz, const Tm &m) {
            // Retrieves the virtual photon from the 'Created' particle group
            return FourVector(react[consts::OrderCreated()][consts::OrderVirtGamma()], px, py, pz, m);
        }
    }

} // namespace rad

//required for physics, use this computation for photo4vector
using rad::electroion::PhotoFourVector;

// #pragma once

// #include "ConfigReaction.h"
// #include "RVecHelpers.h"
// #include "BasicKinematics.h"
// #include "CommonDefines.h"
// #include "ParticleInjector.h" 

// namespace rad {
//   namespace electroion {
//     /// @brief Returns a comma-separated string of beam indices for configuration.
//     inline const std::string BeamIndices() { 
//       return Form("%s,%s", consts::BeamIon().data(), consts::BeamEle().data()); 
//     }
//   }
// }
// using rad::electroion::BeamIndices;

// namespace rad {
    
//     /**
//      * @class ElectroIonReaction
//      * @brief Configuration base class for Electron-Ion scattering experiments.
//      * * This class extends `ConfigReaction` to provide specific utilities for:
//      * 1. **Beam Definitions:** Setting up Electron and Ion beam 4-vectors.
//      * 2. **Scattered Electron:** Helper methods to define the scattered electron candidate.
//      */
//     class ElectroIonReaction : public ConfigReaction {

//     public:
//       // --- Constructors ---
//       ElectroIonReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns ={} );
//       ElectroIonReaction(const std::string_view treeName, const ROOT::RVec<std::string> &filenames, const ROOT::RDF::ColumnNames_t& columns ={} );
//       ElectroIonReaction(ROOT::RDataFrame rdf);

//       // --- Beam Component Definition ---
//       /**
//        * @brief Standardizes the Beam 4-vectors in the output tree.
//        * Creates columns like `BeamEle_px`, `BeamIon_pz`, etc.
//        */
//       void DefineBeamComponents();

//       // --- Scattered Electron Setters ---
      
//       /** @brief Set a fixed single index for the scattered electron (e.g. index 0). */
//       void SetScatElectronIndex(const int idx, const std::string& type = "");
      
//       /** @brief Set a list of potential candidates for the scattered electron. */
//       void SetScatElectronCandidates(const Indices_t& idx, const std::string& type = "");
      
//       /** @brief Define scattered electron candidates using a lambda function. */
//       template<typename Lambda>
//       void SetScatElectronCandidates(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns);
      
//       // Convenience Overloads (Default Type)
//       template<typename Lambda>
//       void SetScatElectronCandidates(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns);

//       template<typename Lambda>
//       void SetScatElectronIndex(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns);
      
//       template<typename Lambda>
//       void SetScatElectronIndex(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns);

//       // --- Beam Setters ---
      
//       void SetBeamElectronIndex(const int idx, const std::string& type = "");
      
//       template<typename Lambda>
//       void SetBeamElectron(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns);
      
//       template<typename Lambda>
//       void SetBeamElectron(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns);

//       void SetBeamIonIndex(const int idx, const std::string& type = "");
      
//       template<typename Lambda>
//       void SetBeamIon(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns);
      
//       template<typename Lambda>
//       void SetBeamIon(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns);


//       // --- Beam Kinematics Storage ---

//       /** @brief Set fixed beam electron momentum (GeV). */
//       void SetBeamElectron(double x, double y, double z);
     
//       /** @brief Set fixed beam ion momentum (GeV). Default mass is Proton. */
//       void SetBeamIon(double x, double y, double z, double m=0.938272);
      
//       void DefineBeamElectron();
//       void DefineBeamIon();

//       /** @brief Force using fixed beam values, ignoring MC truth info. */
//       void FixBeamElectronMomentum(double x,double y,double z);
      
//       /** @brief Force using fixed beam values, ignoring MC truth info. */
//       void FixBeamIonMomentum(double x,double y,double z,double m=0.938272);
   
//       PxPyPzMVector P4BeamIon()const {return _p4ion_beam;}
//       PxPyPzMVector P4BeamEle()const {return _p4el_beam;}

//     protected:
//       PxPyPzMVector _p4el_beam;
//       PxPyPzMVector _p4ion_beam;
//       bool _useBeamsFromMC = false;

//     }; // ElectroIonReaction

//   // =================================================================================
//   // IMPLEMENTATION
//   // =================================================================================

//   // --- Constructors ---
//   inline ElectroIonReaction::ElectroIonReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns) 
//     : ConfigReaction{treeName,fileNameGlob,columns} {}

//   inline ElectroIonReaction::ElectroIonReaction(const std::string_view treeName, const ROOT::RVec<std::string> &filenames, const ROOT::RDF::ColumnNames_t& columns) 
//     : ConfigReaction{treeName,filenames,columns} {}

//   inline ElectroIonReaction::ElectroIonReaction(ROOT::RDataFrame rdf) 
//     : ConfigReaction{rdf} {}


//   // --- Beam Logic ---
//   inline void ElectroIonReaction::DefineBeamComponents() {
//       if(ColumnExists("BeamEle_px")) return;

//       auto def_vec = [&](std::string name, double val) { Define(name, [val](){ return ROOT::RVecD{val}; }, {}); };
//       auto def_pid = [&](std::string name, int val) { Define(name, [val](){ return ROOT::RVecI{val}; }, {}); };

//       // Electron Beam
//       if(_useBeamsFromMC && ColumnExists("MCParticles.momentum.z")) {
//          Define("BeamEle_px", "ROOT::RVecD{MCParticles.momentum.x[0]}"); 
//          Define("BeamEle_py", "ROOT::RVecD{MCParticles.momentum.y[0]}");
//          Define("BeamEle_pz", "ROOT::RVecD{MCParticles.momentum.z[0]}");
//          Define("BeamEle_m",  "ROOT::RVecD{MCParticles.mass[0]}");
//          def_pid("BeamEle_pid", 11);
//       } else {
//          def_vec("BeamEle_px", _p4el_beam.Px());
//          def_vec("BeamEle_py", _p4el_beam.Py());
//          def_vec("BeamEle_pz", _p4el_beam.Pz());
//          def_vec("BeamEle_m",  _p4el_beam.M());
//          def_pid("BeamEle_pid", 11);
//       }

//       // Ion Beam
//       if(_useBeamsFromMC && ColumnExists("MCParticles.momentum.z")) {
//          Define("BeamIon_px", "ROOT::RVecD{MCParticles.momentum.x[1]}"); 
//          Define("BeamIon_py", "ROOT::RVecD{MCParticles.momentum.y[1]}");
//          Define("BeamIon_pz", "ROOT::RVecD{MCParticles.momentum.z[1]}");
//          Define("BeamIon_m",  "ROOT::RVecD{MCParticles.mass[1]}");
//          def_pid("BeamIon_pid", 2212); 
//       } else {
//          def_vec("BeamIon_px", _p4ion_beam.Px());
//          def_vec("BeamIon_py", _p4ion_beam.Py());
//          def_vec("BeamIon_pz", _p4ion_beam.Pz());
//          def_vec("BeamIon_m",  _p4ion_beam.M());
//          def_pid("BeamIon_pid", 2212);
//       }
//   }

//   // --- Scattered Electron ---
//   inline void ElectroIonReaction::SetScatElectronIndex(const int idx, const std::string& type){
//     SetParticleIndex(consts::ScatEle().data(), type.empty() ? GetDefaultType() : type, idx);
//   }

//   inline void ElectroIonReaction::SetScatElectronCandidates(const Indices_t& idx, const std::string& type){
//     SetParticleCandidates(consts::ScatEle().data(), type.empty() ? GetDefaultType() : type, idx);
//   }

//   template<typename Lambda>
//   inline void ElectroIonReaction::SetScatElectronCandidates(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns){
//      SetParticleCandidates(consts::ScatEle().data(), type.empty() ? GetDefaultType() : type, func, columns);
//   }
  
//   template<typename Lambda>
//   inline void ElectroIonReaction::SetScatElectronCandidates(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns){
//      SetParticleCandidates(consts::ScatEle().data(), GetDefaultType(), func, columns);
//   }

//   template<typename Lambda>
//   inline void ElectroIonReaction::SetScatElectronIndex(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns){
//     SetParticleIndex(consts::ScatEle().data(), type.empty() ? GetDefaultType() : type, func, columns);
//   }
  
//   template<typename Lambda>
//   inline void ElectroIonReaction::SetScatElectronIndex(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns){
//     SetParticleIndex(consts::ScatEle().data(), GetDefaultType(), func, columns);
//   }

//   // --- Beam Electron ---
//   inline void ElectroIonReaction::SetBeamElectronIndex(const int idx, const std::string& type){
//     SetParticleIndex(consts::BeamEle().data(), type.empty() ? GetDefaultType() : type, idx);
//   }
//   template<typename Lambda>
//   inline void ElectroIonReaction::SetBeamElectron(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns){
//     SetParticleIndex(consts::BeamEle().data(), type.empty() ? GetDefaultType() : type, func, columns);
//   }
//   template<typename Lambda>
//   inline void ElectroIonReaction::SetBeamElectron(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns){
//     SetParticleIndex(consts::BeamEle().data(), GetDefaultType(), func, columns);
//   }

//   // --- Beam Ion ---
//   inline void ElectroIonReaction::SetBeamIonIndex(const int idx, const std::string& type){
//     SetParticleIndex(consts::BeamIon().data(), type.empty() ? GetDefaultType() : type, idx);
//   }
//   template<typename Lambda>
//   inline void ElectroIonReaction::SetBeamIon(Lambda&& func, const std::string& type, const ROOT::RDF::ColumnNames_t & columns){
//     SetParticleIndex(consts::BeamIon().data(), type.empty() ? GetDefaultType() : type, func, columns);
//   }
//   template<typename Lambda>
//   inline void ElectroIonReaction::SetBeamIon(Lambda&& func, const ROOT::RDF::ColumnNames_t & columns){
//     SetParticleIndex(consts::BeamIon().data(), GetDefaultType(), func, columns);
//   }

//   // --- Fixed Beam Values ---
//   inline void ElectroIonReaction::SetBeamElectron(double x, double y, double z){
//     _p4el_beam = PxPyPzMVector{x,y,z,0.000510999};
//   }
 
//   inline void ElectroIonReaction::SetBeamIon(double x, double y, double z, double m){
//     _p4ion_beam = PxPyPzMVector{x,y,z,m};
//   }
  
//   inline void ElectroIonReaction::DefineBeamElectron(){
//     if( !_useBeamsFromMC ) return;
//     auto p4=_p4el_beam;
//     Define(rad::consts::P4BeamEle(),[p4](){return p4;},{});
//   }
//   inline void ElectroIonReaction::DefineBeamIon(){
//     if( !_useBeamsFromMC ) return;
//     auto p4=_p4ion_beam;
//     Define(rad::consts::P4BeamIon(),[p4](){return p4;},{});
//   }

//   inline void ElectroIonReaction::FixBeamElectronMomentum(double x,double y,double z){
//     SetBeamElectron(x,y,z);
//     _useBeamsFromMC=false; 
//     DefineBeamElectron();
//   }
//   inline void ElectroIonReaction::FixBeamIonMomentum(double x,double y,double z,double m){
//     SetBeamIon(x,y,z,m);
//     _useBeamsFromMC=false; 
//     DefineBeamIon();
//   }


//   namespace electroion{
//     /**
//      * @brief Calculates the Virtual Photon 4-vector (q = k - k').
//      */
//     template<typename Tp, typename Tm>
//     inline PxPyPzMVector PhotoFourVector(const RVecIndexMap& react, const Tp &px, const Tp &py, const Tp &pz, const Tm &m){
//       // OrderVirtGamma() is usually 0 within the Created Particles group.
//       return FourVector(react[consts::OrderCreated()][consts::OrderVirtGamma()],px,py,pz,m);
//     }
//   }
// }
// using rad::electroion::PhotoFourVector;
