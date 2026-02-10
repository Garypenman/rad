/**
 * @file PhysicsSelection.h
 * @brief Manages event selection (Cuts) and mask generation via Lazy Initialization.
 * @details
 * This class provides a centralized interface for defining cuts on kinematic variables.
 * It uses a "Lazy Definition" pattern:
 * 1. Configuration Phase: Store cut parameters via `AddCut...`.
 * 2. Compilation Phase: `Compile()` creates the RDataFrame columns and the master mask.
 * This ensures cuts are defined only AFTER the kinematic variables they depend on exist.
 */

#pragma once

#include "KinematicsProcessor.h"
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>

namespace rad {

    /** * @struct CutDef
     * @brief Configuration storage for deferred definition.
     * @details Stores cut parameters so they can be compiled later.
     */
    struct CutDef {
        enum Type { 
            Range, Min, Max,           // Standard Comparison
            Equal, NotEqual,           // Exact Match (IDs)
            AbsRange, AbsMin, AbsMax   // Magnitude Checks
        };
        std::string name;        ///< Name of the cut (e.g. "MassCut")
        std::string varBaseName; ///< Unresolved variable name (e.g. "MassJ")
        double min;
        double max;
        Type type;
    };

    /**
     * @class PhysicsSelection
     * @brief Aggregates multiple cuts into a single boolean mask.
     */
    class PhysicsSelection {
    public:
        // =====================================================================
        // Constructors
        // =====================================================================

        /**
         * @brief Primary Constructor.
         * @param proc The processor used to resolve variable prefixes/suffixes.
         */
        PhysicsSelection(KinematicsProcessor& proc);

        // =====================================================================
        // Cut Definitions (Standard)
        // =====================================================================

        /** * @brief Define a standard window cut: min < var < max. 
         * @param name Unique name for this cut.
         * @param var Variable name (e.g. "MassJ").
         * @param min Minimum value.
         * @param max Maximum value.
         */
        void AddCutRange(const std::string& name, const std::string& var, double min, double max);

        /** * @brief Define a lower bound cut: var > min. 
         * @param name Unique name for this cut.
         * @param var Variable name.
         * @param min Minimum value.
         */
        void AddCutMin(const std::string& name, const std::string& var, double min);

        /** * @brief Define an upper bound cut: var < max. 
         * @param name Unique name for this cut.
         * @param var Variable name.
         * @param max Maximum value.
         */
        void AddCutMax(const std::string& name, const std::string& var, double max);

        // =====================================================================
        // Cut Definitions (Equality/Identity)
        // =====================================================================

        /** * @brief Define an equality cut: var == val.
         * @details Useful for Detector IDs, Charge, or PID matching.
         * Note: Uses generic double comparison with epsilon tolerance.
         * @param name Unique name for this cut.
         * @param var Variable name.
         * @param val Value to match.
         */
        void AddCutEqual(const std::string& name, const std::string& var, double val);

        /** * @brief Define an inequality cut: var != val. 
         * @param name Unique name for this cut.
         * @param var Variable name.
         * @param val Value to exclude.
         */
        void AddCutNotEqual(const std::string& name, const std::string& var, double val);

        // =====================================================================
        // Cut Definitions (Absolute Value)
        // =====================================================================

        /** * @brief Define an absolute window cut: min < |var| < max. 
         * @details Useful for Vertex Z cuts or mass differences.
         * @param name Unique name for this cut.
         * @param var Variable name.
         * @param min Minimum absolute value.
         * @param max Maximum absolute value.
         */
        void AddCutAbsRange(const std::string& name, const std::string& var, double min, double max);

        /** * @brief Define an absolute upper bound: |var| < max. 
         * @param name Unique name for this cut.
         * @param var Variable name.
         * @param max Maximum absolute value.
         */
        void AddCutAbsMax(const std::string& name, const std::string& var, double max);

        // =====================================================================
        // Execution
        // =====================================================================

        /** * @brief Compiles all added cuts into RDF columns and a master mask. 
         * @details 
         * Loops over the stored configurations and calls Reaction->Define for each.
         * Then creates a master mask column (AND logic).
         * If no cuts are present, creates a "True" vector matching the event size.
         * MUST be called after KinematicsProcessor::Init().
         */
        void Init();

        /** @return The name of the compiled mask column. */
        std::string GetMaskColumn() const;

    private:
        KinematicsProcessor& _proc;
        
        // Configuration Storage (Lazy Definition)
        ROOT::RVec<CutDef> _config;
        
        // Active Column Names (Populated in Compile)
        ROOT::RVec<std::string> _cutNames;
        std::string _finalMask;
    };

    // =========================================================================
    // IMPLEMENTATION
    // =========================================================================

    inline PhysicsSelection::PhysicsSelection(KinematicsProcessor& proc) : _proc(proc) {}

    // --- Standard Cuts (Just Store Config) ---

    inline void PhysicsSelection::AddCutRange(const std::string& name, const std::string& var, double min, double max) {
        _config.push_back({name, var, min, max, CutDef::Range});
    }

    inline void PhysicsSelection::AddCutMin(const std::string& name, const std::string& var, double min) {
        _config.push_back({name, var, min, 0.0, CutDef::Min});
    }

    inline void PhysicsSelection::AddCutMax(const std::string& name, const std::string& var, double max) {
        _config.push_back({name, var, 0.0, max, CutDef::Max});
    }

    inline void PhysicsSelection::AddCutEqual(const std::string& name, const std::string& var, double val) {
        _config.push_back({name, var, val, 0.0, CutDef::Equal});
    }

    inline void PhysicsSelection::AddCutNotEqual(const std::string& name, const std::string& var, double val) {
        _config.push_back({name, var, val, 0.0, CutDef::NotEqual});
    }

    inline void PhysicsSelection::AddCutAbsRange(const std::string& name, const std::string& var, double min, double max) {
        _config.push_back({name, var, min, max, CutDef::AbsRange});
    }

    inline void PhysicsSelection::AddCutAbsMax(const std::string& name, const std::string& var, double max) {
        _config.push_back({name, var, 0.0, max, CutDef::AbsMax});
    }

    // --- Compilation (The Actual Work) ---

    inline void PhysicsSelection::Init() {
        
        _cutNames.clear();
        _finalMask = _proc.GetPrefix() + "Analysis_Mask" + _proc.GetSuffix();
 std::cout << "[PhysicsSelection] Initializing Cuts for Stream: " 
                  << _proc.GetPrefix() << "..." << _proc.GetSuffix() 
                  << " (" << _config.size() << " cuts configured)" << std::endl;

        // 1. Define individual cut columns
        for(const auto& def : _config) {
            std::string col = _proc.FullName(def.varBaseName);
            std::string cutName = _proc.GetPrefix() + def.name + _proc.GetSuffix();
            // DIAGNOSTIC: Print what we are defining
            std::cout << "  -> Defining Cut: " << cutName << " on Variable: " << col << std::endl;

            double min = def.min;
            double max = def.max;

            switch(def.type) {
                // Standard
                case CutDef::Range: 
                    _proc.Reaction()->Define(cutName, [min, max](const RVecResultType& val){ return val > min && val < max; }, {col});
                    break;
                case CutDef::Min:   
                    _proc.Reaction()->Define(cutName, [min](const RVecResultType& val){ return val > min; }, {col});
                    break;
                case CutDef::Max:   
                    _proc.Reaction()->Define(cutName, [max](const RVecResultType& val){ return val < max; }, {col});
                    break;
                
                // Equality
                case CutDef::Equal:    
                    _proc.Reaction()->Define(cutName, [min](const RVecResultType& val){ return ROOT::VecOps::abs(val - min) < 1e-9; }, {col});
                    break;
                case CutDef::NotEqual: 
                    _proc.Reaction()->Define(cutName, [min](const RVecResultType& val){ return ROOT::VecOps::abs(val - min) > 1e-9; }, {col});
                    break;

                // Absolute
                case CutDef::AbsRange: 
                    _proc.Reaction()->Define(cutName, [min, max](const RVecResultType& val){ const RVecResultType& a=ROOT::VecOps::abs(val); return a > min && a < max; }, {col});
                    break;
                case CutDef::AbsMax:   
                    _proc.Reaction()->Define(cutName, [max](const RVecResultType& val){ return ROOT::VecOps::abs(val) < max; }, {col});
                    break;
                default: break;
            }
            _cutNames.push_back(cutName);
        }

        // 2. Define Master Mask (AND logic)
        if (_cutNames.empty()) {
            // CASE: NO CUTS.
            // We cannot use "1" because it creates a scalar boolean. 
            // We need a VECTOR of 1s (all true) with the same length as the combinations.
            // We use the first particle's Px column to determine the size.
            
            auto particle_names = _proc.Creator().GetParticleNames();
            if(particle_names.empty()) {
	      throw std::invalid_argument("PhysicsSelection: particle list cannot be empty.");
                // Fallback for empty/dummy event. "1" is scalar but safe if no data exists.
	      // _proc.Reaction()->Define(_finalMask, "1"); 
            } else {
                std::string refCol = _proc.FullName(particle_names[0] + "_px");
                _proc.Reaction()->Define(_finalMask, [](const ROOT::RVecD& v){
                    // Return vector of same size, all true (1)
                    return ROOT::RVecI(v.size(), 1); 
                }, {refCol});
            }
        } 
        else {
            // CASE: HAS CUTS.
            // "cut1 && cut2" on RVecs produces an RVec<int> (vector mask), which is correct.
            std::stringstream ss;
            for(size_t i=0; i<_cutNames.size(); ++i) {
                ss << _cutNames[i];
                if(i != _cutNames.size() - 1) ss << " && ";
            }
	   
	   // Define a Temporary Boolean Name (The 0/1 result of cuts)
	   std::string boolMaskName = _proc.GetPrefix() + "Analysis_Bool" + _proc.GetSuffix();
          _proc.Reaction()->Define(boolMaskName, ss.str());
	   // B. Convert Boolean -> Indices
	   // Nonzero({1, 0, 1}) -> {0, 2}
	   _proc.Reaction()->Define(_finalMask, 
				    [](const Indices_t& bools) {
				      return ROOT::VecOps::Nonzero(bools);
				    }, 
				    {boolMaskName}
				    );
        }
    }

    inline std::string PhysicsSelection::GetMaskColumn() const {
        return _finalMask;
    }

} // namespace rad
