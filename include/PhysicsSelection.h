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
        std::vector<CutDef> _config;
        
        // Active Column Names (Populated in Compile)
        std::vector<std::string> _cutNames;
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

        // 1. Define individual cut columns
        for(const auto& def : _config) {
            std::string col = _proc.FullName(def.varBaseName);
            std::string cutName = _proc.GetPrefix() + def.name + _proc.GetSuffix();
            
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
	    cout<<"PhysicsSelection "<<ss.str()<<endl;
            _proc.Reaction()->Define(_finalMask, ss.str());
        }
    }

    inline std::string PhysicsSelection::GetMaskColumn() const {
        return _finalMask;
    }

} // namespace rad
// /**
//  * @file PhysicsSelection.h
//  * @brief Helper class to manage "Lazy Masking" of combinatorial candidates.
//  */

// #pragma once

// #include "ConfigReaction.h"
// #include "KinematicsProcessor.h"
// #include "CommonDefines.h"

// #include <string>
// #include <vector>
// #include <sstream>
// #include <cmath>
// #include <iostream>

// namespace rad {

//     /**
//      * @class PhysicsSelection
//      * @brief Builder for defining event selection cuts on combinatorial results.
//      * @details
//      * This class acts as a bridge between the KinematicsProcessor and the Histogrammer.
//      * It manages a list of boolean cut columns and compiles them into a single "Lazy Mask".
//      * * **Workflow:**
//      * 1. Bind to a specific `KinematicsProcessor` (e.g., "rec_").
//      * 2. Define cuts on variables (e.g., "mm2" -> "rec_mm2").
//      * 3. Call `Compile()` to generate the index vector column.
//      * 4. Pass the result to `Histogrammer` or `SnapshotFlat`.
//      */
//     class PhysicsSelection {
//     public:
//         /**
//          * @brief Constructor binding this selection to a specific data stream.
//          * @param processor The processor instance (Rec, Truth, etc.) used to resolve variable names.
//          */
//         explicit PhysicsSelection(KinematicsProcessor& processor) 
//             : _proc(processor), _reaction(processor.Reaction()) {}

//         PhysicsSelection(const PhysicsSelection&) = delete; 
//         PhysicsSelection& operator=(const PhysicsSelection&) = delete;

//       /** @brief Clone Constructor: Copies cuts from 'other' but binds to 'newProc'. */
//       PhysicsSelection(const PhysicsSelection& other, KinematicsProcessor& newProc)
//         : _proc(newProc), _cuts(other._cuts) 
//       {}

//       // =====================================================================
//       // 1. Generic Interface
//       // =====================================================================
      
//         /**
//          * @brief Defines a custom cut using a C++ Lambda or Functor.
//          * @details 
//          * The function is applied via `RDataFrame::Define`. Variable names are automatically
//          * resolved using the processor's prefix/suffix (e.g., "mm2" -> "rec_mm2_miss").
//          * * @tparam Function Type of the lambda/functor.
//          * @param baseName Unique base name for this cut (e.g., "cut_pid").
//          * @param func The function to execute. Must return RVec<int> (0/1) or RVec<bool>.
//          * @param baseCols List of input column names (without prefix/suffix).
//          */
//         template <typename Function>
//         void AddCut(const std::string& baseName, Function&& func, const std::vector<std::string>& baseCols) {
//             std::vector<std::string> fullCols;
//             fullCols.reserve(baseCols.size());
//             for(const auto& c : baseCols) {
//                 fullCols.push_back(_proc.FullName(c));
//             }
            
//             std::string fullName = _proc.FullName(baseName + "_cut");
//             _reaction->Define(fullName, std::forward<Function>(func), fullCols);
//             _cut_names.push_back(fullName);
//         }

//         /**
//          * @brief Defines a custom cut using a String Expression.
//          * @details 
//          * Best for simple math like "x > 0 && y < 5". 
//          * @note Requires the expression to use fully qualified names OR relative names if supported 
//          * by the underlying RDF JIT. For safety, this wrapper defines the output column name 
//          * using the processor's full naming scheme.
//          * * @param baseName Unique base name for this cut.
//          * @param expression The valid RDataFrame string expression.
//          */
//         void AddCut(const std::string& baseName, const std::string& expression) {
//             std::string fullName = _proc.FullName(baseName + "_cut");
//             _reaction->Define(fullName, expression);
//             _cut_names.push_back(fullName);
//         }

//         // =====================================================================
//         // 2. Specific Calculation Filters (Convenience)
//         // =====================================================================

//         /**
//          * @brief Keeps candidates where 'col' is in range [min, max].
//          * @details Automatically handles NaNs (returns false).
//          * @param name Unique name for this cut.
//          * @param col Variable name (without prefix/suffix).
//          * @param min Minimum value (inclusive).
//          * @param max Maximum value (inclusive).
//          */
//         void AddCutRange(const std::string& name, const std::string& col, double min, double max) {
//             auto range_func = [min, max](const ROOT::RVecD& val) {
//                 return val >= min && val <= max;
//             };
//             AddCut(name, range_func, {col});
//         }

//         /** @brief Keeps candidates where 'col' > min. */
//         void AddCutMin(const std::string& name, const std::string& col, double min) {
//             auto min_func = [min](const ROOT::RVecD& val) { return val > min; };
//             AddCut(name, min_func, {col});
//         }

//         /** @brief Keeps candidates where 'col' < max. */
//         void AddCutMax(const std::string& name, const std::string& col, double max) {
//             auto max_func = [max](const ROOT::RVecD& val) { return val < max; };
//             AddCut(name, max_func, {col});
//         }

//         /** @brief Keeps candidates where 'col' == target (Exact Match with epsilon). */
//         void AddCutExact(const std::string& name, const std::string& col, double target) {
//             auto exact_func = [target](const ROOT::RVecD& val) {
//                 return abs(val - target) < 1e-6;
//             };
//             AddCut(name, exact_func, {col});
//         }

//         // =====================================================================
//         // 3. Compilation
//         // =====================================================================

//         /**
//          * @brief Compiles all registered cuts into a single Index Mask.
//          * @details
//          * 1. Creates a Master Boolean Mask (AND logic of all cuts).
//          * 2. Defines a column containing the *indices* of passing candidates.
//          * * @param baseName Base name for the mask column (default: "good_indices").
//          * @return The fully qualified name of the index column.
//          */
//         std::string Compile(const std::string& baseName = "good_indices") {
//             if (_cut_names.empty()) {
//                 std::cerr << "Warning: PhysicsSelection::Compile called with no cuts for processor " 
//                           << _proc.GetPrefix() << ". Returning empty mask string." << std::endl;
//                 return "";
//             }

//             // 1. Build the Master Mask Expression: "cut1 && cut2 && ..."
//             std::stringstream ss;
//             for (size_t i = 0; i < _cut_names.size(); ++i) {
//                 ss << _cut_names[i];
//                 if (i < _cut_names.size() - 1) ss << " && ";
//             }
            
//             std::string maskName = _proc.FullName(baseName + "_bool_mask");
//             _reaction->Define(maskName, ss.str());

//             // 2. Convert Boolean Mask -> Indices using VecOps::Nonzero
//             _final_mask_col = _proc.FullName(baseName);
            
//             // Defines the column that contains {0, 2, 5...} for valid candidates
//             _reaction->Define(_final_mask_col, 
//                 [](const ROOT::RVec<int>& mask) { return ROOT::VecOps::Nonzero(mask); }, 
//                 {maskName}
//             );
            
//             return _final_mask_col;
//         }

//         /** @brief Returns the full name of the compiled index column. */
//         std::string GetMaskColumn() const { return _final_mask_col; }

//     private:
//         KinematicsProcessor& _proc;
//         ConfigReaction* _reaction;
//         std::vector<std::string> _cut_names;
//         std::string _final_mask_col;
//     };

// } // namespace rad
