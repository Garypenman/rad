/**
 * @file PhysicsSelection.h
 * @brief Helper class to manage "Lazy Masking" of combinatorial candidates.
 */

#pragma once

#include "ConfigReaction.h"
#include "KinematicsProcessor.h"
#include "CommonDefines.h"

#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <iostream>

namespace rad {

    /**
     * @class PhysicsSelection
     * @brief Builder for defining event selection cuts on combinatorial results.
     * @details
     * This class acts as a bridge between the KinematicsProcessor and the Histogrammer.
     * It manages a list of boolean cut columns and compiles them into a single "Lazy Mask".
     * * **Workflow:**
     * 1. Bind to a specific `KinematicsProcessor` (e.g., "rec_").
     * 2. Define cuts on variables (e.g., "mm2" -> "rec_mm2").
     * 3. Call `Compile()` to generate the index vector column.
     * 4. Pass the result to `Histogrammer` or `SnapshotFlat`.
     */
    class PhysicsSelection {
    public:
        /**
         * @brief Constructor binding this selection to a specific data stream.
         * @param processor The processor instance (Rec, Truth, etc.) used to resolve variable names.
         */
        explicit PhysicsSelection(KinematicsProcessor& processor) 
            : _proc(processor), _reaction(processor.Reaction()) {}

        PhysicsSelection(const PhysicsSelection&) = delete; 
        PhysicsSelection& operator=(const PhysicsSelection&) = delete;

        // =====================================================================
        // 1. Generic Interface
        // =====================================================================

        /**
         * @brief Defines a custom cut using a C++ Lambda or Functor.
         * @details 
         * The function is applied via `RDataFrame::Define`. Variable names are automatically
         * resolved using the processor's prefix/suffix (e.g., "mm2" -> "rec_mm2_miss").
         * * @tparam Function Type of the lambda/functor.
         * @param baseName Unique base name for this cut (e.g., "cut_pid").
         * @param func The function to execute. Must return RVec<int> (0/1) or RVec<bool>.
         * @param baseCols List of input column names (without prefix/suffix).
         */
        template <typename Function>
        void AddCut(const std::string& baseName, Function&& func, const std::vector<std::string>& baseCols) {
            std::vector<std::string> fullCols;
            fullCols.reserve(baseCols.size());
            for(const auto& c : baseCols) {
                fullCols.push_back(_proc.FullName(c));
            }
            
            std::string fullName = _proc.FullName(baseName + "_cut");
            _reaction->Define(fullName, std::forward<Function>(func), fullCols);
            _cut_names.push_back(fullName);
        }

        /**
         * @brief Defines a custom cut using a String Expression.
         * @details 
         * Best for simple math like "x > 0 && y < 5". 
         * @note Requires the expression to use fully qualified names OR relative names if supported 
         * by the underlying RDF JIT. For safety, this wrapper defines the output column name 
         * using the processor's full naming scheme.
         * * @param baseName Unique base name for this cut.
         * @param expression The valid RDataFrame string expression.
         */
        void AddCut(const std::string& baseName, const std::string& expression) {
            std::string fullName = _proc.FullName(baseName + "_cut");
            _reaction->Define(fullName, expression);
            _cut_names.push_back(fullName);
        }

        // =====================================================================
        // 2. Specific Calculation Filters (Convenience)
        // =====================================================================

        /**
         * @brief Keeps candidates where 'col' is in range [min, max].
         * @details Automatically handles NaNs (returns false).
         * @param name Unique name for this cut.
         * @param col Variable name (without prefix/suffix).
         * @param min Minimum value (inclusive).
         * @param max Maximum value (inclusive).
         */
        void AddCutRange(const std::string& name, const std::string& col, double min, double max) {
            auto range_func = [min, max](const ROOT::RVecD& val) {
                return val >= min && val <= max;
            };
            AddCut(name, range_func, {col});
        }

        /** @brief Keeps candidates where 'col' > min. */
        void AddCutMin(const std::string& name, const std::string& col, double min) {
            auto min_func = [min](const ROOT::RVecD& val) { return val > min; };
            AddCut(name, min_func, {col});
        }

        /** @brief Keeps candidates where 'col' < max. */
        void AddCutMax(const std::string& name, const std::string& col, double max) {
            auto max_func = [max](const ROOT::RVecD& val) { return val < max; };
            AddCut(name, max_func, {col});
        }

        /** @brief Keeps candidates where 'col' == target (Exact Match with epsilon). */
        void AddCutExact(const std::string& name, const std::string& col, double target) {
            auto exact_func = [target](const ROOT::RVecD& val) {
                return abs(val - target) < 1e-6;
            };
            AddCut(name, exact_func, {col});
        }

        // =====================================================================
        // 3. Compilation
        // =====================================================================

        /**
         * @brief Compiles all registered cuts into a single Index Mask.
         * @details
         * 1. Creates a Master Boolean Mask (AND logic of all cuts).
         * 2. Defines a column containing the *indices* of passing candidates.
         * * @param baseName Base name for the mask column (default: "good_indices").
         * @return The fully qualified name of the index column.
         */
        std::string Compile(const std::string& baseName = "good_indices") {
            if (_cut_names.empty()) {
                std::cerr << "Warning: PhysicsSelection::Compile called with no cuts for processor " 
                          << _proc.GetPrefix() << ". Returning empty mask string." << std::endl;
                return "";
            }

            // 1. Build the Master Mask Expression: "cut1 && cut2 && ..."
            std::stringstream ss;
            for (size_t i = 0; i < _cut_names.size(); ++i) {
                ss << _cut_names[i];
                if (i < _cut_names.size() - 1) ss << " && ";
            }
            
            std::string maskName = _proc.FullName(baseName + "_bool_mask");
            _reaction->Define(maskName, ss.str());

            // 2. Convert Boolean Mask -> Indices using VecOps::Nonzero
            _final_mask_col = _proc.FullName(baseName);
            
            // Defines the column that contains {0, 2, 5...} for valid candidates
            _reaction->Define(_final_mask_col, 
                [](const ROOT::RVec<int>& mask) { return ROOT::VecOps::Nonzero(mask); }, 
                {maskName}
            );
            
            return _final_mask_col;
        }

        /** @brief Returns the full name of the compiled index column. */
        std::string GetMaskColumn() const { return _final_mask_col; }

    private:
        KinematicsProcessor& _proc;
        ConfigReaction* _reaction;
        std::vector<std::string> _cut_names;
        std::string _final_mask_col;
    };

} // namespace rad
