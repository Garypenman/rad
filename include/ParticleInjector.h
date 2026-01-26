#pragma once

#include "ConfigReaction.h"
#include "ReactionUtilities.h" // Ensures Concatenate<T> is available
#include "StringUtilities.h"

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

namespace rad {

    /**
     * @class ParticleInjector
     * @brief Manages the merging of disparate data sources into unified RDataFrame columns.
     * * @details
     * In many analyses, particle data comes from different branches or logic streams.
     * For example, "Beams" might come from a fixed 4-vector, while "Tracks" come from 
     * a reconstruction array. The KinematicsProcessor requires a single unified vector 
     * (e.g., `rec_px`) where indices [0,1] are beams and [2...] are tracks.
     * * This class allows you to:
     * 1. Define the schema (types and names) of the data.
     * 2. Add multiple sources (vectors of column names).
     * 3. Automatically generate the RDataFrame code to merge them into a single vector.
     */
    class ParticleInjector {
    public:
        
        /**
         * @brief Construct a new Particle Injector.
         * @param reaction Pointer to the ConfigReaction (interface to RDataFrame).
         */
        explicit ParticleInjector(ConfigReaction* reaction);

        /**
         * @brief Defines the structure and type of the particle data to be injected.
         * * @details
         * Parses strings to determine the C++ type and the variable name.
         * Support multi-word types like "unsigned int" or "long long".
         * Logic: The LAST word is the variable name; everything before it is the type.
         * * @param columns List of definition strings.
         * * **Examples:**
         * - `{"double px", "double py", "int pid"}` -> Explicit types.
         * - `{"unsigned int status"}` -> Multi-word type support.
         * - `{"px", "py"}` -> Defaults to "double" if no type specified.
         */
        void DefineParticleInfo(const std::vector<std::string>& columns);

        /**
         * @brief Adds a source of data to be merged into the unified vectors.
         * * @details
         * Registers a set of columns that correspond to the schema defined in DefineParticleInfo.
         * It creates unique temporary aliases for these columns to prevent name collisions
         * before the final merge.
         * * @param prefix The target prefix for the final vector (e.g., "rec_" or "mc_").
         * @param cols The list of existing RDataFrame column names to read from.
         * Must match the size and order of DefineParticleInfo.
         * @param filter Optional RDataFrame filter string (e.g., "genStatus==1").
         * If provided, the column is read as `colName[filter]`.
         * * @throws std::runtime_error if the number of columns doesn't match the schema.
         */
        void AddSource(const std::string& prefix, const std::vector<std::string>& cols, const std::string& filter = "");

        /**
         * @brief Finalizes the injection by creating the unified vectors.
         * * @details
         * Iterates through all defined properties (px, py, pid...) and all registered sources.
         * Constructs a JIT string for `rad::util::Concatenate<TYPE>(...)` and defines the
         * final column (e.g., "rec_px").
         * * Uses the type specified in DefineParticleInfo directly in the template parameter.
         */
        void CreateUnifiedVectors();

    private:
        ConfigReaction* _reaction;                      ///< Pointer to the reaction interface
        std::vector<std::string> _colNames;             ///< Ordered list of variable names (e.g. px, py)
        std::map<std::string, std::string> _colTypes;   ///< Map of variable name -> C++ type string
        
        // Nested Map Structure:
        // Prefix (rec_) -> List of Sources -> List of Temporary Column Names
        std::map<std::string, std::vector<std::vector<std::string>>> _sources;
    };

    // =========================================================================
    // IMPLEMENTATION
    // =========================================================================

    inline ParticleInjector::ParticleInjector(ConfigReaction* reaction) 
        : _reaction(reaction) {}

    inline void ParticleInjector::DefineParticleInfo(const std::vector<std::string>& columns) {
        _colNames.clear();
        _colTypes.clear();

        for(const auto& entry : columns) {
            std::stringstream ss(entry);
            std::string segment;
            std::vector<std::string> parts;
            
            // Split string by spaces, ignoring empty segments (handles multiple spaces)
            while(std::getline(ss, segment, ' ')) {
                if(!segment.empty()) parts.push_back(segment);
            }

            if(parts.size() >= 2) {
                // Case 1: Explicit Type (e.g., "unsigned int pid")
                // The name is always the last element
                std::string name = parts.back();
                
                // Reassemble the type from all preceding elements
                std::string type = "";
                for(size_t i = 0; i < parts.size() - 1; ++i) {
                    type += parts[i] + " ";
                }
                
                // Remove trailing space from reassembly
                if(!type.empty()) type.pop_back();

                _colNames.push_back(name);
                _colTypes[name] = type; 
            } 
            else if(parts.size() == 1) {
                // Case 2: Name only (e.g., "px")
                // Default to double for compatibility with standard ROOT physics vectors
                _colNames.push_back(parts[0]);
                _colTypes[parts[0]] = "double";
            }
        }
    }

    inline void ParticleInjector::AddSource(const std::string& prefix, const std::vector<std::string>& cols, const std::string& filter) {
        if(cols.size() != _colNames.size()) {
            throw std::runtime_error("ParticleInjector Error: Source column count mismatch. " 
                                     "Expected " + std::to_string(_colNames.size()) + 
                                     ", got " + std::to_string(cols.size()));
        }
        
        // Generate a unique ID for this source to ensure temporary names don't collide.
        // e.g. "rec_px_src0", "rec_px_src1"
        std::string sourceID = "_src" + std::to_string(_sources[prefix].size()); 
        std::vector<std::string> registeredCols;
        
        for(size_t i = 0; i < cols.size(); ++i) {
	  std::string tempName = prefix + _colNames[i] + sourceID + DoNotWriteTag();
            std::string inputCol = cols[i];
            
            // Define the temporary column in RDataFrame
            if(!filter.empty()) {
                // Define filtered version: Define(tempName, "inputCol[filter]")
                _reaction->Define(tempName, inputCol + "[" + filter + "]");
            } else {
                // Simple Alias: Define(tempName, inputCol)
                // Note: We use Define (alias) to expose the column under the temporary name
                _reaction->Define(tempName, inputCol); 
            }
            registeredCols.push_back(tempName);
        }
        
        // Store the temporary names for the final merge step
        _sources[prefix].push_back(registeredCols);
    }

    inline void ParticleInjector::CreateUnifiedVectors() {
        for(auto& [prefix, sources] : _sources) {
            
            // Loop over each property (e.g., "px", then "py", then "pid")
            for(size_t i = 0; i < _colNames.size(); ++i) {
                std::string finalName = prefix + _colNames[i];
                std::string type      = _colTypes[_colNames[i]];
                if(_colNames[i]=="m")cout<<"ParticleInjector::CreateUnifiedVectors() " <<finalName<<" "<<type<<endl;
                // Collect the temporary column names from all sources for this property
                std::vector<std::string> colsToMerge;
                for(const auto& src : sources) {
                    colsToMerge.push_back(src[i]);
                }
                
                // Convert vector of strings to comma-separated string: "col1, col2, col3"
                std::string args = util::ColumnsToStringNoBraces(colsToMerge);

                // DIRECT TEMPLATE INJECTION
                // Generates code: rad::util::Concatenate< type >( col1, col2, ... )
                // This relies on the JIT compiler to instantiate the correct template specialization.
                std::string call = "rad::util::Concatenate<" + type + ">(" + args + ")";
		if(_colNames[i]=="m")cout<<"ParticleInjector::CreateUnifiedVectors() " <<call<<endl;
              
                // Define the final unified column
                _reaction->Define(finalName, call);
            }
        }
    }

} // namespace rad

// #pragma once

// #include "ConfigReaction.h"
// #include "ReactionUtilities.h" // Ensures Concatenate<T> is available
// #include "StringUtilities.h"
// #include <map>
// #include <vector>
// #include <string>
// #include <sstream>

// namespace rad {

//   /**
//    * @class ParticleInjector
//    * @brief Manages the merging of disparate data sources into unified RDataFrame columns.
//    * * @details
//    * In many analyses, particle data comes from different branches or logic streams.
//    * For example, "Beams" might come from a fixed 4-vector, while "Tracks" come from 
//    * a reconstruction array. The KinematicsProcessor requires a single unified vector 
//    * (e.g., `rec_px`) where indices [0,1] are beams and [2...] are tracks.
//    * * This class allows you to:
//    * 1. Define the schema (types and names) of the data.
//    * 2. Add multiple sources (vectors of column names).
//    * 3. Automatically generate the RDataFrame code to merge them into a single vector.
//    */
//   class ParticleInjector {
//   public:
    
//     /**
//      * @brief Construct a new Particle Injector.
//      * @param reaction Pointer to the ConfigReaction (interface to RDataFrame).
//      */
//     ParticleInjector(ConfigReaction* reaction) : _reaction(reaction) {}

//     /**
//      * @brief Defines the structure and type of the particle data to be injected.
//      * * @details
//      * Parses strings to determine the C++ type and the variable name.
//      * Support multi-word types like "unsigned int" or "long long".
//      * Logic: The LAST word is the variable name; everything before it is the type.
//      * * @param columns List of definition strings.
//      * * **Examples:**
//      * - `{"double px", "double py", "int pid"}` -> Explicit types.
//      * - `{"unsigned int status"}` -> Multi-word type support.
//      * - `{"px", "py"}` -> Defaults to "double" if no type specified.
//      */
//     void DefineParticleInfo(const std::vector<std::string>& columns) {
//         _colNames.clear();
//         _colTypes.clear();

//         for(const auto& entry : columns) {
//             std::stringstream ss(entry);
//             std::string segment;
//             std::vector<std::string> parts;
//             // Split string by spaces
//             while(std::getline(ss, segment, ' ')) {
//                 if(!segment.empty()) parts.push_back(segment);
//             }

//             if(parts.size() >= 2) {
//                 // Case 1: Explicit Type (e.g., "unsigned int pid")
//                 // The name is always the last element
//                 std::string name = parts.back();
                
//                 // Reassemble the type from all preceding elements
//                 std::string type = "";
//                 for(size_t i=0; i<parts.size()-1; ++i) type += parts[i] + " ";
                
//                 // Remove trailing space from reassembly
//                 if(!type.empty()) type.pop_back();

//                 _colNames.push_back(name);
//                 _colTypes[name] = type; 
//             } 
//             else if(parts.size() == 1) {
//                 // Case 2: Name only (e.g., "px")
//                 // Default to double for compatibility with standard ROOT physics vectors
//                 _colNames.push_back(parts[0]);
//                 _colTypes[parts[0]] = "double";
//             }
//         }
//     }

//     /**
//      * @brief Adds a source of data to be merged into the unified vectors.
//      * * @details
//      * Registers a set of columns that correspond to the schema defined in DefineParticleInfo.
//      * It creates unique temporary aliases for these columns to prevent name collisions
//      * before the final merge.
//      * * @param prefix The target prefix for the final vector (e.g., "rec_" or "mc_").
//      * @param cols The list of existing RDataFrame column names to read from.
//      * Must match the size and order of DefineParticleInfo.
//      * @param filter Optional RDataFrame filter string (e.g., "genStatus==1").
//      * If provided, the column is read as `colName[filter]`.
//      * * @throws std::runtime_error if the number of columns doesn't match the schema.
//      */
//     void AddSource(const std::string& prefix, const std::vector<std::string>& cols, const std::string& filter = "") {
//         if(cols.size() != _colNames.size()) {
//             throw std::runtime_error("ParticleInjector Error: Source column count mismatch. " 
//                                      "Expected " + std::to_string(_colNames.size()) + 
//                                      ", got " + std::to_string(cols.size()));
//         }
        
//         // Generate a unique ID for this source to ensure temporary names don't collide.
//         // e.g. "rec_px_src0", "rec_px_src1"
//         std::string sourceID = "_src" + std::to_string(_sources[prefix].size()); 
//         std::vector<std::string> registeredCols;
        
//         for(size_t i=0; i<cols.size(); ++i) {
//             std::string tempName = prefix + _colNames[i] + sourceID;
//             std::string inputCol = cols[i];
            
//             // Define the temporary column in RDataFrame
//             if(!filter.empty()) {
//                 // Define filtered version: Define(tempName, "inputCol[filter]")
//                 _reaction->Define(tempName, inputCol + "[" + filter + "]");
//             } else {
//                 // Simple Alias: Define(tempName, inputCol)
//                 // Note: We use Define (alias) to expose the column under the temporary name
//                 _reaction->Define(tempName, inputCol); 
//             }
//             registeredCols.push_back(tempName);
//         }
        
//         // Store the temporary names for the final merge step
//         _sources[prefix].push_back(registeredCols);
//     }

//     /**
//      * @brief Finalizes the injection by creating the unified vectors.
//      * * @details
//      * Iterates through all defined properties (px, py, pid...) and all registered sources.
//      * Constructs a JIT string for `rad::util::Concatenate<TYPE>(...)` and defines the
//      * final column (e.g., "rec_px").
//      * * Uses the type specified in DefineParticleInfo directly in the template parameter.
//      */
//     void CreateUnifiedVectors() {
//         for(auto& [prefix, sources] : _sources) {
            
//             // Loop over each property (e.g., "px", then "py", then "pid")
//             for(size_t i=0; i<_colNames.size(); ++i) {
//                 std::string finalName = prefix + _colNames[i];
//                 std::string type      = _colTypes[_colNames[i]];
                
//                 // Collect the temporary column names from all sources for this property
//                 std::vector<std::string> colsToMerge;
//                 for(const auto& src : sources) {
//                     colsToMerge.push_back(src[i]);
//                 }
                
//                 // Convert vector of strings to comma-separated string: "col1, col2, col3"
//                 std::string args = util::ColumnsToStringNoBraces(colsToMerge);

//                 // DIRECT TEMPLATE INJECTION
//                 // Generates code: rad::util::Concatenate< type >( col1, col2, ... )
//                 // This relies on the JIT compiler to instantiate the correct template specialization.
//                 std::string call = "rad::util::Concatenate<" + type + ">(" + args + ")";
//                 // Define the final unified column
//                 _reaction->Define(finalName, call);
//             }
//         }
//     }

//   private:
//     ConfigReaction* _reaction;                 ///< Pointer to the reaction interface
//     std::vector<std::string> _colNames;        ///< Ordered list of variable names (e.g. px, py)
//     std::map<std::string, std::string> _colTypes; ///< Map of variable name -> C++ type string
    
//     // Nested Map Structure:
//     // Prefix (rec_) -> List of Sources -> List of Temporary Column Names
//     std::map<std::string, std::vector<std::vector<std::string>>> _sources;
//   };

// } // namespace rad

