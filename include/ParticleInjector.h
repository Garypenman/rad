#pragma once

#include "ConfigReaction.h"
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <ROOT/RVec.hxx> // Ensure VecOps is available for Concatenate

namespace rad {
namespace config {

    /**
     * @brief Configuration structure for a single source of particle data.
     */
    struct ParticleSource_t {
      std::string name;                 ///< Unique name for this source (e.g., "Tracks").
      ROOT::RDF::ColumnNames_t columns; ///< The input columns from the tree (e.g. "trk_px", "trk_py").
      std::string filterExpr;           ///< Optional cut string (e.g., "trk_chi2 < 5.0").
    };

    /**
     * @brief Aggregates particle data from multiple detector sources into unified arrays.
     * * This class solves the problem where particle candidates are spread across multiple 
     * branches (e.g. "Tracks" vs "CalorimeterClusters" vs "Neutrals"). It allows the user to:
     * 1. Define a standard set of output columns (e.g., `rec_px`, `rec_py`).
     * 2. Register multiple sources that map to these columns.
     * 3. Apply source-specific filters (e.g. quality cuts) before merging.
     * * The result is a set of unified RVec columns containing all valid particles, 
     * ready for the combinatorial analysis.
     */
    class ParticleInjector {
    public:
      
      /**
       * @brief Constructor.
       * @param cr Pointer to the ConfigReaction object to define columns in.
       */
      explicit ParticleInjector(ConfigReaction* cr) : _reaction{cr} {}

      /**
       * @brief Defines the schema for the unified output vectors.
       * * Example: `{"rec_px", "rec_py", "rec_pz", "rec_m"}`.
       * All added sources must provide exactly this number of columns in this order.
       * @param info List of column names to create.
       */
      void DefineParticleInfo(const ROOT::RDF::ColumnNames_t& info) {
          _info = info;
      }

      /**
       * @brief Adds a source of particles to be merged.
       * * @param name A unique identifier for this source (used for internal column naming).
       * @param cols List of input columns from the RDataFrame tree. Must match the size/order of DefineParticleInfo.
       * @param filterExpr An optional RDataFrame expression string. Only particles satisfying this 
       * condition (mask) will be included in the final merged list. If empty, all entries are taken.
       * * @throws std::runtime_error If the column count does not match `DefineParticleInfo`.
       */
      void AddSource(const std::string& name, 
                     const ROOT::RDF::ColumnNames_t& cols, 
                     const std::string& filterExpr = "") 
      {
        if (cols.size() != _info.size()) {
            throw std::runtime_error("ParticleInjector: Source '" + name + 
                "' column count (" + std::to_string(cols.size()) + 
                ") does not match ParticleInfo (" + std::to_string(_info.size()) + ").");
        }
        _sources.push_back({name, cols, filterExpr});
      }

      /**
       * @brief Executes the definition logic in RDataFrame.
       * * This method:
       * 1. Iterates over all sources.
       * 2. Creates filtered temporary columns if a filter expression was provided.
       * 3. Defines the final Unified Columns (e.g. `rec_px`) by calling `Concatenate(...)` 
       * on the (filtered) source columns.
       */
      void CreateUnifiedVectors() {
        
        // 1. Process each source: Apply filters if necessary
        for(auto& source : _sources) {
            
            // If a filter is provided, we must create new filtered columns
            // so we don't merge rejected particles.
            if(!source.filterExpr.empty()) {
                std::string prefix = source.name + "_";
                
                // Define the mask column
                std::string maskName = prefix + "mask" + config::DoNotWriteTag();
                _reaction->Define(maskName, source.filterExpr);
                
                // Redefine source columns as filtered versions
                // We update source.columns to point to these new "inj_..." names
                ROOT::RDF::ColumnNames_t filtered_cols;
                
                for(const auto& col : source.columns) {
                    std::string newCol = prefix + col; 
                    
                    // Define: newCol = col[mask]
                    // We use a lambda to apply the boolean mask
                    _reaction->Define(newCol, 
                        [](const ROOT::RVecD& val, const ROOT::RVecI& mask) {
                            return val[mask]; 
                        }, {col, maskName});
                        
                    filtered_cols.push_back(newCol);
                }
                
                // Update the source struct to point to the filtered data
                source.columns = filtered_cols; 
            }
        }
        // 2. Merge Sources into Unified Vectors
        // For each info field (e.g. px, then py...), we concatenate all sources.
        for(size_t i=0; i<_info.size(); ++i) {

	  //If only 1 source no need to concatentate
	  if(_sources.size()==1){
	    std::string oldName = _sources[0].columns[i];
	    std::string newName = _sources[0].name + _info[i];

	    cout<<"ParticleInjector::CreateUnifiedVectors() create column "<<newName<<" from "<<oldName<<endl; 
	    _reaction->setBranchAlias(oldName, newName);
	    continue;
 	  }

	  //Multiple sources so need to concatentate them
	  // Build JIT string: "ROOT::VecOps::Concatenate(col1, col2, col3)"
	  std::string expr = "ROOT::VecOps::Concatenate(";
            
            for(size_t k=0; k<_sources.size(); ++k) {
                expr += _sources[k].columns[i];
                if(k < _sources.size()-1) expr += ", ";
            }
            expr += ")";
            
            // Define the Unified Column (e.g. "rec_px")
            // This implicitly handles type promotion (float -> double) if needed.
 	    cout<<"ParticleInjector::CreateUnifiedVectors() create column "<<_info[i]<<" from "<<expr<<endl; 
           _reaction->Define(_sources[0].name + _info[i], expr);
        }
      }

    private:
      ConfigReaction* _reaction = nullptr;
      ROOT::RDF::ColumnNames_t _info;       ///< Target column names (e.g. rec_px)
      std::vector<ParticleSource_t> _sources; ///< Registered sources
    };

} // namespace config
} // namespace rad
