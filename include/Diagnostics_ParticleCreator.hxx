/**
 * @file Diagnostics_ParticleCreator.hxx
 * @brief Diagnostic methods for ParticleCreator.
 * @details Implements PrintReactionMap, PrintGroups, PrintCreationRegistry, etc.
 * This file is included at the end of ParticleCreator.h
 */

#pragma once

#include <iomanip>
#include <algorithm>
#include <vector>
#include <map>

// Assuming Diagnostics.h exists and contains DiagnosticsPrinter.
#include "Diagnostics.h" 

namespace rad {

/**
 * @brief Print the reaction map (particle names to indices).
 * @details Sorts the internal unordered map by index for a readable table output.
 */
inline void ParticleCreator::PrintReactionMap() const {
  diag::DiagnosticsPrinter::PrintSectionHeader("Reaction Map (ParticleCreator)", '=', 80);

  if (_nameIndex.empty()) {
    std::cout << "[empty - Run InitMap() first]" << std::endl;
    diag::DiagnosticsPrinter::PrintBlank();
    return;
  }

  // _nameIndex is unordered, so let's sort it by index for a clean print
  std::vector<std::pair<std::string, int>> sortedMap(_nameIndex.begin(), _nameIndex.end());
  std::sort(sortedMap.begin(), sortedMap.end(), 
    [](const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {
        return a.second < b.second;
    }
  );

  std::vector<int> widths = {25, 10, 15};
  diag::DiagnosticsPrinter::PrintTableRow({"Particle Name", "Index", "Type"}, widths);
  diag::DiagnosticsPrinter::PrintTableSeparator(widths);

  for (const auto& p : sortedMap) {
    // Determine type: check if it exists in the Created list
    bool isCreated = std::find(_p_names.begin(), _p_names.end(), p.first) != _p_names.end();
    std::string type = isCreated ? "Created" : "Input";
    
    diag::DiagnosticsPrinter::PrintTableRow(
        {p.first, std::to_string(p.second), type}, widths);
  }

  diag::DiagnosticsPrinter::PrintBlank();
  diag::DiagnosticsPrinter::PrintKeyValue("Total particles", std::to_string(_nameIndex.size()));
  diag::DiagnosticsPrinter::PrintKeyValue("Input particles", std::to_string(_inputNames.size()));
  diag::DiagnosticsPrinter::PrintKeyValue("Created particles", std::to_string(_p_names.size()));
  diag::DiagnosticsPrinter::PrintBlank();
}

/**
 * @brief Print all registered particle groups (Mesons, Baryons, custom).
 */
inline void ParticleCreator::PrintGroups() const {
  diag::DiagnosticsPrinter::PrintSectionHeader("Particle Groups", '=', 80);

  if (_explicit_groups.empty() && _config_groups.empty()) {
    std::cout << "[no groups defined]" << std::endl;
    diag::DiagnosticsPrinter::PrintBlank();
    return;
  }

  // Print Explicit Groups (Hardcoded lists)
  for (const auto& grp : _explicit_groups) {
    std::cout << "Order [" << grp.first << "] (Explicit): ";
    for (size_t i = 0; i < grp.second.size(); ++i) {
      if (i > 0) std::cout << ", ";
      std::cout << grp.second[i];
    }
    std::cout << std::endl;
  }

  // Print Config Groups (Resolved from ConfigReaction)
  for (const auto& grp : _config_groups) {
    std::cout << "Order [" << grp.first << "] (Config): " << grp.second << std::endl;
  }
  
  diag::DiagnosticsPrinter::PrintBlank();
}

/**
 * @brief Print all registered creation functions and their dependencies.
 */
inline void ParticleCreator::PrintCreationRegistry() const {
  diag::DiagnosticsPrinter::PrintSectionHeader("Creation Registry", '=', 80);

  if (_p_names.empty()) {
    std::cout << "[no created particles registered]" << std::endl;
    diag::DiagnosticsPrinter::PrintBlank();
    return;
  }

  // Iterate through parallel vectors
  for (size_t i = 0; i < _p_names.size(); ++i) {
    std::cout << "[" << i << "] Name: " << std::left << std::setw(15) << _p_names[i];
    
    // Note: _p_creators is a function pointer. We print the dependencies instead.
    
    std::cout << " | Indices: ";
    if (i < _p_dep_indices.size() && !_p_dep_indices[i].empty()) {
        const auto& layer = _p_dep_indices[i]; // This is RVecIndices = vector<vector<int>>
        for(size_t l=0; l < layer.size(); ++l) {
            if(l > 0) std::cout << " + "; 
            std::cout << "{";
            for(size_t k=0; k < layer[l].size(); ++k) {
                if(k>0) std::cout << ",";
                std::cout << layer[l][k];
            }
            std::cout << "}";
        }
    } else {
        // If InitMap/ResolveDependencies hasn't run yet, print the string names
        std::cout << "(Unresolved) ";
        if(i < _p_stru_depends.size()) {
            const auto& layer = _p_stru_depends[i];
             for(size_t l=0; l < layer.size(); ++l) {
                if(l > 0) std::cout << " + "; 
                std::cout << "{";
                for(size_t k=0; k < layer[l].size(); ++k) {
                    if(k>0) std::cout << ",";
                    std::cout << layer[l][k];
                }
                std::cout << "}";
            }
        }
    }
    std::cout << std::endl;
  }
  diag::DiagnosticsPrinter::PrintBlank();
}

/**
 * @brief Print all column aliases for input particles.
 * @details Shows the mapping between the internal alias (Prefix + Name + Suffix) 
 * and the actual column name in the source TTree.
 */
inline void ParticleCreator::PrintAliases() const {
  diag::DiagnosticsPrinter::PrintSectionHeader("Column Aliases (Input Particles)", '=', 80);

  if (_inputNames.empty()) {
    std::cout << "[no inputs defined yet - Run InitMap()]" << std::endl;
    diag::DiagnosticsPrinter::PrintBlank();
    return;
  }

  std::vector<int> widths = {35, 35};
  diag::DiagnosticsPrinter::PrintTableRow({"Alias Name (In Analysis)", "Source Column (In Tree)"}, widths);
  diag::DiagnosticsPrinter::PrintTableSeparator(widths);

  // Replicate the logic from InitMap to show what aliases are generated
  for(const auto& name : _inputNames) {
      std::string colName = _prefix + name + _suffix; // The alias used in the loop
      std::string masterCol = _prefix + name;         // The original data
      
      // Only print if there is actually a difference (e.g. if suffix is used)
      if (colName != masterCol) {
         diag::DiagnosticsPrinter::PrintTableRow({colName, masterCol}, widths);
      } else {
         diag::DiagnosticsPrinter::PrintTableRow({colName, "(Same)"}, widths);
      }
  }
  diag::DiagnosticsPrinter::PrintBlank();
}

/**
 * @brief Comprehensive diagnostic dump of ParticleCreator state.
 */
inline void ParticleCreator::PrintDiagnostics() const {
  PrintReactionMap();
  PrintGroups();
  PrintAliases();
  PrintCreationRegistry();
}

} // namespace rad
