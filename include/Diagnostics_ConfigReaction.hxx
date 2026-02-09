/**
 * @file Diagnostics_ConfigReaction.hxx
 * @brief Diagnostic methods for ConfigReaction.
 * @details Implements PrintParticleCandidates, PrintCombinatoricStructure, etc.
 *          This file is included at the end of ConfigReaction.h
 */

#pragma once

#include "Diagnostics.h"

namespace rad {

///\brief Print registered particle candidates for a given data type.
inline void ConfigReaction::PrintParticleCandidates(const std::string& type) const {
  using namespace diag;

  DiagnosticsPrinter::PrintSubsection("Candidates for type: " + type, '-', 80);

  std::vector<int> widths = {25, 50};
  DiagnosticsPrinter::PrintTableRow(
      {"Particle", "Definition"}, widths);
  DiagnosticsPrinter::PrintTableSeparator(widths);

  bool found = false;

  if (_typeCandidateExpressions.count(type)) {
    for (const auto& cand : _typeCandidateExpressions.at(type)) {
      std::string expr = cand.second;
      if (expr.size() > 47) expr = expr.substr(0, 47) + "...";
      DiagnosticsPrinter::PrintTableRow({cand.first, expr}, widths);
      found = true;
    }
  }

  if (_typeLambdaDependencies.count(type)) {
    for (const auto& cand : _typeLambdaDependencies.at(type)) {
      std::string colStr = "[lambda: ";
      const auto& cols = cand.second;
      for (size_t i = 0; i < cols.size() && i < 2; ++i) {
        if (i > 0) colStr += ", ";
        colStr += cols[i];
      }
      if (cols.size() > 2) colStr += ", ...";
      colStr += "]";
      DiagnosticsPrinter::PrintTableRow({cand.first, colStr}, widths);
      found = true;
    }
  }

  if (!found) {
    std::cout << "[no candidates defined for this type]" << std::endl;
  }
  DiagnosticsPrinter::PrintBlank();
}

///\brief Print the combinatorial structure for all registered data types.
inline void ConfigReaction::PrintCombinatoricStructure() const {
  using namespace diag;

  DiagnosticsPrinter::PrintSectionHeader("Combinatoric Structure", '=', 80);

  if (_types.empty()) {
    std::cout << "[no data types registered]" << std::endl;
    DiagnosticsPrinter::PrintBlank();
    return;
  }

  for (const auto& type : _types) {
    std::cout << "Type: " << type << std::endl;
    std::cout << "  Combo column: " << type << consts::ReactionCombos() << std::endl;
    std::cout << "  Combinatoric generation: ";
    if (_symmetryGroups.empty()) {
      std::cout << "GenerateAllCombinations";
    } else {
      std::cout << "GenerateSymmetricCombinations (with " << _symmetryGroups.size() << " groups)";
    }
    std::cout << std::endl;

    PrintParticleCandidates(type);
  }

  if (!_symmetryGroups.empty()) {
    DiagnosticsPrinter::PrintSubsection("Symmetry Group Details", '-', 80);
    for (size_t i = 0; i < _symmetryGroups.size(); ++i) {
      std::cout << "Group " << i << ": ";
      for (size_t j = 0; j < _symmetryGroups[i].size(); ++j) {
        if (j > 0) std::cout << ", ";
        std::cout << _symmetryGroups[i][j];
      }
      std::cout << std::endl;
    }
    DiagnosticsPrinter::PrintBlank();
  }
}

///\brief Print all registered data types and their P4 component definitions.
inline void ConfigReaction::PrintTypes() const {
  using namespace diag;

  DiagnosticsPrinter::PrintSectionHeader("Registered Data Types", '=', 80);

  if (_types.empty()) {
    std::cout << "[no data types registered]" << std::endl;
    DiagnosticsPrinter::PrintBlank();
    return;
  }

  std::vector<int> widths = {20, 15};
  DiagnosticsPrinter::PrintTableRow(
      {"Type", "Status"}, widths);
  DiagnosticsPrinter::PrintTableSeparator(widths);

  for (const auto& type : _types) {
    std::string status = _type_comps.count(type) ? "configured" : "pending";
    DiagnosticsPrinter::PrintTableRow({type, status}, widths);
  }

  DiagnosticsPrinter::PrintBlank();
  DiagnosticsPrinter::PrintKeyValue("Total types", std::to_string(_types.size()));
  DiagnosticsPrinter::PrintBlank();
}

///\brief Comprehensive diagnostic dump of ConfigReaction state.
inline void ConfigReaction::PrintReactionDiagnostics() const {
  PrintTypes();
  PrintCombinatoricStructure();
}

} // namespace rad
