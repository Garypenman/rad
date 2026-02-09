/**
 * @file Diagnostics_KinematicsProcessor.hxx
 * @brief Diagnostic methods for KinematicsProcessor.
 * @details Implements PrintCalculations, PrintRegisteredVariables, etc.
 * This file is included at the end of KinematicsProcessor.h
 */

#pragma once

#include "Diagnostics.h"

namespace rad {

/** @brief Print all registered calculations and their kernel types. */
inline void KinematicsProcessor::PrintCalculations() const {
  diag::DiagnosticsPrinter::PrintSectionHeader("Registered Calculations", '=', 80);

  if (_calculations.empty()) {
    std::cout << "[no calculations registered]" << std::endl;
    diag::DiagnosticsPrinter::PrintBlank();
    return;
  }

  std::vector<int> widths = {35, 15};
  diag::DiagnosticsPrinter::PrintTableRow(
      {"Calculation Name", "Type"}, widths);
  diag::DiagnosticsPrinter::PrintTableSeparator(widths);

  for (size_t i = 0; i < _calculations.size(); ++i) {
    const auto& calc = _calculations[i];
    std::string fullName = FullName(calc.GetName());
    
    // Assuming KineCalculation has a GetType() returning an Enum or Int
    // and KineCalculation::KernelType enum exists.
    std::string type = (calc.GetType() == KineCalculation::KernelType::Map) ? "MapKernel" : "IndexKernel";
    
    diag::DiagnosticsPrinter::PrintTableRow({fullName, type}, widths);
  }

  diag::DiagnosticsPrinter::PrintBlank();
  diag::DiagnosticsPrinter::PrintKeyValue("Total calculations", std::to_string(_calculations.size()));
  diag::DiagnosticsPrinter::PrintBlank();
}

/** @brief Print the list of variables registered for snapshot/output. */
inline void KinematicsProcessor::PrintRegisteredVariables() const {
  diag::DiagnosticsPrinter::PrintSectionHeader("Variables Registered for Output", '=', 80);

  if (_registered_vars.empty()) {
    std::cout << "[no variables registered]" << std::endl;
    diag::DiagnosticsPrinter::PrintBlank();
    return;
  }

  std::vector<int> widths = {5, 30, 40};
  diag::DiagnosticsPrinter::PrintTableRow({"ID", "Base Name", "Full Column Name"}, widths);
  diag::DiagnosticsPrinter::PrintTableSeparator(widths);

  for (size_t i = 0; i < _registered_vars.size(); ++i) {
    diag::DiagnosticsPrinter::PrintTableRow({
        std::to_string(i), 
        _registered_vars[i], 
        FullName(_registered_vars[i])
    }, widths);
  }

  diag::DiagnosticsPrinter::PrintBlank();
  diag::DiagnosticsPrinter::PrintKeyValue("Total variables", std::to_string(_registered_vars.size()));
  diag::DiagnosticsPrinter::PrintKeyValue("Prefix", _prefix);
  diag::DiagnosticsPrinter::PrintKeyValue("Suffix", _suffix.empty() ? "[none]" : _suffix);
  diag::DiagnosticsPrinter::PrintKeyValue("Initialized", _isInitialized ? "yes" : "no");
  diag::DiagnosticsPrinter::PrintBlank();
}

/** @brief Print group overrides (custom group definitions). */
inline void KinematicsProcessor::PrintGroupOverrides() const {
  diag::DiagnosticsPrinter::PrintSectionHeader("Group Overrides", '=', 80);

  if (_groupOverrides.empty()) {
    std::cout << "[no group overrides]" << std::endl;
    diag::DiagnosticsPrinter::PrintBlank();
    return;
  }

  for (const auto& ovr : _groupOverrides) {
    std::cout << ovr.name << ": ";
    for (size_t i = 0; i < ovr.particles.size(); ++i) {
      if (i > 0) std::cout << ", ";
      std::cout << ovr.particles[i];
    }
    std::cout << std::endl;
  }

  diag::DiagnosticsPrinter::PrintBlank();
}

/** @brief Comprehensive diagnostic dump of KinematicsProcessor state. */
inline void KinematicsProcessor::PrintProcessorDiagnostics() const {
  diag::DiagnosticsPrinter::PrintSectionHeader("KINEMATICS PROCESSOR DIAGNOSTICS", '=', 90);

  diag::DiagnosticsPrinter::PrintKeyValue("Processor type (prefix)", _prefix);
  diag::DiagnosticsPrinter::PrintKeyValue("Suffix", _suffix.empty() ? "[none]" : _suffix);
  diag::DiagnosticsPrinter::PrintKeyValue("Initialized", _isInitialized ? "yes" : "no");

  diag::DiagnosticsPrinter::PrintBlank();

  Creator().PrintDiagnostics();

  PrintGroupOverrides();
  PrintCalculations();
  PrintRegisteredVariables();
}

} // namespace rad
