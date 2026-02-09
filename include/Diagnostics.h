/**
 * @file Diagnostics.h
 * @brief Centralized diagnostics and introspection utilities for the rad framework.
 * @details Provides pretty-printing and state inspection for debugging and onboarding.
 */

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

namespace rad {
namespace diag {

  ///\brief Utility class for formatting diagnostic output
  class DiagnosticsPrinter {
  public:
    ///\brief Print a section header with separators
    static void PrintSectionHeader(const std::string& title, char sep = '=', int width = 80) {
      std::cout << std::string(width, sep) << std::endl;
      std::cout << "  " << title << std::endl;
      std::cout << std::string(width, sep) << std::endl;
    }

    ///\brief Print a subsection header
    static void PrintSubsection(const std::string& title, char sep = '-', int width = 80) {
      std::cout << std::string(width, sep) << std::endl;
      std::cout << title << std::endl;
      std::cout << std::string(width, sep) << std::endl;
    }

    ///\brief Print a key-value pair with aligned columns
    static void PrintKeyValue(const std::string& key, const std::string& value, int keyWidth = 30) {
      std::cout << std::left << std::setw(keyWidth) << key << " : " << value << std::endl;
    }

    ///\brief Print a table row
    static void PrintTableRow(const std::vector<std::string>& cols, const std::vector<int>& widths) {
      for (size_t i = 0; i < cols.size(); ++i) {
        if (i < widths.size()) {
          std::cout << std::left << std::setw(widths[i]) << cols[i] << " ";
        } else {
          std::cout << cols[i] << " ";
        }
      }
      std::cout << std::endl;
    }

    ///\brief Print a horizontal table separator
    static void PrintTableSeparator(const std::vector<int>& widths) {
      for (size_t i = 0; i < widths.size(); ++i) {
        std::cout << std::string(widths[i], '-') << "-";
      }
      std::cout << std::endl;
    }

    ///\brief Print empty line
    static void PrintBlank() {
      std::cout << std::endl;
    }
  };

} // namespace diag
} // namespace rad
