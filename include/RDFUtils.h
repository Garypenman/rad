/**
 * @file RDFUtils.h
 * @brief Common types and utilities for RDataFrame interactions.
 */
#pragma once

#include <string>
#include <stdexcept>
#include <TString.h> // For TString operations if preferred, or std::string

namespace rad {
  
  using ROOT::RDF::RNode;

  // =========================================================================
    // Shared Column Type Enum
    // =========================================================================
  /**
   * @enum ColType
   * @brief internal enum to manage type casting for SnapshotCombi and Tree writing.
   */
  enum class ColType { 
    Undef, 
    Int,      // int, Int_t
    UInt,     // unsigned int, UInt_t
    Float,    // float, Float_t
    Double,   // double, Double_t
    Short,    // short, Short_t, char, UChar_t (mapped to smallest bucket)
    Bool,     // bool, Bool_t
    Long,     // long, Long_t
    Long64    // long long, Long64_t, ULong64_t
  };
      /**
     * @brief deduced simple enum types from RDataFrame type strings.
     * @param col_type The type string returned by RDF (e.g. "vector<double>", "Int_t")
     * @details Order matters! Specific types (Long64, UInt) are checked before generic ones (Long, Int).
     */
    inline ColType DeduceTypeFromString(const std::string& typeStr) {
        // Use TString for easy case-insensitive checks if needed, or string find
        TString col_type = typeStr.c_str(); // Adapter to existing logic

	// 1. Floating Point (High Priority)
        if (col_type.Contains("Double", TString::kIgnoreCase)) return ColType::Double;
        if (col_type.Contains("Float",  TString::kIgnoreCase)) return ColType::Float;

        // 2. 64-bit Integers (Check BEFORE "Long" or "Int")
        // Matches "Long64_t", "long long", "ULong64_t"
        if (col_type.Contains("Long64",    TString::kIgnoreCase) || 
            col_type.Contains("long long", TString::kIgnoreCase)) return ColType::Long64;

        // 3. Long Integers (Check AFTER "Long64")
        if (col_type.Contains("Long", TString::kIgnoreCase)) return ColType::Long;

        // 4. Unsigned Integers (Check BEFORE "Int")
        // Matches "UInt_t", "unsigned int", "unsigned"
        if (col_type.Contains("UInt",         TString::kIgnoreCase) || 
            col_type.Contains("unsigned int", TString::kIgnoreCase)) return ColType::UInt;

        // 5. Standard Integers (Check AFTER "UInt")
        if (col_type.Contains("Int", TString::kIgnoreCase)) return ColType::Int;

        // 6. Shorts (Matches "Short_t", "short", "unsigned short")
        // We map both signed and unsigned short to Short to preserve 16-bit size.
        if (col_type.Contains("Short", TString::kIgnoreCase)) return ColType::Short;

        // 7. Booleans
        if (col_type.Contains("Bool", TString::kIgnoreCase)) return ColType::Bool;

   
        throw std::runtime_error("RDFUtils: Cannot deduce a type for " + typeStr); 
        return ColType::Undef;
    }
  /**
   * @brief Helper to deduce simple enum types from RDataFrame type strings.
   */
  template<typename T>
  inline ColType DeduceColumnVectorType(T* const radf, const string& name) {
    // Get the type string (e.g., "std::vector<float>", "ROOT::VecOps::RVec<Long64_t>")
    std::string col_type = radf->ColObjTypeString(name);
    return DeduceTypeFromString(col_type);
        
  }

    // =========================================================================
    // Type Deduction Helper
    // =========================================================================


  void PrintDefinedColumnNames(RNode  df){
    std::cout<<"Print Column Names : ";
    auto cols =  df.GetDefinedColumnNames();
    for(auto& col:cols){
      std::cout<<col<<", ";
    }
    cout<<"\n";
  }
  void PrintAllColumnNames(RNode  df){
    std::cout<<"Print Column Names : ";
    auto cols =  df.GetColumnNames();
    for(auto& col:cols){
      std::cout<<col<<", ";
    }
    cout<<"\n";
  }

  
} // namespace rad
