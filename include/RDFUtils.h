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
    enum class ColType { Undef, Int, UInt, Float, Double, Short, Bool, Long };

   
    /**
     * @brief Helper to deduce simple enum types from RDataFrame type strings.
     */
    template<typename T>
      inline ColType DeduceColumnVectorType(T* const radf, const string& name) {
        TString col_type = radf->ColObjTypeString(name);
        
        if (col_type.Contains("UInt_t")  || col_type.Contains("uint"))   return ColType::UInt;
        if (col_type.Contains("Float_t") || col_type.Contains("float"))  return ColType::Float;
        if (col_type.Contains("Double_t")|| col_type.Contains("double")) return ColType::Double;
        if (col_type.Contains("Short_t") || col_type.Contains("short"))  return ColType::Short;
        if (col_type.Contains("Bool_t")  || col_type.Contains("bool"))   return ColType::Bool;
        if (col_type.Contains("Long_t")  || col_type.Contains("long"))   return ColType::Long;
        if (col_type.Contains("Int_t")   || col_type.Contains("int"))    return ColType::Int;

        throw std::runtime_error(std::string("RDFInterface : DeduceColumnVectorType cannot deduce a type for ") + 
                                 name + " which is " + col_type.Data()); 
        return ColType::Undef;
    }

    // =========================================================================
    // Type Deduction Helper
    // =========================================================================
    /**
     * @brief deduced simple enum types from RDataFrame type strings.
     * @param col_type The type string returned by RDF (e.g. "vector<double>", "Int_t")
     */
    inline ColType DeduceTypeFromString(const std::string& typeStr) {
        // Use TString for easy case-insensitive checks if needed, or string find
        TString col_type = typeStr.c_str(); // Adapter to existing logic
        
        if (col_type.Contains("UInt_t")  || col_type.Contains("uint"))   return ColType::UInt;
        if (col_type.Contains("Float_t") || col_type.Contains("float"))  return ColType::Float;
        if (col_type.Contains("Double_t")|| col_type.Contains("double")) return ColType::Double;
        if (col_type.Contains("Short_t") || col_type.Contains("short"))  return ColType::Short;
        if (col_type.Contains("Bool_t")  || col_type.Contains("bool"))   return ColType::Bool;
        if (col_type.Contains("Long_t")  || col_type.Contains("long"))   return ColType::Long;
        if (col_type.Contains("Int_t")   || col_type.Contains("int"))    return ColType::Int;

        throw std::runtime_error("RDFUtils: Cannot deduce a type for " + typeStr); 
        return ColType::Undef;
    }

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
