#pragma once

#include "DefineNames.h"

#include <TString.h>
#include <ROOT/RDFHelpers.hxx>

namespace rad{
    namespace util{

      using rad::consts::data_type::Rec;
      using rad::consts::data_type::Truth;

      /**
       * @brief Defines columns counting the occurrences of standard PIDs (gamma, pi, K, p, n, e).
       * * Creates columns like `rec_Ngamma`, `rec_Npip`, etc.
       * @tparam T The RDataFrame type.
       * @param rdf Pointer to the RDataFrame interface.
       * @param type Data prefix (e.g. "rec_").
       */
      template<typename T>
      void CountParticles(T* rdf, const std::string& type);

      /**
       * @brief Calculates simple resolution: (Truth - Rec).
       * * Requires inputs to be aliased/matched first.
       * Creates column `res_{var}`.
       */
      template<typename T> 
      void Resolution(T* const rdf, const std::string& var);
      
      /**
       * @brief Calculates fractional resolution: (Truth - Rec) / Truth.
       * Creates column `res_{var}`.
       */
      template<typename T> 
      void ResolutionFraction(T* const rdf, const std::string& var);
 
    } // namespace util
} // namespace rad

// =================================================================================
// IMPLEMENTATION
// =================================================================================

namespace rad {
  namespace util {

      template<typename T>
      void CountParticles(T* rdf, const std::string& type){
        rdf->Define(type+"Ngamma",   Form("rad::util::Count(%spid,22)",   type.data()) );
        rdf->Define(type+"Npip",     Form("rad::util::Count(%spid,211)",  type.data()) );
        rdf->Define(type+"Npim",     Form("rad::util::Count(%spid,-211)", type.data()) );
        rdf->Define(type+"NKp",      Form("rad::util::Count(%spid,321)",  type.data()) );
        rdf->Define(type+"NKm",      Form("rad::util::Count(%spid,-321)", type.data()) );
        rdf->Define(type+"Nele",     Form("rad::util::Count(%spid,11)",   type.data()) );
        rdf->Define(type+"Npos",     Form("rad::util::Count(%spid,-11)",  type.data()) );
        rdf->Define(type+"Npro",     Form("rad::util::Count(%spid,2212)", type.data()) );
        rdf->Define(type+"Nneutron", Form("rad::util::Count(%spid,2112)", type.data()) );
      }

      template<typename T> 
      void Resolution(T* const rdf, const std::string& var){
        rdf->Define(string("res_")+var, Form("%s-%s", (Truth()+var).data(), (Rec()+var).data() ));
      }

      template<typename T> 
      void ResolutionFraction(T* const rdf, const std::string& var){
        rdf->Define(string("res_")+var, Form("(%s-%s)/%s", (Truth()+var).data(), (Rec()+var).data(), (Truth()+var).data() ));
      }

  }
}
