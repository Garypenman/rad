#pragma once

#include "RDFInterface.h" 
// Note: We do NOT include KinematicsProcessor.h here to avoid circular dependencies.
// The template <typename T> handles the processor type.

namespace rad {

  // Import helpers for column type deduction
  using rad::config::DeduceColumnVectorType;
  using rad::config::ColType;

  // Aliases for the Packed Aux Columns (Must match KinematicsProcessor)
  using RVecRVecD = ROOT::RVec<ROOT::RVecD>;
  using RVecRVecI = ROOT::RVec<ROOT::RVecI>;

  /**
   * @brief Registers the KinematicsProcessor execution loop with RDataFrame.
   * * This function:
   * 1. Constructs the column names (Indices, P4, and Packed Aux Data).
   * 2. Handles type deduction (Float vs Double) for the P4 columns.
   * 3. Calls RDataFrame::Define with a lambda that invokes processor.operator().
   */
  template <typename T>
  void DefineKinematicsProcessor(config::ConfigReaction& cr,
                                 T& processor, 
                                 const std::string& type) 
  {
    // 1. Get Suffix (e.g., "_miss")
    std::string suffix = processor.GetSuffix();

    // 2. Construct Column Names
    // Main Indices and Output Name
    auto indicesName = type + names::KineIndices() + suffix;
    auto name = type + names::KineComponents() + suffix;

    // Auxiliary Data Packs (Created by processor.DefineAux)
    auto auxPreD  = type + "aux_pre_d" + suffix+config::DoNotWriteTag();
    auto auxPreI  = type + "aux_pre_i" + suffix+config::DoNotWriteTag();
    auto auxPostD = type + "aux_post_d" + suffix+config::DoNotWriteTag();
    auto auxPostI = type + "aux_post_i" + suffix+config::DoNotWriteTag();

    // 3. Build Argument List for RDataFrame
    // Order: [Indices, Px, Py, Pz, M, PreD, PreI, PostD, PostI]
    ROOT::RDF::ColumnNames_t columns = {
        indicesName, 
        type + "px", type + "py", type + "pz", type + "m",
        auxPreD, auxPreI, auxPostD, auxPostI
    };

    // 4. Deduce Input Types (Float vs Double)
    auto vecTypeP = DeduceColumnVectorType(&cr, type + "px");
    auto vecTypeM = DeduceColumnVectorType(&cr, type + "m");
    
    // Helper to generate the define lambda to reduce code duplication
    // T_Px, T_M: The types of the momentum and mass vectors (RVecF or RVecD)
    auto do_define = [&](auto dummy_px, auto dummy_m) {
        
        using TypePx = decltype(dummy_px);
        using TypeM  = decltype(dummy_m);

        cr.Define(name, 
            [processor](const RVecIndices& indices, 
                        const TypePx& px, const TypePx& py, 
                        const TypePx& pz, const TypeM& m,
                        const RVecRVecD& apd, const RVecRVecI& api,
                        const RVecRVecD& bpd, const RVecRVecI& bpi) 
            {
                // Invoke the template operator() on the processor
                return processor.template operator()<TypePx, TypeM>(
                    indices, px, py, pz, m, apd, api, bpd, bpi);
            }, 
            columns
        );
    };

    // 5. Dispatch based on detected types
    if (vecTypeP == ColType::Float && vecTypeM == ColType::Double) {
        do_define(ROOT::RVecF{}, ROOT::RVecD{});
    }
    else if (vecTypeP == ColType::Double && vecTypeM == ColType::Double) {
        do_define(ROOT::RVecD{}, ROOT::RVecD{});
    }
    else if (vecTypeP == ColType::Double && vecTypeM == ColType::Float) {
        do_define(ROOT::RVecD{}, ROOT::RVecF{});
    }
    else if (vecTypeP == ColType::Float && vecTypeM == ColType::Float) {
        do_define(ROOT::RVecF{}, ROOT::RVecF{});
    }
    else {
        throw std::runtime_error("KinematicsDispatch: Cannot deduce template types for " + type + " " + name);
    }

    // 6. Optional Hook: Unpack components if needed
    processor.DefineNewComponentVecs();
     
    return;
  }

} // namespace rad


// #pragma once

// #include "RDFInterface.h" 
// // We don't include KinematicsProcessor.h directly to avoid circular dependencies
// // if KinematicsProcessor includes this file. The template handles the type.

// namespace rad {

//   // Import helpers for column type deduction
//   using rad::config::DeduceColumnVectorType;
//   using rad::config::ColType;

//   /**
//    * @brief Registers the KinematicsProcessor execution loop with RDataFrame.
//    * * This function handles the type deduction (Float vs Double) of the input columns
//    * and defines the main output column ("kine_comps") which stores the transposed data.
//    * * @tparam T The Processor Type (KinematicsProcessor or derived).
//    * @param cr The ConfigReaction interface.
//    * @param processor The processor instance (captured by value for thread safety).
//    * @param type The data type prefix (e.g. "rec_", "tru_").
//    */
//   template <typename T>
//   void DefineKinematicsProcessor(config::ConfigReaction& cr,
//                                  T& processor, 
//                                  const std::string& type) 
//   {
//     // 1. Get Suffix (e.g., "_miss") to ensure unique column names
//     std::string suffix = processor.GetSuffix();

//     // 2. Construct Column Names
//     // Input: The indices generated by ParticleCreator::InitMap (suffixed)
//     auto indicesName = type + names::KineIndices() + suffix;
    
//     // Output: The massive component array (suffixed)
//     auto name = type + names::KineComponents() + suffix;

//     // Columns required by the operator(): Indices + Raw P4 Components
//     ROOT::RDF::ColumnNames_t columns = {
//         indicesName, 
//         type + "px", 
//         type + "py", 
//         type + "pz", 
//         type + "m"
//     };

//     // 3. Deduce Input Types (Float vs Double)
//     // The Input arrays (px, py...) usually define the template instantiation.
//     auto vecTypeP = DeduceColumnVectorType(&cr, type + "px");
//     auto vecTypeM = DeduceColumnVectorType(&cr, type + "m");
    
//     // 4. Register the Define based on detected types
//     if (vecTypeP == ColType::Float && vecTypeM == ColType::Double) {
//       cr.Define(name, 
//         [processor](const RVecIndices& indices, 
//                     const ROOT::RVecF& px, const ROOT::RVecF& py, 
//                     const ROOT::RVecF& pz, const ROOT::RVecD& m) 
//         {
//           return processor.template operator()<ROOT::RVecF, ROOT::RVecD>(indices, px, py, pz, m);
//         }, 
//         columns
//       );
//     }
//     else if (vecTypeP == ColType::Double && vecTypeM == ColType::Double) {
//       cr.Define(name, 
//         [processor](const RVecIndices& indices, 
//                     const ROOT::RVecD& px, const ROOT::RVecD& py, 
//                     const ROOT::RVecD& pz, const ROOT::RVecD& m) 
//         {
//           return processor.template operator()<ROOT::RVecD, ROOT::RVecD>(indices, px, py, pz, m);
//         }, 
//         columns
//       );
//     }
//     else if (vecTypeP == ColType::Double && vecTypeM == ColType::Float) {
//       cr.Define(name, 
//         [processor](const RVecIndices& indices, 
//                     const ROOT::RVecD& px, const ROOT::RVecD& py, 
//                     const ROOT::RVecD& pz, const ROOT::RVecF& m) 
//         {
//           return processor.template operator()<ROOT::RVecD, ROOT::RVecF>(indices, px, py, pz, m);
//         }, 
//         columns
//       );
//     }
//     else if (vecTypeP == ColType::Float && vecTypeM == ColType::Float) {
//       cr.Define(name, 
//         [processor](const RVecIndices& indices, 
//                     const ROOT::RVecF& px, const ROOT::RVecF& py, 
//                     const ROOT::RVecF& pz, const ROOT::RVecF& m) 
//         {
//           return processor.template operator()<ROOT::RVecF, ROOT::RVecF>(indices, px, py, pz, m);
//         }, 
//         columns
//       );
//     }
//     else {
//       throw std::runtime_error("KinematicsDispatch: Cannot deduce template types for " + type + " " + name);
//     }

//     // 5. Trigger unpacking (Optional hook)
//     processor.DefineNewComponentVecs();
     
//     return;
//   }

// } // namespace rad

