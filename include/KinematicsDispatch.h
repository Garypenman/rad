#pragma once

#include "RDFInterface.h" 
// Note: We do NOT include KinematicsProcessor.h here to avoid circular dependencies.
// The template <typename T> handles the processor type.

namespace rad {

  // Import helpers for column type deduction
  using rad::DeduceColumnVectorType;
  using rad::ColType;

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
  void DefineKinematicsProcessor(ConfigReaction& cr,
                                 T& processor, 
                                 const std::string& type) 
  {
    // 1. Get Suffix (e.g., "_miss")
    std::string suffix = processor.GetSuffix();

    // 2. Construct Column Names
    // Main Indices and Output Name
    auto indicesName = type + consts::KineIndices() + suffix;
    auto name = type + consts::KineComponents() + suffix;

    // Auxiliary Data Packs (Created by processor.DefineAux)
    auto auxPreD  = type + "aux_pre_d" + suffix+DoNotWriteTag();
    auto auxPreI  = type + "aux_pre_i" + suffix+DoNotWriteTag();
    auto auxPostD = type + "aux_post_d" + suffix+DoNotWriteTag();
    auto auxPostI = type + "aux_post_i" + suffix+DoNotWriteTag();

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
