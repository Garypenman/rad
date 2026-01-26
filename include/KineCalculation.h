/**
 * @file KineCalculation.h
 * @brief Wrapper class for configuring and registering custom kinematic calculations.
 */

#pragma once

#include "CommonDefines.h"
#include <string>
#include <vector>
#include <functional>

namespace rad {

  class KinematicsProcessor; // Forward declaration

  /**
   * @class KineCalculation
   * @brief Helper class to register custom kinematic calculations with the KinematicsProcessor.
   * @details
   * This class acts as a configuration store for physics calculations (Strategy Pattern). 
   * It bridges the gap between standalone C++ functions (kernels) and the RDataFrame execution loop.
   * * It supports two types of kernels:
   * 1. **Map Kernels:** Functions that have access to the full reaction topology (all particle indices).
   * 2. **Index Kernels:** Functions that operate on specific, pre-resolved groups of particles.
   * * Instead of executing immediately upon construction, this class stores the configuration
   * and applies it when `Define(processor)` is called. This allows the calculation to be
   * copied/cloned to other processors (e.g. Truth stream) where indices resolve differently.
   */
  class KineCalculation {
  public:
    
    // =========================================================================
    // Type Definitions
    // =========================================================================

    /**
     * @brief Signature for a Standard Map Kernel.
     * @details Receives the full `RVecIndexMap` containing all particle indices.
     */
    using MapKernel = ResultType_t(*)(const RVecIndexMap&, 
                                      const RVecResultType&, const RVecResultType&, 
                                      const RVecResultType&, const RVecResultType&);

    /**
     * @brief Signature for an Indexed Kernel.
     * @details Receives a specific list of pre-resolved indices (`RVecIndices`).
     */
    using IndexKernel = ResultType_t(*)(const RVecIndices&, 
                                        const RVecResultType&, const RVecResultType&, 
                                        const RVecResultType&, const RVecResultType&);


    // =========================================================================
    // Constructors
    // =========================================================================

    /**
     * @brief Constructor for a Standard Map Kernel (Global Topology).
     * @param name The name of the output column.
     * @param func The C++ function pointer to the MapKernel.
     */
    KineCalculation(std::string name, MapKernel func);

    /**
     * @brief Constructor for an Indexed Kernel (Generic Adapter).
     * @param name The name of the output column.
     * @param func The generic C++ function pointer (IndexKernel).
     * @param particles Vector of particle name lists to bind to the kernel.
     */
    KineCalculation(std::string name, IndexKernel func, std::vector<ParticleNames_t> particles);

    // =========================================================================
    // Execution
    // =========================================================================

    /**
     * @brief Resolves indices against the provided processor and registers the kernel.
     * @details
     * This is called inside `KinematicsProcessor::Init()`. It performs the resolution
     * of string particle names to integer indices specific to the processor instance
     * (e.g. Rec vs Truth) and registers the lambda with RDataFrame.
     * @param processor Pointer to the KinematicsProcessor instance.
     */
    void Define(KinematicsProcessor* processor);

  private:
    enum class KernelType { Map, Index };
    
    std::string _name;
    KernelType _kern_type;
    
    // Storage (Only one function pointer is active based on _type)
    MapKernel   _mapFunc = nullptr;
    IndexKernel _indexFunc = nullptr;
    std::vector<ParticleNames_t> _particles;
  };

  // =========================================================================
  // Implementation
  // =========================================================================

  inline KineCalculation::KineCalculation(std::string name, MapKernel func) 
      : _name(name), _mapFunc(func), _kern_type(KernelType::Map) {}

  inline KineCalculation::KineCalculation(std::string name, IndexKernel func, std::vector<ParticleNames_t> particles)
      : _name(name), _indexFunc(func), _particles(particles), _kern_type(KernelType::Index) {}

} // namespace rad

// #pragma once
// #include <string>
// #include <vector>
// #include <functional>
// #include "CommonDefines.h" 

// namespace rad {

//   /**
//    * @brief Helper class to register custom kinematic calculations with the KinematicsProcessor.
//    * * This class acts as a bridge between standalone C++ functions (kernels) and the 
//    * RDataFrame-based event loop. It simplifies the registration of custom physics calculations.
//    * * It supports two types of kernels:
//    * 1. **Map Kernels:** Functions that have access to the full reaction topology (all particle indices).
//    * 2. **Index Kernels:** Functions that operate on specific, pre-resolved groups of particles (e.g., calculating the mass of a specific pair).
//    * * Key features:
//    * - **Name Resolution:** Automatically converts string particle names to efficient integer indices at configuration time.
//    * - **Adaptation:** Wraps generic functions into the signature required by RDataFrame's `Define`.
//    * - **Registration:** Handles the call to `processor.DefineKernel`.
//    */
//   class KineCalculation {
//   public:
    
//     // =========================================================================
//     // Type Definitions
//     // =========================================================================

//     /**
//      * @brief Signature for a Standard Map Kernel.
//      * * This function signature receives the full `RVecIndexMap` containing all particle indices defined in the reaction.
//      * The function logic is responsible for looking up specific particles from this map.
//      * * **Prototype:**
//      * `ResultType_t Func(const RVecIndexMap& map, const RVecResultType& px, const RVecResultType& py, const RVecResultType& pz, const RVecResultType& m)`
//      */
//     using MapKernel = ResultType_t(*)(const RVecIndexMap&, 
//                                       const RVecResultType&, const RVecResultType&, 
//                                       const RVecResultType&, const RVecResultType&);

//     /**
//      * @brief Signature for an Indexed Kernel.
//      * * This function signature receives a specific list of pre-resolved indices (`RVecIndices`).
//      * This is ideal for generic algorithms (like invariant mass or missing mass) where the target particles are specified by the user.
//      * * **Prototype:**
//      * `ResultType_t Func(const RVecIndices& indices, const RVecResultType& px, const RVecResultType& py, const RVecResultType& pz, const RVecResultType& m)`
//      */
//     using IndexKernel = ResultType_t(*)(const RVecIndices&, 
//                                         const RVecResultType&, const RVecResultType&, 
//                                         const RVecResultType&, const RVecResultType&);


//     // =========================================================================
//     // Constructors
//     // =========================================================================

//     /**
//      * @brief Constructor for registering a Standard Map Kernel.
//      * * Use this constructor when your calculation requires knowledge of the full event topology
//      * or when particle lookups are embedded within the function logic itself.
//      * * @tparam ProcessorT The type of the KinematicsProcessor (e.g., `KinematicsProcElectro`).
//      * @param kine Reference to the KinematicsProcessor instance.
//      * @param name The name of the output column. The processor will automatically append the appropriate suffix (e.g., "_miss").
//      * @param func The C++ function pointer to the MapKernel.
//      */
//     template <typename ProcessorT>
//     KineCalculation(ProcessorT& kine, std::string name, MapKernel func) {
//        // The processor's DefineKernel handles the suffix generation internally.
//        // We simply pass the function pointer directly.
//        kine.DefineKernel(name, func);
//     }

//     /**
//      * @brief Constructor for registering an Indexed Kernel (Generic Adapter).
//      * * Use this constructor for generic functions where you want to specify the target particles at runtime
//      * (e.g., "Calculate Mass of {ele, pos}" or "Missing Mass of {Beam} - {Scat}").
//      * * **Mechanism:**
//      * 1. **Resolve:** Converts the provided `particle_groups` (strings) into integer indices using the processor's `ParticleCreator`.
//      * 2. **Adapt:** Creates a lambda function that captures these indices and calls your `func` with them.
//      * 3. **Register:** Registers this lambda as a new column in the RDataFrame.
//      * * @tparam ProcessorT The type of the KinematicsProcessor.
//      * @param kine Reference to the KinematicsProcessor instance.
//      * @param name The name of the output column.
//      * @param func The generic C++ function pointer (IndexKernel).
//      * @param particle_groups A vector of particle name lists. 
//      * Example: `{{ "ele", "pos" }, { "Jpsi" }}` creates a `RVecIndices` structure with two groups.
//      */
//     template <typename ProcessorT>
//     KineCalculation(ProcessorT& kine, std::string name, IndexKernel func, 
//                     std::vector<ParticleNames_t> particle_groups) 
//     {
//  	// 1. Resolve String Names -> Integer Indices
// 	// This happens once at configuration time.
//        RVecIndices resolved_indices;
//        for(const auto& group : particle_groups) {
//            Indices_t idxs;
//            for(const auto& pname : group) {
//              // Look up the fixed reaction index for the given particle name
// 	     idxs.push_back(kine.Creator().GetReactionIndex(pname));
//            }
//            resolved_indices.push_back(idxs);
//        }
//        // 2. Capture and Adapt (Lambda Generation)
//        // We create a lambda that matches the signature expected by DefineKernel (MapKernel style),
//        // but internally calls the IndexKernel using our pre-calculated 'resolved_indices'.
//        // The 'map' argument is unused here because indices are already resolved.
//        auto adapter_lambda = [resolved_indices, func](const RVecIndexMap&, 
//                                                       const RVecResultType& px, const RVecResultType& py, 
//                                                       const RVecResultType& pz, const RVecResultType& m) 
//        {
//            return func(resolved_indices, px, py, pz, m);
//        };

//        // 3. Register with RDataFrame
//        // The DefineKernel method will compile this lambda into the analysis loop.
//        kine.DefineKernel(name, adapter_lambda);
//     }

//   };

// } // namespace rad
