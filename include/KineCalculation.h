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
    // Public Enums
    // =========================================================================

    /** * @brief Discriminator for the type of kernel stored.
     * @note Must be public for Diagnostic tools to inspect the type.
     */
    enum class KernelType { Map, Index };

    // =========================================================================
    // Type Definitions
    // =========================================================================

    /**
     * @brief Signature for a Standard Map Kernel.
     * @details Receives the full `RVecIndexMap` containing all particle indices.
     * The kernel is expected to return a single scalar result (e.g. Mass) for *one* combination.
     * @return ResultType_t (Scalar double).
     */
    using MapKernel = ResultType_t(*)(const RVecIndexMap&, 
                                      const RVecResultType&, const RVecResultType&, 
                                      const RVecResultType&, const RVecResultType&);

    /**
     * @brief Signature for an Indexed Kernel.
     * @details Receives a specific list of pre-resolved indices (`RVecIndices`).
     * The kernel is expected to return a single scalar result (e.g. Mass) for *one* combination.
     * @return ResultType_t (Scalar double).
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
    KineCalculation(std::string name, IndexKernel func, ROOT::RVec<ParticleNames_t> particles);

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

    // =========================================================================
    // Accessors
    // =========================================================================

    /** @return The base name of the calculation column. */
    std::string GetName() const;

    /** @return The type of kernel (Map or Index). */
    KernelType GetType() const;

  private:
    std::string _name;
    KernelType _kern_type;
    
    // Storage (Only one function pointer is active based on _type)
    MapKernel   _mapFunc = nullptr;
    IndexKernel _indexFunc = nullptr;
    
    ROOT::RVec<ParticleNames_t> _particles;
  };

  // =========================================================================
  // Implementation
  // =========================================================================

  inline KineCalculation::KineCalculation(std::string name, MapKernel func) 
      : _name(name), _mapFunc(func), _kern_type(KernelType::Map) {}

  inline KineCalculation::KineCalculation(std::string name, IndexKernel func, ROOT::RVec<ParticleNames_t> particles)
      : _name(name), _indexFunc(func), _particles(particles), _kern_type(KernelType::Index) {}

  inline std::string KineCalculation::GetName() const { return _name; }
  
  inline KineCalculation::KernelType KineCalculation::GetType() const { return _kern_type; }

} // namespace rad
