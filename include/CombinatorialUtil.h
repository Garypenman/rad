#pragma once

#include <ROOT/RVec.hxx>
#include <type_traits> // For std::invoke_result_t
#include <utility>     // For std::forward
#include "CommonDefines.h" 

/**
 * @file CombinatorialUtil.h
 * @brief Core execution templates for vectorized combinatorial analysis.
 * * @details
 * This file contains the generic "Engine" functions that drive the analysis.
 * They bridge the gap between the "Structure of Arrays" (SoA) layout of RDataFrame columns
 * and the "Array of Structures" (AoS) logic often needed by physics kernels.
 * * **Key Concept: Vectorization**
 * Instead of users writing `for(int i=0; i<n_combos; ++i)`, these templates 
 * handle the looping internally. This allows the compiler (and ROOT's JIT) to 
 * see the full loop structure, potentially enabling SIMD optimizations.
 */

namespace rad {
  namespace util {

    //---------------------------------------------------------
    // Generic Combinatorial Helper Pattern
    //---------------------------------------------------------
 
    /**
     * @brief Vectorizes a generic function across all unique combinations in an event.
     * * @details
     * This is the master loop for general logic (e.g. non-kinematic functions).
     * It reconstructs the exact particle arguments for a given combination 
     * from the structured input index maps.
     * * **Optimization Note:** Memory allocations for the index buffers are hoisted 
     * out of the loop to minimize overhead when processing events with many combinations.
     * * @tparam F Type of the function/lambda.
     * @tparam Args Types of additional arguments (e.g. momentum vectors).
     * * @param singleComboFunc The function to apply. Signature: `(const RVecIndices&, Args...)`.
     * @param combo_indices Structure holding candidate indices: `RVec<RVecCombis>[set][particle][combo]`.
     * @param args Additional data columns to pass to the function.
     * @return `ROOT::RVec<T>` Vector of results, one per combination.
     */ 
    template <typename F, typename... Args>
    auto ApplyCombinations(
			   F&& singleComboFunc,
			   const ROOT::RVec<RVecCombis>& combo_indices,
			   Args&&... args) {
      
      // Deduce return type T from the function signature
      using T = std::invoke_result_t<F, const RVecIndices&, Args...>;       
 
      // 1. Validation & Sizing
      if (combo_indices.empty()) return ROOT::RVec<T>{};

      std::vector<UInt_t> n_particles;
      n_particles.reserve(combo_indices.size());

      UInt_t n_combos = 0;
      bool first_set = true;

      for(const auto& indices : combo_indices){
          n_particles.push_back(indices.size()); // Number of particles in this set
          
          if(indices.empty()) continue; 
          
          auto n = indices[0].size(); // Number of combinations
          
          if(!first_set && n != n_combos){
              throw std::runtime_error("ApplyCombinations : Index sets have inconsistent combination counts: " 
                                       + std::to_string(n_combos) + " vs " + std::to_string(n));
          }
          if(first_set) {
              n_combos = n;
              first_set = false;
          }
      }
      
      if (n_combos == 0) return ROOT::RVec<T>{};

      // 2. Execution Loop
      ROOT::RVec<T> results(n_combos);
      const auto n_idx_sets = n_particles.size();
     
      // -- PERF: Hoist allocation outside the loop --
      // We reuse these buffers for every combination to save alloc/free cycles.
      RVecIndices indices_for_F(n_idx_sets);
      // Pre-allocate the inner vectors
      for(size_t i=0; i<n_idx_sets; ++i) {
          indices_for_F[i].resize(n_particles[i]);
      }

      for (size_t i_combo = 0; i_combo < n_combos; ++i_combo) {
          
        // Reconstruct the index mapping for this specific combination
        for (size_t i_idx = 0; i_idx < n_idx_sets; ++i_idx) {
            // Note: We access the pre-allocated vector directly
            auto& current_indices = indices_for_F[i_idx]; 
            
            for(size_t i_particle = 0; i_particle < n_particles[i_idx]; ++i_particle){
                current_indices[i_particle] = combo_indices[i_idx][i_particle][i_combo];
            }
        }
        
        // Execute User Function
        results[i_combo] = std::forward<F>(singleComboFunc)(
                                indices_for_F, 
                                std::forward<Args>(args)...
                            );
      }
      
      return results;
    }

    /**
     * @brief Specialized vectorization loop for Kinematic Kernels.
     * * @details
     * This is the high-performance loop for physics calculations (Mass, Q2, etc.).
     * Unlike `ApplyCombinations`, it is specialized to extract Momentum (Px,Py,Pz,M)
     * efficiently from the transposed cache provided by `KinematicsProcessor`.
     * * @tparam F The kernel function type.
     * * @param singleComboKernel The kernel to call. Signature: `(map, px, py, pz, m)`.
     * @param fixed_map The static ReactionMap (fixed indices for particles).
     * @param components The Transposed Kinematic Cache: `[Combo][Component][Particle]`.
     * @return `ROOT::RVec<ResultType_t>` Vector of scalar results.
     */
    template <typename F>
    ROOT::RVec<ResultType_t> ApplyKinematics(
        F&& singleComboKernel,
        const RVecIndexMap& fixed_map,
        const ROOT::RVec<ROOT::RVec<RVecResultType>>& components) 
    {
      // Safety check: Ensure components exist
      if(components.empty()) {
          return {};
      }

      const size_t N_combis = components.size(); 
      // const size_t N_components = components[0].size(); // Expected to be 4 (Px,Py,Pz,M)
               
      ROOT::RVec<ResultType_t> results(N_combis);
      
      // --- Core Execution Loop ---
      for (size_t i_combo = 0; i_combo < N_combis; ++i_combo) {
          
          // We pass the slices directly to the kernel.
          // components[i_combo] is RVec<RVecResultType> of size 4.
          // Index 0=Px, 1=Py, 2=Pz, 3=M (See CommonDefines.h: OrderX...)
          
          results[i_combo] = singleComboKernel(fixed_map, 
                         components[i_combo][OrderX()],
                         components[i_combo][OrderY()],
                         components[i_combo][OrderZ()],
                         components[i_combo][OrderM()]);
      }

      return results;
    }
    
  } // namespace util
} // namespace rad

// #pragma once

// #include <ROOT/RVec.hxx>
// #include <type_traits> // For std::invoke_result_t
// #include <utility>     // For std::forward
// //#include "ConfigReaction.h" // For config::RVecIndexMap definition
// #include "CommonDefines.h" 

// namespace rad {
//   namespace util {

//     //---------------------------------------------------------
//     // Generic Combinatorial Helper Pattern
//     //---------------------------------------------------------
 
//     /**
//      * @brief Vectorizes a single-combination kinematic function across all unique combinations in an event.
//      * * This template acts as the central execution engine for combinatorial analysis. It reconstructs 
//      * the exact particle index arguments for a given combination from the structured input data, 
//      * and calls the user-provided kinematic function once for each combination, returning a vector
//      * of results.
//      * * @tparam F The type of the single-combination kinematic function (must be explicitly instantiated 
//      * with C++ types, e.g., FourVectorMassCalc<double, double>).
//      * @tparam Args The types of any additional, non-indexed kinematic arguments (e.g., RVec<double> px, py, etc.).
//      * * @param singleComboFunc The function to apply (must accept const RVecIndices& as its first argument).
//      * @param combo_indices The structured input holding all candidate indices, for each set of indices (e.g. addition particles, subtraction particles): 
//      * RVec<RVecCombis>[set][particle][combo].
//      * @param args Any additional columns/data needed by the kinematic function (e.g., momentum components).
//      * * @return ROOT::RVec<T> A vector where T is the return type of singleComboFunc, containing one result per combination.
//      * * @note This function ensures that the number of combinations is consistent across all index sets.
//      * @see RVecCombis, RVecIndices
//      */ 
//     template <typename F, typename... Args>
//     auto ApplyCombinations(
// 			   F&& singleComboFunc,
// 			   const ROOT::RVec<RVecCombis>& combo_indices,
// 			   Args&&... args) {
//       // Determine the return type of the single-combo function
//       using T = std::invoke_result_t<F, const RVecIndices&, Args...>;       
 
//       //check how many sets of combi indices this functions requires
//       std::vector<UInt_t> n_particles;
//       UInt_t n_combos = 0;
//       for(const auto& indices:combo_indices){
// 	n_particles.push_back(indices.size()); //indices.size = Number of particle indices
// 	if(n_particles.back()==0 ) continue;
       
// 	auto n = indices[0].size();
// 	if(n_combos>0 && n!=n_combos){
// 	  throw std::runtime_error("ApplyCombinations : we have indices with different numbers of combos ! :" + std::to_string(n)+ ","+std::to_string(n_combos));
// 	}
       
// 	n_combos = n;  //indices[0].size = Number of combinations,
// 	//must be same for every element of combo_indices
//       }
     
//       // Initialize the result vector with the size of the combinations vector
//       // and deduced type T
//       ROOT::RVec<T> results(n_combos);

//       //should check combos1[0].size()==combos2[0].size()
//       const auto n_idx_sets=n_particles.size();
     
//       // Loop over all combinations 
//       for (size_t i_combo = 0; i_combo < n_combos; ++i_combo) {
// 	//Loop over sets of indices and apply the function

// 	// Build the consolidated RVec<RVecIndices> argument for F
// 	RVecIndices indices_for_F(n_idx_sets);
    
// 	for (size_t i_idx = 0; i_idx < n_idx_sets; ++i_idx) {
// 	  //sets refer to different categories of indices
// 	  //e.g. positive particles and negatice particles
// 	  Indices_t indices(n_particles[i_idx]);
	 
// 	  for(size_t i_particle=0;i_particle<n_particles[i_idx];++i_particle){
// 	    indices[i_particle]=combo_indices[i_idx][i_particle][i_combo];
// 	  }//particle

// 	  indices_for_F[i_idx]=indices;
// 	}//indices
// 	//now have all indices sorted for this combi
// 	results[i_combo] = std::forward<F>(singleComboFunc)(
// 							    indices_for_F, 
// 							    std::forward<Args>(args)...
// 							    );
	
//       }//combi
      
//       return results;
//     }

//     /**
//      * @brief Vectorizes a single-combination kinematic kernel over all calculated combinations.
//      * * This function loops over the combination index (i_combo) and extracts the scalar momentum 
//      * components from the consolidated arrays (kine_px, kine_py, etc.) before calling the 
//      * specific kinematic calculation kernel.
//      * * @tparam F The single-combination kernel function (e.g., Q2_CombiKernel).
//      * @tparam Args A pack of the consolidated component RVecs (kine_px, kine_py, etc.).
//      * * @param singleComboKernel The calculation kernel (must be explicitly instantiated/non-templated).
//      * @param fixed_map The static ReactionMap containing the fixed array indices for lookup.
//      * @param components The pack of consolidated momentum component RVecs (kine_px, kine_py, ...)
//      * -> ComponentRVecs[icombi][icomp][iparti]
//      * * @return ROOT::RVec<ResultType_t> A vector of scalar results, one per combination.
//      */
    
//     template <typename F>
//     ROOT::RVec<ResultType_t> ApplyKinematics(
//         F&& singleComboKernel,
//         const RVecIndexMap& fixed_map,
//         const ROOT::RVec<ROOT::RVec<RVecResultType>>& components) 
//     {
//         // Safety check: Ensure at least one component vector is present to determine size.
//       if(components.empty()) {
// 	return {};
//       }

//       const size_t N_combis = components.size(); //should be 4
//       const size_t N_components = components[0].size(); //should be 4
               
//       ROOT::RVec<ResultType_t> results(N_combis);
      
//       // --- Core Execution Loop ---
//       for (size_t i_combo = 0; i_combo < N_combis; ++i_combo) {
// 	// 1. Execute the single-combination kernel.
// 	// The kernel expects the four component RVecs (Px, Py, Pz, M) as separate arguments.
// 	// We pass the inner RVecs from the slice by index (0, 1, 2, 3).
        
// 	// NOTE: The kernel receives the particle-indexed RVecs (size N_particles) for the current combo.
        
// 	results[i_combo] = singleComboKernel(fixed_map, 
// 					     components[i_combo][OrderX()],components[i_combo][OrderY()],
// 					     components[i_combo][OrderZ()],components[i_combo][OrderM()]);
// 	////////////////////
//       }

//       return results;
//     }
    
//     /**
//      * @brief Vectorizes a single-combination kinematic kernel over all calculated combinations.
//      * * This function loops over the combination index (i_combo) and extracts the scalar momentum 
//      * components from the consolidated arrays (kine_px, kine_py, etc.) before calling the 
//      * specific kinematic calculation kernel.
//      * * @tparam F The single-combination kernel function (e.g., Q2_CombiKernel).
//      * @tparam Args A pack of the consolidated component RVecs (kine_px, kine_py, etc.).
//      * * @param singleComboKernel The calculation kernel (must be explicitly instantiated/non-templated).
//      * @param fixed_map The static ReactionMap containing the fixed array indices for lookup.
//      * @param components The pack of consolidated momentum component RVecs (kine_px, kine_py, ...)
//      * -> ComponentRVecs[icombi][icomp][iparti]
//      * * @return ROOT::RVec<ResultType_t> A vector of scalar results, one per combination.
//      */
    
//     // template <typename F>
//     // ROOT::RVec<ResultType_t> ApplyToParticles(
//     //     F&& singleComboKernel,
//     //     const ROOT::RVec<RVecIndices>& indices,// [group][particle][combi]
//     //     const ROOT::RVec<ROOT::RVec<RVecResultType>>& components) 
//     // {
//     //     // Safety check: Ensure at least one component vector is present to determine size.
//     //   if(components.empty()) {
//     // 	return {};
//     //   }

//     //   const size_t N_combis = components.size(); //should be 4
//     //   const size_t N_components = components[0].size(); //should be 4
               
//     //   ROOT::RVec<ResultType_t> results(N_combis);
      
//     //   // --- Core Execution Loop ---
//     //   for (size_t i_combo = 0; i_combo < N_combis; ++i_combo) {
//     // 	// 1. Execute the single-combination kernel.
//     // 	// The kernel expects the four component RVecs (Px, Py, Pz, M) as separate arguments.
//     // 	// We pass the inner RVecs from the slice by index (0, 1, 2, 3).
        
//     // 	// NOTE: The kernel receives the particle-indexed RVecs (size N_particles) for the current combo.
        
//     // 	results[i_combo] = singleComboKernel(fixed_map, 
//     // 					     components[i_combo][OrderX()],components[i_combo][OrderY()],
//     // 					     components[i_combo][OrderZ()],components[i_combo][OrderM()]);
//     // 	////////////////////
//     //   }

//     //   return results;
//     // }
    
//   } // namespace util
// } // namespace rad
