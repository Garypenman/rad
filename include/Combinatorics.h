#pragma once

#include <ROOT/RVec.hxx>
#include <vector>
#include <string>
#include <map>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include "DefineNames.h" // For access to reaction component names
//#include "ConfigReaction.h" // For config::RVecIndexMap definition
#include "Constants.h" // For InvalidEntry
#include "CommonDefines.h"

namespace rad {
namespace combinatorics {

    using ROOT::RVecI;
    using ROOT::RVec;
    using std::vector;
    using std::string;
    using rad::constant::InvalidEntry;
    
    //---------------------------------------------------------
    // Core Combinatorial Logic Helper
    //---------------------------------------------------------
/**
     * @brief Generates unique combinations by enforcing that no index is repeated within a single combination.
     * * This filters out combinations resulting from overlapping candidate lists using an unordered_set 
     * for O(1) average-time uniqueness checking.
     * * @param candidates_vec RVec of candidate lists (RVecIndices[particle]).
     * @return RVecIndices where RVecIndices[particle] holds the particle index for every valid combination.
     */
    RVecIndices GenerateAllCombinations(const RVecIndices& candidates_vec) {
        
        const size_t n_particles = candidates_vec.size();
        if (n_particles == 0) {
            return RVecIndices{};
        }

        // --- 1. Initial Setup and Pre-calculation ---
        std::vector<size_t> max_indices;
        max_indices.reserve(n_particles);
        
        size_t total_combos = 1;
        for (const auto& vec : candidates_vec) {
            if (vec.empty()) return RVecIndices{};
            max_indices.push_back(vec.size());
            total_combos *= vec.size();
        }

        std::vector<size_t> current_indices(n_particles, 0);

        // --- 2. Initialize Transposed Output Structure ---
        RVecIndices result_by_particle(n_particles);
        // We cannot reserve accurately due to skipped combinations, but reserve total_combos
        // is the best estimate to minimize reallocations.
        for (size_t i = 0; i < n_particles; ++i) {
            result_by_particle[i].reserve(total_combos / n_particles); // Conservative estimate
        }

        // --- 3. Iterative Combinatorial Loop with Uniqueness Filter ---
        
        for (size_t i_combo = 0; i_combo < total_combos; ++i_combo) {
            
            // Temporary structure to hold the indices of the current combination
            std::vector<int> current_combination_indices(n_particles); 
            // The set performs the uniqueness check (O(1) average complexity)
            std::unordered_set<int> unique_indices;
            bool is_unique = true;

            // --- A. Build Combination and Check Uniqueness ---
            for (size_t i_particle = 0; i_particle < n_particles; ++i_particle) {
                
                int particle_index = candidates_vec[i_particle][current_indices[i_particle]];
                
                // Try to insert the index. If insert() returns {iterator, false}, the index already exists.
                if (!unique_indices.insert(particle_index).second) {
                    is_unique = false;
                    break; // Combination is invalid, stop building it.
                }
                current_combination_indices[i_particle] = particle_index;
            }

            // --- B. Record Valid Combination and Transpose ---
	    if (is_unique) {
                // If unique, transpose the indices into the result structure
                for (size_t i_particle = 0; i_particle < n_particles; ++i_particle) {
                    result_by_particle[i_particle].push_back(current_combination_indices[i_particle]);
                }
            }
            
            // --- C. Increment the N-dimensional counter (for the next loop iteration) ---
            size_t k = 0;
            while (k < n_particles) {
                current_indices[k]++;
                if (current_indices[k] < max_indices[k]) {
                    break;
                }
                current_indices[k] = 0;
                k++;
            }
        }
        return result_by_particle;
    }
  /**
     * @brief Generates all unique combinations (Cartesian Product) from the candidate lists
     * and transposes the output to be indexed by particle particle.
     * * * @param candidates_vec RVec of candidate lists (RVecIndices[particle]).
     * @return RVecIndices where RVecIndices[particle] returns an Indices_t containing 
     * the index of that particle for every combination. (size = N_combos).
     */
    // RVecIndices GenerateAllCombinations(const RVecIndices& candidates_vec) {
        
    //     const size_t n_particles = candidates_vec.size();
    //     if (n_particles == 0) {
    //         return RVecIndices{};
    //     }

    //     // --- 1. Initial Setup and Validation ---
    //     std::vector<size_t> max_indices;
    //     max_indices.reserve(n_particles);

    //     size_t total_combos = 1;
    //     for (const auto& vec : candidates_vec) {
    //         if (vec.empty()) {
    //             return RVecIndices{};
    //         }
    //         max_indices.push_back(vec.size());
    //         total_combos *= vec.size();
    //     }

    //     // Vector to track the current index for each particle particle (the counter)
    //     std::vector<size_t> current_indices(n_particles, 0);

    //     // --- 2. Initialize Transposed Output Structure ---
    //     // The output is RVecIndices[particle], where each element is an Indices_t of size total_combos
    //     RVecIndices result_by_particle(n_particles);
    //     for (size_t i = 0; i < n_particles; ++i) {
    //         result_by_particle[i].reserve(total_combos);
    //     }

    //     // --- 3. Iterative Combinatorial Loop (N-Dimensional Counter) ---
        
    //     for (size_t i_combo = 0; i_combo < total_combos; ++i_combo) {
            
    //         // --- A. Extract and Transpose Indices ---
    //         for (size_t i_particle = 0; i_particle < n_particles; ++i_particle) {
                
    //             // Get the single index from the candidate RVecI using the counter value
    //             // candidates_vec[i_particle] is the RVecI candidate list.
    //             int particle_index = candidates_vec[i_particle][current_indices[i_particle]];
                
    //             // Transpose: Append this index to the RVec dedicated to this specific particle.
    //             result_by_particle[i_particle].push_back(particle_index); 
    //         }
            
    //         // --- B. Increment the N-dimensional counter ---
    //         size_t k = 0;
    //         while (k < n_particles) {
    //             current_indices[k]++;
    //             if (current_indices[k] < max_indices[k]) {
    //                 break;
    //             }
    //             // Overflow: reset this counter and move to the next dimension
    //             current_indices[k] = 0;
    //             k++;
    //         }
    //     }
        
    //     // The output RVecIndices is now correctly structured: [particle][index_for_combo]
    //     return result_by_particle;
    // }

  
    /**
     * @brief Internal helper to generate all combinations from collected RVecI lists.
     * @param candidates_vec A vector of pointers to the Indices_t candidate for each particle.
     * @return RVecCombis combi indices for each particle, size = number of combos
     */
    //   RVecCombis GenerateAllCombinations(const RVecIndices& candidates_vec) {
        
    //     size_t n_particles = candidates_vec.size();
    // 	RVecCombis all_combos(n_particles);
    // 	std::cout<<"GenerateAllCombinations " <<n_particles<<" "<<candidates_vec<<std::endl;
    //     // Check for empty candidate lists
    //     for (const auto vec : candidates_vec) {
    //         if (vec.empty()) {
    //             return all_combos;
    //         }
    //     }

    //     // --- 1. Initial setup for indexing ---
    //     vector<size_t> current_indices(n_particles, InvalidEntry<int>());
    //     /* vector<size_t> max_indices; */
    //     /* for (const auto* vec : candidates_vec) { */
    //     /*     max_indices.push_back(vec->size()); */
    //     /* } */

    //     // --- 2. Iterative Combinatorial Loop (N-Dimensional Counter) ---
    //     bool done = false;

    //     while (!done) {
            
    // 	  // --- A. Build the RVecIndexMap for the current combination ---
    // 	  // RVecI current_map(n_candidates); 
            
    //         for (size_t iparticle = 0; iparticle < n_particles; ++iparticle) {
    //             // Get the index from the candidate RVecI at the position specified by current_indices[i]
    // 	      int particle_index = (candidates_vec[iparticle])[0];
              
    // 	      // Store the particle index for this combination.
    // 	      // We use the RVecI{index} wrapper because RVecIndexMap is RVec<RVecI>.
    // 	      //current_map[i] = RVecI{particle_index};
	      
    // 	      //set the particle index for this combo
    // 	      all_combos[iparticle].push_back(particle_index);
    //         }
            
    // 	    //all_combos.push_back(current_map);
    // 	    done=true;
    // 	    /*
    //         // --- B. Increment the N-dimensional counter ---
    //         size_t k = 0;
    //         while (k < n_candidates) {
    //             current_indices[k]++;
    //             if (current_indices[k] < max_indices[k]) {
    //                 break;
    //             }
    //             current_indices[k] = 0;
    //             k++;
    //         }

    //         if (k == n_candidates) {
    //             done = true;
    //         }
    // 	    */
    //     }
    // 	std::cout<<"GenerateAllCombinations " <<all_combos<<std::endl;

    //     return all_combos;
    // }
    
   /* RVec<RVecIndexMap> GenerateAllCombinations_impl(const std::vector<const RVecI*>& candidates_vec) { */
        
   /*      RVec<RVecIndexMap> all_combos; */
   /*      size_t n_candidates = candidates_vec.size(); */

   /*      // Check for empty candidate lists */
   /*      for (const auto* vec : candidates_vec) { */
   /*          if (vec->empty()) { */
   /*              return all_combos; */
   /*          } */
   /*      } */

   /*      // --- 1. Initial setup for indexing --- */
   /*      vector<size_t> current_indices(n_candidates, InvalidEntry<int>()); */
   /*      /\* vector<size_t> max_indices; *\/ */
   /*      /\* for (const auto* vec : candidates_vec) { *\/ */
   /*      /\*     max_indices.push_back(vec->size()); *\/ */
   /*      /\* } *\/ */

   /*      // --- 2. Iterative Combinatorial Loop (N-Dimensional Counter) --- */
   /*      bool done = false; */

   /*      while (!done) { */
            
   /*          // --- A. Build the RVecIndexMap for the current combination --- */
   /*          // RVecIndexMap is defined as ROOT::RVec<ROOT::RVecI> */
   /*          RVecIndexMap current_map(n_candidates);  */
            
   /*          for (size_t i = 0; i < n_candidates; ++i) { */
   /*              // Get the index from the candidate RVecI at the position specified by current_indices[i] */
   /*              int particle_index = (*candidates_vec[i])[current_indices[i]]; */
                
   /*              // Store the particle index for this combination. */
   /*              // We use the RVecI{index} wrapper because RVecIndexMap is RVec<RVecI>. */
   /*              current_map[i] = RVecI{particle_index};  */
   /*          } */
            
   /*          all_combos.push_back(current_map); */
   /* 	    done=true; */
   /* 	    /\* */
   /*          // --- B. Increment the N-dimensional counter --- */
   /*          size_t k = 0; */
   /*          while (k < n_candidates) { */
   /*              current_indices[k]++; */
   /*              if (current_indices[k] < max_indices[k]) { */
   /*                  break; */
   /*              } */
   /*              current_indices[k] = 0; */
   /*              k++; */
   /*          } */

   /*          if (k == n_candidates) { */
   /*              done = true; */
   /*          } */
   /* 	    *\/ */
   /*      } */
        
   /*      return all_combos; */
   /*  } */

    //---------------------------------------------------------
    // RDataFrame Interface (Variadic Template)
    //---------------------------------------------------------

    /**
     * @brief Variadic template wrapper to accept multiple RVecI columns from RDataFrame.
     * @tparam Args A pack of types (all expected to be const ROOT::RVecI&).
     * @param candidates A variadic list of RVecI candidate index lists.
     * @return RVec<RVecIndexMap> containing all unique index combinations.
     */
    //  template<typename... Args>
    //RVec<RVecIndexMap> GenerateAllCombinations(const Args&... candidates) {
    /* RVec<RVecIndexMap> GenerateAllCombinations(const RVec<RVecI>& candidates) { */
        
    /*     // Use a vector of pointers to collect all RVecI inputs */
    /*   // std::vector<const RVecI*> candidates_vec; */

    /*     // Fold expression (C++17) to push all arguments into the vector */
    /*     // This effectively collects the RVecI columns passed by RDataFrame */
    /*     (candidates_vec.push_back(&candidates), ...); */
        
    /*     // Delegate the actual heavy lifting to the internal helper */
    /*     return GenerateAllCombinations_impl(candidates_vec); */
    /* } */

} // namespace combinatorics
} // namespace rad
