#pragma once

#include <ROOT/RVec.hxx>
#include <vector>
#include <string>
#include <map>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <unordered_set> 
#include "DefineNames.h" 
#include "Constants.h" 
#include "CommonDefines.h"

/**
 * @file Combinatorics.h
 * @brief Utilities for generating particle combinations.
 * * @details
 * This file contains the logic to generate the "Cartesian Product" of all particle candidates.
 * It converts a list of candidates (e.g., 2 electrons, 1 positron) into a set of 
 * unique combinatorial events.
 */

namespace rad {
  namespace combinatorics {

    using ROOT::RVecI;
    using ROOT::RVec;
    using std::vector;
    using std::string;
    using rad::consts::InvalidEntry;
    
    //---------------------------------------------------------
    // Core Combinatorial Logic Helper
    //---------------------------------------------------------
    
    /**
     * @brief Generates all valid, unique combinations of candidates.
     * * @details
     * This algorithm performs two main tasks:
     * 1. **Cartesian Product:** It iterates through every possible permutation of the input candidates.
     * (e.g., if Ele has 2 candidates and Pos has 3, it tests 2*3=6 combinations).
     * 2. **Uniqueness Filter:** It discards any combination where the same underlying detector object 
     * (index) is used for multiple roles (e.g., the same track used as both an electron and a pion).
     * * **Output Format:**
     * The output is a "Structure of Arrays" (SoA).
     * `result[particle_role_index]` is a vector containing the candidate index for that role 
     * for every valid combination.
     * * @param candidates_vec A vector of candidate lists. `candidates_vec[particle_role]` contains the 
     * list of valid indices for that role (e.g., {0, 2, 5} for electrons).
     * @return RVecIndices The structure of valid combinations.
     */
    inline RVecIndices GenerateAllCombinations(const RVecIndices& candidates_vec) {
        
      const size_t n_particles = candidates_vec.size();
      if (n_particles == 0) {
        return RVecIndices(n_particles);
      }

      // --- 1. Initial Setup and Pre-calculation ---
      // Determine the bounds for the N-dimensional counter
      std::vector<size_t> max_indices;
      max_indices.reserve(n_particles);
        
      size_t total_combos = 1;
      for (const auto& vec : candidates_vec) {
        if (vec.empty()) return RVecIndices(n_particles); // If any required particle has 0 candidates, 0 combos exist.
        max_indices.push_back(vec.size());
        total_combos *= vec.size();
      }

      // The N-dimensional counter (state of the current permutation)
      std::vector<size_t> current_indices(n_particles, 0);

      // --- 2. Initialize Transposed Output Structure ---
      RVecIndices result_by_particle(n_particles);
      // We cannot reserve accurately due to skipped combinations (overlaps), 
      // but reserving 'total_combos' avoids reallocations in the worst case.
      for (size_t i = 0; i < n_particles; ++i) {
        result_by_particle[i].reserve(total_combos); 
      }

      // --- 3. Iterative Combinatorial Loop ---
      
      // We reuse these buffers to avoid malloc/free overhead for every combination.
      std::vector<int> current_combination_buffer(n_particles);
      std::unordered_set<int> unique_checker;
      unique_checker.reserve(n_particles * 2); 

      for (size_t i_combo = 0; i_combo < total_combos; ++i_combo) {
            
        unique_checker.clear();
        bool is_unique = true;

        // --- A. Build Combination and Check Uniqueness ---
        for (size_t i_particle = 0; i_particle < n_particles; ++i_particle) {
                
          // Lookup the actual Data Index for this role in the current permutation
          int particle_index = candidates_vec[i_particle][current_indices[i_particle]];
                
          // Try to insert the index. If insert() returns {iterator, false}, the index is a duplicate.
          if (!unique_checker.insert(particle_index).second) {
            is_unique = false;
            break; // Stop checking this combination immediately
          }
          current_combination_buffer[i_particle] = particle_index;
        }

        // --- B. Record Valid Combination ---
        if (is_unique) {
          // Transpose: Append this combo's indices to the columnar output
          for (size_t i_particle = 0; i_particle < n_particles; ++i_particle) {
            result_by_particle[i_particle].push_back(current_combination_buffer[i_particle]);
          }
        }
            
        // --- C. Increment the N-dimensional counter ---
        // Simulates a nested loop of depth 'n_particles'
        size_t k = 0;
        while (k < n_particles) {
          current_indices[k]++;
          if (current_indices[k] < max_indices[k]) {
            break; // Carry propagation done
          }
          current_indices[k] = 0; // Reset this digit and carry over
          k++;
        }
      }
      return result_by_particle;
    }
 
  } // namespace combinatorics
} // namespace rad
