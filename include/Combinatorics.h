/**
 * @file Combinatorics.h
 * @brief Combinatorial engines for generating particle pairs, triplets, etc.
 */

#pragma once

#include <ROOT/RVec.hxx>
#include <vector>
#include <tuple>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace rad {
namespace combinatorics {

  // Output Type: [ColumnIndex][CombinationIndex]
  using Combis_t = ROOT::RVec<ROOT::RVec<int>>;

  // =================================================================================
  // 1. STANDARD CARTESIAN PRODUCT (GenerateAllCombinations)
  // =================================================================================

  namespace impl {
    /**
     * @brief Recursive helper to build the Cartesian product.
     * @tparam I Current column index being processed.
     * @tparam ColTypes Variadic types of the input index vectors.
     */
    template <std::size_t I, typename... ColTypes>
    struct CartesianProduct {
      static void Generate(const std::tuple<ColTypes...>& cols, 
                           ROOT::RVec<int>& current_indices, 
                           Combis_t& result) {
        
        const auto& current_col = std::get<I>(cols);
        
        // Loop over the candidates in the current column
        for (const auto& idx : current_col) {
          current_indices[I] = idx;
          
          // Recurse to the next column
          CartesianProduct<I + 1, ColTypes...>::Generate(cols, current_indices, result);
        }
      }
    };

    // Base case: All columns processed, push to result
    template <typename... ColTypes>
    struct CartesianProduct<sizeof...(ColTypes), ColTypes...> {
      static void Generate(const std::tuple<ColTypes...>&, 
                           ROOT::RVec<int>& current_indices, 
                           Combis_t& result) {
        for (std::size_t i = 0; i < sizeof...(ColTypes); ++i) {
          result[i].push_back(current_indices[i]);
        }
      }
    };
  }

  /**
   * @brief Generates all possible combinations (Cartesian Product) of the input indices.
   * @param args Variadic list of RVec<int> columns (one per candidate).
   * @return RVec<RVec<int>> The combinations in SoA format.
   */
  template <typename... Args>
  Combis_t GenerateAllCombinations(const Args&... args) {
    // 1. Setup Output
    const std::size_t N_cols = sizeof...(Args);
    Combis_t result(N_cols);
    
    // 2. Prepare Recursion
    auto cols_tuple = std::make_tuple(args...);
    ROOT::RVec<int> current_indices(N_cols);
    
    impl::CartesianProduct<0, Args...>::Generate(cols_tuple, current_indices, result);
    
    return result;
  }

  // =================================================================================
  // 2. SYMMETRIC COMBINATIONS (GenerateSymmetricCombinations)
  // =================================================================================

  namespace impl {
    
    /** * @brief Checks if the current set of indices satisfies the symmetry constraints.
     * @details Ensures that for any defined symmetry group {A, B}, the values verify val(A) < val(B).
     */
    inline bool CheckSymmetry(const ROOT::RVec<int>& indices, 
                              const ROOT::RVec<ROOT::RVec<std::string>>& /*raw_groups_unused*/, 
                              // Note: Logic below assumes groups are passed as vector<vector<int>> indices
                              const ROOT::RVec<ROOT::RVec<int>>& groups) {
      
      for (const auto& group : groups) {
        // Group contains column indices (e.g. {0, 1} for gamma1, gamma2)
        // We strictly enforce indices[0] < indices[1] < ... to prevent swaps
        for (size_t k = 0; k < group.size() - 1; ++k) {
          int col_idx_a = group[k];
          int col_idx_b = group[k+1];
          if (indices[col_idx_a] >= indices[col_idx_b]) return false;
        }
      }
      return true;
    }

    /** @brief Recursive generator with symmetry check hook. */
    template <std::size_t I, typename... ColTypes>
    struct SymmetricProduct {
      static void Generate(const std::tuple<ColTypes...>& cols, 
                           const ROOT::RVec<ROOT::RVec<int>>& groups,
                           ROOT::RVec<int>& current_indices, 
                           Combis_t& result) {
        
        const auto& current_col = std::get<I>(cols);
        for (const auto& idx : current_col) {
          current_indices[I] = idx;
          SymmetricProduct<I + 1, ColTypes...>::Generate(cols, groups, current_indices, result);
        }
      }
    };

    template <typename... ColTypes>
    struct SymmetricProduct<sizeof...(ColTypes), ColTypes...> {
      static void Generate(const std::tuple<ColTypes...>&, 
                           const ROOT::RVec<ROOT::RVec<int>>& groups,
                           ROOT::RVec<int>& current_indices, 
                           Combis_t& result) {
        // Here is the "Symmetry Filter"
        if (CheckSymmetry(current_indices, {}, groups)) {
           for (std::size_t i = 0; i < sizeof...(ColTypes); ++i) {
             result[i].push_back(current_indices[i]);
           }
        }
      }
    };

    // --- Variadic Argument Peeler ---
    // This machinery separates the parameter pack (columns) from the last argument (groups).

    template <typename... T> struct SymHelper;

    // Recursive Step: Still have more than 2 args (Head + Tail...)
    // We keep building the tuple of columns.
    template <typename Head, typename... Tail>
    struct SymHelper<Head, Tail...> {
       static Combis_t Run(std::tuple<Head> head_tuple, const Tail&... tail) {
           return SymHelper<Tail...>::RecurseAndRun(head_tuple, tail...);
       }
       
       template <typename... PreviousCols>
       static Combis_t RecurseAndRun(std::tuple<PreviousCols...> prev, const Head& current, const Tail&... rest) {
           auto next_tuple = std::tuple_cat(prev, std::make_tuple(current));
           return SymHelper<Tail...>::RecurseAndRun(next_tuple, rest...);
       }
    };

    // Base Case: Only 1 type left in the pack -> This MUST be the Groups vector
    template <typename Last>
    struct SymHelper<Last> {
        template <typename... CollectedCols>
        static Combis_t RecurseAndRun(std::tuple<CollectedCols...> cols, const Last& groups) {
            // "Last" is expected to be ROOT::RVec<ROOT::RVec<int>> (or convertible)
            // "cols" is the tuple of all RVec columns.
            
            const std::size_t N_cols = sizeof...(CollectedCols);
            Combis_t result(N_cols);
            ROOT::RVec<int> current_indices(N_cols);
            
            SymmetricProduct<0, CollectedCols...>::Generate(cols, groups, current_indices, result);
            return result;
        }
    };

    // Helper to extract the first argument and kickstart the peeling
    // MOVED INSIDE IMPL NAMESPACE
    template <typename First, typename... Rest>
    Combis_t GenerateSymmetricCombinations_Launcher(const First& first, const Rest&... rest) {
        return SymHelper<Rest...>::RecurseAndRun(std::make_tuple(first), rest...);
    }

  } // namespace impl

  /**
   * @brief Generates combinations enforcing index ordering for symmetric groups.
   * * usage: GenerateSymmetricCombinations(col1, col2, col3, {{0,1}, {2,3}})
   * @tparam Args... The columns followed by one vector<vector<int>>.
   */
  template <typename... Args>
  Combis_t GenerateSymmetricCombinations(const Args&... args) {
     // Start the recursive peeling process
     // We assume at least 1 column + 1 group arg, so Args size >= 2
     static_assert(sizeof...(Args) >= 2, "GenerateSymmetricCombinations requires at least 1 column and 1 group vector.");
     
     // We need to unpack the first argument to start the tuple
     return impl::GenerateSymmetricCombinations_Launcher(args...);
  }

} // namespace combinatorics
} // namespace rad


// #pragma once

// #include <ROOT/RVec.hxx>
// #include <vector>
// #include <string>
// #include <map>
// #include <numeric>
// #include <algorithm>
// #include <stdexcept>
// #include <unordered_set> 
// #include "DefineNames.h" 
// #include "Constants.h" 
// #include "CommonDefines.h"

// /**
//  * @file Combinatorics.h
//  * @brief Utilities for generating particle combinations.
//  * * @details
//  * This file contains the logic to generate the "Cartesian Product" of all particle candidates.
//  * It converts a list of candidates (e.g., 2 electrons, 1 positron) into a set of 
//  * unique combinatorial events.
//  */

// namespace rad {
//   namespace combinatorics {

//     using ROOT::RVecI;
//     using ROOT::RVec;
//     using ROOT::RVec;
//     using std::string;
//     using rad::consts::InvalidEntry;
    
//     //---------------------------------------------------------
//     // Core Combinatorial Logic Helper
//     //---------------------------------------------------------
    
//     /**
//      * @brief Generates all valid, unique combinations of candidates.
//      * * @details
//      * This algorithm performs two main tasks:
//      * 1. **Cartesian Product:** It iterates through every possible permutation of the input candidates.
//      * (e.g., if Ele has 2 candidates and Pos has 3, it tests 2*3=6 combinations).
//      * 2. **Uniqueness Filter:** It discards any combination where the same underlying detector object 
//      * (index) is used for multiple roles (e.g., the same track used as both an electron and a pion).
//      * * **Output Format:**
//      * The output is a "Structure of Arrays" (SoA).
//      * `result[particle_role_index]` is a vector containing the candidate index for that role 
//      * for every valid combination.
//      * * @param candidates_vec A vector of candidate lists. `candidates_vec[particle_role]` contains the 
//      * list of valid indices for that role (e.g., {0, 2, 5} for electrons).
//      * @return RVecIndices The structure of valid combinations.
//      */
//     inline RVecIndices GenerateAllCombinations(const RVecIndices& candidates_vec) {
        
//       const size_t n_particles = candidates_vec.size();
//       if (n_particles == 0) {
//         return RVecIndices(n_particles);
//       }

//       // --- 1. Initial Setup and Pre-calculation ---
//       // Determine the bounds for the N-dimensional counter
//       ROOT::RVec<size_t> max_indices;
//       max_indices.reserve(n_particles);
        
//       size_t total_combos = 1;
//       for (const auto& vec : candidates_vec) {
//         if (vec.empty()) return RVecIndices(n_particles); // If any required particle has 0 candidates, 0 combos exist.
//         max_indices.push_back(vec.size());
//         total_combos *= vec.size();
//       }

//       // The N-dimensional counter (state of the current permutation)
//       ROOT::RVec<size_t> current_indices(n_particles, 0);

//       // --- 2. Initialize Transposed Output Structure ---
//       RVecIndices result_by_particle(n_particles);
//       // We cannot reserve accurately due to skipped combinations (overlaps), 
//       // but reserving 'total_combos' avoids reallocations in the worst case.
//       for (size_t i = 0; i < n_particles; ++i) {
//         result_by_particle[i].reserve(total_combos); 
//       }

//       // --- 3. Iterative Combinatorial Loop ---
      
//       // We reuse these buffers to avoid malloc/free overhead for every combination.
//       ROOT::RVec<int> current_combination_buffer(n_particles);
//       std::unordered_set<int> unique_checker;
//       unique_checker.reserve(n_particles * 2); 

//       for (size_t i_combo = 0; i_combo < total_combos; ++i_combo) {
            
//         unique_checker.clear();
//         bool is_unique = true;

//         // --- A. Build Combination and Check Uniqueness ---
//         for (size_t i_particle = 0; i_particle < n_particles; ++i_particle) {
                
//           // Lookup the actual Data Index for this role in the current permutation
//           int particle_index = candidates_vec[i_particle][current_indices[i_particle]];
                
//           // Try to insert the index. If insert() returns {iterator, false}, the index is a duplicate.
//           if (!unique_checker.insert(particle_index).second) {
//             is_unique = false;
//             break; // Stop checking this combination immediately
//           }
//           current_combination_buffer[i_particle] = particle_index;
//         }

//         // --- B. Record Valid Combination ---
//         if (is_unique) {
//           // Transpose: Append this combo's indices to the columnar output
//           for (size_t i_particle = 0; i_particle < n_particles; ++i_particle) {
//             result_by_particle[i_particle].push_back(current_combination_buffer[i_particle]);
//           }
//         }
            
//         // --- C. Increment the N-dimensional counter ---
//         // Simulates a nested loop of depth 'n_particles'
//         size_t k = 0;
//         while (k < n_particles) {
//           current_indices[k]++;
//           if (current_indices[k] < max_indices[k]) {
//             break; // Carry propagation done
//           }
//           current_indices[k] = 0; // Reset this digit and carry over
//           k++;
//         }
//       }
//       return result_by_particle;
//     }
 
//   } // namespace combinatorics
// } // namespace rad
