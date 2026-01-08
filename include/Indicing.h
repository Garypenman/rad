#pragma once

#include "RVecHelpers.h"
#include "Constants.h"
#include <ROOT/RVec.hxx>

/**
 * @file Indicing.h
 * @brief Utilities for finding and filtering particle indices in RVecs.
 * * @details
 * This file provides "Function Factories" (Lambda Generators).
 * These functions do not return indices directly; instead, they return **Lambda Functions**
 * that are registered with `ParticleCreator` and executed later inside the RDataFrame event loop.
 * * **Why Factories?**
 * Because the particle search criteria (e.g., "Sort electrons by Energy and take the 2nd one")
 * must be defined at configuration time (C++ setup) but executed at runtime (on the data arrays).
 */

namespace rad {
  namespace index {

    // =========================================================================
    //  Particle Candidate Selection (Function Factories)
    // =========================================================================

    /**
     * @brief Creates a search strategy to find the N-th occurrence of a specific value.
     * * **Usage:** `useNthOccurance(2, 11)` -> Returns a lambda that finds the index of the 2nd electron (PDG 11).
     * * @param n_occurance Which instance to find (0-based: 0 = 1st, 1 = 2nd).
     * @param value The value (e.g., PDG code) to search for in the target column.
     * @return A lambda with signature: `(const RVecI& values) -> int`
     */ 
    inline auto useNthOccurance(const int n_occurance, const int value) {
      return [n_occurance, value](const ROOT::RVecI& values) {
	return util::findNthIndex(values, n_occurance, value);
      };
    }

    /**
     * @brief Creates a search strategy to find the N-th occurrence of a value, 
     * but prioritizes particles based on a sorting criteria (e.g., Energy, pT).
     * * **Usage:** `useNthOccuranceSortedBy<double>(0, 11)` 
     * -> Finds the highest-ranked electron according to the provided sorter vector.
     * * **Sorting Logic:** The default behavior of `StableArgsort` is usually ascending.
     * If you need descending order (e.g., highest energy first), the user must pass the 
     * sorter vector as negative (e.g., `-energy`).
     * * @tparam T Type of the sorting vector (float, double, int).
     * @param n_occurance Which instance to find (0-based) after sorting.
     * @param value The value (e.g. PDG code) to filter by.
     * @return A lambda with signature: `(const RVecI& values, const RVec<T>& sorter) -> int`
     */
    template <typename T>
    inline auto useNthOccuranceSortedBy(const int n_occurance, const int value) {
      return [n_occurance, value](const ROOT::RVecI& values, const ROOT::RVec<T>& sorter) {
	// 1. Get indices that would sort the 'sorter' vector
	auto sortIndices = ROOT::VecOps::StableArgsort(sorter);
            
	// 2. Reorder the 'values' (e.g. PDG codes) using that sort order
	auto reorderedValues = ROOT::VecOps::Take(values, sortIndices);
            
	// 3. Find the Nth occurrence in this new sorted list
	// Note: 'k' is the index in the REORDERED array.
	int k = util::findNthIndex(reorderedValues, n_occurance, value);
            
	if(k == -1) return rad::consts::InvalidIndex();
            
	// 4. Return the ORIGINAL index
	// We must map 'k' back to the original event index using sortIndices.
	return (int)sortIndices[k];
      };
    }

    /**
     * @brief Creates a simple strategy to treat a specific entry in an ID vector as the particle index.
     * * **Use Case:** When an upstream algorithm (e.g., a trigger or external finder) has already 
     * determined which particle index to use, and stored it in a branch.
     * * @param entry The position in the ID vector to read (usually 0 if the vector holds just the best candidate).
     * @param offset Optional offset to subtract (e.g., to account for removed beam particles or array shifts).
     * @return A lambda with signature: `(const RVecI& ids) -> int`
     */
    inline auto UseAsID(size_t entry, int offset = 0) {
      return [entry, offset](const ROOT::RVecI& id) -> int {
	if (entry >= id.size()) return rad::consts::InvalidIndex();
	return id[entry] - offset;
      };
    }

    /**
     * @brief Creates a strategy to filter candidate indices by value.
     * * **Usage:** `FilterIndices(11)` -> Returns list of indices where `vec[i] == 11`.
     * This is useful for creating a list of **Candidates** (e.g. "All Electrons") rather than a single particle.
     * * @tparam T Type of value to compare (int, double, etc).
     * @param val The value to match.
     * @return A lambda with signature: `(const RVec<T>& vec) -> RVecI`
     */
    template <typename T>
    inline auto FilterIndices(const T val) {
      return [val](const ROOT::RVec<T>& vec) -> ROOT::RVecI {
	// Nonzero returns the indices where the condition is true (vec[i] == val)
	return ROOT::VecOps::Nonzero(vec == val);
      };
    }


    // =========================================================================
    //  Index Utilities (Helper Functions)
    // =========================================================================

    /**
     * @brief Checks if a vector of indices contains any invalid values (-1).
     * * @tparam T Vector type.
     * @param indices Vector of indices to check.
     * @return `true` if any index in the vector is `InvalidIndex` (-1).
     */
    template<typename T>
    inline bool InvalidIndices(const ROOT::RVec<T>& indices) {
      // Efficient vectorized check using ROOT::VecOps::Any
      return ROOT::VecOps::Any(indices == rad::consts::InvalidIndex());
    }

    /**
     * @brief Validates a single index against a whitelist of allowed indices.
     * If 'index' is NOT found in the 'valid' list, it is reset to `InvalidIndex` (-1).
     * * @tparam T Type of the valid list (e.g., RVecI).
     * @tparam Ti Type of the index to check (e.g., int).
     * @param valid List of acceptable/valid indices.
     * @param index Reference to the index to check. Will be modified in-place if invalid.
     */
    template<typename T, typename Ti>
    inline void InvalidateIndices(const ROOT::RVec<T>& valid, Ti& index) {
      if (ROOT::VecOps::Any(valid == index)) {
	return; // Found, it's valid
      }
      index = rad::consts::InvalidIndex();
    }

  } // namespace index
} // namespace rad
