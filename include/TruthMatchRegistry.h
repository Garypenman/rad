/**
 * @file TruthMatchRegistry.h
 * @brief Registry for assigning Truth Roles to reconstructed particles.
 */

#pragma once
#include <string>
#include <map>
#include <vector>

namespace rad {

  /**
   * @class TruthMatchRegistry
   * @brief Manages the storage of Truth Roles (MCParticles indices) for analysis candidates.
   */
  class TruthMatchRegistry {
  public:
    /**
     * @brief Registers a Truth Role for a flat candidate.
     * @param name The name of the candidate (e.g. "scat_ele").
     * @param truthIdx The index of the corresponding particle in the MCParticles collection.
     */
    void AddParticleMatch(const std::string& name, int truthIdx) {
      _particleMatches[name] = truthIdx;
    }

    /**
     * @brief Registers a Truth Role for a topological parent.
     * @param parentName The name of the parent particle (e.g. "pi0").
     * @param truthIdx The index of the parent in the MCParticles collection.
     * * The validation logic checks if the set of children matches the truth children.
     */
    void AddParentMatch(const std::string& parentName, int truthIdx) {
      _parentMatches[parentName] = truthIdx;
    }

    /**
     * @brief Retrieves all registered particle matches.
     * @return A map of candidate names to truth indices.
     */
    const std::map<std::string, int>& GetParticleMatches() const { return _particleMatches; }

    /**
     * @brief Retrieves all registered parent matches.
     * @return A map of parent names to truth indices.
     */
    const std::map<std::string, int>& GetParentMatches() const { return _parentMatches; }

  private:
    std::map<std::string, int> _particleMatches; ///< Map for flat candidates
    std::map<std::string, int> _parentMatches;   ///< Map for topological parents
  };

} // namespace rad
