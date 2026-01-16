/**
 * @file ParticleCreator.h
 * @brief Manages the definition, indexing, and creation of intermediate particles.
 */

#pragma once

#include "ConfigReaction.h"
#include "StringUtilities.h"
#include "ParticleCreatorMethods.h"
#include "CommonDefines.h"
#include <map>
#include <set>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm> 
#include <iostream>

namespace rad {

    using StructuredNames_t = std::vector<ParticleNames_t>;
    using IndexMap_t = std::unordered_map<std::string, int>;
    
    using ParticleCreatorFunc_t = void (*)(
        const Indice_t position, 
        const RVecIndices&, 
        ROOT::RVecD&, ROOT::RVecD&, ROOT::RVecD&, ROOT::RVecD&
    );
    
    /**
     * @class ParticleCreator
     * @brief Manages the creation of intermediate particles and the resolution of reaction topology.
     * @details 
     * This class acts as the topology builder. It handles:
     * - Registration of input particles from the analysis tree.
     * - Definition of composite particles via Sum/Diff algorithms.
     * - Mapping of string names to efficient integer indices.
     * - Management of particle groups (e.g., Mesons, Baryons) ensuring they update 
     * correctly when composite particles are added.
     */
    class ParticleCreator {

    public:
      
      /** @brief Default constructor. */
      ParticleCreator() = default;

      /** * @brief Constructor for a new topology configuration.
       * @param cr Pointer to the ConfigReaction object managing RDF nodes.
       * @param prefix Prefix for branch names (e.g., "rec_", "mc_").
       * @param suffix Suffix for branch names (e.g., "", "_sysUp").
       */
      explicit ParticleCreator(ConfigReaction* cr, const std::string& prefix, const std::string& suffix = "");

    
      /** * @brief Copy Constructor for cloning topologies (e.g. Rec -> Truth).
       * @param other The source ParticleCreator to copy configuration from.
       * @param new_suffix The new suffix to apply to the cloned topology.
       */
      ParticleCreator(const ParticleCreator& other, const std::string& new_suffix);
      
      /** @brief Sets the ConfigReaction pointer. */
      void SetReaction(ConfigReaction* reaction);
      
      /** @return Pointer to the associated ConfigReaction. */
      ConfigReaction* Reaction() const;
      
      /** @return The prefix string used for this topology. */
      std::string GetPrefix() const;
      
      /** @return The suffix string used for this topology. */
      std::string GetSuffix() const;

      // =======================================================================
      // Configuration Interface
      // =======================================================================

      /** * @brief Associates a standard group ID (e.g., Baryons) with a group name.
       * @param order The integer ID of the group (from rad::consts).
       * @param groupName The string name of the group (e.g., "Baryons").
       */
      void DefineGroup(int order, const std::string& groupName);

      /** * @brief Overrides a standard group with an explicit list of particles.
       * @param order The integer ID of the group.
       * @param particles Vector of particle names to include in this group.
       */
      void OverrideGroup(int order, const std::vector<std::string>& particles);
      
      /** * @brief Helper to override a group by its abstract name.
       * @param abstractName The abstract name (e.g., "Baryons", "Mesons").
       * @param particles Vector of particle names to include.
       */
      void OverrideGroup(const std::string& abstractName, const std::vector<std::string>& particles);

      /** * @brief Sets the names of the beam particles.
       * @param beams Vector of beam particle names (usually {Ion, Electron}).
       */
      void SetBeamNames(const std::vector<std::string>& beams);
      
      /** * @brief Forces an input particle to be registered in the map.
       * @details Essential for Master processors when a Linked Clone needs a particle 
       * (e.g., "recoil") that the Master itself does not explicitly use.
       * @param name The name of the input particle to require.
       */
      void RequireParticle(const std::string& name);

      /** * @brief Copies configuration (groups, beams, requirements) from another creator.
       * @param other The ParticleCreator to copy from.
       */
      void CopyConfigurationFrom(const ParticleCreator& other);

      /** * @brief Ports group definitions (like ScatEle) from a parent Creator to this one.
       * @details Crucial for Clones: reads parent groups (e.g. "rec_scat_ele_group"), 
       * swaps the prefix (to "mc_"), and registers the new group.
       */
      void PortGroupsFrom(const ParticleCreator& parent);

      // =======================================================================
      // Particle Definition Logic
      // =======================================================================

      /** * @brief Registers a generic particle with a custom creation function.
       * @param name The unique name of the new particle.
       * @param func Function pointer to the creation logic.
       * @param depends List of dependencies required by the creation function.
       */
      void AddParticle(const std::string& name, ParticleCreatorFunc_t func, const StructuredNames_t& depends={{}});

      /** * @brief Defines a particle as the 4-vector Sum of others.
       * @param name The name of the composite particle (e.g., "Jpsi").
       * @param depends List of daughter particles (e.g., {{"ele", "pos"}}).
       */
      void Sum(const std::string& name, const StructuredNames_t& depends={{}});

      /** * @brief Defines a particle as Sum and registers it for Truth Matching.
       * @param name The name of the composite particle.
       * @param truthRole The integer Role ID for truth matching comparison.
       * @param depends List of daughter particles.
       */
      void SumTruthMatch(const std::string& name, int truthRole, const StructuredNames_t& depends={{}});

      /** * @brief Defines a particle via subtraction (Missing Mass).
       * @details P_miss = Sum(Initial) - Sum(Final).
       * @param name The name of the missing particle (e.g., "n_miss").
       * @param depends List containing {{Initial State}, {Final State}}.
       */
      void Diff(const std::string& name, const StructuredNames_t& depends={{}});

      // =======================================================================
      // Indexing & Initialization 
      // =======================================================================

      /** * @brief Resolves indices, registers columns, and builds the reaction map.
       * @details Includes logic for idempotent registration (safe for shared configs)
       * and handles the definition of helper columns and groups.
       */
      void InitMap();

      /** * @brief Checks if a particle exists in the current map.
       * @param name The particle name to check.
       * @return True if the particle is registered.
       */
      bool HasParticle(const std::string& name) const;

      /** @return Reference to the internal name-to-index map. */
      const IndexMap_t& GetIndexMap() const;

      /** @brief Resolves internal dependency indices for all created particles. */
      void ResolveDependencies();
      
      /** * @brief Adopts indices and configuration from a Master processor.
       * @details Used by Linked Processors to ensure index consistency.
       * @param master The master ParticleCreator to copy data from.
       */
      void AdoptIndices(const ParticleCreator& master);
      
      /** * @brief Rebuilds the map for a Linked Processor.
       * @details Ensures groups and helper columns are defined correctly for the 
       * new suffix, dealing with any overrides specific to the clone.
       */
      void RebuildReactionMap();

      // =======================================================================
      // Execution
      // =======================================================================

      /** * @brief Executes the particle creation functions to generate 4-vectors.
       * @param px Reference to Px vector to write to.
       * @param py Reference to Py vector to write to.
       * @param pz Reference to Pz vector to write to.
       * @param m Reference to Mass vector to write to.
       */
      void ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, 
                         ROOT::RVecD& pz, ROOT::RVecD& m) const;

      // =======================================================================
      // Accessors
      // =======================================================================

      /** @return The number of created (non-input) particles. */
      size_t GetNCreated() const;
      
      /** * @brief Gets the reaction index for a named particle.
       * @param name The particle name.
       * @return The integer index. Throws runtime_error if not found.
       */
      Int_t GetReactionIndex(const std::string& name) const;
      
      /** * @brief Translates an input index to a reaction index.
       * @param input The raw input index.
       * @return The reaction index.
       */
      Int_t GetReactionIndex(size_t input) const;
      
      /** @return The name of the RDataFrame column containing the map indices. */
      std::string GetMapName() const;

      /** @return A list of input columns required by this creator. */
      std::vector<std::string> GetPriorDependencies();
      
      /** * @brief Helper to assign indices to a list of names.
       * @param nameIndex The map to populate.
       * @param names The list of names to index.
       * @param idx Reference to the current index counter.
       * @return A vector of assigned indices.
       */
      Indices_t CreateIndices(IndexMap_t& nameIndex, const ParticleNames_t& names, size_t& idx);

      /** * @brief Returns a list of all particle names currently registered in the index map.
       * @details Used by KinematicsProcessor to iterate over particles and define output columns.
       */
      ROOT::RVec<std::string> GetParticleNames() const {
	ROOT::RVec<std::string> names;
          names.reserve(_nameIndex.size());
          for(const auto& pair : _nameIndex) {
              names.push_back(pair.first);
          }
          return names;
      }
      
    private:
      ConfigReaction* _reaction = nullptr;
      std::string _prefix; 
      std::string _suffix; 
      
      std::vector<std::string> _beam_names;
      std::vector<std::string> _forced_inputs; 

      std::map<int, std::string> _config_groups; 
      std::map<int, std::vector<std::string>> _explicit_groups; 

      ParticleNames_t _p_names;
      std::vector<ParticleCreatorFunc_t> _p_creators;
      std::vector<ParticleNames_t> _p_required;
      std::vector<StructuredNames_t> _p_stru_depends;
      ROOT::RVec<RVecIndices> _p_dep_indices; 
      std::set<std::string> _dependencies;
      
      ParticleNames_t _inputNames;
      RVecIndexMap _mapIndices;
      IndexMap_t _nameIndex;        
      IndexMap_t _nameInputIndex;   
      IndexMap_t _createIndex;
      Indices_t _input2ReactionIndex; 

      /** @brief Helper to redefine group columns safely handling RVec types. */
      void RedefineGroupColumn(const std::string& name, const std::vector<std::string>& cols);
      
      /** @brief Safety wrapper for index lookup to prevent confusing out_of_range errors. */
      int GetIndexSafe(const std::string& name) const;
    };

    // ========================================================================
    // Implementation
    // ========================================================================

    inline ParticleCreator::ParticleCreator(ConfigReaction* cr, const std::string& prefix, const std::string& suffix) 
        : _reaction{cr}, _prefix{prefix}, _suffix{suffix} {
          _beam_names = {consts::BeamIon(), consts::BeamEle()};
          DefineGroup(consts::OrderBaryons(), consts::Baryons());
          DefineGroup(consts::OrderMesons(), consts::Mesons());
    }

    inline ParticleCreator::ParticleCreator(const ParticleCreator& other, const std::string& new_suffix) 
        : _reaction(other._reaction), _prefix(other._prefix), _suffix(new_suffix),
          _beam_names(other._beam_names), _forced_inputs(other._forced_inputs),
          _config_groups(other._config_groups), _explicit_groups(other._explicit_groups),
          _p_names(other._p_names), _p_creators(other._p_creators),
          _p_required(other._p_required), _p_stru_depends(other._p_stru_depends),
          _dependencies(other._dependencies), _inputNames(other._inputNames),
          _createIndex(other._createIndex)
    {}

    inline void ParticleCreator::SetReaction(ConfigReaction* reaction) { _reaction = reaction; }
    inline ConfigReaction* ParticleCreator::Reaction() const { return _reaction; }
    inline std::string ParticleCreator::GetPrefix() const { return _prefix; }
    inline std::string ParticleCreator::GetSuffix() const { return _suffix; }

    inline void ParticleCreator::DefineGroup(int order, const std::string& groupName) {
          _config_groups[order] = groupName;
          _explicit_groups.erase(order); 
    }
    inline void ParticleCreator::OverrideGroup(int order, const std::vector<std::string>& particles) {
          _explicit_groups[order] = particles;
          _config_groups.erase(order);    
    }
    inline void ParticleCreator::OverrideGroup(const std::string& abstractName, const std::vector<std::string>& particles) {
          if(abstractName == consts::Baryons()) OverrideGroup(consts::OrderBaryons(), particles);
          else if(abstractName == consts::Mesons()) OverrideGroup(consts::OrderMesons(), particles);
    }
    inline void ParticleCreator::SetBeamNames(const std::vector<std::string>& beams) { _beam_names = beams; }
    
    inline void ParticleCreator::RequireParticle(const std::string& name) {
          if(std::find(_forced_inputs.begin(), _forced_inputs.end(), name) == _forced_inputs.end())
              _forced_inputs.push_back(name);
    }

    inline void ParticleCreator::CopyConfigurationFrom(const ParticleCreator& other) {
          _config_groups = other._config_groups; _explicit_groups = other._explicit_groups;
          _beam_names = other._beam_names; _forced_inputs = other._forced_inputs;
          _p_names = other._p_names; _p_creators = other._p_creators;
          _p_required = other._p_required; _p_stru_depends = other._p_stru_depends;
          _dependencies = other._dependencies; _createIndex = other._createIndex;
    }

    inline void ParticleCreator::PortGroupsFrom(const ParticleCreator& parent) {
        if (parent.GetPrefix() == _prefix) return; 

        // Iterate over the parent's abstract groups (e.g. ScatEle, Baryons)
        for(const auto& [order, abstractName] : parent._config_groups) {
            // Reconstruct the Parent's specific group name (e.g. "rec_scat_ele_group")
            std::string parentGroupFull = parent.GetPrefix() + abstractName + parent.GetSuffix();
            
            try {
                // Fetch the list of particles from the Parent's RDF definition
                auto particles = _reaction->GetGroup(parentGroupFull);
                std::vector<std::string> newParticles;
                
                // For each particle, strip the parent prefix and prep for the new prefix
                // e.g. "rec_scat_ele" -> "scat_ele" (which will later become "mc_scat_ele")
                for(const auto& p : particles) {
                    std::string base = p;
                    if(p.find(parent.GetPrefix()) == 0) base = p.substr(parent.GetPrefix().length());
                    newParticles.push_back(base);
                }
                
                // Register the New Group under the current processor's prefix
                // SetGroupParticles automatically prepends _prefix to the column names
                _reaction->SetGroupParticles(abstractName + _suffix, _prefix, newParticles);
            } catch (...) {
                // If the parent group isn't defined, we simply skip it.
            }
        }
    }

    inline void ParticleCreator::AddParticle(const std::string& name, ParticleCreatorFunc_t func, const StructuredNames_t& depends) {
        _p_names.push_back(name);
        auto flat_depends = util::flattenColumnNames(depends);
        _p_required.push_back(flat_depends);
        _dependencies.insert(flat_depends.begin(), flat_depends.end());
        _p_stru_depends.push_back(depends);
        _p_creators.push_back(func);
        _createIndex[name] = GetNCreated() - 1;
    }
    inline void ParticleCreator::Sum(const std::string& name, const StructuredNames_t& depends) { AddParticle(name, ParticleCreateBySum, depends); }
    inline void ParticleCreator::SumTruthMatch(const std::string& name, int truthRole, const StructuredNames_t& depends) {
        Sum(name, depends);
        if(_reaction) _reaction->GetTruthMatchRegistry().AddParentMatch(name, truthRole);
    }
    inline void ParticleCreator::Diff(const std::string& name, const StructuredNames_t& depends) { AddParticle(name, ParticleCreateByDiff, depends); }

    inline std::vector<std::string> ParticleCreator::GetPriorDependencies() {
        std::vector<std::string> vec_deps(_dependencies.begin(), _dependencies.end());
        util::removeExistingStrings(vec_deps, _p_names);
        return vec_deps;
    }

    inline int ParticleCreator::GetIndexSafe(const std::string& name) const {
        auto it = _nameIndex.find(name);
        if (it == _nameIndex.end()) {
            std::cerr << "\n[ParticleCreator FATAL] Missing Particle Index: " << name << "\n";
            std::cerr << "  Processor: " << _prefix << " (Suffix: '" << _suffix << "')\n";
            throw std::runtime_error("Particle '" + name + "' missing in map. Check inputs.");
        }
        return it->second;
    }

    inline Indices_t ParticleCreator::CreateIndices(IndexMap_t& nameIndex, const ParticleNames_t& names, size_t& idx) {
        Indices_t indices{};
        for(const auto& name : names) {
          if(nameIndex.count(name) == 0) { indices.push_back(idx); nameIndex[name] = idx; ++idx; } 
          else { indices.push_back(nameIndex[name]); }
        }
        return indices;
    }

    inline void ParticleCreator::RedefineGroupColumn(const std::string& name, const std::vector<std::string>& cols) {
        size_t n = cols.size();
        using RVecI = ROOT::VecOps::RVec<int>;
        if(n == 1) _reaction->Redefine(name, [](const RVecI& a){ return a; }, cols);
        else if(n == 2) _reaction->Redefine(name, [](const RVecI& a, const RVecI& b){ RVecI v=a; v.insert(v.end(),b.begin(),b.end()); return v; }, cols);
        else if(n == 3) _reaction->Redefine(name, [](const RVecI& a, const RVecI& b, const RVecI& c){ RVecI v=a; v.insert(v.end(),b.begin(),b.end()); v.insert(v.end(),c.begin(),c.end()); return v; }, cols);
        else if(n == 4) _reaction->Redefine(name, [](const RVecI& a, const RVecI& b, const RVecI& c, const RVecI& d){ RVecI v=a; v.insert(v.end(),b.begin(),b.end()); v.insert(v.end(),c.begin(),c.end()); v.insert(v.end(),d.begin(),d.end()); return v; }, cols);
        else if(n == 5) _reaction->Redefine(name, [](const RVecI& a, const RVecI& b, const RVecI& c, const RVecI& d, const RVecI& e){ RVecI v=a; v.insert(v.end(),b.begin(),b.end()); v.insert(v.end(),c.begin(),c.end()); v.insert(v.end(),d.begin(),d.end()); v.insert(v.end(),e.begin(),e.end()); return v; }, cols);
        else _reaction->RedefineExpr(name, rad::util::createPackVectorString(cols));
    }

    /**
     * @brief Builds the Reaction Map, defining particles and resolving indices.
     * @details This function performs the critical setup:
     * 1. Gathers all Input Particles (Beams, Manually Forced inputs, and Configured Groups).
     * 2. Checks for namespace collisions to prevent overwriting existing processors.
     * 3. Assigns integer indices to all particles (String -> Int map).
     * 4. Defines RDF helper columns for single particles (e.g. converting int -> RVec<int>).
     * 5. Refreshes/Redefines Group columns to ensure they point to the correct suffix-specific columns.
     */
    inline void ParticleCreator::InitMap() {
      // 1. GATHER INPUT PARTICLES
      // Start with beams and forced inputs (manually required by user)
      ParticleNames_t logicalInputNames;
      logicalInputNames.insert(logicalInputNames.end(), _beam_names.begin(), _beam_names.end());
      logicalInputNames.insert(logicalInputNames.end(), _forced_inputs.begin(), _forced_inputs.end());

      // Collect particles from defined Groups (e.g. "Baryons", "ScatEle")
      // We check for both explicitly overridden groups and abstract config groups.
      std::vector<int> orders = {consts::OrderBaryons(), consts::OrderMesons(), consts::OrderScatEle()};
      for(int order : orders) {
          if(_explicit_groups.count(order)) {
              const auto& list = _explicit_groups[order];
              logicalInputNames.insert(logicalInputNames.end(), list.begin(), list.end());
          } 
          else if (_config_groups.count(order)) {
              try {
                 std::string uniqueGroupName = _config_groups[order] + _suffix;
                 std::vector<std::string> particles;
                 // Try getting group with full typed name (e.g. "rec_Baryons_miss")
                 try { particles = _reaction->GetGroup(_prefix + uniqueGroupName); } 
                 catch(...) { particles = _reaction->GetGroup(_prefix + _config_groups[order]); }
                 
                 // Normalize names by stripping prefix for internal storage
                 for(const auto& p : particles) {
                     if(p.find(_prefix) == 0) logicalInputNames.push_back(p.substr(_prefix.length()));
                     else logicalInputNames.push_back(p); 
                 }
              } catch (...) { }
          }
      }

      // Add dependencies required by created particles (e.g. "Jpsi" depends on "ele" and "pos")
      auto dep_names = GetPriorDependencies();
      util::removeExistingStrings(dep_names, logicalInputNames);
      
      _inputNames = logicalInputNames;
      _inputNames.insert(_inputNames.end(), dep_names.begin(), dep_names.end());
      util::removeExistingStrings(_inputNames, _p_names); // Ensure inputs are not created particles

      // 2. COLLISION CHECK & MAP REGISTRATION
      // Ensure we aren't defining a processor that already exists (same prefix + suffix)
      std::string kineCol = consts::KineIndices() + _suffix;
      if (_reaction->ColumnExists(_prefix + kineCol)) {
          throw std::runtime_error("Topology Collision: Processor [" + _prefix + "] suffix [" + _suffix + "] already exists.");
      }
      _reaction->SetGroupParticles(kineCol, _prefix, _inputNames);

      // 3. BUILD INDEX MAPS
      // Map string names to integer indices [0, N] for efficient internal processing
      size_t in_idx = 0;
      CreateIndices(_nameInputIndex, _inputNames, in_idx);
        
      size_t idx = 0;
      Indices_t idxBeam = CreateIndices(_nameIndex, _beam_names, idx);
      
      // Helper to resolve indices for a specific group order
      auto get_indices_for_group = [&](int order) {
          if (_explicit_groups.count(order)) return CreateIndices(_nameIndex, _explicit_groups[order], idx);
          if (_config_groups.count(order)) {
              try {
                  std::string uniqueName = _config_groups[order] + _suffix;
                  std::vector<std::string> particles;
                  try { particles = _reaction->GetGroup(_prefix + uniqueName); } 
                  catch(...) { particles = _reaction->GetGroup(_prefix + _config_groups[order]); }
                  
                  ParticleNames_t logicalNames;
                  for(auto p : particles) {
                      if(p.find(_prefix) == 0) logicalNames.push_back(p.substr(_prefix.length()));
                      else logicalNames.push_back(p);
                  }
                  return CreateIndices(_nameIndex, logicalNames, idx);
              } catch (...) {}
          }
          return Indices_t{};
      };

      Indices_t idxBaryons  = get_indices_for_group(consts::OrderBaryons());
      Indices_t idxMesons   = get_indices_for_group(consts::OrderMesons());
      Indices_t idxScat_ele = get_indices_for_group(consts::OrderScatEle());
      Indices_t idxDeps     = CreateIndices(_nameIndex, dep_names, idx);
      Indices_t idxCreate   = CreateIndices(_nameIndex, _p_names, idx); // Created particles get highest indices

      // Register the complete index map into RDF for use by kernels
      RVecIndexMap retIndices{idxBeam, idxBaryons, idxMesons, idxScat_ele, idxDeps, idxCreate};
      _mapIndices = retIndices;
      _reaction->Define(GetMapName(), [retIndices]() { return retIndices; }, {});

      // DEFINE INPUT ALIASES 
      // For a Linked Clone (suffix="_miss"), we alias "tru_scat_ele" -> "tru_scat_ele_miss"
      // This allows Loop 5 to uniformly add suffixes to everything without breaking.
      for(const auto& name : _inputNames) {
           std::string colName = _prefix + name + _suffix;
           std::string masterCol = _prefix + name;
           
           // Only define if we have a suffix, the alias doesn't exist, but the master does
           if(!_suffix.empty() && !_reaction->ColumnExists(colName) && _reaction->ColumnExists(masterCol)) {
               _reaction->Define(colName, masterCol); 
           }
      }
      
      // 4. DEFINE HELPER COLUMNS
      // For every single particle, define an RDF column returning RVec<int> (for uniformity with groups)
      for(const auto& particle : _nameIndex) {
        bool isInput = std::find(_inputNames.begin(), _inputNames.end(), particle.first) != _inputNames.end();

        // Register particle index if new
        if (!isInput) {
            try { _reaction->SetParticleIndex(particle.first, _prefix, particle.second); } 
            catch (const std::invalid_argument&) {}
        }

        if (isInput) continue; // Inputs are already defined by SetParticleCandidates

        // Define column "prefix_name_suffix" -> returns [index]
        std::string colName = _prefix  + particle.first + _suffix;
        auto func = [idx = particle.second](){ return ROOT::VecOps::RVec<int>{idx}; };
        
        if(_reaction->ColumnExists(colName)) _reaction->Redefine(colName, func, {});
        else _reaction->Define(colName, func, {});
      }

      // 5. REFRESH GROUPS (CRITICAL for Clones)
      // Ensure group columns (e.g. "rec_Baryons_miss") point to the correct suffix-specific particle columns.
      for(auto const& [order, groupName] : _config_groups) {
           std::string uniqueGroupName = groupName + _suffix;
           std::string typedUniqueName = _prefix + uniqueGroupName;
           
           std::vector<std::string> abstractParticles;
           bool groupFound = false;

           // A. Try fetching exact group from RDF (e.g. "tru_scat_ele_group")
           try { 
               abstractParticles = _reaction->GetGroup(typedUniqueName); 
               groupFound = true;
           } 
           catch(...) { 
               // B. Try fetching base group (e.g. "scat_ele_group") - mostly for Rec
               try {
                   abstractParticles = _reaction->GetGroup(_prefix + groupName); 
                   groupFound = true;
               } catch (...) {}
           }

           // C. Fallback: Auto-define from Particle Name (Handles Truth ScatEle case)
           // If "scat_ele_group" is missing but "scat_ele" particle exists, make a group of 1.
           if (!groupFound) {
               std::string candidateName = groupName;
               if (HasParticle(candidateName)) {
                   abstractParticles = {candidateName};
               } else {
                   // Try stripping "_group" suffix
                   std::string suffixStr = "_group";
                   if (candidateName.length() > suffixStr.length() && 
                       candidateName.compare(candidateName.length() - suffixStr.length(), suffixStr.length(), suffixStr) == 0) {
                       candidateName = candidateName.substr(0, candidateName.length() - suffixStr.length());
                       if (HasParticle(candidateName)) abstractParticles = {candidateName};
                   }
               }

               if (!abstractParticles.empty()) {
                   // Auto-register this found particle as the group content
                   _reaction->SetGroupParticles(uniqueGroupName, _prefix, abstractParticles);
               } else {
                   std::cerr << "\n[ParticleCreator ERROR] Missing Group: " << _prefix + groupName << "\n";
                   throw std::runtime_error("Group definition failed. Ensure particle candidates are defined.");
               }
           }
           
           // Build dependency list for the group column
           std::vector<std::string> colDependencies;
           for(const auto& p : abstractParticles) {
               std::string pName = p;
               if(p.find(_prefix) == 0) pName = p.substr(_prefix.length());
               colDependencies.push_back(_prefix + pName + _suffix);
           }
           
           // Redefine the group column to use the correct dependencies
           if(_reaction->ColumnExists(typedUniqueName)) {
               RedefineGroupColumn(typedUniqueName, colDependencies);
           } else {
               std::vector<std::string> clean_particles;
               for(const auto& col : colDependencies) clean_particles.push_back(col.substr(_prefix.length()));
               _reaction->SetGroupParticles(uniqueGroupName, _prefix, clean_particles);
           }
      }
      
      // Finalize Index Lookups
      _input2ReactionIndex.resize(_nameInputIndex.size());
      for(const auto& particle : _nameInputIndex) {
        _input2ReactionIndex[particle.second] = GetIndexSafe(particle.first);
      }
      ResolveDependencies();
    }
  
    inline void ParticleCreator::ResolveDependencies(){
        for (size_t i = 0; i < GetNCreated(); ++i) {
          RVecIndices vec_indices;
          for(const auto& type_index : _p_stru_depends[i]) {
            Indices_t indices;
            for(const auto& particle : type_index) indices.push_back(GetIndexSafe(particle));
            vec_indices.push_back(indices);
          }
          _p_dep_indices.push_back(vec_indices);
        }
    }
  
    inline bool ParticleCreator::HasParticle(const std::string& name) const { return _nameIndex.find(name) != _nameIndex.end(); }
    inline const IndexMap_t& ParticleCreator::GetIndexMap() const { return _nameIndex; }
    inline size_t ParticleCreator::GetNCreated() const { return _p_names.size(); }
    inline Int_t ParticleCreator::GetReactionIndex(const std::string& name) const { return GetIndexSafe(name); }
    inline Int_t ParticleCreator::GetReactionIndex(size_t input) const { return _input2ReactionIndex[input]; }
    inline std::string ParticleCreator::GetMapName() const { return _prefix + consts::ReactionMap() + _suffix + DoNotWriteTag(); }

    inline void ParticleCreator::AdoptIndices(const ParticleCreator& master) {
        for(const auto& name : _p_names) {
            if(!master.HasParticle(name)) throw std::runtime_error("Linked Topology Error: Particle '" + name + "' missing in Master.");
        }
        _nameIndex = master.GetIndexMap();
        _input2ReactionIndex = master._input2ReactionIndex; 
        _inputNames = master._inputNames;
        _nameInputIndex = master._nameInputIndex;
        _beam_names = master._beam_names;
    }

    inline void ParticleCreator::RebuildReactionMap() {
        auto get_indices = [&](const ParticleNames_t& names) {
          Indices_t ret;
          for(const auto& n : names) ret.push_back(GetIndexSafe(n));
          return ret;
        };

        auto get_group_indices = [&](int order) {
          if (_explicit_groups.count(order)) return get_indices(_explicit_groups[order]);
          if (_config_groups.count(order)) {
            try {
              std::string uniqueName = _config_groups[order] + _suffix;
              std::vector<std::string> particles;
              try { particles = _reaction->GetGroup(_prefix + uniqueName); } 
              catch(...) { particles = _reaction->GetGroup(_prefix + _config_groups[order]); }
              ParticleNames_t logicalNames;
              for(auto p : particles) {
                if(p.find(_prefix) == 0) logicalNames.push_back(p.substr(_prefix.length()));
                else logicalNames.push_back(p);
              }
              return get_indices(logicalNames);
            } catch (...) {}
          }
          return Indices_t{};
        };

        Indices_t idxBeam     = get_indices(_beam_names);
        Indices_t idxBaryons  = get_group_indices(consts::OrderBaryons());
        Indices_t idxMesons   = get_group_indices(consts::OrderMesons());
        Indices_t idxScat_ele = get_group_indices(consts::OrderScatEle());
        Indices_t idxDeps     = get_indices(GetPriorDependencies());
        Indices_t idxCreate   = get_indices(_p_names);

        RVecIndexMap retIndices{idxBeam, idxBaryons, idxMesons, idxScat_ele, idxDeps, idxCreate};
        _mapIndices = retIndices;
        _reaction->Define(GetMapName(), [retIndices]() { return retIndices; }, {});
        _reaction->SetGroupParticles(consts::KineIndices() + _suffix, _prefix, _inputNames);

        for(const auto& name : _inputNames) {
             std::string colName = _prefix + name + _suffix;
             std::string masterCol = _prefix + name;
             if(!_suffix.empty() && !_reaction->ColumnExists(colName) && _reaction->ColumnExists(masterCol)) {
                 _reaction->Define(colName, masterCol); 
             }
        }

        for(const auto& particle : _nameIndex) {
            bool isInput = std::find(_inputNames.begin(), _inputNames.end(), particle.first) != _inputNames.end();
            if (!isInput) {
                try { _reaction->SetParticleIndex(particle.first, _prefix, particle.second); } 
                catch (const std::invalid_argument&) {}
            }
            if (isInput) continue;
            std::string colName = _prefix  + particle.first + _suffix;
            auto func = [idx = particle.second](){ return ROOT::VecOps::RVec<int>{idx}; };
            if(_reaction->ColumnExists(colName)) _reaction->Redefine(colName, func, {});
            else _reaction->Define(colName, func, {});
        }

        for(auto const& [order, groupName] : _config_groups) {
               std::string uniqueGroupName = groupName + _suffix;
               std::string typedUniqueName = _prefix + uniqueGroupName;
               std::vector<std::string> abstractParticles;
               try { abstractParticles = _reaction->GetGroup(typedUniqueName); } 
               catch(...) { abstractParticles = _reaction->GetGroup(_prefix + groupName); }
               std::vector<std::string> colDependencies;
               for(const auto& p : abstractParticles) {
                   std::string pName = p;
                   if(p.find(_prefix) == 0) pName = p.substr(_prefix.length());
                   colDependencies.push_back(_prefix + pName + _suffix);
               }
               if(_reaction->ColumnExists(typedUniqueName)) RedefineGroupColumn(typedUniqueName, colDependencies);
               else {
                   std::vector<std::string> clean_particles;
                   for(const auto& col : colDependencies) clean_particles.push_back(col.substr(_prefix.length()));
                   _reaction->SetGroupParticles(uniqueGroupName, _prefix, clean_particles);
               }
        }
        
        _input2ReactionIndex.resize(_nameInputIndex.size());
        for(const auto& particle : _nameInputIndex) {
           _input2ReactionIndex[particle.second] = GetIndexSafe(particle.first);
        }
        ResolveDependencies();
    }  

    inline void ParticleCreator::ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, 
                       ROOT::RVecD& pz, ROOT::RVecD& m) const 
    {
       for (size_t i = 0; i < GetNCreated(); ++i) {
         _p_creators[i](GetIndexSafe(_p_names[i]), _p_dep_indices[i], px, py, pz, m);
      }
    }
} // end rad
