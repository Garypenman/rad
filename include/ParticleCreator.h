/**
 * @file ParticleCreator.h
 * @brief Manages the definition, indexing, and creation of intermediate particles.
 * @details
 * This class acts as the "Topology Builder" for the analysis. It is responsible for:
 * 1. Registering Input Particles (from the source TTree).
 * 2. Defining Composite Particles (Resonances, Missing Mass, etc.).
 * 3. Resolving string names to efficient integer indices for the kernels.
 * 4. Managing Particle Groups (Mesons, Baryons) for combinatorics.
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

    using StructuredNames_t = ROOT::RVec<ParticleNames_t>;
    using IndexMap_t = std::unordered_map<std::string, int>;
    
    using ParticleCreatorFunc_t = void (*)(
        const Indice_t position, 
        const RVecIndices&, 
        ROOT::RVecD&, ROOT::RVecD&, ROOT::RVecD&, ROOT::RVecD&
    );
    
    /**
     * @class ParticleCreator
     * @brief Manages the creation of intermediate particles and the resolution of reaction topology.
     */
    class ParticleCreator {

    public:
      
      /** @brief Default constructor. */
      ParticleCreator() = default;

      /**
       * @brief Constructor for a new topology configuration.
       * @param cr Pointer to the ConfigReaction object managing RDF nodes.
       * @param prefix Prefix for branch names (e.g., "rec_", "mc_").
       * @param suffix Suffix for branch names (e.g., "", "_sysUp").
       */
      explicit ParticleCreator(ConfigReaction* cr, const std::string& prefix, const std::string& suffix = "");

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

      /**
       * @brief Associates a standard group ID (e.g., Baryons) with a group name.
       * @param order The integer ID of the group (from rad::consts).
       * @param groupName The string name of the group (e.g., "Baryons").
       */
      void DefineGroup(int order, const std::string& groupName);

      /**
       * @brief Overrides a standard group with an explicit list of particles.
       * @param order The integer ID of the group.
       * @param particles Vector of particle names to include in this group.
       */
      void OverrideGroup(int order, const ROOT::RVec<std::string>& particles);

      /**
       * @brief Helper to override a group by its abstract name.
       * @param abstractName The abstract name (e.g., "Baryons", "Mesons").
       * @param particles Vector of particle names to include.
       */
      void OverrideGroup(const std::string& abstractName, const ROOT::RVec<std::string>& particles);

      /**
       * @brief Sets the names of the beam particles.
       * @param beams Vector of beam particle names (usually {Ion, Electron}).
       */
      void SetBeamNames(const ROOT::RVec<std::string>& beams);
      
      /**
       * @brief Forces an input particle to be registered in the map.
       * @details Essential if a particle exists in the input tree but isn't part of any default group.
       * @param name The name of the input particle to require.
       */
      void RequireParticle(const std::string& name);

      // =======================================================================
      // Particle Definition Logic
      // =======================================================================

      /**
       * @brief Registers a generic particle with a custom creation function.
       * @param name The unique name of the new particle.
       * @param func Function pointer to the creation logic.
       * @param depends List of dependencies required by the creation function.
       */
      void AddParticle(const std::string& name, ParticleCreatorFunc_t func, const StructuredNames_t& depends={{}});

      /**
       * @brief Defines a particle as the 4-vector Sum of others.
       * @param name The name of the composite particle (e.g., "Jpsi").
       * @param depends List of daughter particles (e.g., {{"ele", "pos"}}).
       */
      void Sum(const std::string& name, const StructuredNames_t& depends={{}});

      /**
       * @brief Defines a particle as Sum and registers it for Truth Matching.
       * @param name The name of the composite particle.
       * @param truthRole The integer Role ID for truth matching comparison.
       * @param depends List of daughter particles.
       */
      void SumTruthMatch(const std::string& name, int truthRole, const StructuredNames_t& depends={{}});

      /**
       * @brief Defines a particle via subtraction (Missing Mass).
       * @details P_miss = Sum(Initial) - Sum(Final).
       * @param name The name of the missing particle (e.g., "n_miss").
       * @param depends List containing {{Initial State}, {Final State}}.
       */
      void Diff(const std::string& name, const StructuredNames_t& depends={{}});

      // =======================================================================
      // Indexing & Initialization 
      // =======================================================================

      /**
       * @brief Resolves indices, registers columns, and builds the reaction map.
       * @details 
       * 1. Gathers all inputs and created particles.
       * 2. Checks for collisions.
       * 3. Assigns integer indices.
       * 4. Registers the map with RDataFrame.
       */
      void InitMap();

      /** @brief Checks if a particle exists in the current map. */
      bool HasParticle(const std::string& name) const;

      /** @return Reference to the internal name-to-index map. */
      const IndexMap_t& GetIndexMap() const;

      /** @brief Resolves internal dependency indices for all created particles. */
      void ResolveDependencies();
      
      // =======================================================================
      // Execution
      // =======================================================================

      /**
       * @brief Executes the particle creation functions to generate 4-vectors.
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
      
      /** @return The integer index for a named particle. Throws if not found. */
      Int_t GetReactionIndex(const std::string& name) const;
      
      /** @return The reaction index corresponding to an input index. */
      Int_t GetReactionIndex(size_t input) const;
      
      /** @return The name of the RDataFrame column containing the map indices. */
      std::string GetMapName() const;

      /** @return A list of input columns required by this creator (calculated dependencies). */
      ROOT::RVec<std::string> GetPriorDependencies();
      
      /** @brief Assigns indices to a list of names. */
      Indices_t CreateIndices(IndexMap_t& nameIndex, const ParticleNames_t& names, size_t& idx);

      /** @return A list of all particle names currently registered in the index map. */
      ROOT::RVec<std::string> GetParticleNames() const {
        ROOT::RVec<std::string> names;
          names.reserve(_nameIndex.size());
          for(const auto& pair : _nameIndex) {
              names.push_back(pair.first);
          }
          return names;
      }

      //////////////////////
      ///Diagnostics
      ///\brief Print the reaction map.
      void PrintReactionMap() const;

      ///\brief Print all particle groups.
      void PrintGroups() const;

      ///\brief Print creation function registry.
      void PrintCreationRegistry() const;

      ///\brief Print column aliases.
      void PrintAliases() const;

      ///\brief Print comprehensive diagnostics.
      void PrintDiagnostics() const;

    private:
      ConfigReaction* _reaction = nullptr;
      std::string _prefix; 
      std::string _suffix; 
      
      ROOT::RVec<std::string> _beam_names;
      ROOT::RVec<std::string> _forced_inputs; 

      std::map<int, std::string> _config_groups; 
      std::map<int, ROOT::RVec<std::string>> _explicit_groups; 

      ParticleNames_t _p_names;
      ROOT::RVec<ParticleCreatorFunc_t> _p_creators;
      ROOT::RVec<ParticleNames_t> _p_required;
      ROOT::RVec<StructuredNames_t> _p_stru_depends;
      ROOT::RVec<RVecIndices> _p_dep_indices; 
      std::set<std::string> _dependencies;
      
      ParticleNames_t _inputNames;
      RVecIndexMap _mapIndices;
      IndexMap_t _nameIndex;        
      IndexMap_t _nameInputIndex;   
      IndexMap_t _createIndex;
      Indices_t _input2ReactionIndex; 

      void RedefineGroupColumn(const std::string& name, const ROOT::RVec<std::string>& cols);
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

    inline void ParticleCreator::SetReaction(ConfigReaction* reaction) { _reaction = reaction; }
    inline ConfigReaction* ParticleCreator::Reaction() const { return _reaction; }
    inline std::string ParticleCreator::GetPrefix() const { return _prefix; }
    inline std::string ParticleCreator::GetSuffix() const { return _suffix; }

    inline void ParticleCreator::DefineGroup(int order, const std::string& groupName) {
          _config_groups[order] = groupName;
          _explicit_groups.erase(order); 
    }
    inline void ParticleCreator::OverrideGroup(int order, const ROOT::RVec<std::string>& particles) {
          _explicit_groups[order] = particles;
          _config_groups.erase(order);     
    }
    inline void ParticleCreator::OverrideGroup(const std::string& abstractName, const ROOT::RVec<std::string>& particles) {
          if(abstractName == consts::Baryons()) OverrideGroup(consts::OrderBaryons(), particles);
          else if(abstractName == consts::Mesons()) OverrideGroup(consts::OrderMesons(), particles);
    }
    inline void ParticleCreator::SetBeamNames(const ROOT::RVec<std::string>& beams) { _beam_names = beams; }
    
    inline void ParticleCreator::RequireParticle(const std::string& name) {
          if(std::find(_forced_inputs.begin(), _forced_inputs.end(), name) == _forced_inputs.end())
              _forced_inputs.push_back(name);
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

    inline ROOT::RVec<std::string> ParticleCreator::GetPriorDependencies() {
        ROOT::RVec<std::string> vec_deps(_dependencies.begin(), _dependencies.end());
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

    inline void ParticleCreator::InitMap() {
      // 1. GATHER INPUT PARTICLES
      ParticleNames_t logicalInputNames;
      logicalInputNames.insert(logicalInputNames.end(), _beam_names.begin(), _beam_names.end());
      logicalInputNames.insert(logicalInputNames.end(), _forced_inputs.begin(), _forced_inputs.end());

      ROOT::RVec<int> orders = {consts::OrderBaryons(), consts::OrderMesons(), consts::OrderScatEle()};
      for(int order : orders) {
          if(_explicit_groups.count(order)) {
              const auto& list = _explicit_groups[order];
              logicalInputNames.insert(logicalInputNames.end(), list.begin(), list.end());
          } 
          else if (_config_groups.count(order)) {
              try {
                 std::string uniqueGroupName = _config_groups[order] + _suffix;
                 ROOT::RVec<std::string> particles;
                 try { particles = _reaction->GetGroup(_prefix + uniqueGroupName); } 
                 catch(...) { particles = _reaction->GetGroup(_prefix + _config_groups[order]); }
                 
                 for(const auto& p : particles) {
                     if(p.find(_prefix) == 0) logicalInputNames.push_back(p.substr(_prefix.length()));
                     else logicalInputNames.push_back(p); 
                 }
              } catch (...) { }
          }
      }

      auto dep_names = GetPriorDependencies();
      util::removeExistingStrings(dep_names, logicalInputNames);
      
      _inputNames = logicalInputNames;
      _inputNames.insert(_inputNames.end(), dep_names.begin(), dep_names.end());
      util::removeExistingStrings(_inputNames, _p_names); 

      // 2. COLLISION CHECK & MAP REGISTRATION
      std::string kineCol = consts::KineIndices() + _suffix;
      if (_reaction->ColumnExists(_prefix + kineCol)) {
          throw std::runtime_error("Topology Collision: Processor [" + _prefix + "] suffix [" + _suffix + "] already exists.");
      }
      _reaction->SetGroupParticles(kineCol, _prefix, utils::as_stdvector(_inputNames));

      // 3. BUILD INDEX MAPS
      size_t in_idx = 0;
      CreateIndices(_nameInputIndex, _inputNames, in_idx);
        
      size_t idx = 0;
      Indices_t idxBeam = CreateIndices(_nameIndex, _beam_names, idx);
      
      auto get_indices_for_group = [&](int order) {
          if (_explicit_groups.count(order)) return CreateIndices(_nameIndex, _explicit_groups[order], idx);
          if (_config_groups.count(order)) {
              try {
                  std::string uniqueName = _config_groups[order] + _suffix;
                  ROOT::RVec<std::string> particles;
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
      Indices_t idxCreate   = CreateIndices(_nameIndex, _p_names, idx); 

      RVecIndexMap retIndices{idxBeam, idxBaryons, idxMesons, idxScat_ele, idxDeps, idxCreate};
      _mapIndices = retIndices;
      _reaction->Define(GetMapName(), [retIndices]() { return retIndices; }, {});

      // DEFINE INPUT ALIASES 
      for(const auto& name : _inputNames) {
            std::string colName = _prefix + name + _suffix;
            std::string masterCol = _prefix + name;
            
            if(!_suffix.empty() && !_reaction->ColumnExists(colName) && _reaction->ColumnExists(masterCol)) {
                _reaction->Define(colName, masterCol); 
            }
      }
      
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
    inline Indice_t ParticleCreator::GetReactionIndex(const std::string& name) const { return GetIndexSafe(name); }
    inline Indice_t ParticleCreator::GetReactionIndex(size_t input) const { return _input2ReactionIndex[input]; }
    inline std::string ParticleCreator::GetMapName() const { return _prefix + consts::ReactionMap() + _suffix + DoNotWriteTag(); }

    inline void ParticleCreator::ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, 
                        ROOT::RVecD& pz, ROOT::RVecD& m) const 
    {
       for (size_t i = 0; i < GetNCreated(); ++i) {
         _p_creators[i](GetIndexSafe(_p_names[i]), _p_dep_indices[i], px, py, pz, m);
      }
    }
} // end rad

#include "Diagnostics_ParticleCreator.hxx"
