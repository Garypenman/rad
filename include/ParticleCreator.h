#pragma once

#include "ConfigReaction.h"
#include "StringUtilities.h"
#include "ParticleCreatorMethods.h"
#include <map>
#include <set>
#include <vector>
#include <string>
#include <stdexcept>

namespace rad {

    using StructuredNames_t = std::vector<ParticleNames_t>;
    using IndexMap_t = std::unordered_map<std::string, int>;
    
    /**
     * @brief Function pointer type for concrete, non-templated particle creation logic.
     * Arguments: (Target Index, Dependency Indices, Px, Py, Pz, M).
     */
    using ParticleCreatorFunc_t = void (*)(
        const Indice_t position, 
        const RVecIndices&, 
        ROOT::RVecD&, ROOT::RVecD&, ROOT::RVecD&, ROOT::RVecD&
    );
    
    /**
     * @class ParticleCreator
     * @brief Manages the definition, indexing, and creation of intermediate particles for a SINGLE data type.
     * * @details
     * This class acts as the "Architect" of the reaction topology for a specific prefix (e.g., "rec_" or "tru_").
     * It maps human-readable names ("Jpsi") to efficient integer indices in the structure-of-arrays data.
     * * **Architectural Modes:**
     * 1. **Master Mode:**
     * - Calculates new indices for all particles.
     * - Defines the primary `ReactionMap` column in RDataFrame.
     * - Used for the primary reconstruction or truth processor.
     * 2. **Linked Mode:**
     * - "Adopts" the index structure from a Master creator via `AdoptIndices`.
     * - Ensures that the "tru_" topology matches the "rec_" topology exactly 1:1.
     * - Reuses the indices for zero-cost topological alignment.
     */
    class ParticleCreator {

    public:
      
      ParticleCreator() = default;

      /**
       * @brief Standard Constructor.
       * @param cr Pointer to the shared ConfigReaction interface.
       * @param prefix The data type prefix this instance manages (e.g., "rec_", "tru_").
       * @param suffix Optional suffix for output columns (e.g. "_miss" for background studies).
       */
      explicit ParticleCreator(ConfigReaction* cr, const std::string& prefix, const std::string& suffix = "") 
        : _reaction{cr}, _prefix{prefix}, _suffix{suffix} {
          // Default Beam Names (Abstract)
          _beam_names = {consts::BeamIon(), consts::BeamEle()};
          
          // Default Group Mappings
          // These are abstract IDs; InitMap will prepend _prefix to look them up.
          DefineGroup(consts::OrderBaryons(), consts::Baryons());
          DefineGroup(consts::OrderMesons(), consts::Mesons());
      }

      /**
       * @brief Copy Constructor (The Fork).
       * Copies configuration but allows divergence via a new suffix.
       * Does NOT copy runtime indices (must call Init/InitLinked again).
       */
      ParticleCreator(const ParticleCreator& other, const std::string& new_suffix) 
        : _reaction(other._reaction),
          _prefix(other._prefix), // Keeps the same type (e.g. forking rec -> rec_miss)
          _suffix(new_suffix),
          _beam_names(other._beam_names),
          _config_groups(other._config_groups),
          _explicit_groups(other._explicit_groups),
          _p_names(other._p_names),
          _p_creators(other._p_creators),
          _p_required(other._p_required),
          _p_stru_depends(other._p_stru_depends),
          _dependencies(other._dependencies),
          _inputNames(other._inputNames),
          _createIndex(other._createIndex)
      {}
      
      void SetReaction(ConfigReaction* reaction) { _reaction = reaction; }
      ConfigReaction* Reaction() const { return _reaction; }
      std::string GetPrefix() const { return _prefix; }

      // =======================================================================
      // Configuration Interface
      // =======================================================================

      /**
       * @brief Maps a fixed ParticleGroupOrder slot to a ConfigReaction group name.
       * @param order The group order ID (e.g. OrderMesons).
       * @param groupName The abstract group name (e.g. "Mesons"). InitMap will look up "_prefix + Mesons".
       */
      void DefineGroup(int order, const std::string& groupName) {
          _config_groups[order] = groupName;
          _explicit_groups.erase(order); 
      }

      /**
       * @brief Manually defines the particles in a group (Override).
       * Useful for defining "Missing Mass" topologies (e.g. Baryon = "n_miss").
       * @param order The group order ID.
       * @param particles List of abstract particle names (e.g. "pi+"). 
       */
      void OverrideGroup(int order, const std::vector<std::string>& particles) {
          _explicit_groups[order] = particles;
          _config_groups.erase(order);   
      }

      void SetBeamNames(const std::vector<std::string>& beams) {
          _beam_names = beams;
      }

      // =======================================================================
      // Particle Definition Logic
      // =======================================================================

      /**
       * @brief Registers a new intermediate particle calculation.
       * @param name Name of the new particle (e.g., "Miss").
       * @param func Function pointer to the calculation logic.
       * @param depends List of dependencies (other particle names).
       */
      void AddParticle(const std::string& name, ParticleCreatorFunc_t func, const StructuredNames_t& depends={{}}) {
        _p_names.push_back(name);
        
        auto flat_depends = util::flattenColumnNames(depends);
        _p_required.push_back(flat_depends);
        _dependencies.insert(flat_depends.begin(), flat_depends.end());
        _p_stru_depends.push_back(depends);
        _p_creators.push_back(func);
        _createIndex[name] = GetNCreated() - 1;
        
        // Register the column for THIS prefix only
       _reaction->setParticleIndex(name, _prefix, consts::InvalidEntry<int>());
      }

      /** @brief Create a particle by Summing 4-vectors (e.g. J/psi -> e+ e-). */
      void Sum(const std::string& name, const StructuredNames_t& depends={{}}) {
        AddParticle(name, ParticleCreateBySum, depends);
      }

      /** @brief Create a particle by Difference (e.g. Missing Neutron). */
      void Diff(const std::string& name, const StructuredNames_t& depends={{}}) {
        AddParticle(name, ParticleCreateByDiff, depends);
      }

      // =======================================================================
      // Indexing & Initialization 
      // =======================================================================

      /**
       * @brief MASTER MODE: Calculates indices and defines the ReactionMap.
       * Call this only for the primary processor of a type.
       */
      void InitMap();

      bool HasParticle(const std::string& name) const {
          return _nameIndex.find(name) != _nameIndex.end();
      }

      const IndexMap_t& GetIndexMap() const { return _nameIndex; }

      /** @brief Translates dependency strings into integer indices for runtime creation. */
      void ResolveDependecies();
      
      /**
       * @brief LINKED MODE: Adopts indices from a Master creator.
       * Ensures that this creator uses the same memory slots as the Master.
       */
      void AdoptIndices(const ParticleCreator& master);

      /**
       * @brief LINKED MODE: Builds a local ReactionMap using adopted indices.
       * Allows for different group definitions (e.g. Baryon="n_miss") using the same data pool.
       */
      void RebuildReactionMap();

      // =======================================================================
      // Execution
      // =======================================================================

      /**
       * @brief Runs the particle creation logic (the actual math).
       * Called inside the Event Loop.
       */
      void ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, 
                         ROOT::RVecD& pz, ROOT::RVecD& m) const;

      // =======================================================================
      // Accessors
      // =======================================================================

      size_t GetNCreated() const { return _p_names.size(); }
      
      /** @brief Get the fixed array index for a particle name. */
      Int_t GetReactionIndex(const std::string& name) const { return _nameIndex.at(name); }
      
      Int_t GetReactionIndex(size_t input) const { return _input2ReactionIndex[input]; }
      
      /** @brief Returns RDataFrame column name: e.g., "rec_ReactionMap" */
      std::string GetMapName() const { return _prefix + consts::ReactionMap() + _suffix + DoNotWriteTag(); }

      std::vector<std::string> GetPriorDependencies();
      Indices_t CreateIndices(IndexMap_t& nameIndex, const ParticleNames_t& names, size_t& idx);

    private:
      ConfigReaction* _reaction = nullptr;
      std::string _prefix; // e.g. "rec_"
      std::string _suffix; // e.g. "_miss"
      
      // Configuration
      std::vector<std::string> _beam_names;
      std::map<int, std::string> _config_groups; 
      std::map<int, std::vector<std::string>> _explicit_groups; 

      // Particle Definitions
      ParticleNames_t _p_names;
      std::vector<ParticleCreatorFunc_t> _p_creators;
      std::vector<ParticleNames_t> _p_required;
      std::vector<StructuredNames_t> _p_stru_depends;
      ROOT::RVec<RVecIndices> _p_dep_indices; // Resolved indices for creators
      std::set<std::string> _dependencies;
      
      // Index Maps
      ParticleNames_t _inputNames;
      RVecIndexMap _mapIndices;
      IndexMap_t _nameIndex;       // Name -> ReactionIndex
      IndexMap_t _nameInputIndex;  // Name -> InputIndex
      IndexMap_t _createIndex;
      Indices_t _input2ReactionIndex; // InputIndex -> ReactionIndex
    };

    // ========================================================================
    // Implementation
    // ========================================================================

    inline std::vector<std::string> ParticleCreator::GetPriorDependencies() {
        std::vector<std::string> vec_deps(_dependencies.begin(), _dependencies.end());
        util::removeExistingStrings(vec_deps, _p_names);
        return vec_deps;
    }

    inline Indices_t ParticleCreator::CreateIndices(IndexMap_t& nameIndex, const ParticleNames_t& names, size_t& idx) {
        Indices_t indices{};
        for(const auto& name : names) {
          if(nameIndex.count(name) == 0) {
            indices.push_back(idx);
            nameIndex[name] = idx;
            ++idx;
          } else {
             indices.push_back(nameIndex[name]);
          }
        }
        return indices;
    }

    inline void ParticleCreator::InitMap() {
      // 1. Gather Input Particles for THIS prefix
      // Logic: Abstract names (pi+) -> lookup Typed Columns (rec_pi+) -> Store Abstract Name
      ParticleNames_t logicalInputNames;
      logicalInputNames.insert(logicalInputNames.end(), _beam_names.begin(), _beam_names.end());

      std::vector<int> orders = {consts::OrderBaryons(), consts::OrderMesons(), consts::OrderScatEle()};
      
      for(int order : orders) {
          if(_explicit_groups.count(order)) {
              const auto& list = _explicit_groups[order];
              logicalInputNames.insert(logicalInputNames.end(), list.begin(), list.end());
          } 
          else if (_config_groups.count(order)) {
              try {
                  // Lookup e.g., "rec_Mesons"
                 std::string typedGroupName = _prefix + _config_groups[order];
                  auto particles = _reaction->getGroup(typedGroupName);
                  
                  // Strip prefix (rec_pi+ -> pi+) for internal mapping
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

      // 2. Register KineIndices for THIS prefix
      // "rec_KineIndices" -> { "rec_pip", "rec_Beam", ... }
      _reaction->setGroupParticles(consts::KineIndices() + _suffix, _prefix, _inputNames);

      // 3. Build Index Maps
      size_t in_idx = 0;
      CreateIndices(_nameInputIndex, _inputNames, in_idx);
        
      size_t idx = 0;
      Indices_t idxBeam = CreateIndices(_nameIndex, _beam_names, idx);
      
      auto get_indices_for_group = [&](int order) {
          if (_explicit_groups.count(order)) return CreateIndices(_nameIndex, _explicit_groups[order], idx);
          if (_config_groups.count(order)) {
              try {
                  std::string typedGroupName = _prefix + _config_groups[order];
                  auto particles = _reaction->getGroup(typedGroupName);
                  
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
      
      // Define helper scalar columns: "rec_Jpsi" -> 6
      for(const auto& particle : _nameIndex) {
	std::string colName = _prefix  + particle.first + _suffix;
	auto func = [idx = particle.second](){ return idx; };
	
	if(_reaction->ColumnExists(colName)) {
	  _reaction->Redefine(colName, func, {});
	} else {
	  _reaction->Define(colName, func, {});
	}
      }
      
      _input2ReactionIndex.resize(_nameInputIndex.size());
      for(const auto& particle : _nameInputIndex) {
        _input2ReactionIndex[particle.second] = GetReactionIndex(particle.first);
      }
      
      // 4. Resolve Dependencies
      ResolveDependecies();
    }
  
  inline void ParticleCreator::ResolveDependecies(){
    for (size_t i = 0; i < GetNCreated(); ++i) {
      RVecIndices vec_indices;
      for(const auto& type_index : _p_stru_depends[i]) {
	Indices_t indices;
	for(const auto& particle : type_index) indices.push_back(GetReactionIndex(particle));
	vec_indices.push_back(indices);
      }
      _p_dep_indices.push_back(vec_indices);
    }
  }
  
  inline void ParticleCreator::AdoptIndices(const ParticleCreator& master) {
        for(const auto& name : _p_names) {
            if(!master.HasParticle(name)) {
                throw std::runtime_error("Linked Topology Error: Particle '" + name + 
                    "' required by Linked processor is missing in Master.");
            }
        }
        _nameIndex = master.GetIndexMap();
        // Crucial: Copy the input mapping so we point to the correct slots in the vector
        _input2ReactionIndex = master._input2ReactionIndex; 
    }

  inline void ParticleCreator::RebuildReactionMap() {
    // 1. Build the local group lists, but use the ADOPTED indices
    auto get_indices = [&](const ParticleNames_t& names) {
      Indices_t ret;
      for(const auto& n : names) ret.push_back(_nameIndex.at(n));
      return ret;
    };

    auto get_group_indices = [&](int order) {
      if (_explicit_groups.count(order)) return get_indices(_explicit_groups[order]);
      if (_config_groups.count(order)) {
	try {
	  std::string typedGroupName = _prefix + _config_groups[order];
	  auto particles = _reaction->getGroup(typedGroupName);
                
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
    
    // 2. Define the Map (The data itself)
    _reaction->Define(GetMapName(), [retIndices]() { return retIndices; }, {});

    // 3. Re-Register Group Particles for KinematicsDispatch
    _reaction->setGroupParticles(consts::KineIndices() + _suffix, _prefix, _inputNames);

    ResolveDependecies();
  }  

  inline void ParticleCreator::ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, 
					     ROOT::RVecD& pz, ROOT::RVecD& m) const 
  {
     for (size_t i = 0; i < GetNCreated(); ++i) {
       _p_creators[i](_nameIndex.at(_p_names[i]), _p_dep_indices[i], px, py, pz, m);
    }
  }
} // end rad
