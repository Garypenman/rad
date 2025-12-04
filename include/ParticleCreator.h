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
namespace config {

    using StructuredNames_t = std::vector<ParticleNames_t>;
    using IndexMap_t = std::unordered_map<std::string, int>;
    
    /**
     * @brief Function pointer type for concrete, non-templated particle creation logic.
     */
    using ParticleCreatorFunc_t = void (*)(
        const Indice_t position, 
        const RVecIndices&, 
        ROOT::RVecD&, ROOT::RVecD&, ROOT::RVecD&, ROOT::RVecD&
    );
    
    /**
     * @brief Manages the definition, indexing, and creation of intermediate particles.
     * * This class acts as the "Architect" of the reaction topology.
     * * Features:
     * 1. Topology Agnostic: Configures groups via the controlling Kinematics class.
     * 2. Twin Topology Support: Can copy configuration and add suffixes to avoid collisions.
     * 3. Index Adoption: Can "borrow" indices from a Master creator for zero-cost re-analysis.
     */
    class ParticleCreator {

    public:
      
      ParticleCreator() = default;

      /**
       * @brief Standard Constructor.
       * Sets up default groups (Beams, Baryons, Mesons).
       * @param cr The ConfigReaction interface.
       * @param suffix Optional suffix for all output columns (e.g. "_miss").
       */
      explicit ParticleCreator(ConfigReaction* cr, const std::string& suffix = "") 
        : _reaction{cr}, _suffix{suffix} {
          // Default Beam Names
          _beam_names = {names::BeamIon(), names::BeamEle()};
          
          // Default Group Mappings (Look up in ConfigReaction)
          DefineGroup(names::OrderBaryons(), names::Baryons());
          DefineGroup(names::OrderMesons(), names::Mesons());
      }

      /**
       * @brief Copy Constructor (The Fork).
       * Copies all definitions (Jpsi, Z, etc.) from another creator but sets a new suffix.
       * Used to create a divergent topology from a common base.
       */
      ParticleCreator(const ParticleCreator& other, const std::string& new_suffix) 
        : _reaction(other._reaction),
          _suffix(new_suffix),
          _beam_names(other._beam_names),
          _config_groups(other._config_groups),
          _explicit_groups(other._explicit_groups),
          // Copy internal definitions
          _p_names(other._p_names),
          _p_creators(other._p_creators),
          _p_required(other._p_required),
          _p_stru_depends(other._p_stru_depends),
          _dependencies(other._dependencies),
          _createIndex(other._createIndex)
          // Runtime maps (_nameIndex, etc.) are NOT copied; they are rebuilt on Init.
      {}
      
      void SetReaction(ConfigReaction* reaction) { _reaction = reaction; }
      ConfigReaction* Reaction() const { return _reaction; }

      // =======================================================================
      // Configuration Interface
      // =======================================================================

      /**
       * @brief Maps a fixed ParticleGroupOrder slot to a ConfigReaction group name.
       * Example: DefineGroup(names::OrderScatEle(), names::ScatGroup())
       */
      void DefineGroup(int order, const std::string& groupName) {
          _config_groups[order] = groupName;
          _explicit_groups.erase(order); // Priority rule: Specific definition clears override
      }

      /**
       * @brief Manually defines the particles in a group (Override).
       * Bypasses ConfigReaction lookup. Useful for "n_miss" topologies.
       * Example: OverrideGroup(names::OrderBaryons(), {"n_miss"})
       */
      void OverrideGroup(int order, const std::vector<std::string>& particles) {
          _explicit_groups[order] = particles;
          _config_groups.erase(order);   // Clear config lookup for this slot
      }

      void SetBeamNames(const std::vector<std::string>& beams) {
          _beam_names = beams;
      }

      // =======================================================================
      // Particle Definition Logic
      // =======================================================================

      void AddParticle(const std::string& name, ParticleCreatorFunc_t func, const StructuredNames_t& depends={{}}) {
        _p_names.push_back(name);
        
        auto flat_depends = utils::flattenColumnNames(depends);
        _p_required.push_back(flat_depends);
        _dependencies.insert(flat_depends.begin(), flat_depends.end());
        _p_stru_depends.push_back(depends);
        _p_creators.push_back(func);
        _createIndex[name] = GetNCreated() - 1;
        
        _reaction->setParticleIndex(name, constant::InvalidEntry<int>());
      }

      void Sum(const std::string& name, const StructuredNames_t& depends={{}}) {
        AddParticle(name, ParticleCreateBySum, depends);
      }

      void Diff(const std::string& name, const StructuredNames_t& depends={{}}) {
        AddParticle(name, ParticleCreateByDiff, depends);
      }

      // =======================================================================
      // Indexing & Initialization (Master Mode)
      // =======================================================================

      /**
       * @brief Calculates the static ReactionMap and defines RDataFrame columns.
       * Used by the MASTER processor.
       */
      void InitMap(const std::vector<std::string>& types);

      // =======================================================================
      // Indexing & Initialization (Linked/Adoption Mode)
      // =======================================================================

      /**
       * @brief Check if a particle is defined in this creator.
       */
      bool HasParticle(const std::string& name) const {
          return _nameIndex.find(name) != _nameIndex.end();
      }

      const IndexMap_t& GetIndexMap() const { return _nameIndex; }

      /**
       * @brief Adoption Mode: Copy indices from Master instead of calculating new ones.
       * This ensures the Linked processor reads the correct slots in the shared array.
       */
      void AdoptIndices(const ParticleCreator& master);

      /**
       * @brief Rebuilds the ReactionMap using Adopted indices but Local groups.
       * Used by the LINKED processor.
       */
      void RebuildReactionMap();

      // =======================================================================
      // Execution
      // =======================================================================

      /**
       * @brief Executes creation logic for an event.
       */
      void ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, 
                         ROOT::RVecD& pz, ROOT::RVecD& m) const;

      // =======================================================================
      // Accessors & Helpers
      // =======================================================================

      size_t GetNCreated() const { return _p_names.size(); }
      
      // Index Accessors
      Int_t GetReactionIndex(const std::string& name) const { return _nameIndex.at(name); }
      Int_t GetReactionIndex(size_t input) const { return _input2ReactionIndex[input]; }
      
      // Get the suffixed name of the ReactionMap column
      std::string GetMapName() const { return names::ReactionMap() + _suffix + config::DoNotWriteTag(); }

      std::vector<std::string> GetPriorDependencies();
      Indices_t CreateIndices(IndexMap_t& nameIndex, const ParticleNames_t& names, size_t& idx);

    private:
      ConfigReaction* _reaction = nullptr;
      std::string _suffix;
      
      // Configuration Storage
      std::vector<std::string> _beam_names;
      std::map<int, std::string> _config_groups; // Map<OrderEnum, ConfigGroupName>
      std::map<int, std::vector<std::string>> _explicit_groups; // Map<OrderEnum, ExplicitList>

      // Particle Definitions
      ParticleNames_t _p_names;
      std::vector<ParticleCreatorFunc_t> _p_creators;
      std::vector<ParticleNames_t> _p_required;
      std::vector<StructuredNames_t> _p_stru_depends;
      ROOT::RVec<RVecIndices> _p_dep_indices;
      std::set<std::string> _dependencies;
      
      // Runtime Maps
      ParticleNames_t _inputNames;
      RVecIndexMap _mapIndices;
      IndexMap_t _nameIndex;
      IndexMap_t _nameInputIndex;
      IndexMap_t _createIndex;
      Indices_t _input2ReactionIndex;
    };

    // ========================================================================
    // Implementation
    // ========================================================================

    inline std::vector<std::string> ParticleCreator::GetPriorDependencies() {
        std::vector<std::string> vec_deps(_dependencies.begin(), _dependencies.end());
        utils::removeExistingStrings(vec_deps, _p_names);
        return vec_deps;
    }

    inline Indices_t ParticleCreator::CreateIndices(IndexMap_t& nameIndex, const ParticleNames_t& names, size_t& idx) {
        Indices_t indices{};
        for(const auto& name : names) {
          // If already mapped, skip (unless we want to enforce re-mapping?)
          // Standard logic: First come, first served for index assignment.
          if(nameIndex.count(name) == 0) {
            indices.push_back(idx);
            nameIndex[name] = idx;
            ++idx;
          } else {
             // If existing, just reuse the index (don't increment idx)
             indices.push_back(nameIndex[name]);
          }
        }
        return indices;
    }

    inline void ParticleCreator::InitMap(const std::vector<std::string>& types) {
      
      // 1. Gather Input Particles
      ParticleNames_t all_group_particles;
      
      // Beams
      all_group_particles.insert(all_group_particles.end(), _beam_names.begin(), _beam_names.end());

      // Iterate over ALL possible group slots (Baryons, Mesons, ScatEle...)
      // We check Explicit overrides first, then Config lookups.
      // We scan a reasonable range of enums or specific known ones.
      std::vector<int> orders = {names::OrderBaryons(), names::OrderMesons(), names::OrderScatEle()};
      
      for(int order : orders) {
          if(_explicit_groups.count(order)) {
              // Use Explicit List
              const auto& list = _explicit_groups[order];
              all_group_particles.insert(all_group_particles.end(), list.begin(), list.end());
          } else if (_config_groups.count(order)) {
              // Use Config Lookup
              try {
                 auto particles = _reaction->getGroup(_config_groups[order]);
                 all_group_particles.insert(all_group_particles.end(), particles.begin(), particles.end());
              } catch (...) { /* Group undefined in config, skip */ }
          }
      }

      // Dependencies
      auto dep_names = GetPriorDependencies();
      utils::removeExistingStrings(dep_names, all_group_particles);
      
      // Master List
      _inputNames = all_group_particles;
      _inputNames.insert(_inputNames.end(), dep_names.begin(), dep_names.end());
      utils::removeExistingStrings(_inputNames, _p_names);

      // Define KineIndices in RDF (with Suffix if needed, usually Master defines these)
      // Note: KineIndices are generally shared if topologies align, but safest to suffix them 
      // if this is a standalone init.
      for(const auto& type : types) {
        auto typeNames = utils::prependToAll(_inputNames, type);
        _reaction->setGroupParticles(type + names::KineIndices() + _suffix, typeNames); 
      }

      size_t in_idx = 0;
      CreateIndices(_nameInputIndex, _inputNames, in_idx);
       
      // 2. Calculate Fixed Positions (ReactionMap)
      size_t idx = 0;
      
      // A. Beams
      Indices_t idxBeam = CreateIndices(_nameIndex, _beam_names, idx);
      
      // Helper to fetch indices based on configuration state
      auto get_indices_for_group = [&](int order) {
          if (_explicit_groups.count(order)) {
              return CreateIndices(_nameIndex, _explicit_groups[order], idx);
          }
          if (_config_groups.count(order)) {
              try {
                  return CreateIndices(_nameIndex, _reaction->getGroup(_config_groups[order]), idx);
              } catch (...) {}
          }
          return Indices_t{};
      };

      Indices_t idxBaryons  = get_indices_for_group(names::OrderBaryons());
      Indices_t idxMesons   = get_indices_for_group(names::OrderMesons());
      Indices_t idxScat_ele = get_indices_for_group(names::OrderScatEle());

      Indices_t idxDeps = CreateIndices(_nameIndex, dep_names, idx);
      Indices_t idxCreate = CreateIndices(_nameIndex, _p_names, idx);

      // Build Map
      RVecIndexMap retIndices{idxBeam, idxBaryons, idxMesons, idxScat_ele, idxDeps, idxCreate};
      
      _mapIndices = retIndices;
      _reaction->Define(GetMapName(), [retIndices]() { return retIndices; }, {});
      
      // Helper Columns
      for(const auto& particle : _nameIndex) {
        _reaction->Define(names::data_type::Kine() + particle.first + _suffix, [idx = particle.second](){ return idx; }, {});
      }

      // Conversion Map
      _input2ReactionIndex.resize(_nameInputIndex.size());
      for(const auto& particle : _nameInputIndex) {
        _input2ReactionIndex[particle.second] = GetReactionIndex(particle.first);
      }
      
      // 3. Resolve Dependencies (Map Names -> Indices)
      for (size_t i = 0; i < GetNCreated(); ++i) {
        RVecIndices vec_indices;
        for(const auto& type_index : _p_stru_depends[i]) {
          Indices_t indices;
          for(const auto& particle : type_index) {
            indices.push_back(GetReactionIndex(particle));
          }
          vec_indices.push_back(indices);
        }
        _p_dep_indices.push_back(vec_indices);
      }
    }

    inline void ParticleCreator::AdoptIndices(const ParticleCreator& master) {
        // 1. Validation
        for(const auto& name : _p_names) {
            if(!master.HasParticle(name)) {
                throw std::runtime_error("Linked Topology Error: Particle '" + name + 
                    "' required by Linked processor is missing in Master.");
            }
        }
        
        // 2. Copy Maps
        _nameIndex = master.GetIndexMap();
        _input2ReactionIndex = master._input2ReactionIndex; 
    }

    inline void ParticleCreator::RebuildReactionMap() {
        // We use the ADOPTED indices (from _nameIndex) to build a NEW map based on LOCAL groups.
        
        // Helper to look up adopted indices
        auto get_indices = [&](const ParticleNames_t& names) {
            Indices_t ret;
            for(const auto& n : names) ret.push_back(_nameIndex.at(n));
            return ret;
        };

        // Helper to resolve groups based on local config (Explicit or Config)
        auto get_group_indices = [&](int order) {
            if (_explicit_groups.count(order)) {
                return get_indices(_explicit_groups[order]);
            }
            if (_config_groups.count(order)) {
                try {
                    return get_indices(_reaction->getGroup(_config_groups[order]));
                } catch (...) {}
            }
            return Indices_t{};
        };

        Indices_t idxBeam     = get_indices(_beam_names);
        Indices_t idxBaryons  = get_group_indices(names::OrderBaryons());
        Indices_t idxMesons   = get_group_indices(names::OrderMesons());
        Indices_t idxScat_ele = get_group_indices(names::OrderScatEle());
        Indices_t idxDeps     = get_indices(GetPriorDependencies());
        Indices_t idxCreate   = get_indices(_p_names);

        RVecIndexMap retIndices{idxBeam, idxBaryons, idxMesons, idxScat_ele, idxDeps, idxCreate};
        
        // Define ONLY the Map (Helper columns and KineIndices are aliased by Processor)
        _mapIndices = retIndices;
        _reaction->Define(GetMapName(), [retIndices]() { return retIndices; }, {});
    }

    inline void ParticleCreator::ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, 
                                               ROOT::RVecD& pz, ROOT::RVecD& m) const 
    {
      for (size_t i = 0; i < GetNCreated(); ++i) {
        _p_creators[i](_nameIndex.at(_p_names[i]), _p_dep_indices[i], px, py, pz, m);
      }
    }

  } // end config
} // end rad


// #pragma once

// #include "ConfigReaction.h"
// #include "StringUtilities.h"
// #include "ParticleCreatorMethods.h"
// #include <map>

// namespace rad {
// namespace config {

//     using StructuredNames_t = std::vector<ParticleNames_t>;
//     using IndexMap_t = std::unordered_map<std::string, int>;
//     using ParticleCreatorFunc_t = void (*)(const Indice_t, const RVecIndices&, ROOT::RVecD&, ROOT::RVecD&, ROOT::RVecD&, ROOT::RVecD&);
    
//     /**
//      * @brief Manages the definition, indexing, and creation of intermediate particles.
//      * * This class is Topology Agnostic. It relies on the controlling Kinematics class
//      * to define which particle groups (Scattered Electrons, etc.) exist in the reaction.
//      */
//     class ParticleCreator {

//     public:
//       ParticleCreator() = default;
//       explicit ParticleCreator(ConfigReaction* cr) : _reaction{cr} {
//           // Set Defaults for standard topologies (can be overridden)
//           _beam_names = {names::BeamIon(), names::BeamEle()};
          
//           // Default Groups (Baryons and Mesons are standard across almost all fixed-target physics)
//           DefineGroup(names::OrderBaryons(), names::Baryons());
//           DefineGroup(names::OrderMesons(), names::Mesons());
          
//           // NOTE: ScatEle is NOT defined by default. 
//           // ElectroKinematics must explicitly register it.
//       }
      
//       // --- Configuration Interface ---

//       /**
//        * @brief Maps a ConfigReaction group name to a fixed ParticleGroupOrder slot.
//        * * Example: DefineGroup(names::OrderScatEle(), "scat_ele")
//        * @param order The fixed slot index (from names::ParticleGroupOrder).
//        * @param groupName The name of the group in ConfigReaction (e.g., "baryons").
//        */
//       void DefineGroup(int order, const std::string& groupName) {
//           _topology_groups[order] = groupName;
//       }

//       /**
//        * @brief Sets the specific names of the beam particles.
//        * Default is {names::BeamIon(), names::BeamEle()}.
//        */
//       void SetBeamNames(const std::vector<std::string>& beams) {
//           _beam_names = beams;
//       }

//       // --- Core Functionality ---

//       void AddParticle(const std::string& name, ParticleCreatorFunc_t func, const StructuredNames_t& depends={{}}) {
//         _p_names.push_back(name);
//         auto flat_depends = utils::flattenColumnNames(depends);
//         _p_required.push_back(flat_depends);
//         _dependencies.insert(flat_depends.begin(), flat_depends.end());
//         _p_stru_depends.push_back(depends);
//         _p_creators.push_back(func);
//         _createIndex[name] = GetNCreated() - 1;
//         _reaction->setParticleIndex(name, constant::InvalidEntry<int>());
//       }

//       void Sum(const std::string& name, const StructuredNames_t& depends={{}}) {
//         AddParticle(name, ParticleCreateBySum, depends);
//       }

//       void Diff(const std::string& name, const StructuredNames_t& depends={{}}) {
//         AddParticle(name, ParticleCreateByDiff, depends);
//       }

//       // --- Indexing & Execution ---

//       void InitMap(const std::vector<std::string>& types);
      
//       void ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, 
//                          ROOT::RVecD& pz, ROOT::RVecD& m) const;

//       // Accessors
//       size_t GetNCreated() const { return _p_names.size(); }
//       Int_t GetReactionIndex(const std::string& name) const { return _nameIndex.at(name); }
//       Int_t GetReactionIndex(size_t input) const { return _input2ReactionIndex[input]; }

//       // Helper
//       std::vector<std::string> GetPriorDependencies();
//       Indices_t CreateIndices(IndexMap_t& nameIndex, const ParticleNames_t& names, size_t& idx);

//     private:
//       ConfigReaction* _reaction = nullptr;
      
//       // Configuration
//       std::vector<std::string> _beam_names;
//       std::map<int, std::string> _topology_groups; // Maps Order ID -> Group Name

//       // Internal State
//       ParticleNames_t _p_names;
//       std::vector<ParticleCreatorFunc_t> _p_creators;
//       std::vector<ParticleNames_t> _p_required;
//       std::vector<StructuredNames_t> _p_stru_depends;
//       ROOT::RVec<RVecIndices> _p_dep_indices;
//       std::set<std::string> _dependencies;
      
//       ParticleNames_t _inputNames;
//       RVecIndexMap _mapIndices;
//       IndexMap_t _nameIndex;
//       IndexMap_t _nameInputIndex;
//       IndexMap_t _createIndex;
//       Indices_t _input2ReactionIndex;
//     };

//     // ========================================================================
//     // Implementation
//     // ========================================================================

//     inline std::vector<std::string> ParticleCreator::GetPriorDependencies() {
//         std::vector<std::string> vec_deps(_dependencies.begin(), _dependencies.end());
//         utils::removeExistingStrings(vec_deps, _p_names);
//         return vec_deps;
//     }

//     inline Indices_t ParticleCreator::CreateIndices(IndexMap_t& nameIndex, const ParticleNames_t& names, size_t& idx) {
//         Indices_t indices{};
//         for(const auto& name : names) {
//           if(nameIndex.count(name) == 0) {
//             indices.push_back(idx);
//             nameIndex[name] = idx;
//             ++idx;
//           }
//         }
//         return indices;
//     }

//     inline void ParticleCreator::InitMap(const std::vector<std::string>& types) {
      
//       // 1. Resolve Configured Groups
//       ParticleNames_t all_group_particles;
      
//       // Add Beams
//       all_group_particles.insert(all_group_particles.end(), _beam_names.begin(), _beam_names.end());

//       // Add Registered Groups (Baryons, Mesons, ScatEle if defined)
//       cout<<"Add Registered Groups number = "<<_topology_groups.size()<<endl;
//       for(auto const& [order, groupName] : _topology_groups) {
// 	  cout<<order<<" "<<groupName <<endl;
// 	  auto particles = _reaction->getGroup(groupName);
//           all_group_particles.insert(all_group_particles.end(), particles.begin(), particles.end());
//      }

//       // Add Dependencies
//       auto dep_names = GetPriorDependencies();
//       utils::removeExistingStrings(dep_names, all_group_particles);
      
//       // 2. Define Master Input List
//       _inputNames = all_group_particles;
//       _inputNames.insert(_inputNames.end(), dep_names.begin(), dep_names.end());
//       utils::removeExistingStrings(_inputNames, _p_names); // Safety check

//       // Define KineIndices in RDF
//       cout<<"Define KineIndices in RDF"<<endl;
//      for(const auto& type : types) {
//         auto typeNames = utils::prependToAll(_inputNames, type);
//         _reaction->setGroupParticles(type + names::KineIndices(), typeNames);
// 	cout<<type + names::KineIndices()<<" "<<typeNames.size()<<endl;
//       }

//       // Store input map
//       size_t in_idx = 0;
//       CreateIndices(_nameInputIndex, _inputNames, in_idx);
       
//       // 3. Calculate Fixed Positions (ReactionMap)
//       // We must construct the RVecIndexMap in the strict order of names::ParticleGroupOrder
//       // enum: {Beams, Baryons, Mesons, ScatEle, Deps, Createds}
      
//       size_t idx = 0;
      
//       // A. Beams
//       Indices_t idxBeam = CreateIndices(_nameIndex, _beam_names, idx);
      
//       // B. Dynamic Groups (Baryons, Mesons, ScatEle)
//       // We assume the user has registered these in the correct enum order, 
//       // or we explicitly fetch them by enum ID.
//       Indices_t idxBaryons = CreateIndices(_nameIndex, _reaction->getGroup(_topology_groups[names::OrderBaryons()]), idx);
//       Indices_t idxMesons = CreateIndices(_nameIndex, _reaction->getGroup(_topology_groups[names::OrderMesons()]), idx);
      
//       // C. Scattered Electron (Optional)
//       Indices_t idxScat_ele;
//       if(_topology_groups.count(names::OrderScatEle())) {
//           idxScat_ele = CreateIndices(_nameIndex, _reaction->getGroup(_topology_groups[names::OrderScatEle()]), idx);
//       }
//       // If not defined, idxScat_ele remains empty {} which is correct.

//       // D. Dependencies & Created
//       Indices_t idxDeps = CreateIndices(_nameIndex, dep_names, idx);
//       Indices_t idxCreate = CreateIndices(_nameIndex, _p_names, idx);

//       // Build Map
//       RVecIndexMap retIndices{idxBeam, idxBaryons, idxMesons, idxScat_ele, idxDeps, idxCreate};
      
//       _mapIndices = retIndices;
//       _reaction->Define(names::ReactionMap(), [retIndices]() { return retIndices; }, {});
      
//       // Helper Columns
//       cout<<" Helper Columns "<<endl;
//       for(const auto& particle : _nameIndex) {
//         _reaction->Define(names::data_type::Kine() + particle.first, [idx = particle.second](){ return idx; }, {});
// 	cout<< names::data_type::Kine() + particle.first <<" "<<idx<<endl;
//       }

//       // Conversion Map
//       _input2ReactionIndex.resize(_nameInputIndex.size());
//       for(const auto& particle : _nameInputIndex) {
//         _input2ReactionIndex[particle.second] = GetReactionIndex(particle.first);
//       }
      
//       // 4. Resolve Dependencies
//       for (size_t i = 0; i < GetNCreated(); ++i) {
//         RVecIndices vec_indices;
//         for(const auto& type_index : _p_stru_depends[i]) {
//           Indices_t indices;
//           for(const auto& particle : type_index) {
//             indices.push_back(GetReactionIndex(particle));
//           }
//           vec_indices.push_back(indices);
//         }
//         _p_dep_indices.push_back(vec_indices);
//       }
//     }

//     inline void ParticleCreator::ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, ROOT::RVecD& pz, ROOT::RVecD& m) const {
//       for (size_t i = 0; i < GetNCreated(); ++i) {
//         _p_creators[i](_nameIndex.at(_p_names[i]), _p_dep_indices[i], px, py, pz, m);
//       }
//     }

//   } // end config
// } // end rad
