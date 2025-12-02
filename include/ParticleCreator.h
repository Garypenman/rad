#pragma once

#include "ConfigReaction.h"
#include "StringUtilities.h"
#include "ParticleCreatorMethods.h"
#include <map>

namespace rad {
namespace config {

    using StructuredNames_t = std::vector<ParticleNames_t>;
    using IndexMap_t = std::unordered_map<std::string, int>;
    using ParticleCreatorFunc_t = void (*)(const Indice_t, const RVecIndices&, ROOT::RVecD&, ROOT::RVecD&, ROOT::RVecD&, ROOT::RVecD&);
    
    /**
     * @brief Manages the definition, indexing, and creation of intermediate particles.
     * * This class is Topology Agnostic. It relies on the controlling Kinematics class
     * to define which particle groups (Scattered Electrons, etc.) exist in the reaction.
     */
    class ParticleCreator {

    public:
      ParticleCreator() = default;
      explicit ParticleCreator(ConfigReaction* cr) : _reaction{cr} {
          // Set Defaults for standard topologies (can be overridden)
          _beam_names = {names::BeamIon(), names::BeamEle()};
          
          // Default Groups (Baryons and Mesons are standard)
          DefineGroup(names::OrderBaryons(), names::Baryons());
          DefineGroup(names::OrderMesons(), names::Mesons());
          
          // NOTE: ScatEle is NOT defined by default. 
          // ElectroKinematics must explicitly register it.
      }
      
      // --- Configuration Interface ---

      /**
       * @brief Maps a ConfigReaction group name to a fixed ParticleGroupOrder slot.
       * * Example: DefineGroup(names::OrderScatEle(), "scat_ele")
       * @param order The fixed slot index (from names::ParticleGroupOrder).
       * @param groupName The name of the group in ConfigReaction (e.g., "baryons").
       */
      void DefineGroup(int order, const std::string& groupName) {
          _topology_groups[order] = groupName;
      }

      /**
       * @brief Sets the specific names of the beam particles.
       * Default is {names::BeamIon(), names::BeamEle()}.
       */
      void SetBeamNames(const std::vector<std::string>& beams) {
          _beam_names = beams;
      }

      // --- Core Functionality ---

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

      // --- Indexing & Execution ---

      void InitMap(const std::vector<std::string>& types);
      
      void ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, 
                         ROOT::RVecD& pz, ROOT::RVecD& m) const;

      // Accessors
      size_t GetNCreated() const { return _p_names.size(); }
      Int_t GetReactionIndex(const std::string& name) const { return _nameIndex.at(name); }
      Int_t GetReactionIndex(size_t input) const { return _input2ReactionIndex[input]; }

      // Helper
      std::vector<std::string> GetPriorDependencies();
      Indices_t CreateIndices(IndexMap_t& nameIndex, const ParticleNames_t& names, size_t& idx);

    private:
      ConfigReaction* _reaction = nullptr;
      
      // Configuration
      std::vector<std::string> _beam_names;
      std::map<int, std::string> _topology_groups; // Maps Order ID -> Group Name

      // Internal State
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
          if(nameIndex.count(name) == 0) {
            indices.push_back(idx);
            nameIndex[name] = idx;
            ++idx;
          }
        }
        return indices;
    }

    inline void ParticleCreator::InitMap(const std::vector<std::string>& types) {
      
      // 1. Resolve Configured Groups
      ParticleNames_t all_group_particles;
      
      // Add Beams
      all_group_particles.insert(all_group_particles.end(), _beam_names.begin(), _beam_names.end());

      // Add Registered Groups (Baryons, Mesons, ScatEle if defined)
      cout<<"Add Registered Groups number = "<<_topology_groups.size()<<endl;
      for(auto const& [order, groupName] : _topology_groups) {
	  cout<<order<<" "<<groupName <<endl;
	  auto particles = _reaction->getGroup(groupName);
          all_group_particles.insert(all_group_particles.end(), particles.begin(), particles.end());
     }

      // Add Dependencies
      auto dep_names = GetPriorDependencies();
      utils::removeExistingStrings(dep_names, all_group_particles);
      
      // 2. Define Master Input List
      _inputNames = all_group_particles;
      _inputNames.insert(_inputNames.end(), dep_names.begin(), dep_names.end());
      utils::removeExistingStrings(_inputNames, _p_names); // Safety check

      // Define KineIndices in RDF
      cout<<"Define KineIndices in RDF"<<endl;
     for(const auto& type : types) {
        auto typeNames = utils::prependToAll(_inputNames, type);
        _reaction->setGroupParticles(type + names::KineIndices(), typeNames);
	cout<<type + names::KineIndices()<<" "<<typeNames.size()<<endl;
      }

      // Store input map
      size_t in_idx = 0;
      CreateIndices(_nameInputIndex, _inputNames, in_idx);
       
      // 3. Calculate Fixed Positions (ReactionMap)
      // We must construct the RVecIndexMap in the strict order of names::ParticleGroupOrder
      // enum: {Beams, Baryons, Mesons, ScatEle, Deps, Createds}
      
      size_t idx = 0;
      
      // A. Beams
      Indices_t idxBeam = CreateIndices(_nameIndex, _beam_names, idx);
      
      // B. Dynamic Groups (Baryons, Mesons, ScatEle)
      // We assume the user has registered these in the correct enum order, 
      // or we explicitly fetch them by enum ID.
      Indices_t idxBaryons = CreateIndices(_nameIndex, _reaction->getGroup(_topology_groups[names::OrderBaryons()]), idx);
      Indices_t idxMesons = CreateIndices(_nameIndex, _reaction->getGroup(_topology_groups[names::OrderMesons()]), idx);
      
      // C. Scattered Electron (Optional)
      Indices_t idxScat_ele;
      if(_topology_groups.count(names::OrderScatEle())) {
          idxScat_ele = CreateIndices(_nameIndex, _reaction->getGroup(_topology_groups[names::OrderScatEle()]), idx);
      }
      // If not defined, idxScat_ele remains empty {} which is correct.

      // D. Dependencies & Created
      Indices_t idxDeps = CreateIndices(_nameIndex, dep_names, idx);
      Indices_t idxCreate = CreateIndices(_nameIndex, _p_names, idx);

      // Build Map
      RVecIndexMap retIndices{idxBeam, idxBaryons, idxMesons, idxScat_ele, idxDeps, idxCreate};
      
      _mapIndices = retIndices;
      _reaction->Define(names::ReactionMap(), [retIndices]() { return retIndices; }, {});
      
      // Helper Columns
      cout<<" Helper Columns "<<endl;
      for(const auto& particle : _nameIndex) {
        _reaction->Define(names::data_type::Kine() + particle.first, [idx = particle.second](){ return idx; }, {});
	cout<< names::data_type::Kine() + particle.first <<" "<<idx<<endl;
      }

      // Conversion Map
      _input2ReactionIndex.resize(_nameInputIndex.size());
      for(const auto& particle : _nameInputIndex) {
        _input2ReactionIndex[particle.second] = GetReactionIndex(particle.first);
      }
      
      // 4. Resolve Dependencies
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

    inline void ParticleCreator::ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, ROOT::RVecD& pz, ROOT::RVecD& m) const {
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

// namespace rad{
//   namespace config{
//     using StructuredNames_t = std::vector<ParticleNames_t>;
//     using IndexMap_t = std::unordered_map<std::string,int>;
    
//     /**
//      * @brief Type alias for the concrete, non-templated particle creation function.
//      * * All creator functions must match this signature for storage in std::vector.
//      */
//     using ParticleCreatorFunc_t = void (*)(
// 					   const Indice_t position, //position to add at
// 					   const RVecIndices&,  // indices
// 					   ROOT::RVecD&,   // px
// 					   ROOT::RVecD&,   // py
// 					   ROOT::RVecD&,   // pz
// 					   ROOT::RVecD&    // m
// 					   );
    
//     /////////////////////////////////////////////////////////////////////////////////
//     class ParticleCreator {

//     public:
      
//       ParticleCreator() = default;
//       ParticleCreator(ConfigReaction* cr):_reaction{cr}{};
      
    
//       void SetReaction(ConfigReaction* reaction){_reaction=reaction;}
//       ConfigReaction* Reaction() const { return _reaction;}

   
//       void AddParticle(const std::string& name,ParticleCreatorFunc_t func,const StructuredNames_t& depends={{}}){
// 	//all _p_ datamembers must be in synch
// 	_p_names.push_back(name);
// 	//keep particle required to creat this one
// 	auto flat_depends = utils::flattenColumnNames(depends);
// 	_p_required.push_back(flat_depends);
// 	//add to unique dependencies container
// 	_dependencies.insert(flat_depends.begin(),flat_depends.end());
// 	//save origin depends to recreate the indice structure
// 	_p_stru_depends.push_back(depends);
// 	//save the function
// 	_p_creators.push_back(func);
// 	//save the index
// 	_createIndex[name]=GetNCreated()-1;
// 	//tell reaction about this particle
// 	_reaction->setParticleIndex(name,constant::InvalidEntry<int>());
//       }
//       size_t GetNCreated() const {return _p_names.size();}

//       std::vector<std::string> GetPriorDependencies(){

// 	//convert set to vector
// 	std::vector<std::string> vec_deps(_dependencies.begin(), _dependencies.end());
// 	//get dependencies which are created particles
// 	auto depCreated = utils::getCommonStrings(_p_names,vec_deps);
	
// 	//must remove deps that are themselves created particles from vec_deps
// 	utils::removeExistingStrings(vec_deps,_p_names);
	
// 	return vec_deps;
//       }
      
//       Int_t GetInputIndex(const std::string& name) const{return _nameInputIndex.at(name);}
//       Int_t GetReactionIndex(const std::string& name) const {return _nameIndex.at(name);}
//       Int_t GetReactionIndex(size_t input) const {return _input2ReactionIndex[input];}
      
//       ParticleCreatorFunc_t GetMethod(const std::string& name) const{ return _p_creators[_createIndex.at(name)];}
//       ///given list of names create indices for them
//       ///and map them to the particle name
//       Indices_t CreateIndices(IndexMap_t& nameIndex ,const ParticleNames_t& names, size_t& idx){
       
// 	Indices_t indices{};
// 	for(const auto& name:names){
// 	  if(nameIndex.count(name) > 0 ){
// 	    std::cout<<" Warning : ParticleCreator::SortIndices I have multiple "<<name<<" for index map"<<endl;
// 	  }
// 	  else{
// 	    indices.push_back(idx);
// 	    nameIndex[name] = idx;
// 	    ++idx;
// 	  }
// 	}
// 	return indices;
//       }
      
//       void InitMap(const std::vector<std::string>& types);
//       void ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, 
// 			 ROOT::RVecD& pz, ROOT::RVecD& m) const ;

//       //User Shortcuts
//       //e.g. Creator().Sum("interm",{{"p1","p2"}})  interm = p1 + p2
//       void Sum(const std::string& name,const  StructuredNames_t& depends={{}}){
// 	AddParticle(name,ParticleCreateBySum,depends);
//       }
//       //e.g. Creator().Sum("pdiff",{{"pplus1","pplus2"},{"pminus1"}})
//       //=> pdiff = pplus1+ pplus2 - pminus
//       void Diff(const std::string& name,const  StructuredNames_t& depends={{}}){
// 	AddParticle(name,ParticleCreateByDiff,depends);
//       }
       
//     private:
      
//       ConfigReaction* _reaction=nullptr;
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

  
//     void ParticleCreator::InitMap(const std::vector<std::string>& types){
      
//       // --- 1. Create initial kinematic groups ---
//       auto beam_names = {names::BeamIon(), names::BeamEle()};
//       auto baryon_names = _reaction->getGroup(names::Baryons());
//       auto meson_names = _reaction->getGroup(names::Mesons());
//       auto scat_names = {names::ScatEle()}; 
//       auto dep_names=GetPriorDependencies();
//       //can remove particles already existing in other groups from input dep_names
//       utils::removeExistingStrings(dep_names,
// 				   utils::concatenateStringVectors(beam_names, baryon_names, meson_names, scat_names));
//       // --- 2. Define the input particles already Defined ---
//       // Consolidate all particle names into one list to define the order of the final momentum arrays.
//       _inputNames = utils::concatenateStringVectors(beam_names, baryon_names, meson_names, scat_names, dep_names);
      
//       //remove created particles from list of inputs
//       utils::removeExistingStrings(_inputNames,_p_names);

//       //Need to define a set of KineIndices for each type
//       for(const auto& type:types){
// 	//append type to indices
// 	auto typeNames = utils::prependToAll(_inputNames, type);
// 	//Define input indices in RDF
// 	_reaction->setGroupParticles(type+names::KineIndices(), typeNames); // Defines the ordered column
// 	cout<<" input indices "<< utils::combineVectorToQuotedString(typeNames) <<endl;
//       }

//       //store index map for input particles
//       size_t in_idx=0;
//       CreateIndices(_nameInputIndex,_inputNames,in_idx);
       
//       // --- 3. Calculate Fixed Positions for ReactionMap ---
//       size_t idx = 0;
//       // The position indices are calculated sequentially based on group sizes.
//       Indices_t idxBeam=CreateIndices(_nameIndex,beam_names,idx);
//       Indices_t idxBaryons =CreateIndices(_nameIndex,baryon_names,idx);
//       Indices_t idxMesons = CreateIndices(_nameIndex,meson_names,idx);            
//       Indices_t idxScat_ele = CreateIndices(_nameIndex,scat_names,idx);

//       //now indices for particle which are required by created particles
//       //this makes sure their indices are passed to KinematicProcessor and 4-vectors constructed
//       Indices_t idxDeps = CreateIndices(_nameIndex,dep_names,idx);

//       //Now create the indices for the created particles
//       //in cases they exist in existing groups, they will keep that index
//       //only particles not included in prior groups will get a new index
//       Indices_t idxCreate = CreateIndices(_nameIndex,_p_names,idx);

//       // Stores the fixed array index (0, 1, 2...) for each particle role in the consolidated array.
//       // same ordering as reactionNames, ParticleGroupOrder
//       // enum class names::ParticleGroupOrder{BeamIon, BeamEle, Baryons, Mesons, ScatEle, VirtGamma, Deps};
//       RVecIndexMap retIndices{idxBeam, idxBaryons, idxMesons,
// 			      idxScat_ele, idxDeps, idxCreate};
//       //cout<<" RVecIndexMap retIndices "<<retIndices<<endl;
//       //retIndices now specifies the fixed reaction map
//       //Define it in RDF
//       _mapIndices=retIndices;
//       _reaction->Define(names::ReactionMap(),
// 			[retIndices]() {
// 			  return retIndices;
// 			},
// 			{});
      
//       //Define index for each particle in components array
//       for(const auto& particle:_nameIndex){
// 	auto pindex  = particle.second;
// 	//dont use setParticleIndex as want integer indices not Indices_t
// 	//as kinematic components are in fixed order
// 	_reaction->Define(names::data_type::Kine()+particle.first,[pindex](){return pindex;},{});
//       }
//       //create conversion indices from input to reaction
//       _input2ReactionIndex.resize(_nameInputIndex.size());
//       for(const auto& particle: _nameInputIndex){
// 	cout<<"input2ReactionIndex "<<particle.first<<" "<<particle.second<<" "<<GetReactionIndex(particle.first)<<endl;
// 	_input2ReactionIndex[particle.second]=GetReactionIndex(particle.first);
//       }
      
//       //Finally convert dependency names to RVecIndices of
//       //different indice types for each created particle
//       for (size_t i = 0; i < GetNCreated(); ++i) {
// 	RVecIndices vec_indices;
// 	//get different types from structured dependency names vector
// 	for(const auto& type_index:_p_stru_depends[i]){
// 	  Indices_t indices;
// 	  for(const auto& particle:type_index){
// 	    indices.push_back(GetReactionIndex(particle));
// 	  }
// 	  vec_indices.push_back(indices);
// 	}
// 	_p_dep_indices.push_back(vec_indices);

//       }

//     }
//     /**
//      * @brief Executes all registered particle creation methods for the current combination.
//      * * It iterates through the stored creation functions and applies the calculated 
//      * momentum to the corresponding position in the component arrays.
//      * * @param px, py, pz, m References to the component RVecs (modified in place).
//      */
//     void ParticleCreator::ApplyCreation(ROOT::RVecD& px, ROOT::RVecD& py, 
// 		       ROOT::RVecD& pz, ROOT::RVecD& m) const {
    
//       // Loop through all registered creation methods and their dependencies.
//       for (size_t i = 0; i < GetNCreated(); ++i) {
        
//         // 1. Get the creation function pointer.
//         ParticleCreatorFunc_t func = _p_creators[i];
        
//         // 2. Get the target position (fixed array index) for this created particle.
//         // The fixed index is found by looking up the particle name in the ReactionMap's structure.
//         const std::string& created_name = _p_names[i];
//         const int position = _nameIndex.at(created_name); // Get fixed position of this particle in Reaction
//          // 3. Get the input particle names/roles required by this specific creator.
//         const auto& indices= _p_dep_indices[i];

// 	//	cout<<"ParticleCreator::ApplyCreation "<<position<<" "<<indices<<" "<<pz<<m<<endl;
// 	// 5. Execute the creation function.
//         // NOTE: The function requires the single, fixed position and the resolved indices.
//         func(position, indices, px, py, pz, m);
//       }
//     }

  
  
//   }//end config
// }
