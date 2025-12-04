#pragma once

// Assuming RDFInterface.h contains the base class definition
#include "RDFInterface.h" 
#include "CommonDefines.h" 
#include "ConfigUtils.h" 
#include "Combinatorics.h" 
#include "CombinatorialUtil.h" 

#include "Constants.h" // For InvalidEntry
#include "DefineNames.h"
#include "RVecHelpers.h"
#include "ReactionUtilities.h"
#include "StringUtilities.h"
#include "Random.h"

// Forward declarations for external combinatorial helpers
namespace rad {

  
  namespace combinatorics {
    // Assuming this is the C++ kernel that takes RVec<RVecI>... and returns RVec<RVecIndexMap>
    ROOT::RVec<RVecIndexMap> GenerateAllCombinations(...); 
  }
  namespace util {
    // The generic wrapper from Refactor.md
    template <typename F, typename... Args>
    auto ApplyCombinationsGeneric(F&& singleComboFunc, const ROOT::RVec<RVecIndexMap>& combos, Args&&... args);
  }
}

namespace rad {
  namespace config {

    class ConfigReaction : public RDFInterface {

    public:
      //------------------ Constructors ------------------

      // Calls the appropriate base class constructor
      ConfigReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns)
        : RDFInterface(treeName, fileNameGlob, columns) {}

      ConfigReaction(const std::string_view treeName, const std::vector<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns)
        : RDFInterface(treeName, filenames, columns) {}

      ConfigReaction(ROOT::RDataFrame rdf)
        : RDFInterface(rdf) {}

      //------------------ RDFInterface Overrides (Pure Virtuals) ------------------

  
      void Snapshot(const string& filename) override {
        try {
          RDFstep final_df = CurrFrame();
          
          // Flatten the frame if in combinatorial mode (Stage III Logic)
          if (_isCombinatorialMode) {
              // NOTE: The actual Unroll target/logic must be implemented here.
              // For now, we assume a representative column, but this needs refinement.
              // Example: final_df = final_df.Unroll("rec_" + _firstDefinedCombiColumn); 
          }

          auto cols = final_df.GetDefinedColumnNames();
          RemoveSnapshotColumns(cols);
          final_df.Snapshot("rad_tree", filename, cols);
        } catch (const std::exception& ex) {
          std::cerr << "Snapshot failed: " << ex.what() << std::endl;
          throw;
        }
      }

      void BookLazySnapshot(const string& filename) override {
        try {
          RDFstep final_df = CurrFrame();

          // Flatten the frame if in combinatorial mode (Stage III Logic)
          if (_isCombinatorialMode) {
              // NOTE: Implementation of the Unroll logic goes here.
          }
          
          ROOT::RDF::RSnapshotOptions opts;
          opts.fLazy = true;
          auto cols = final_df.GetDefinedColumnNames();
          RemoveSnapshotColumns(cols);
          auto snapshot_result = final_df.Snapshot("rad_tree", filename, cols, opts);
          _triggerSnapshots.emplace_back([snapshot = std::move(snapshot_result)]() mutable {});
        } catch (const std::exception& ex) {
          std::cerr << "BookLazySnapshot failed: " << ex.what() << std::endl;
          throw;
        }
      }

      // Hook to remove internal columns
      void RemoveSnapshotColumns(std::vector<string>& cols) override {
        // Remove internal columns specific to ConfigReaction (e.g., the map and alias data)
        cols.erase(std::remove(cols.begin(), cols.end(), names::ReactionMap()), cols.end());

        // Remove any columns with the DoNotWriteTag
        auto tag = DoNotWriteTag();
        cols.erase(std::remove_if(cols.begin(), cols.end(),
				  [&tag](const string& col) -> bool { return col.find(tag) != std::string::npos; }),cols.end());

        // Call base class cleanup (though RDFInterface currently has no columns to remove)
        RDFInterface::RemoveSnapshotColumns(cols);
      }

      //------------------ Combinatorial API ------------------

      /**
       * @brief Define a collection of particle candidates using a string expression.
       */
      void setParticleCandidatesExpr(const string& name, const string& expression) {
	//void setParticleIndex(const string& name, const string& expression) {
        if (_candidateExpressions.count(name) || _lambdaCandidateDependencies.count(name)) {
          throw std::invalid_argument("Candidate '" + name + "' already defined.");
        }
        _candidateExpressions[name] = expression;
	AddParticleName(name);
	AddFinalParticleName(name);
      }

      /**
       * @brief Define a collection of particle candidates using a C++ function (lambda/functor).
       */
      template<typename Lambda>
      void setParticleCandidates(const string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
	//void setParticleIndex(const string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
        if (_candidateExpressions.count(name) || _lambdaCandidateDependencies.count(name)) {
          throw std::invalid_argument("Candidate '" + name + "' already defined.");
        }

        // Define the internal column using the Lambda right away
        //string indices_name = DoNotWriteTag() + name + "_combi";
        Define(name, std::forward<Lambda>(func), columns);

        // Store the dependency names for tracking
        _lambdaCandidateDependencies[name] = columns; 

	AddParticleName(name);
	AddFinalParticleName(name);
      }
      
      /** 
       * Set constant index in collection for particle
       * This assumes constant position in collection (e.g in some HepMC3 files)
       * and update the current frame to the aliased one
       */
      void setParticleIndex(const string& name, const int idx, int pdg=0 ){
	// Check if particle index column already exists to avoid accidental overwrite
	if (ColumnExists(name, CurrFrame())) {
	  throw std::invalid_argument("setParticleIndex: Index column for particle '" + name + "' already exists!");
	}
	//call combi version with single combination
	setParticleCandidates(name,[idx](){return RVecI{idx};},{});
	AddParticleName(name);
 	AddFinalParticleName(name);
      }
      /** 
       * Set constant index in collection for particle
       * This assumes constant position in collection (e.g in some HepMC3 files)
       * and update the current frame to the aliased one
       */
      void setParticleCandidates(const string& name, const Indices_t idx  ){
	// Check if particle index column already exists to avoid accidental overwrite
	if (ColumnExists(name, CurrFrame())) {
	  throw std::invalid_argument("setParticleIndex: Index column for particle '" + name + "' already exists!");
	}
	//call combi version with single combination
	setParticleCandidates(name,[idx](){return idx;},{});
	AddParticleName(name);
 	AddFinalParticleName(name);
      }
      /**
       * @brief Groups particle candidate indices for the Meson component of the final state.
       * * * This function consolidates the candidate index lists (Indices_t) defined for individual 
       * final-state particles into a single container representing the Meson group candidates.
       * * * If the input list is empty, a column is still defined, containing only the 
       * invalid index as a placeholder.
       * * @param particles List of column names (Indices_t) representing the particles 
       * belonging to the Meson component (e.g., {"pi+", "pi-"}).
       * @see setGroupParticles
       */
      void setMesonParticles(const ROOT::RDF::ColumnNames_t& particles) {
	if (particles.empty() == true) {
	  std::cout << "setMesonParticles " << as_string(names::Mesons()) << " no particles, " << "give default index " << std::endl;
	  // The default return type must match the output of setGroupParticles (RVecIndices)
	  Define(as_string(names::Mesons()), [](){return RVecIndices{{constant::InvalidIndex()}}; }, {}); 
	  return;
	}
	setGroupParticles(as_string(names::Mesons()), particles);
      }

      //-----------------------------------------------------------------------------

      /**
       * @brief Groups particle candidate indices for the Baryon component of the final state.
       * * * This function consolidates the candidate index lists (Indices_t) defined for individual 
       * final-state particles into a single container representing the Baryon group candidates.
       * * * If the input list is empty, a column is still defined, containing only the 
       * invalid index as a placeholder.
       * * @param particles List of column names (Indices_t) representing the particles 
       * belonging to the Baryon component (e.g., {"p", "n"}).
       * @see setGroupParticles
       */
      void setBaryonParticles(const ROOT::RDF::ColumnNames_t& particles) {
	if (particles.empty() == true) {
	  std::cout << "setBaryonParticles " << as_string(names::Baryons()) << " no particles, " << "give default index " << std::endl;
	  // The default return type must match the output of setGroupParticles (RVecIndices)
	  Define(as_string(names::Baryons()), [](){return RVecIndices{{constant::InvalidIndex()}}; }, {}); 
	  return;
	}
	setGroupParticles(as_string(names::Baryons()), particles);
      }

      //-----------------------------------------------------------------------------

      /**
       * @brief Defines a new RDataFrame column by aggregating multiple existing particle index vectors.
       * * * This function uses the generic rad::helpers::Group template to take multiple 
       * RVecI (Indices_t) columns and group them into a single RVecIndices column.
       * * @param name The name of the new column to be defined (e.g., "Mesons" or "Baryons").
       * @param particles List of existing column names (Indices_t) to be grouped.
       * @note The resulting column type is RVecIndices (ROOT::RVec<Indices_t>).
       */
      void setGroupParticles(const string& name, const ROOT::RDF::ColumnNames_t& particles) {
	
	auto pstring = reaction::util::ColumnsToString(particles); // Creates string like "{p1,p2,p3,p4,...}"
	pstring = pstring.substr(1, pstring.size() - 2); // Removes outer braces to get: "p1,p2,p3,p4,..."
	
	// Defines the new column using the Group<Indices_t> template, which returns RVec<Indices_t>
	Define(name, Form("rad::helpers::Group<rad::Indices_t>(%s)", pstring.data()));

	_groupMap[name]=particles;
	
	return;
      }
      
      
      /**
       * @brief Generates all possible combinations of the defined particle candidates.
       */
      void makeCombinations() {
        if (_candidateExpressions.empty() && _lambdaCandidateDependencies.empty()) {
          throw std::runtime_error("makeCombinations: No particle candidates defined.");
        }
        _isCombinatorialMode = true;

        ROOT::RDF::ColumnNames_t candidateCols;

        // 1. Process String-Defined Candidates: Define and collect name
        for (const auto& pair : _candidateExpressions) {
          const string& name = pair.first;
          string indices_name = name;// + "_combi";
          Define(indices_name, pair.second);
          candidateCols.push_back(indices_name);
        }

        // 2. Process Lambda-Defined Candidates: Just collect the name (already defined)
        for (const auto& pair : _lambdaCandidateDependencies) {
          const string& name = pair.first;
          string indices_name = name;// + "_combi";
          candidateCols.push_back(indices_name);
        }

        // 3. Generate the RVecCombis column named ReactionCombos(), each element is a set particle
        Define(names::ReactionCombos(),
	       utils::createFunctionCallStringFromVec("rad::combinatorics::GenerateAllCombinations",
						      {rad::reaction::util::ColumnsToString(candidateCols)}));
	       //  Form("rad::combinatorics::GenerateAllCombinations(%s)", rad::reaction::util::ColumnsToString(candidateCols).data()));

	// 4. Split ReactionCombos into individal particle vectors, reuse particle indice names
	// generate combitorial indices for each data type
	for(auto &atype:_type_comps){
	  for(size_t ip=0;ip<candidateCols.size();++ip){
	    Define(atype.first + candidateCols[ip],[ip](const RVecIndices& part_combos){return part_combos[ip];},{names::ReactionCombos().data()});
	  }
	}
	for(size_t ip=0;ip<candidateCols.size();++ip){
	  Redefine(candidateCols[ip],[ip](const RVecIndices& part_combos){return part_combos[ip];},{names::ReactionCombos().data()});
	}
	
	
      }


      //------------------ Custom Define Overrides ------------------
      /**
       * Call Define for predefined data types, pprepend type string
       */
      void DefineForAllTypes(const string& name,const string& expression){
	for(auto &atype:_type_comps){
	  if (atype.second.find(names::P4Components()) == atype.second.end() ||
	      atype.second.find(names::P3Components()) == atype.second.end()) {
            throw std::runtime_error("DefineForAllTypes: Missing 'components_p4' or 'components_p3' for type: " + atype.first);
	  }
	  TString type_expr = expression.data();
	  type_expr.ReplaceAll(names::P4Components(),atype.second[names::P4Components()]);
	  type_expr.ReplaceAll(names::P3Components(),atype.second[names::P3Components()]);
	  Define(atype.first + name.data(),type_expr.Data());
	}
      }
      /**
       * @brief Defines a new column in the RDataFrame by applying a vectorized combinatorial function
       * for every registered physics object type (e.g., "rec_", "tru_").
       *
       * This function acts as an interface to RDataFrame::Define() and implements the combinatorial
       * vectorization logic using `rad::util::ApplyCombinations`. It automatically handles the
       * template instantiation of the core single-combination function based on the
       * component types (float, double, etc.) available for each object prefix.
       *
       * @note This method assumes the existence of `_type_comps` which maps a type prefix (e.g., "rec_")
       * to a map of component names and their types/values. It specifically requires
       * `names::P4Components()` and `names::P3Components()` to be defined for type substitution.
       *
       * @param name The base name for the new column (e.g., "Mass"), which will be prefixed by the
       * object type (e.g., "rec_Mass", "tru_Mass").
       * @param sfunc The string name of the single-combination function to be applied (e.g., "MyMassFunc").
       * This function will be automatically templated based on the underlying object types.
       * @param indices The name of the column containing the combination indices (e.g., "combos" or user-defined).
       * @param arguments A string containing the arguments passed to the combination function.
       * It must use placeholders for the component names (e.g., "components_p4" or "components_p3")
       * which will be automatically substituted with the actual column names
       * (e.g., "Pt,Eta,Phi,Mass") registered for the current type.
       *
       * @exception std::runtime_error if either "components_p4" or "components_p3" are missing for a type
       * in the internal `_type_comps` structure, as these are required for type substitution.
       *
       * @par Generated Expression Structure:
       * @code
       * rad::util::ApplyCombinations(sfunc<obj_types>, indices, arguments_with_substitution)
       * @endcode
       */
      void DefineForAllTypes(const string& name, const string& sfunc, const string& indices,const string& arguments ) {
        // NOTE: A column naming convention is needed to reliably track the first defined column
        // for the final Snapshot::Unroll() call.

	
        for (auto& atype : _type_comps) {
          if (atype.second.find(names::P4Components()) == atype.second.end() ||
            atype.second.find(names::P3Components()) == atype.second.end()) {
            throw std::runtime_error("DefineForAllTypes: Missing 'components_p4' or 'components_p3' for type: " + atype.first);
          }
	/// We have to account for different types (rec_, tru_,..)
	/// having different object types (float, double,...)
	/// we can get the type of a Column from the dataframe
	/// for arguments which are 4vector or 3vectors assume
	/// the 3-momentum components are same type, but mass can be different 
	  std::cout<<"DefineForAllTypes " << atype.first <<" "<<atype.second.size()<<" "<<atype.second[names::P4Components()]<<" "<<atype.second[names::P3Components()]<<" "<<(atype.second.find(names::P4Components())== atype.second.end())<<std::endl;
	  
	  TString args = arguments.data();
	  string obj_types;

	  
	  if(args.Contains(names::P4Components())){
	    obj_types =TypeComponentsTypeString(atype.first,names::P4Components());
	    args.ReplaceAll(names::P4Components(), atype.second[names::P4Components()]);
	  }
	  else if(args.Contains(names::P3Components())){
	    obj_types =TypeComponentsTypeString(atype.first,names::P3Components());
	    args.ReplaceAll(names::P3Components(), atype.second[names::P3Components()]);
	  }
	  
 	  auto function_expr = sfunc + "<" + obj_types + ">";
 
	  cout<<"DefineForAllTypes "<<function_expr<<" ; "<<indices<<" ; "<<args<<endl;
	  //if function takes combi indices
	  if(indices.empty()==false){
	    auto defString = Form("rad::util::ApplyCombinations(%s, %s, %s)",
				  function_expr.data(),
				  indices.data(), // "combos"
				  args.Data()); 
	    Define(atype.first + name.data(), defString);
	  }
	  //if function does not take combi indices
	  else{
	    auto defString = Form("%s(%s)",
				  function_expr.data(),
				  args.Data()); 
	    Define(atype.first + name.data(), defString);
	  }
	  
	}
      }
      
      void AddParticleName(const std::string& particle){_particleNames.push_back(particle);}
      void AddFinalParticleName(const std::string& particle){_finalNames.push_back(particle);}
      const ROOT::RDF::ColumnNames_t& ParticleNames() const {return _particleNames;}
      const ROOT::RDF::ColumnNames_t& FinalParticleNames() const {return _finalNames;}

      // NOTE: Other existing methods (setParticleIndex, makeParticleMap, setGroupParticles, etc.) remain in ConfigReaction
      // but should be reviewed for compatibility with the new combinatorial mode and removed if made redundant.

      /**
       * Make map that links particle names to indices in user functions
       * in C++ functions you can use the RVecIndexMap object indexed by 
       * name of the reaction component you need
       */
      virtual void makeParticleMap() {
	// std::string particle_func("1E6+");
	// for(auto& part : _particleNames){
	//   particle_func+=part+"+";
	// }
	// particle_func.pop_back(); //remove last +

	// Filter(particle_func.data(),"particle_list");

	PostParticles();
      }
      /**
       *Any additional stuff to be done after all particles have been indiced
       */
      virtual void PostParticles(){

      }

      /**
       * create shortcut string for 3 and 4 momentum components
       */
      void AddType(const string& atype){
	if(_primary_type.empty()==true) _primary_type=atype;
	_type_comps[atype][names::P4Components()] = Form("%spx,%spy,%spz,%sm",atype.data(),atype.data(),atype.data(),atype.data());
	_type_comps[atype][names::P3Components()] = Form("%spx,%spy,%spz",atype.data(),atype.data(),atype.data());

	_types.push_back(atype);
     }

      /** 
       * Return columns types as strings to include in expression
       * for templated functions
       */
      string TypeComponentsTypeString(const string& type,const string& var){
	/// for 4-vector p and m may have different types
	/// so must return 2 types
	if(var==names::P4Components()){
	  return ColObjTypeString(type+"px")+","+ColObjTypeString(type+"m");
	}
	else if (var==names::P3Components()){
	  return ColObjTypeString(type+"px");
	}
	else{
	  cout<<"TypeColTypeString no valid vars"<<endl;exit(0);
	}
      }
      /**
       * Get names of all types of data
       */
      std::vector<std::string> GetTypes() const {return _types;}
      
      /**
       *  change name of this first type column to same name without type prefix
       */
      void AliasToPrimaryType(const string& name){
	if(_primary_type.empty()==true) return;
	std::string fullName = _primary_type + name;
	if (!OriginalColumnExists(fullName)) {
	  throw std::invalid_argument("AliasToPrimaryType: Column '" + fullName + "' does not exist.");
	}
	setBranchAlias(_primary_type+name,name);
      }

      const ROOT::RDF::ColumnNames_t getGroup(const string& name) const {return _groupMap.at(name);}
      
    protected:

      bool _useBeamsFromMC=false; 
       //------------------ Private Members ------------------
    
    private:
      // Combinatorial Members
      std::map<string, std::string> _candidateExpressions;
      std::map<string, ROOT::RDF::ColumnNames_t> _lambdaCandidateDependencies; 
      std::map<string, ROOT::RDF::ColumnNames_t> _groupMap; 
      bool _isCombinatorialMode = false;
      // std::string _firstDefinedCombiColumn; // Needed for robust Unroll in Snapshot

      std::map<string, std::map<string, string>> _type_comps;
      std::vector<std::string> _types;
      std::string _primary_type;
      ROOT::RDF::ColumnNames_t _particleNames;
      ROOT::RDF::ColumnNames_t _finalNames;
      // ... other existing members ...
      
    }; // class ConfigReaction

  } // namespace config
} // namespace rad
