#pragma once

#include "RDFInterface.h" 
#include "CommonDefines.h" 
#include "ConfigUtils.h" 
#include "Combinatorics.h" 
#include "CombinatorialUtil.h" 
#include "Constants.h" 
#include "DefineNames.h"
#include "RVecHelpers.h"
#include "ReactionUtilities.h"
#include "StringUtilities.h"
#include "Random.h"
#include <algorithm> 


namespace rad {

    /**
     * @class ConfigReaction
     * @brief The Central Configuration Manager for Reaction Analysis.
     * * @details
     * This class extends `RDFInterface` to add "Physics Awareness". 
     * It manages Data Types, Candidate Selection, and Combinatorial Generation.
     */
    class ConfigReaction : public RDFInterface {

    public:
      // --- Constructors ---
      ConfigReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns);
      ConfigReaction(const std::string_view treeName, const std::vector<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns);
      ConfigReaction(ROOT::RDataFrame rdf);

      // --- Output Management ---
      void Snapshot(const string& filename) override;
      void BookLazySnapshot(const string& filename) override;
      void RemoveSnapshotColumns(std::vector<string>& cols) override;

      // --- Type System ---
      /** @brief Registers a data type (e.g. "rec_") for analysis. */
      void AddType(const string& atype);

      /** @brief Validates that a type has been registered. */
      void ValidateType(const string& type) const;

      std::vector<std::string> GetTypes() const;

      // --- Candidate Definition (Input Selection) ---

      /** @brief Define candidates using a string expression (JIT). */
      void setParticleCandidatesExpr(const string& name, const string& type, const string& expression);

      /** @brief Define candidates using a C++ Lambda. */
      template<typename Lambda>
      void setParticleCandidates(const string& name, const string& type, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns);
      
      /** @brief Shortcut: Set a single fixed index as the candidate. */
      void setParticleIndex(const string& name, const string& type, const int idx);

      /** @brief Shortcut: Set a fixed list of indices as candidates. */
      void setParticleCandidates(const string& name, const string& type, const Indices_t idx);

      // --- Overloads relying on Default Type ---
      void setParticleCandidatesExpr(const string& name, const string& expression);
      
      template<typename Lambda>
      void setParticleCandidates(const string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns);
      
      void setParticleIndex(const string& name, const int idx);
      void setParticleCandidates(const string& name, const Indices_t idx);

      // --- Grouping Logic ---

      void setGroupParticles(const string& name, const string& type, const ROOT::RDF::ColumnNames_t& particles);

      void setMesonParticles(const string& type, const ROOT::RDF::ColumnNames_t& particles);
      void setBaryonParticles(const string& type, const ROOT::RDF::ColumnNames_t& particles);

      // Global Group Setters (Apply to ALL registered types)
      void setMesonParticles(const ROOT::RDF::ColumnNames_t& particles);
      void setBaryonParticles(const ROOT::RDF::ColumnNames_t& particles);
      void setGroupParticles(const string& name, const ROOT::RDF::ColumnNames_t& particles);

      // --- Combinatorial Engine ---

      /** @brief Triggers the generation of combinatorial events. */
      void makeCombinations();

      // --- Generic Definition Interface ---

      /** @brief Defines a column for ALL registered types using string substitution. */
      void DefineForAllTypes(const string& name, const string& expression);

      /** @brief Defines a generic function call for ALL registered types (Templated). */
      void DefineForAllTypes(const string& name, const string& sfunc, const string& indices, const string& arguments);
      
      // --- Utilities ---

      void AddParticleName(const std::string& particle) { _particleNames.push_back(particle); }
      void AddFinalParticleName(const std::string& particle) { _finalNames.push_back(particle); }
      
      const ROOT::RDF::ColumnNames_t& ParticleNames() const { return _particleNames; }
      const ROOT::RDF::ColumnNames_t& FinalParticleNames() const { return _finalNames; }

      virtual void makeParticleMap();
      virtual void PostParticles() {}

      const ROOT::RDF::ColumnNames_t getGroup(const string& name) const;
      
      string TypeComponentsTypeString(const string& type, const string& var);

      void AliasToPrimaryType(const string& name);

    protected:
      bool _useBeamsFromMC = false; 
      const string& GetDefaultType() const;

    private:
      void RegisterParticleName(const string& name);

      std::map<string, std::map<string, std::string>> _typeCandidateExpressions;
      std::map<string, std::map<string, ROOT::RDF::ColumnNames_t>> _typeLambdaDependencies; 

      std::map<string, ROOT::RDF::ColumnNames_t> _groupMap; 
      bool _isCombinatorialMode = false;

      std::map<string, std::map<string, string>> _type_comps;
      std::vector<std::string> _types;
      std::string _primary_type;
      
      ROOT::RDF::ColumnNames_t _particleNames;
      ROOT::RDF::ColumnNames_t _finalNames;
      
    }; 

    // =======================================================================
    // IMPLEMENTATION
    // =======================================================================

    inline ConfigReaction::ConfigReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns)
      : RDFInterface(treeName, fileNameGlob, columns) {}

    inline ConfigReaction::ConfigReaction(const std::string_view treeName, const std::vector<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns)
      : RDFInterface(treeName, filenames, columns) {}

    inline ConfigReaction::ConfigReaction(ROOT::RDataFrame rdf)
      : RDFInterface(rdf) {}

    // --- Output Management ---
    inline void ConfigReaction::Snapshot(const string& filename) {
      try {
        RDFstep final_df = CurrFrame();
        auto cols = final_df.GetDefinedColumnNames();
        RemoveSnapshotColumns(cols);
        final_df.Snapshot("rad_tree", filename, cols);
      } catch (const std::exception& ex) {
        std::cerr << "Snapshot failed: " << ex.what() << std::endl;
        throw;
      }
    }

    inline void ConfigReaction::BookLazySnapshot(const string& filename) {
      try {
        RDFstep final_df = CurrFrame();
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

    inline void ConfigReaction::RemoveSnapshotColumns(std::vector<string>& cols) {
      cols.erase(std::remove(cols.begin(), cols.end(), consts::ReactionMap()), cols.end());
      auto tag = DoNotWriteTag();
      cols.erase(std::remove_if(cols.begin(), cols.end(),
               [&tag](const string& col) -> bool { return col.find(tag) != std::string::npos; }),cols.end());
      RDFInterface::RemoveSnapshotColumns(cols);
    }

    // --- Type System ---
    inline void ConfigReaction::AddType(const string& atype) {
      if(_primary_type.empty()) _primary_type = atype;
      _type_comps[atype][consts::P4Components()] = Form("%spx,%spy,%spz,%sm", atype.data(), atype.data(), atype.data(), atype.data());
      _type_comps[atype][consts::P3Components()] = Form("%spx,%spy,%spz", atype.data(), atype.data(), atype.data());
      _types.push_back(atype);
    }

    inline void ConfigReaction::ValidateType(const string& type) const {
      if (std::find(_types.begin(), _types.end(), type) == _types.end()) {
        throw std::invalid_argument("Error: Data type '" + type + "' is not registered.");
      }
    }

    inline std::vector<std::string> ConfigReaction::GetTypes() const { return _types; }

    // --- Candidate Definition ---
    inline void ConfigReaction::setParticleCandidatesExpr(const string& name, const string& type, const string& expression) {
      ValidateType(type);
      if (_typeCandidateExpressions[type].count(name) || _typeLambdaDependencies[type].count(name)) {
        throw std::invalid_argument("Candidate '" + name + "' already defined for type '" + type + "'.");
      }
      _typeCandidateExpressions[type][name] = expression;
      RegisterParticleName(name);
    }

    template<typename Lambda>
    inline void ConfigReaction::setParticleCandidates(const string& name, const string& type, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
      ValidateType(type);
      if (_typeCandidateExpressions[type].count(name) || _typeLambdaDependencies[type].count(name)) {
        throw std::invalid_argument("Candidate '" + name + "' already defined for type '" + type + "'.");
      }
      string colName = type + name; 
      Define(colName, std::forward<Lambda>(func), columns);
      _typeLambdaDependencies[type][name] = columns; 
      RegisterParticleName(name);
    }

    inline void ConfigReaction::setParticleIndex(const string& name, const string& type, const int idx) {
       setParticleCandidates(name, type, [idx](){ return RVecI{idx}; }, {});
    }

    inline void ConfigReaction::setParticleCandidates(const string& name, const string& type, const Indices_t idx) {
       setParticleCandidates(name, type, [idx](){ return idx; }, {});
    }

    // --- Overloads ---
    inline void ConfigReaction::setParticleCandidatesExpr(const string& name, const string& expression) {
        setParticleCandidatesExpr(name, GetDefaultType(), expression);
    }
    
    template<typename Lambda>
    inline void ConfigReaction::setParticleCandidates(const string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
        setParticleCandidates(name, GetDefaultType(), std::forward<Lambda>(func), columns);
    }
    
    inline void ConfigReaction::setParticleIndex(const string& name, const int idx) {
        setParticleIndex(name, GetDefaultType(), idx);
    }
    
    inline void ConfigReaction::setParticleCandidates(const string& name, const Indices_t idx) {
        setParticleCandidates(name, GetDefaultType(), idx);
    }

    // --- Grouping Logic ---
    inline void ConfigReaction::setGroupParticles(const string& name, const string& type, const ROOT::RDF::ColumnNames_t& particles) {
      ValidateType(type);
      ROOT::RDF::ColumnNames_t typedParticles;
      for(const auto& p : particles) {
          typedParticles.push_back(type + p);
      }
      auto pstring = util::ColumnsToString(typedParticles); 
      pstring = pstring.substr(1, pstring.size() - 2); 

      string groupColName = type + name;
      Define(groupColName, Form("rad::util::Group<rad::Indices_t>(%s)", pstring.data()));
      _groupMap[groupColName] = typedParticles;
    }

    inline void ConfigReaction::setMesonParticles(const string& type, const ROOT::RDF::ColumnNames_t& particles) {
      ValidateType(type);
      if (particles.empty()) {
        Define(type + as_string(consts::Mesons()), [](){return RVecIndices{{consts::InvalidIndex()}}; }, {}); 
        return;
      }
      setGroupParticles(as_string(consts::Mesons()), type, particles);
    }

    inline void ConfigReaction::setBaryonParticles(const string& type, const ROOT::RDF::ColumnNames_t& particles) {
      ValidateType(type);
      if (particles.empty()) {
        Define(type + as_string(consts::Baryons()), [](){return RVecIndices{{consts::InvalidIndex()}}; }, {}); 
        return;
      }
      setGroupParticles(as_string(consts::Baryons()), type, particles);
    }

    inline void ConfigReaction::setMesonParticles(const ROOT::RDF::ColumnNames_t& particles) {
       if(_types.empty()) throw std::runtime_error("setMesonParticles: No types registered.");
       for(const auto& type : _types) setMesonParticles(type, particles);
    }
    
    inline void ConfigReaction::setBaryonParticles(const ROOT::RDF::ColumnNames_t& particles) {
       if(_types.empty()) throw std::runtime_error("setBaryonParticles: No types registered.");
       for(const auto& type : _types) setBaryonParticles(type, particles);
    }
    
    inline void ConfigReaction::setGroupParticles(const string& name, const ROOT::RDF::ColumnNames_t& particles) {
       if(_types.empty()) throw std::runtime_error("setGroupParticles: No types registered.");
       for(const auto& type : _types) setGroupParticles(name, type, particles);
    }

    // --- Combinatorial Engine ---
    inline void ConfigReaction::makeCombinations() {
      if (_types.empty()) {
          throw std::runtime_error("makeCombinations: No types (rec_, tru_) registered via AddType.");
      }
      _isCombinatorialMode = true;

      for (const auto& type : _types) {
          
          ROOT::RDF::ColumnNames_t currentTypeCandidateCols;
          
          // 1. Gather Candidates (Strings)
          if (_typeCandidateExpressions.count(type)) {
              for (const auto& pair : _typeCandidateExpressions[type]) {
                  const string& name = pair.first;
                  string colName = type + name; 
                  Define(colName, pair.second);
                  currentTypeCandidateCols.push_back(colName);
              }
          }

          // 2. Gather Candidates (Lambdas)
          if (_typeLambdaDependencies.count(type)) {
              for (const auto& pair : _typeLambdaDependencies[type]) {
                  currentTypeCandidateCols.push_back(type + pair.first);
              }
          }

          if (currentTypeCandidateCols.empty()) {
              std::cerr << "WARNING: Type '" << type << "' is registered but has no particle candidates set." << std::endl;
              continue; 
          }

          // 3. Generate Combinations (The Cartesian Product)
          string comboColName = type + consts::ReactionCombos();
          Define(comboColName,
             util::createFunctionCallStringFromVec("rad::combinatorics::GenerateAllCombinations",
                    {rad::util::ColumnsToString(currentTypeCandidateCols)}));

          // 4. Overwrite Particle Columns (Redefine as single index for current combo)
          for(size_t ip=0; ip < currentTypeCandidateCols.size(); ++ip) {
              Redefine(currentTypeCandidateCols[ip], 
                      [ip](const RVecIndices& part_combos){ return part_combos[ip]; },
                      {comboColName});
          }
      }
    }

    // --- Generic Definition Interface ---
    inline void ConfigReaction::DefineForAllTypes(const string& name, const string& expression) {
      for(auto &atype : _type_comps){
        if (atype.second.find(consts::P4Components()) == atype.second.end() ||
            atype.second.find(consts::P3Components()) == atype.second.end()) {
           throw std::runtime_error("DefineForAllTypes: Missing components for type: " + atype.first);
        }
        TString type_expr = expression.data();
        type_expr.ReplaceAll(consts::P4Components(), atype.second[consts::P4Components()]);
        type_expr.ReplaceAll(consts::P3Components(), atype.second[consts::P3Components()]);
        Define(atype.first + name.data(), type_expr.Data());
      }
    }

    inline void ConfigReaction::DefineForAllTypes(const string& name, const string& sfunc, const string& indices, const string& arguments) {
      for (auto& atype : _type_comps) {
        if (atype.second.find(consts::P4Components()) == atype.second.end()) {
          throw std::runtime_error("DefineForAllTypes: Missing 'components_p4' for type: " + atype.first);
        }
        
        TString args = arguments.data();
        string obj_types;
        
        if(args.Contains(consts::P4Components())){
          obj_types = TypeComponentsTypeString(atype.first, consts::P4Components());
          args.ReplaceAll(consts::P4Components(), atype.second[consts::P4Components()]);
        }
        else if(args.Contains(consts::P3Components())){
          obj_types = TypeComponentsTypeString(atype.first, consts::P3Components());
          args.ReplaceAll(consts::P3Components(), atype.second[consts::P3Components()]);
        }
        
        auto function_expr = sfunc + "<" + obj_types + ">";
        
        if(!indices.empty()){
          auto defString = Form("rad::util::ApplyCombinations(%s, %s, %s)",
                function_expr.data(), indices.data(), args.Data()); 
          Define(atype.first + name.data(), defString);
        }
        else{
          auto defString = Form("%s(%s)", function_expr.data(), args.Data()); 
          Define(atype.first + name.data(), defString);
        }
      }
    }

    // --- Utilities ---
    inline void ConfigReaction::makeParticleMap() { 
        PostParticles(); 
    }

    inline const ROOT::RDF::ColumnNames_t ConfigReaction::getGroup(const string& name) const {
      return _groupMap.at(name);
    }

    inline string ConfigReaction::TypeComponentsTypeString(const string& type, const string& var) {
      if(var == consts::P4Components()){
        return ColObjTypeString(type + "px") + "," + ColObjTypeString(type + "m");
      }
      else if (var == consts::P3Components()){
        return ColObjTypeString(type + "px");
      }
      throw std::runtime_error("TypeComponentsTypeString: Invalid variable placeholder.");
    }

    inline void ConfigReaction::AliasToPrimaryType(const string& name) {
      if(_primary_type.empty()) return;
      std::string fullName = _primary_type + name;
      if (!OriginalColumnExists(fullName)) {
        throw std::invalid_argument("AliasToPrimaryType: Column '" + fullName + "' does not exist.");
      }
      setBranchAlias(_primary_type + name, name);
    }

    inline const string& ConfigReaction::GetDefaultType() const {
        if (_types.empty()) {
            throw std::runtime_error("Reaction Class Error: No types registered. Call AddType() first.");
        }
        if (_types.size() > 1) {
            std::cerr << "WARNING: Defaulting to first type: '" << _types[0] << "'." << std::endl;
        }
        return _types[0];
    }

    inline void ConfigReaction::RegisterParticleName(const string& name) {
        if(std::find(_particleNames.begin(), _particleNames.end(), name) == _particleNames.end()) {
           AddParticleName(name);
           AddFinalParticleName(name);
      }
    }

} // namespace rad
