/**
 * @file ConfigReaction.h
 * @brief Central configuration manager for RAD analysis.
 */

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
#include "TruthMatchRegistry.h" 
#include "Random.h"
#include <algorithm> 
#include <set>

namespace rad {

    class ConfigReaction : public RDFInterface {

    public:
      // --- Constructors ---
      ConfigReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns);
      ConfigReaction(const std::string_view treeName, const std::vector<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns);
      ConfigReaction(ROOT::RDataFrame rdf);

      // =======================================================================
      // Truth Matching Interface
      // =======================================================================
      
      /** * @brief Define a candidate and immediately register its Truth Role. */
      void SetParticleTruthMatch(const std::string& name, int truthRole, const std::string& type, const Indices_t idx);

      template<typename Lambda>
      void SetParticleTruthMatch(const std::string& name, int truthRole, const std::string& type, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
          SetParticleCandidates(name, type, std::forward<Lambda>(func), columns);
          _truthMatchRegistry.AddParticleMatch(name, truthRole);
          SetParticleIndex(name, consts::data_type::Truth(), truthRole);
      }

      void AddTruthMatch(const std::string& name, int truthRole) {
          _truthMatchRegistry.AddParticleMatch(name, truthRole);
      }
      
      TruthMatchRegistry& GetTruthMatchRegistry() { return _truthMatchRegistry; }

      /**
       * @brief Generates the combinatorial truth flag "is_signal_combi".
       * @param matchIdCol The column containing Truth IDs (e.g. "rec_match_id").
       * @param type The prefix of the candidates to check (Default: "rec_").
       */
      void DefineTrueMatchedCombi(const std::string& matchIdCol, const std::string& type = rad::consts::data_type::Rec());

      // --- Symmetry Interface ---
      template<typename... Args>
      void SetSymmetryParticles(Args... args) {
          std::vector<std::string> group = {args...};
          if(group.size() > 1) _symmetryGroups.push_back(group);
      }

      // --- Type System ---
      void AddType(const std::string& atype);
      void ValidateType(const std::string& type) const;
      std::vector<std::string> GetTypes() const;

      // --- Candidate Definition ---
      void SetParticleCandidatesExpr(const std::string& name, const std::string& type, const std::string& expression);

      template<typename Lambda>
      void SetParticleCandidates(const std::string& name, const std::string& type, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns);

       void SetParticleIndex(const std::string& name, const std::string& type, const int idx);
      void SetParticleCandidates(const std::string& name, const std::string& type, const Indices_t idx);

      // Overloads
      void SetParticleCandidatesExpr(const std::string& name, const std::string& expression);
      
      template<typename Lambda>
      void SetParticleCandidates(const std::string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns);
      
      void SetParticleIndex(const std::string& name, const int idx);
      void SetParticleCandidates(const std::string& name, const Indices_t idx);

      void SetParticleCandidates(const std::string& name, int truthRole,  const Indices_t idx){//use default type and mcmatch
	SetParticleTruthMatch(name, truthRole, GetDefaultType(), idx);
      }
      
      template<typename Lambda>
      void SetParticleCandidates(const std::string& name, int truthRole, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {//use default type and mcmatch
	SetParticleTruthMatch(name,truthRole, GetDefaultType(),func, columns);
      }
      
      
      // --- Grouping Logic ---
      void SetGroupParticles(const std::string& name, const std::string& type, const ROOT::RDF::ColumnNames_t& particles);
      void SetMesonParticles(const std::string& type, const ROOT::RDF::ColumnNames_t& particles);
      void SetBaryonParticles(const std::string& type, const ROOT::RDF::ColumnNames_t& particles);

      void SetMesonParticles(const ROOT::RDF::ColumnNames_t& particles);
      void SetBaryonParticles(const ROOT::RDF::ColumnNames_t& particles);
      void SetGroupParticles(const std::string& name, const ROOT::RDF::ColumnNames_t& particles);

      // --- Combinatorial Engine ---
      void MakeCombinations();

      // --- Generic Definition ---
      void DefineForAllTypes(const std::string& name, const std::string& expression);
      void DefineForAllTypes(const std::string& name, const std::string& sfunc, const std::string& indices, const std::string& arguments);
      
      // --- Utilities ---
      void AddParticleName(const std::string& particle) { _particleNames.push_back(particle); }
      void AddFinalParticleName(const std::string& particle) { _finalNames.push_back(particle); }
      const ROOT::RDF::ColumnNames_t& ParticleNames() const { return _particleNames; }
      const ROOT::RDF::ColumnNames_t& FinalParticleNames() const { return _finalNames; }

      virtual void MakeParticleMap();
      virtual void PostParticles() {}
      const ROOT::RDF::ColumnNames_t GetGroup(const std::string& name) const;
      
      void Snapshot(const std::string& filename) override;
      void BookLazySnapshot(const std::string& filename) override;
      void RemoveSnapshotColumns(std::vector<std::string>& cols) override;

    protected:
      bool _useBeamsFromMC = false; 
      const std::string& GetDefaultType() const;
      TruthMatchRegistry _truthMatchRegistry; 
      std::vector<std::vector<std::string>> _symmetryGroups;

    private:
      void RegisterParticleName(const std::string& name);
      std::string TypeComponentsTypeString(const std::string& type, const std::string& var);

      std::map<std::string, std::map<std::string, std::string>> _typeCandidateExpressions;
      std::map<std::string, std::map<std::string, ROOT::RDF::ColumnNames_t>> _typeLambdaDependencies; 

      std::map<std::string, ROOT::RDF::ColumnNames_t> _groupMap; 
      bool _isCombinatorialMode = false;

      std::map<std::string, std::map<std::string, std::string>> _type_comps;
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
    inline ConfigReaction::ConfigReaction(ROOT::RDataFrame rdf) : RDFInterface(rdf) {}

    // --- Candidate Definition & Truth ---

    inline void ConfigReaction::SetParticleTruthMatch(const std::string& name, int truthRole, const std::string& type, const Indices_t idx) {
        SetParticleCandidates(name, type, idx); 
        _truthMatchRegistry.AddParticleMatch(name, truthRole);
        SetParticleIndex(name, consts::data_type::Truth(), truthRole);
    }

    // Lambda now iterates over the RVecI (pIndices) which represents the combinations
    inline void ConfigReaction::DefineTrueMatchedCombi(const std::string& matchIdCol, const std::string& type) {
        std::string logic = "";
        
        for (auto const& match : _truthMatchRegistry.GetParticleMatches()) {
            const std::string& name = match.first; 
            int role = match.second;
            
            std::string flagName = name + "_is_true" + DoNotWriteTag();
            std::string colName = type + name; // e.g. "rec_ele"
            
            // Check: pIndices[i] is the candidate index for the i-th combination.
            // matchIds is the array of truth IDs for ALL tracks.
            Define(flagName, 
                [role](const Indices_t& pIndices, const Indices_t& matchIds) {
                    Indices_t result(pIndices.size());
                     for(size_t i = 0; i < pIndices.size(); ++i) {
                         int pIdx = pIndices[i];
                        if (pIdx >= 0  && matchIds[pIdx] == role) {
                             result[i] = 1;
                         } else {
                             result[i] = 0;
                         }
                     }
                    return result;
                }, 
                {colName, matchIdCol});

            if (!logic.empty()) logic += " && ";
            logic += flagName;
        }

        if(logic.empty()) logic = "1";
        Define("is_signal_combi", logic);
    }

    // ... [Rest of implementation remains unchanged] ...
    inline void ConfigReaction::SetParticleCandidatesExpr(const std::string& name, const std::string& type, const std::string& expression) {
      ValidateType(type);
      if (_typeCandidateExpressions[type].count(name) || _typeLambdaDependencies[type].count(name)) {
        throw std::invalid_argument("Candidate '" + name + "' already defined for type '" + type + "'.");
      }
      _typeCandidateExpressions[type][name] = expression;
      RegisterParticleName(name);
    }
    template<typename Lambda>
    inline void ConfigReaction::SetParticleCandidates(const std::string& name, const std::string& type, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
      
      ValidateType(type);
      if (_typeCandidateExpressions[type].count(name) || _typeLambdaDependencies[type].count(name)) {
        throw std::invalid_argument("Candidate '" + name + "' already defined for type '" + type + "'.");
      }
      std::string colName = type + name; 
      Define(colName, std::forward<Lambda>(func), columns);
      _typeLambdaDependencies[type][name] = columns; 
      RegisterParticleName(name);
    }
    inline void ConfigReaction::SetParticleIndex(const std::string& name, const std::string& type, const int idx) {
       SetParticleCandidates(name, type, [idx](){ return RVecI{idx}; }, {});
    }
    inline void ConfigReaction::SetParticleCandidates(const std::string& name, const std::string& type, const Indices_t idx) {
       SetParticleCandidates(name, type, [idx](){ return idx; }, {});
    }
    inline void ConfigReaction::SetParticleCandidatesExpr(const std::string& name, const std::string& expression) {
        SetParticleCandidatesExpr(name, GetDefaultType(), expression);
    }
    template<typename Lambda>
    inline void ConfigReaction::SetParticleCandidates(const std::string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
        SetParticleCandidates(name, GetDefaultType(), std::forward<Lambda>(func), columns);
    }
    inline void ConfigReaction::SetParticleIndex(const std::string& name, const int idx) {
        SetParticleIndex(name, GetDefaultType(), idx);
    }
    inline void ConfigReaction::SetParticleCandidates(const std::string& name, const Indices_t idx) {
        SetParticleCandidates(name, GetDefaultType(), idx);
    }
    
    // --- Combinatorial Engine ---
    inline void ConfigReaction::MakeCombinations() {
      if (_types.empty()) throw std::runtime_error("MakeCombinations: No types registered via AddType.");
      _isCombinatorialMode = true;

      for (const auto& type : _types) {
          ROOT::RDF::ColumnNames_t candidateCols;
          std::vector<std::string> rawNames; 

          auto collect = [&](const std::string& name) {
             rawNames.push_back(name);
             if (_typeCandidateExpressions[type].count(name)) {
                 std::string colName = type + name; 
                 Define(colName, _typeCandidateExpressions[type][name]);
                 candidateCols.push_back(colName);
             } else {
                 candidateCols.push_back(type + name);
             }
          };

          for(const auto& name : _particleNames) {
              if(_typeCandidateExpressions[type].count(name) || _typeLambdaDependencies[type].count(name)) {
                  collect(name);
              }
          }
          if (candidateCols.empty()) continue; 

          std::string comboColName = type + consts::ReactionCombos();
          std::string colList = rad::util::ColumnsToStringNoBraces(candidateCols);
          if(colList.size() >= 2 && colList.front() == '{' && colList.back() == '}') colList = colList.substr(1, colList.size() - 2);

          if (_symmetryGroups.empty()) {
              Define(comboColName, Form("rad::combinatorics::GenerateAllCombinations(%s)", colList.c_str()));
          } else {
              std::vector<std::string> groupStrs;
              for(const auto& group : _symmetryGroups) {
                  std::vector<std::string> idxStrs;
                  for(const auto& pName : group) {
                      auto it = std::find(rawNames.begin(), rawNames.end(), pName);
                      if(it != rawNames.end()) idxStrs.push_back(std::to_string(std::distance(rawNames.begin(), it)));
                  }
                  if(idxStrs.size() > 1) groupStrs.push_back("{" + util::combineVectorToString(idxStrs) + "}");
              }
              std::string symString = "{" + util::combineVectorToString(groupStrs) + "}";
              Define(comboColName, Form("rad::combinatorics::GenerateSymmetricCombinations(%s, %s)", colList.c_str(), symString.c_str()));
          }

          for(size_t ip=0; ip < candidateCols.size(); ++ip) {
              Redefine(candidateCols[ip], [ip](const RVecIndices& part_combos){ return part_combos[ip]; }, {comboColName});
          }
      }
    }

    // --- Grouping Logic ---
    inline void ConfigReaction::SetGroupParticles(const std::string& name, const std::string& type, const ROOT::RDF::ColumnNames_t& particles) {
      ValidateType(type);
      ROOT::RDF::ColumnNames_t typedParticles;
      for(const auto& p : particles) typedParticles.push_back(type + p);
      auto pstring = util::ColumnsToStringNoBraces(typedParticles); 
      if(pstring.size() >= 2 && pstring.front() == '{' && pstring.back() == '}') pstring = pstring.substr(1, pstring.size() - 2);
      std::string groupColName = type + name;
      Define(groupColName, Form("rad::util::Group<rad::Indices_t>(%s)", pstring.data()));
      _groupMap[groupColName] = typedParticles;
    }
    inline void ConfigReaction::SetMesonParticles(const std::string& type, const ROOT::RDF::ColumnNames_t& particles) {
      ValidateType(type);
      if (particles.empty()) { Define(type + as_string(consts::Mesons()), [](){return RVecIndices{{consts::InvalidIndex()}}; }, {}); return; }
      SetGroupParticles(as_string(consts::Mesons()), type, particles);
    }
    inline void ConfigReaction::SetBaryonParticles(const std::string& type, const ROOT::RDF::ColumnNames_t& particles) {
      ValidateType(type);
      if (particles.empty()) { Define(type + as_string(consts::Baryons()), [](){return RVecIndices{{consts::InvalidIndex()}}; }, {}); return; }
      SetGroupParticles(as_string(consts::Baryons()), type, particles);
    }
    inline void ConfigReaction::SetMesonParticles(const ROOT::RDF::ColumnNames_t& particles) { for(const auto& type : _types) SetMesonParticles(type, particles); }
    inline void ConfigReaction::SetBaryonParticles(const ROOT::RDF::ColumnNames_t& particles) { for(const auto& type : _types) SetBaryonParticles(type, particles); }
    inline void ConfigReaction::SetGroupParticles(const std::string& name, const ROOT::RDF::ColumnNames_t& particles) { for(const auto& type : _types) SetGroupParticles(name, type, particles); }

    // --- Utilities & Generic Definitions ---
    inline void ConfigReaction::DefineForAllTypes(const std::string& name, const std::string& expression) {
      for(auto &atype : _type_comps){
        TString type_expr = expression.data();
        type_expr.ReplaceAll(consts::P4Components(), atype.second[consts::P4Components()]);
        type_expr.ReplaceAll(consts::P3Components(), atype.second[consts::P3Components()]);
        Define(atype.first + name.data(), type_expr.Data());
      }
    }
    inline void ConfigReaction::DefineForAllTypes(const std::string& name, const std::string& sfunc, const std::string& indices, const std::string& arguments) {
      for (auto& atype : _type_comps) {
        TString args = arguments.data();
        std::string obj_types;
        if(args.Contains(consts::P4Components())){ obj_types = TypeComponentsTypeString(atype.first, consts::P4Components()); args.ReplaceAll(consts::P4Components(), atype.second[consts::P4Components()]); }
        else if(args.Contains(consts::P3Components())){ obj_types = TypeComponentsTypeString(atype.first, consts::P3Components()); args.ReplaceAll(consts::P3Components(), atype.second[consts::P3Components()]); }
        auto function_expr = sfunc + "<" + obj_types + ">";
        if(!indices.empty()) Define(atype.first + name.data(), Form("rad::util::ApplyCombinations(%s, %s, %s)", function_expr.data(), indices.data(), args.Data())); 
        else Define(atype.first + name.data(), Form("%s(%s)", function_expr.data(), args.Data())); 
      }
    }
    inline void ConfigReaction::MakeParticleMap() { PostParticles(); }
    inline const ROOT::RDF::ColumnNames_t ConfigReaction::GetGroup(const std::string& name) const { return _groupMap.at(name); }
    inline std::string ConfigReaction::TypeComponentsTypeString(const std::string& type, const std::string& var) {
      if(var == consts::P4Components()) return ColObjTypeString(type + "px") + "," + ColObjTypeString(type + "m");
      else if (var == consts::P3Components()) return ColObjTypeString(type + "px");
      return "";
    }
    inline const std::string& ConfigReaction::GetDefaultType() const {
        if (_types.empty()) throw std::runtime_error("Reaction Class Error: No types registered. Call AddType() first.");
        return _types[0];
    }
    inline void ConfigReaction::RegisterParticleName(const std::string& name) {
        if(std::find(_particleNames.begin(), _particleNames.end(), name) == _particleNames.end()) {
           AddParticleName(name); AddFinalParticleName(name);
      }
    }
    inline void ConfigReaction::Snapshot(const std::string& filename) {
        RDFstep final_df = CurrFrame();
        auto cols = final_df.GetDefinedColumnNames();
        RemoveSnapshotColumns(cols);
        final_df.Snapshot("rad_tree", filename, cols);
    }
    inline void ConfigReaction::BookLazySnapshot(const std::string& filename) {
        RDFstep final_df = CurrFrame();
        ROOT::RDF::RSnapshotOptions opts; opts.fLazy = true;
        auto cols = final_df.GetDefinedColumnNames();
        RemoveSnapshotColumns(cols);
        auto snapshot_result = final_df.Snapshot("rad_tree", filename, cols, opts);
        _triggerSnapshots.emplace_back([snapshot = std::move(snapshot_result)]() mutable {});
    }
    inline void ConfigReaction::RemoveSnapshotColumns(std::vector<std::string>& cols) {
      cols.erase(std::remove(cols.begin(), cols.end(), consts::ReactionMap()), cols.end());
      auto tag = DoNotWriteTag();
      cols.erase(std::remove_if(cols.begin(), cols.end(), [&tag](const std::string& col) -> bool { return col.find(tag) != std::string::npos; }),cols.end());
      RDFInterface::RemoveSnapshotColumns(cols);
    }
    inline void ConfigReaction::AddType(const std::string& atype) {
      if(_primary_type.empty()) _primary_type = atype;
      _type_comps[atype][consts::P4Components()] = Form("%spx,%spy,%spz,%sm", atype.data(), atype.data(), atype.data(), atype.data());
      _type_comps[atype][consts::P3Components()] = Form("%spx,%spy,%spz", atype.data(), atype.data(), atype.data());
      _types.push_back(atype);
    }
    inline void ConfigReaction::ValidateType(const std::string& type) const {
      if (std::find(_types.begin(), _types.end(), type) == _types.end()) throw std::invalid_argument("Error: Data type '" + type + "' is not registered.");
    }
    inline std::vector<std::string> ConfigReaction::GetTypes() const { return _types; }

} // namespace rad// #pragma once

// #include "RDFInterface.h" 
// #include "CommonDefines.h" 
// #include "ConfigUtils.h" 
// #include "Combinatorics.h" 
// #include "CombinatorialUtil.h" 
// #include "Constants.h" 
// #include "DefineNames.h"
// #include "RVecHelpers.h"
// #include "ReactionUtilities.h"
// #include "StringUtilities.h"
// #include "Random.h"
// #include <algorithm> 


// namespace rad {

//     /**
//      * @class ConfigReaction
//      * @brief The Central Configuration Manager for Reaction Analysis.
//      * * @details
//      * This class extends `RDFInterface` to add "Physics Awareness". 
//      * It manages Data Types, Candidate Selection, and Combinatorial Generation.
//      */
//     class ConfigReaction : public RDFInterface {

//     public:
//       // --- Constructors ---
//       ConfigReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns);
//       ConfigReaction(const std::string_view treeName, const std::vector<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns);
//       ConfigReaction(ROOT::RDataFrame rdf);

//       // --- Output Management ---
//       void Snapshot(const string& filename) override;
//       void BookLazySnapshot(const string& filename) override;
//       void RemoveSnapshotColumns(std::vector<string>& cols) override;

//       // --- Type System ---
//       /** @brief Registers a data type (e.g. "rec_") for analysis. */
//       void AddType(const string& atype);

//       /** @brief Validates that a type has been registered. */
//       void ValidateType(const string& type) const;

//       std::vector<std::string> GetTypes() const;

//       // --- Candidate Definition (Input Selection) ---

//       /** @brief Define candidates using a string expression (JIT). */
//       void setParticleCandidatesExpr(const string& name, const string& type, const string& expression);

//       /** @brief Define candidates using a C++ Lambda. */
//       template<typename Lambda>
//       void setParticleCandidates(const string& name, const string& type, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns);
      
//       /** @brief Shortcut: Set a single fixed index as the candidate. */
//       void setParticleIndex(const string& name, const string& type, const int idx);

//       /** @brief Shortcut: Set a fixed list of indices as candidates. */
//       void setParticleCandidates(const string& name, const string& type, const Indices_t idx);

//       // --- Overloads relying on Default Type ---
//       void setParticleCandidatesExpr(const string& name, const string& expression);
      
//       template<typename Lambda>
//       void setParticleCandidates(const string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns);
      
//       void setParticleIndex(const string& name, const int idx);
//       void setParticleCandidates(const string& name, const Indices_t idx);

//       // --- Grouping Logic ---

//       void setGroupParticles(const string& name, const string& type, const ROOT::RDF::ColumnNames_t& particles);

//       void setMesonParticles(const string& type, const ROOT::RDF::ColumnNames_t& particles);
//       void setBaryonParticles(const string& type, const ROOT::RDF::ColumnNames_t& particles);

//       // Global Group Setters (Apply to ALL registered types)
//       void setMesonParticles(const ROOT::RDF::ColumnNames_t& particles);
//       void setBaryonParticles(const ROOT::RDF::ColumnNames_t& particles);
//       void setGroupParticles(const string& name, const ROOT::RDF::ColumnNames_t& particles);

//       // --- Combinatorial Engine ---

//       /** @brief Triggers the generation of combinatorial events. */
//       void makeCombinations();

//       // --- Generic Definition Interface ---

//       /** @brief Defines a column for ALL registered types using string substitution. */
//       void DefineForAllTypes(const string& name, const string& expression);

//       /** @brief Defines a generic function call for ALL registered types (Templated). */
//       void DefineForAllTypes(const string& name, const string& sfunc, const string& indices, const string& arguments);
      
//       // --- Utilities ---

//       void AddParticleName(const std::string& particle) { _particleNames.push_back(particle); }
//       void AddFinalParticleName(const std::string& particle) { _finalNames.push_back(particle); }
      
//       const ROOT::RDF::ColumnNames_t& ParticleNames() const { return _particleNames; }
//       const ROOT::RDF::ColumnNames_t& FinalParticleNames() const { return _finalNames; }

//       virtual void makeParticleMap();
//       virtual void PostParticles() {}

//       const ROOT::RDF::ColumnNames_t getGroup(const string& name) const;
      
//       string TypeComponentsTypeString(const string& type, const string& var);

 
//     protected:
//       bool _useBeamsFromMC = false; 
//       const string& GetDefaultType() const;

//     private:
//       void RegisterParticleName(const string& name);

//       std::map<string, std::map<string, std::string>> _typeCandidateExpressions;
//       std::map<string, std::map<string, ROOT::RDF::ColumnNames_t>> _typeLambdaDependencies; 

//       std::map<string, ROOT::RDF::ColumnNames_t> _groupMap; 
//       bool _isCombinatorialMode = false;

//       std::map<string, std::map<string, string>> _type_comps;
//       std::vector<std::string> _types;
//       std::string _primary_type;
      
//       ROOT::RDF::ColumnNames_t _particleNames;
//       ROOT::RDF::ColumnNames_t _finalNames;
      
//     }; 

//     // =======================================================================
//     // IMPLEMENTATION
//     // =======================================================================

//     inline ConfigReaction::ConfigReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns)
//       : RDFInterface(treeName, fileNameGlob, columns) {}

//     inline ConfigReaction::ConfigReaction(const std::string_view treeName, const std::vector<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns)
//       : RDFInterface(treeName, filenames, columns) {}

//     inline ConfigReaction::ConfigReaction(ROOT::RDataFrame rdf)
//       : RDFInterface(rdf) {}

//     // --- Output Management ---
//     inline void ConfigReaction::Snapshot(const string& filename) {
//       try {
//         RDFstep final_df = CurrFrame();
//         auto cols = final_df.GetDefinedColumnNames();
//         RemoveSnapshotColumns(cols);
//         final_df.Snapshot("rad_tree", filename, cols);
//       } catch (const std::exception& ex) {
//         std::cerr << "Snapshot failed: " << ex.what() << std::endl;
//         throw;
//       }
//     }

//     inline void ConfigReaction::BookLazySnapshot(const string& filename) {
//       try {
//         RDFstep final_df = CurrFrame();
//         ROOT::RDF::RSnapshotOptions opts;
//         opts.fLazy = true;
//         auto cols = final_df.GetDefinedColumnNames();
//         RemoveSnapshotColumns(cols);
//         auto snapshot_result = final_df.Snapshot("rad_tree", filename, cols, opts);
//         _triggerSnapshots.emplace_back([snapshot = std::move(snapshot_result)]() mutable {});
//       } catch (const std::exception& ex) {
//         std::cerr << "BookLazySnapshot failed: " << ex.what() << std::endl;
//         throw;
//       }
//     }

//     inline void ConfigReaction::RemoveSnapshotColumns(std::vector<string>& cols) {
//       cols.erase(std::remove(cols.begin(), cols.end(), consts::ReactionMap()), cols.end());
//       auto tag = DoNotWriteTag();
//       cols.erase(std::remove_if(cols.begin(), cols.end(),
//                [&tag](const string& col) -> bool { return col.find(tag) != std::string::npos; }),cols.end());
//       RDFInterface::RemoveSnapshotColumns(cols);
//     }

//     // --- Type System ---
//     inline void ConfigReaction::AddType(const string& atype) {
//       if(_primary_type.empty()) _primary_type = atype;
//       _type_comps[atype][consts::P4Components()] = Form("%spx,%spy,%spz,%sm", atype.data(), atype.data(), atype.data(), atype.data());
//       _type_comps[atype][consts::P3Components()] = Form("%spx,%spy,%spz", atype.data(), atype.data(), atype.data());
//       _types.push_back(atype);
//     }

//     inline void ConfigReaction::ValidateType(const string& type) const {
//       if (std::find(_types.begin(), _types.end(), type) == _types.end()) {
//         throw std::invalid_argument("Error: Data type '" + type + "' is not registered.");
//       }
//     }

//     inline std::vector<std::string> ConfigReaction::GetTypes() const { return _types; }

//     // --- Candidate Definition ---
//     inline void ConfigReaction::setParticleCandidatesExpr(const string& name, const string& type, const string& expression) {
//       ValidateType(type);
//       if (_typeCandidateExpressions[type].count(name) || _typeLambdaDependencies[type].count(name)) {
//         throw std::invalid_argument("Candidate '" + name + "' already defined for type '" + type + "'.");
//       }
//       _typeCandidateExpressions[type][name] = expression;
//       RegisterParticleName(name);
//     }

//     template<typename Lambda>
//     inline void ConfigReaction::setParticleCandidates(const string& name, const string& type, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
//       ValidateType(type);
//       if (_typeCandidateExpressions[type].count(name) || _typeLambdaDependencies[type].count(name)) {
//         throw std::invalid_argument("Candidate '" + name + "' already defined for type '" + type + "'.");
//       }
//       string colName = type + name; 
//       Define(colName, std::forward<Lambda>(func), columns);
//       _typeLambdaDependencies[type][name] = columns; 
//       RegisterParticleName(name);
//     }

//     inline void ConfigReaction::setParticleIndex(const string& name, const string& type, const int idx) {
//        setParticleCandidates(name, type, [idx](){ return RVecI{idx}; }, {});
//     }

//     inline void ConfigReaction::setParticleCandidates(const string& name, const string& type, const Indices_t idx) {
//        setParticleCandidates(name, type, [idx](){ return idx; }, {});
//     }

//     // --- Overloads ---
//     inline void ConfigReaction::setParticleCandidatesExpr(const string& name, const string& expression) {
//         setParticleCandidatesExpr(name, GetDefaultType(), expression);
//     }
    
//     template<typename Lambda>
//     inline void ConfigReaction::setParticleCandidates(const string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
//         setParticleCandidates(name, GetDefaultType(), std::forward<Lambda>(func), columns);
//     }
    
//     inline void ConfigReaction::setParticleIndex(const string& name, const int idx) {
//         setParticleIndex(name, GetDefaultType(), idx);
//     }
    
//     inline void ConfigReaction::setParticleCandidates(const string& name, const Indices_t idx) {
//         setParticleCandidates(name, GetDefaultType(), idx);
//     }

//     // --- Grouping Logic ---
//     inline void ConfigReaction::setGroupParticles(const string& name, const string& type, const ROOT::RDF::ColumnNames_t& particles) {
//       ValidateType(type);
//       ROOT::RDF::ColumnNames_t typedParticles;
//       for(const auto& p : particles) {
//           typedParticles.push_back(type + p);
//       }
//       auto pstring = util::ColumnsToString(typedParticles); 
//       pstring = pstring.substr(1, pstring.size() - 2); 

//       string groupColName = type + name;
//       Define(groupColName, Form("rad::util::Group<rad::Indices_t>(%s)", pstring.data()));
//       _groupMap[groupColName] = typedParticles;
//     }

//     inline void ConfigReaction::setMesonParticles(const string& type, const ROOT::RDF::ColumnNames_t& particles) {
//       ValidateType(type);
//       if (particles.empty()) {
//         Define(type + as_string(consts::Mesons()), [](){return RVecIndices{{consts::InvalidIndex()}}; }, {}); 
//         return;
//       }
//       setGroupParticles(as_string(consts::Mesons()), type, particles);
//     }

//     inline void ConfigReaction::setBaryonParticles(const string& type, const ROOT::RDF::ColumnNames_t& particles) {
//       ValidateType(type);
//       if (particles.empty()) {
//         Define(type + as_string(consts::Baryons()), [](){return RVecIndices{{consts::InvalidIndex()}}; }, {}); 
//         return;
//       }
//       setGroupParticles(as_string(consts::Baryons()), type, particles);
//     }

//     inline void ConfigReaction::setMesonParticles(const ROOT::RDF::ColumnNames_t& particles) {
//        if(_types.empty()) throw std::runtime_error("setMesonParticles: No types registered.");
//        for(const auto& type : _types) setMesonParticles(type, particles);
//     }
    
//     inline void ConfigReaction::setBaryonParticles(const ROOT::RDF::ColumnNames_t& particles) {
//        if(_types.empty()) throw std::runtime_error("setBaryonParticles: No types registered.");
//        for(const auto& type : _types) setBaryonParticles(type, particles);
//     }
    
//     inline void ConfigReaction::setGroupParticles(const string& name, const ROOT::RDF::ColumnNames_t& particles) {
//        if(_types.empty()) throw std::runtime_error("setGroupParticles: No types registered.");
//        for(const auto& type : _types) setGroupParticles(name, type, particles);
//     }

//     // --- Combinatorial Engine ---
//     inline void ConfigReaction::makeCombinations() {
//       if (_types.empty()) {
//           throw std::runtime_error("makeCombinations: No types (rec_, tru_) registered via AddType.");
//       }
//       _isCombinatorialMode = true;

//       for (const auto& type : _types) {
          
//           ROOT::RDF::ColumnNames_t currentTypeCandidateCols;
          
//           // 1. Gather Candidates (Strings)
//           if (_typeCandidateExpressions.count(type)) {
//               for (const auto& pair : _typeCandidateExpressions[type]) {
//                   const string& name = pair.first;
//                   string colName = type + name; 
//                   Define(colName, pair.second);
//                   currentTypeCandidateCols.push_back(colName);
//               }
//           }

//           // 2. Gather Candidates (Lambdas)
//           if (_typeLambdaDependencies.count(type)) {
//               for (const auto& pair : _typeLambdaDependencies[type]) {
//                   currentTypeCandidateCols.push_back(type + pair.first);
//               }
//           }

//           if (currentTypeCandidateCols.empty()) {
//               std::cerr << "WARNING: Type '" << type << "' is registered but has no particle candidates set." << std::endl;
//               continue; 
//           }

//           // 3. Generate Combinations (The Cartesian Product)
//           string comboColName = type + consts::ReactionCombos();
//           Define(comboColName,
//              util::createFunctionCallStringFromVec("rad::combinatorics::GenerateAllCombinations",
//                     {rad::util::ColumnsToString(currentTypeCandidateCols)}));

//           // 4. Overwrite Particle Columns (Redefine as single index for current combo)
//           for(size_t ip=0; ip < currentTypeCandidateCols.size(); ++ip) {
//               Redefine(currentTypeCandidateCols[ip], 
//                       [ip](const RVecIndices& part_combos){ return part_combos[ip]; },
//                       {comboColName});
//           }
//       }
//     }

//     // --- Generic Definition Interface ---
//     inline void ConfigReaction::DefineForAllTypes(const string& name, const string& expression) {
//       for(auto &atype : _type_comps){
//         if (atype.second.find(consts::P4Components()) == atype.second.end() ||
//             atype.second.find(consts::P3Components()) == atype.second.end()) {
//            throw std::runtime_error("DefineForAllTypes: Missing components for type: " + atype.first);
//         }
//         TString type_expr = expression.data();
//         type_expr.ReplaceAll(consts::P4Components(), atype.second[consts::P4Components()]);
//         type_expr.ReplaceAll(consts::P3Components(), atype.second[consts::P3Components()]);
//         Define(atype.first + name.data(), type_expr.Data());
//       }
//     }

//     inline void ConfigReaction::DefineForAllTypes(const string& name, const string& sfunc, const string& indices, const string& arguments) {
//       for (auto& atype : _type_comps) {
//         if (atype.second.find(consts::P4Components()) == atype.second.end()) {
//           throw std::runtime_error("DefineForAllTypes: Missing 'components_p4' for type: " + atype.first);
//         }
        
//         TString args = arguments.data();
//         string obj_types;
        
//         if(args.Contains(consts::P4Components())){
//           obj_types = TypeComponentsTypeString(atype.first, consts::P4Components());
//           args.ReplaceAll(consts::P4Components(), atype.second[consts::P4Components()]);
//         }
//         else if(args.Contains(consts::P3Components())){
//           obj_types = TypeComponentsTypeString(atype.first, consts::P3Components());
//           args.ReplaceAll(consts::P3Components(), atype.second[consts::P3Components()]);
//         }
        
//         auto function_expr = sfunc + "<" + obj_types + ">";
        
//         if(!indices.empty()){
//           auto defString = Form("rad::util::ApplyCombinations(%s, %s, %s)",
//                 function_expr.data(), indices.data(), args.Data()); 
//           Define(atype.first + name.data(), defString);
//         }
//         else{
//           auto defString = Form("%s(%s)", function_expr.data(), args.Data()); 
//           Define(atype.first + name.data(), defString);
//         }
//       }
//     }

//     // --- Utilities ---
//     inline void ConfigReaction::makeParticleMap() { 
//         PostParticles(); 
//     }

//     inline const ROOT::RDF::ColumnNames_t ConfigReaction::getGroup(const string& name) const {
//       return _groupMap.at(name);
//     }

//     inline string ConfigReaction::TypeComponentsTypeString(const string& type, const string& var) {
//       if(var == consts::P4Components()){
//         return ColObjTypeString(type + "px") + "," + ColObjTypeString(type + "m");
//       }
//       else if (var == consts::P3Components()){
//         return ColObjTypeString(type + "px");
//       }
//       throw std::runtime_error("TypeComponentsTypeString: Invalid variable placeholder.");
//     }

 
//     inline const string& ConfigReaction::GetDefaultType() const {
//         if (_types.empty()) {
//             throw std::runtime_error("Reaction Class Error: No types registered. Call AddType() first.");
//         }
//         if (_types.size() > 1) {
//             std::cerr << "WARNING: Defaulting to first type: '" << _types[0] << "'." << std::endl;
//         }
//         return _types[0];
//     }

//     inline void ConfigReaction::RegisterParticleName(const string& name) {
//         if(std::find(_particleNames.begin(), _particleNames.end(), name) == _particleNames.end()) {
//            AddParticleName(name);
//            AddFinalParticleName(name);
//       }
//     }

// } // namespace rad
