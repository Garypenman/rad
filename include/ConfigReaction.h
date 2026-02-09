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
      ConfigReaction(const std::string_view treeName, const ROOT::RVec<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns);
      ConfigReaction(ROOT::RDataFrame rdf);
      ConfigReaction(ROOT::RDF::RNode rdf);


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
       * @brief Generates the combinatorial truth flag TruthMatchedCombi().
       * @param matchIdCol The column containing Truth IDs (e.g. "rec_match_id").
       * @param type The prefix of the candidates to check (Default: "rec_").
       */
      void DefineTrueMatchedCombi(const std::string& matchIdCol, const std::string& type);

      /** @brief Get the standard name for the matching column (e.g. "rec_match_id"). */
      virtual std::string GetMatchName(const std::string& type) const {
        // Standard convention: prefix + "match_id"
        return type + "match_id";
      }

      /** @brief Define Signal Combination flag using standard naming. 
       * @details Requires 'SetParticleCandidates' to have been called with truth roles.
       */
      void DefineTrueMatchedCombi(const std::string& type) {
        std::string col = GetMatchName(type);
        if(ColumnExists(col)) {
	  DefineTrueMatchedCombi(col, type);
        } else {
	  // If no matching column exists (e.g. Data or Truth stream), 
	  // define a default "Always True" signal flag so snapshots don't break.
	  std::cout<<"[ConfigReaction] DefineTrueMatchCombi " <<col<<" does not exist all set TruthMatchedCombi()=1 for all."<< std::endl;
	  Define(consts::TruthMatchedCombi(), "1"); 
        }
      }
      
      // --- Symmetry Interface ---
      template<typename... Args>
      void SetSymmetryParticles(Args... args) {
	ROOT::RVec<std::string> group = {args...};
          if(group.size() > 1) _symmetryGroups.push_back(group);
      }

      // --- Type System ---
      void AddType(const std::string& atype);
      void ValidateType(const std::string& type) const;
      ROOT::RVec<std::string> GetTypes() const;
      bool TypeExists(const std::string& type) const;
      
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
	cout<< "SetParticleCandidates "<< name<<endl;
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
      void RemoveSnapshotColumns(ROOT::RVec<std::string>& cols) override;

      //////////////////////
      ///Diagnostics
      ///\brief Print particle candidates for a given type.
      void PrintParticleCandidates(const std::string& type) const;

      ///\brief Print combinatoric structure.
      void PrintCombinatoricStructure() const;

      ///\brief Print all registered types.
      void PrintTypes() const;

      ///\brief Print comprehensive reaction diagnostics.
      void PrintReactionDiagnostics() const;

    protected:
      bool _useBeamsFromMC = false; 
      const std::string& GetDefaultType() const;
      TruthMatchRegistry _truthMatchRegistry; 
      ROOT::RVec<ROOT::RVec<std::string>> _symmetryGroups;

    private:
      void RegisterParticleName(const std::string& name);
      std::string TypeComponentsTypeString(const std::string& type, const std::string& var);

      std::map<std::string, std::map<std::string, std::string>> _typeCandidateExpressions;
      std::map<std::string, std::map<std::string, ROOT::RDF::ColumnNames_t>> _typeLambdaDependencies; 

      std::map<std::string, ROOT::RDF::ColumnNames_t> _groupMap; 
      bool _isCombinatorialMode = false;

      std::map<std::string, std::map<std::string, std::string>> _type_comps;
      ROOT::RVec<std::string> _types;
      std::string _primary_type;
      
      ROOT::RDF::ColumnNames_t _particleNames;
      ROOT::RDF::ColumnNames_t _finalNames;
      
    }; 

    // =======================================================================
    // IMPLEMENTATION
    // =======================================================================

    inline ConfigReaction::ConfigReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns)
      : RDFInterface(treeName, fileNameGlob, columns) {}
    inline ConfigReaction::ConfigReaction(const std::string_view treeName, const ROOT::RVec<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns)
      : RDFInterface(treeName, filenames, columns) {}
    inline ConfigReaction::ConfigReaction(ROOT::RDataFrame rdf) : RDFInterface(rdf) {}
    inline ConfigReaction::ConfigReaction(ROOT::RDF::RNode rdf) : RDFInterface(rdf) {}

    // --- Candidate Definition & Truth ---

    inline void ConfigReaction::SetParticleTruthMatch(const std::string& name, int truthRole, const std::string& type, const Indices_t idx) {
        SetParticleCandidates(name, type, idx); 
        _truthMatchRegistry.AddParticleMatch(name, truthRole);
        SetParticleIndex(name, consts::data_type::Truth(), truthRole);
    }

    // Lambda now iterates over the RVecI (pIndices) which represents the combinations
    inline void ConfigReaction::DefineTrueMatchedCombi(const std::string& matchIdCol, const std::string& type) {
        std::string logic = "";
	cout<<"ConfigReaction::DefineTrueMatchedCombi  "<<type<<endl;
	
        for (auto const& match : _truthMatchRegistry.GetParticleMatches()) {
            const std::string& name = match.first; 
            int role = match.second;
            
            std::string flagName = name + "_is_true";// + DoNotWriteTag();
            std::string colName = type + name; // e.g. "rec_ele"
            
            // Check: pIndices[i] is the candidate index for the i-th combination.
            // matchIds is the array of truth IDs for ALL tracks.
            Define(flagName, 
                [role](const Indices_t& pIndices, const Indices_t& matchIds) {
		  using namespace ROOT::VecOps;

		  // Create a mask of valid indices (filters out -1)
		  auto valid = pIndices >= 0;
		  auto currentMatches = Take(matchIds, pIndices);
		  
		  // Logic: (Is match?) AND (Was index valid?)
		  //    RVec logic operations automatically return vectors of 1s and 0s.
		  return (currentMatches == role) && valid;
    
              
                }, 
                {colName, matchIdCol});

            if (!logic.empty()) logic += " && ";
            logic += flagName;
        }

        if(logic.empty()) logic = "1";
	Define(consts::TruthMatchedCombi(), logic);
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

      cout<<"ConfigReaction::MakeCombinations "<<_types.size()<<endl;
      for (const auto& type : _types) {
          ROOT::RDF::ColumnNames_t candidateCols;
          ROOT::RVec<std::string> rawNames; 

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

	  cout<<"ConfigReaction::MakeCombinations "<<candidateCols.size()<<endl;
          if (candidateCols.empty()) continue; 

          std::string comboColName = type + consts::ReactionCombos();
          std::string colList = rad::util::ColumnsToStringNoBraces(candidateCols);
          if(colList.size() >= 2 && colList.front() == '{' && colList.back() == '}') colList = colList.substr(1, colList.size() - 2);

          if (_symmetryGroups.empty()) {
              Define(comboColName, Form("rad::combinatorics::GenerateAllCombinations(%s)", colList.c_str()));
          } else {
              ROOT::RVec<std::string> groupStrs;
              for(const auto& group : _symmetryGroups) {
                  ROOT::RVec<std::string> idxStrs;
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
      //Now have defined particle indices for all types in terms of combis
      //Apply truth matching if configured
      if( _truthMatchRegistry.GetParticleMatches().empty()==false ){
	if(TypeExists(consts::data_type::Rec())){
	  DefineTrueMatchedCombi(consts::data_type::Rec());
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
      if (particles.empty()) { Define(type + utils::as_string(consts::Mesons()), [](){return RVecIndices{{consts::InvalidIndex()}}; }, {}); return; }
      SetGroupParticles(utils::as_string(consts::Mesons()), type, particles);
    }
    inline void ConfigReaction::SetBaryonParticles(const std::string& type, const ROOT::RDF::ColumnNames_t& particles) {
      ValidateType(type);
      if (particles.empty()) { Define(type + utils::as_string(consts::Baryons()), [](){return RVecIndices{{consts::InvalidIndex()}}; }, {}); return; }
      SetGroupParticles(utils::as_string(consts::Baryons()), type, particles);
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
        auto cols = utils::as_rvec(final_df.GetDefinedColumnNames());
        RemoveSnapshotColumns(cols);
        final_df.Snapshot("rad_tree", filename, utils::as_stdvector(cols));
    }
    inline void ConfigReaction::BookLazySnapshot(const std::string& filename) {
        RDFstep final_df = CurrFrame();
        ROOT::RDF::RSnapshotOptions opts; opts.fLazy = true;
        auto cols = utils::as_rvec(final_df.GetDefinedColumnNames());
        RemoveSnapshotColumns(cols);
        auto snapshot_result = final_df.Snapshot("rad_tree", filename, utils::as_stdvector(cols), opts);
        _triggerSnapshots.emplace_back([snapshot = std::move(snapshot_result)]() mutable {});
    }
    inline void ConfigReaction::RemoveSnapshotColumns(ROOT::RVec<std::string>& cols) {
      cols.erase(std::remove(cols.begin(), cols.end(), consts::ReactionMap()), cols.end());
      auto tag = DoNotWriteTag();
      cols.erase(std::remove_if(cols.begin(), cols.end(), [&tag](const std::string& col) -> bool { return col.find(tag) != std::string::npos; }),cols.end());
      RDFInterface::RemoveSnapshotColumns(cols);
    }
    inline void ConfigReaction::AddType(const std::string& atype) {
      if(_primary_type.empty()) _primary_type = atype;
      _type_comps[atype][consts::P4Components()] = Form("%s%s,%s%s,%s%s,%s%s", atype.data(),consts::NamePx().data(), atype.data(),consts::NamePy().data(), atype.data(), consts::NamePz().data(),atype.data(),consts::NameM().data());
      _type_comps[atype][consts::P3Components()] = Form("%s%s,%s%s,%s%s", atype.data(),consts::NamePx().data(), atype.data(),consts::NamePy().data(), atype.data(), consts::NamePz().data());
      _types.push_back(atype);
    }
    inline void ConfigReaction::ValidateType(const std::string& type) const {
      if (std::find(_types.begin(), _types.end(), type) == _types.end()) throw std::invalid_argument("Error: Data type '" + type + "' is not registered. "+_types.size());
    }
  inline bool ConfigReaction::TypeExists(const std::string& type) const {
      if (std::find(_types.begin(), _types.end(), type) == _types.end()) return false;
      return true;
    }
    inline ROOT::RVec<std::string> ConfigReaction::GetTypes() const { return _types; }

} // namespace rad

#include "Diagnostics_ConfigReaction.hxx"
