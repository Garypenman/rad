#pragma once

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

// Forward declaration of utility functions, if needed
// using rad::names::data_type::Rec; 
// using rad::names::data_type::Tru;

namespace rad {
  namespace config {

    using ROOT::RVecI;
    using RDFstep = ROOT::RDF::RNode;
    using std::string;
    using std::string_view;

    // Helper to convert string_view to string for RDataFrame
    std::string as_string(std::string_view v) { 
      return {v.data(), v.size()}; 
    }

    //string to append to columns to stop them being snapshotted
    const std::string DoNotWriteTag(){return "__dnwtag";};

    //---------------------------------------------------------
    // RDataFrame Interface Base Class
    //---------------------------------------------------------

    class RDFInterface {

    public:
      // Constructors
      RDFInterface(const string_view treeName, const string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns) 
	: _orig_df{treeName, {fileNameGlob.data()}, columns}, _curr_df{_orig_df}, _base_df{_orig_df}, _treeName{as_string(treeName)}, _fileName{as_string(fileNameGlob)} {
	if (fileNameGlob.empty()) {
	  throw std::invalid_argument("RDFInterface: fileNameGlob cannot be empty.");
	}
	_orig_col_names = _orig_df.GetColumnNames();
      }

      RDFInterface(const string_view treeName, const std::vector<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns) 
	: _orig_df{treeName, filenames, columns}, _curr_df{_orig_df}, _base_df{_orig_df}, _treeName{as_string(treeName)}, _fileNames{filenames} {
	if (filenames.empty()) {
	  throw std::invalid_argument("RDFInterface: filenames list cannot be empty.");
	}
	_orig_col_names = _orig_df.GetColumnNames();
      }

      // if creating from alternative data source
      RDFInterface(ROOT::RDataFrame rdf) : _orig_df{rdf}, _curr_df{rdf}, _base_df{rdf} {
	_orig_col_names = _orig_df.GetColumnNames();
      }

      // Destructor handles lazy snapshots
      virtual ~RDFInterface() { 
	for (auto& trigger : _triggerSnapshots) {
	  if (trigger) trigger();
	}
      }

      //------------------ RDataFrame State Access ------------------

      RDFstep CurrFrame() { return _curr_df; }
      void setCurrFrame(RDFstep df) { _curr_df = df; }
      RDFstep getBaseFrame() const { return _base_df; }
      void setBaseFrame(RDFstep step) { _base_df = step; }
      void setMyBaseFrame() { _base_df = CurrFrame(); }
      RDFstep getOrigFrame() const { return _orig_df; }

      //------------------ RDataFrame Actions ------------------
        
      // Define (string expression)
      void Define(const string_view name, const string& expression) {
	setCurrFrame(CurrFrame().Define(name, expression));
      }
        
      // Define (lambda/functor)
      template<typename Lambda>
      void Define(const string_view name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
	setCurrFrame(CurrFrame().Define(name, func, columns));
      }

      // Redefine (string expression)
      void Redefine(const string& name, const string& expression) {
	setCurrFrame(CurrFrame().Redefine(name, expression));
      }

      // Redefine (lambda/functor)
      template<typename Lambda>
      void Redefine(const string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns = {}) {
	setCurrFrame(CurrFrame().Redefine(name, func, columns));
      }
	
      /**
       * Interface to RDataFrame Redefine via any aliases that may be used
       */
      //The redefine with alias does not work as the branchname
      //is not a mathematical expression. This is not checked for
      //the Lambda version below and works.
      // void RedefineViaAlias(const string& alias,const string& expression){
      // 	RedefineExpr(_aliasMap[alias],expression);
      // }
      template<typename Lambda>
      void RedefineViaAlias(const string& alias,Lambda&& func,const ROOT::RDF::ColumnNames_t& columns ){
	//	Redefine(_aliasMap[alias],func,columns);
	  
	auto it = _aliasMap.find(alias);
	if (it == _aliasMap.end()) {
	  throw std::invalid_argument("RedefineViaAlias: alias '" + alias + "' does not exist in _aliasMap.");
	}
	Redefine(it->second, std::forward<Lambda>(func), columns);
      }
 
      // Filter (string expression)
      void Filter(const std::string& expression, const std::string& name = "") {
	setCurrFrame(CurrFrame().Filter(expression, name));
      }

      // Filter (lambda/functor)
      template<typename Lambda>
      void Filter(Lambda&& func, const ROOT::RDF::ColumnNames_t& columns = {}, std::string name = "") {
	setCurrFrame(CurrFrame().Filter(func, columns, name));
      }
        
      // Alias
      //  virtual void setBranchAlias(const string& old_name, const string& new_name) = 0; // MUST be implemented by derived class

      //------------------ Snapshot ------------------
        
      virtual void Snapshot(const string& filename) = 0; // MUST be implemented by derived class (to handle aliases/maps)
      virtual void BookLazySnapshot(const string& filename) = 0; // MUST be implemented by derived class

      // Virtual hook for column cleanup before snapshot
      virtual void RemoveSnapshotColumns(std::vector<string>& cols) {
	// Only removes internal RDFInterface columns if any were introduced here.
	// Derived class (ConfigReaction) will handle its own cleanup.
      }

      //------------------ Getters ------------------

      std::string GetTreeName() const { return _treeName; }
      std::string GetFileName() const { return _fileName; }
      std::vector<std::string> GetFileNames() const { return _fileNames; }
      bool OriginalColumnExists(const string& col) {
	return std::find(_orig_col_names.begin(), _orig_col_names.end(), col) != _orig_col_names.end();
      }
        
      // Utility method used by ConfigReaction
      bool ColumnExists(const string& col, RDFstep df) {
	auto cols = df.GetDefinedColumnNames();
	return std::find(cols.begin(), cols.end(), col) != cols.end();
      }
	
      void setBranchAlias(const string& old_name, const string& new_name) {
	if (!OriginalColumnExists(old_name)) {
	  throw std::invalid_argument("setBranchAlias: Source column '" + old_name + "' does not exist in the DataFrame.");
	}
	_aliasMap[new_name] = old_name;
	setCurrFrame(CurrFrame().Alias(new_name, old_name));
      }
      /**
       * check if alias is used
       */
      bool CheckAlias(const string& alias){
	if(_aliasMap.find(alias) != _aliasMap.end()) return true;
	else return false;
      }

      const std::map<string,string>& AliasMap() const {return _aliasMap;}

      string ColObjTypeString(const string& name){return CurrFrame().GetColumnType(name);}
	
    protected:
      ROOT::RDataFrame _orig_df;
      RDFstep _curr_df;
      RDFstep _base_df;
        
      std::vector<std::function<void()>> _triggerSnapshots; // for BookLazySnapshot
        
      std::vector<std::string> _fileNames;
      std::string _fileName;
      std::string _treeName;
      ROOT::RDF::ColumnNames_t _orig_col_names;

      std::map<string, string> _aliasMap;
 
    }; // class RDFInterface


    //Helpers for column vector run-time type deduction
    enum class ColType{Undef,Int,UInt,Float,Double,Short,Bool,Long};
  
     ColType DeduceColumnVectorType(RDFInterface* const radf,const string& name){

      TString col_type = radf->ColObjTypeString(name);
       
      if(col_type.Contains("UInt_t")||col_type.Contains("uint")) return ColType::UInt;
      if(col_type.Contains("Float_t")||col_type.Contains("float")) return ColType::Float;
      if(col_type.Contains("Double_t")||col_type.Contains("double")) return ColType::Double;
      if(col_type.Contains("Short_t")||col_type.Contains("short")) return ColType::Short;
      if(col_type.Contains("Bool_t")||col_type.Contains("bool")) return ColType::Bool;
      if(col_type.Contains("Long_t")||col_type.Contains("long")) return ColType::Long;
      if(col_type.Contains("Int_t")||col_type.Contains("int")) return ColType::Int;

      throw std::runtime_error(std::string("RDFInterface : DeduceColumnVectorType cannot deduce a type for ") + 
			       name + 
			       " which is " + 
			       col_type.Data()); 
      return ColType::Undef;
    }
  
  } // namespace config
} // namespace rad
