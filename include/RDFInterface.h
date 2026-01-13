#pragma once

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <vector>
#include <string>
#include <map>
#include <functional>
#include <stdexcept>
#include <iostream>

namespace rad {

    using ROOT::RVecI;
    using RDFstep = ROOT::RDF::RNode;
    using std::string;
    using std::string_view;

    // --- Helpers ---
    inline std::string as_string(std::string_view v) { 
      return {v.data(), v.size()}; 
    }

    inline const std::string DoNotWriteTag(){ return "__dnwtag"; };

    // =========================================================================
    // RDataFrame Interface Base Class
    // =========================================================================

    /**
     * @class RDFInterface
     * @brief A stateful wrapper around ROOT::RDataFrame.
     * * @details
     * This class manages the chain of RDataFrame operations. Unlike raw RDataFrame, which 
     * returns a new node for every operation, RDFInterface maintains the state of the 
     * "Current Frame" (_curr_df).
     */
    class RDFInterface {

    public:
      // --- Constructors & Destructor ---
      RDFInterface(const string_view treeName, const string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns);
      RDFInterface(const string_view treeName, const std::vector<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns);
      RDFInterface(ROOT::RDataFrame rdf);
      
      virtual ~RDFInterface();

      // --- State Management ---
      RDFstep CurrFrame();
      void setCurrFrame(RDFstep df);
      
      RDFstep getBaseFrame() const;
      void setBaseFrame(RDFstep step);
      void setMyBaseFrame();
      
      RDFstep getOrigFrame() const;

      // --- RDataFrame Actions (Define, Filter, Redefine) ---
      
      void Define(const string_view name, const string& expression);
      
      template<typename Lambda>
      void Define(const string_view name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns);

      void RedefineExpr(const string& name, const string& expression);

      template<typename Lambda>
      void Redefine(const string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns = {});
    
  
      void Filter(const std::string& expression, const std::string& name = "");

      template<typename Lambda>
      void Filter(Lambda&& func, const ROOT::RDF::ColumnNames_t& columns = {}, std::string name = "");
        
      // --- Output & Snapshot ---
      
      virtual void Snapshot(const string& filename) = 0;
      virtual void BookLazySnapshot(const string& filename) = 0;
      virtual void RemoveSnapshotColumns(std::vector<string>& cols);

      // --- Metadata & Utilities ---

      std::string GetTreeName() const;
      std::string GetFileName() const;
      std::vector<std::string> GetFileNames() const;

      bool OriginalColumnExists(const string& col);
      bool ColumnExists(const string& col);
    
      string ColObjTypeString(const string& name);
    
    protected:
      ROOT::RDataFrame _orig_df;
      RDFstep _curr_df;
      RDFstep _base_df;
        
      std::vector<std::function<void()>> _triggerSnapshots;
        
      std::vector<std::string> _fileNames;
      std::string _fileName;
      std::string _treeName;
      ROOT::RDF::ColumnNames_t _orig_col_names;

      std::map<string, string> _aliasMap;
 
    }; // class RDFInterface


    // =========================================================================
    // Type Deduction Utilities
    // =========================================================================

    enum class ColType{ Undef, Int, UInt, Float, Double, Short, Bool, Long };
  
    inline ColType DeduceColumnVectorType(RDFInterface* const radf, const string& name){
      TString col_type = radf->ColObjTypeString(name);
       
      if(col_type.Contains("UInt_t")  || col_type.Contains("uint"))   return ColType::UInt;
      if(col_type.Contains("Float_t") || col_type.Contains("float"))  return ColType::Float;
      if(col_type.Contains("Double_t")|| col_type.Contains("double")) return ColType::Double;
      if(col_type.Contains("Short_t") || col_type.Contains("short"))  return ColType::Short;
      if(col_type.Contains("Bool_t")  || col_type.Contains("bool"))   return ColType::Bool;
      if(col_type.Contains("Long_t")  || col_type.Contains("long"))   return ColType::Long;
      if(col_type.Contains("Int_t")   || col_type.Contains("int"))    return ColType::Int;

      throw std::runtime_error(std::string("RDFInterface : DeduceColumnVectorType cannot deduce a type for ") + 
                   name + " which is " + col_type.Data()); 
      return ColType::Undef;
    }


    // =========================================================================
    // IMPLEMENTATION
    // =========================================================================

    // --- Constructors ---
    inline RDFInterface::RDFInterface(const string_view treeName, const string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns) 
    : _orig_df{treeName, {fileNameGlob.data()}, columns}, 
          _curr_df{_orig_df}, 
          _base_df{_orig_df}, 
          _treeName{as_string(treeName)}, 
          _fileName{as_string(fileNameGlob)} 
    {
      if (fileNameGlob.empty()) throw std::invalid_argument("RDFInterface: fileNameGlob cannot be empty.");
      _orig_col_names = _orig_df.GetColumnNames();
    }

    inline RDFInterface::RDFInterface(const string_view treeName, const std::vector<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns) 
    : _orig_df{treeName, filenames, columns}, 
          _curr_df{_orig_df}, 
          _base_df{_orig_df}, 
          _treeName{as_string(treeName)}, 
          _fileNames{filenames} 
    {
      if (filenames.empty()) throw std::invalid_argument("RDFInterface: filenames list cannot be empty.");
      _orig_col_names = _orig_df.GetColumnNames();
    }

    inline RDFInterface::RDFInterface(ROOT::RDataFrame rdf) 
          : _orig_df{rdf}, _curr_df{rdf}, _base_df{rdf} 
    {
      _orig_col_names = _orig_df.GetColumnNames();
    }

    inline RDFInterface::~RDFInterface() { 
      for (auto& trigger : _triggerSnapshots) {
        if (trigger) trigger();
      }
    }

    // --- State Management ---
    inline RDFstep RDFInterface::CurrFrame() { return _curr_df; }
    inline void RDFInterface::setCurrFrame(RDFstep df) { _curr_df = df; }
    
    inline RDFstep RDFInterface::getBaseFrame() const { return _base_df; }
    inline void RDFInterface::setBaseFrame(RDFstep step) { _base_df = step; }
    inline void RDFInterface::setMyBaseFrame() { _base_df = CurrFrame(); }
    inline RDFstep RDFInterface::getOrigFrame() const { return _orig_df; }

    // --- Actions ---
    inline void RDFInterface::Define(const string_view name, const string& expression) {
        setCurrFrame(CurrFrame().Define(name, expression));
    }
        
    template<typename Lambda>
    inline void RDFInterface::Define(const string_view name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
        setCurrFrame(CurrFrame().Define(name, func, columns));
    }

    inline void RDFInterface::RedefineExpr(const string& name, const string& expression) {
        setCurrFrame(CurrFrame().Redefine(name, expression));
    }

    template<typename Lambda>
    inline void RDFInterface::Redefine(const string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
        setCurrFrame(CurrFrame().Redefine(name, func, columns));
    }
    
   
    inline void RDFInterface::Filter(const std::string& expression, const std::string& name) {
        setCurrFrame(CurrFrame().Filter(expression, name));
    }

    template<typename Lambda>
    inline void RDFInterface::Filter(Lambda&& func, const ROOT::RDF::ColumnNames_t& columns, std::string name) {
        setCurrFrame(CurrFrame().Filter(func, columns, name));
    }

    // --- Utilities ---
    inline void RDFInterface::RemoveSnapshotColumns(std::vector<string>& cols) {
        // Base implementation does nothing, but virtual for override
    }

    inline std::string RDFInterface::GetTreeName() const { return _treeName; }
    inline std::string RDFInterface::GetFileName() const { return _fileName; }
    inline std::vector<std::string> RDFInterface::GetFileNames() const { return _fileNames; }

    inline bool RDFInterface::OriginalColumnExists(const string& col) {
        return std::find(_orig_col_names.begin(), _orig_col_names.end(), col) != _orig_col_names.end();
    }
        
    inline bool RDFInterface::ColumnExists(const string& col) {
        auto cols = CurrFrame().GetColumnNames();
        return std::find(cols.begin(), cols.end(), col) != cols.end();
    }
    
   
    inline string RDFInterface::ColObjTypeString(const string& name){ return CurrFrame().GetColumnType(name); }
  
} // namespace rad
