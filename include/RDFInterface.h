/**
 * @file RDFInterface.h
 * @brief Stateful wrapper around ROOT::RDataFrame to manage analysis chains.
 * @details
 * This class encapsulates a ROOT::RDataFrame workflow. Unlike the raw RDataFrame,
 * which returns a new node for every transformation, RDFInterface maintains the
 * state of the "Current Frame". This allows users to chain operations (Define, Filter)
 * naturally without manually managing intermediate node objects.
 * * Key Features:
 * - **State Management:** Tracks the current node (`_curr_df`) and the original base node.
 * - **Lazy Snapshots:** Supports deferred triggering of snapshots (execution actions).
 * - **SnapshotCombi:** Advanced flattening of combinatorial data into N-Tuples.
 */

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
#include <algorithm>
#include <type_traits>

// Include the Action for SnapshotCombi
#include "RDFUtils.h" 
#include "SnapshotCombi.h" 

namespace rad {

    // Shorthand aliases
    using RVecI = ROOT::RVec<int>;
    using RDFstep = ROOT::RDF::RNode;
    using std::string;
    using std::string_view;

    /**
     * @brief Helper: Convert string_view to string safely.
     */
    inline std::string as_string(std::string_view v) { 
        return {v.data(), v.size()}; 
    }

    /**
     * @brief Tag used to mark columns that should be excluded from automatic writing.
     */
    inline const std::string DoNotWriteTag() { return "__dnwtag"; }


    // =========================================================================
    // CLASS DECLARATION
    // =========================================================================

    class RDFInterface {

    public:
        // --- Constructors & Destructor ---
        
        /**
         * @brief Constructs a workflow from a file glob (e.g. "*.root").
         * @param treeName Name of the TTree in the files.
         * @param fileNameGlob Glob pattern or single filename.
         * @param columns List of columns to read (optional).
         */
        RDFInterface(const string_view treeName, const string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns = {});

        /**
         * @brief Constructs a workflow from a list of specific filenames.
         * @param treeName Name of the TTree in the files.
         * @param filenames Vector of file paths.
         * @param columns List of columns to read (optional).
         */
        RDFInterface(const string_view treeName, const std::vector<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns = {});

        /**
         * @brief Wraps an existing RDataFrame object.
         */
        RDFInterface(ROOT::RDataFrame rdf);
        
         /**
         * @brief Wraps an existing RDataFrame object.
         */
        RDFInterface(ROOT::RDF::RNode rdf);
        
        /**
         * @brief Destructor. Executes any pending deferred actions (Lazy Snapshots).
         */
        virtual ~RDFInterface();

        // --- State Management ---
        
        /** @brief Returns the current RDataFrame node. */
        RDFstep CurrFrame();

        /** @brief Manually updates the current node. */
        void SetCurrFrame(RDFstep df);
        
        /** @brief Returns the base node (before any user modifications). */
        RDFstep GetBaseFrame() const;

        /** @brief Sets a specific node as the "Base". */
        void SetBaseFrame(RDFstep step);

        /** @brief Sets the *Current* frame as the new Base. */
        void SetMyBaseFrame();
        
        /** @brief Returns the absolutely original node (file loader). */
        RDFstep GetOrigFrame() const;

        // --- Transformations (Define, Filter) ---
        
        /**
         * @brief Defines a new column using a string expression (JIT).
         * @param name Name of the new column.
         * @param expression C++ expression string (e.g. "x*y").
         */
        void Define(const string_view name, const string& expression);
        
        /**
         * @brief Defines a new column using a C++ function/lambda.
         * @param name Name of the new column.
         * @param func Function to execute.
         * @param columns List of input column names.
         */
        template<typename Lambda>
        void Define(const string_view name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns);

        /**
         * @brief Redefines an existing column (String JIT).
         */
        void RedefineExpr(const string& name, const string& expression);

        /**
         * @brief Redefines an existing column (Lambda).
         */
        template<typename Lambda>
        void Redefine(const string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns = {});
    
        /**
         * @brief Filters the dataframe (String JIT).
         * @param expression Boolean expression string.
         * @param name Optional name for the filter node (for stats).
         */
        void Filter(const std::string& expression, const std::string& name = "");

        /**
         * @brief Filters the dataframe (Lambda).
         */
        template<typename Lambda>
        void Filter(Lambda&& func, const ROOT::RDF::ColumnNames_t& columns = {}, std::string name = "");
        
        // --- Output & Snapshots ---
        
        /** @brief Immediate Snapshot (Virtual, implementation specific). */
        virtual void Snapshot(const string& filename) = 0;
        
        /** @brief Lazy Snapshot (Virtual, implementation specific). */
        virtual void BookLazySnapshot(const string& filename) = 0;
        
        /** @brief Helper to remove specific columns from snapshot list. */
        virtual void RemoveSnapshotColumns(std::vector<string>& cols);

        /**
         * @brief Schedules a "Combi Snapshot" (Flat Tree).
         * @details 
         * Creates a TTree where 1 Entry = 1 Good Combination.
         * - **Vectors** are flattened (element at combination index).
         * - **Scalars** are broadcast (repeated).
         * Uses JIT compilation to deduce column types automatically.
         * Execution is deferred to the destructor (Lazy).
         * * @param filename Output file name.
         * @param treename Output tree name.
         * @param columns List of columns to save.
         * @param maskCol Name of the index mask column (from PhysicsSelection).
         */
        void BookSnapshotCombi(const std::string& filename, const std::string& treename,
                               const std::vector<std::string>& columns,
                               const std::string& maskCol);

      /** * @brief Manually triggers the event loop for all booked snapshots.
       * @details Call this in your macro after booking snapshots to ensure they write to disk immediately.
       */
      void TriggerSnapshots();
      
      /** * @brief Clears pending snapshot triggers. 
       * @details Essential when cloning a reaction to avoid duplicating the parent's output actions.
       */
      void ClearTriggers() { 
        _triggerSnapshots.clear(); 
      }

      // --- Metadata & Utilities ---

        std::string GetTreeName() const;
        std::string GetFileName() const;
        std::vector<std::string> GetFileNames() const;

        bool OriginalColumnExists(const string& col);
        bool ColumnExists(const string& col);
        
        /** @brief Returns the C++ type string of a column (e.g. "double", "vector<int>"). */
        string ColObjTypeString(const string& name);
        
    protected:
        ROOT::RDataFrame _orig_df;
        RDFstep _curr_df;
        RDFstep _base_df;
        
        // List of functions to execute on destruction (deferred snapshots)
        std::vector<std::function<void()>> _triggerSnapshots;
        
        std::vector<std::string> _fileNames;
        std::string _fileName;
        std::string _treeName;
        ROOT::RDF::ColumnNames_t _orig_col_names;

        std::map<string, string> _aliasMap;

    }; // class RDFInterface


   


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
   inline RDFInterface::RDFInterface(ROOT::RDF::RNode rdf) 
     : _orig_df(0), _curr_df{rdf}, _base_df{rdf} 
    {
        _orig_col_names = _orig_df.GetColumnNames();
    }

    inline RDFInterface::~RDFInterface() { 
      TriggerSnapshots();
    }

    // --- State Management ---

    inline RDFstep RDFInterface::CurrFrame() { return _curr_df; }
    inline void RDFInterface::SetCurrFrame(RDFstep df) { _curr_df = df; }
    
    inline RDFstep RDFInterface::GetBaseFrame() const { return _base_df; }
    inline void RDFInterface::SetBaseFrame(RDFstep step) { _base_df = step; }
    inline void RDFInterface::SetMyBaseFrame() { _base_df = CurrFrame(); }
    inline RDFstep RDFInterface::GetOrigFrame() const { return _orig_df; }

    // --- Transformations ---

    inline void RDFInterface::Define(const string_view name, const string& expression) {
        SetCurrFrame(CurrFrame().Define(name, expression));
    }
        
    template<typename Lambda>
    inline void RDFInterface::Define(const string_view name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
        SetCurrFrame(CurrFrame().Define(name, func, columns));
    }

    inline void RDFInterface::RedefineExpr(const string& name, const string& expression) {
        SetCurrFrame(CurrFrame().Redefine(name, expression));
    }

    template<typename Lambda>
    inline void RDFInterface::Redefine(const string& name, Lambda&& func, const ROOT::RDF::ColumnNames_t& columns) {
        SetCurrFrame(CurrFrame().Redefine(name, func, columns));
    }
    
    inline void RDFInterface::Filter(const std::string& expression, const std::string& name) {
        SetCurrFrame(CurrFrame().Filter(expression, name));
    }

    template<typename Lambda>
    inline void RDFInterface::Filter(Lambda&& func, const ROOT::RDF::ColumnNames_t& columns, std::string name) {
        SetCurrFrame(CurrFrame().Filter(func, columns, name));
    }

    // --- Output & Snapshots ---

    inline void RDFInterface::RemoveSnapshotColumns(std::vector<string>& cols) {
        // Virtual hook: default implementation does nothing.
    }

  inline void RDFInterface::BookSnapshotCombi(const std::string& filename, const std::string& treename,
                                                const std::vector<std::string>& columns,
                                                const std::string& maskCol) 
    {
        // 1. Deduce Types
        std::vector<ColType> types;
        types.reserve(columns.size());
        
        for(const auto& col : columns) {
            types.push_back(DeduceColumnVectorType(this, col));
        }

        // 2. Prepare Cols
        std::vector<std::string> all_cols = columns;
        all_cols.push_back(maskCol);

        // 3. Book Action
        rad::io::SnapshotCombi action(filename, treename, columns, types);
        auto action_node = CurrFrame().Book(std::move(action), all_cols);

        // 4. Register Trigger
        // Dereferencing the RResultPtr (now valid!) triggers the event loop.
        _triggerSnapshots.push_back([snapshot = std::move(action_node)]()  mutable { auto trigger = *snapshot;});

	//auto snapshot_result = final_df.Snapshot("rad_tree", filename, cols, opts);
	//_triggerSnapshots.push_back([snapshot = std::move(snapshot_result)]() mutable {});

    }
  /** * @brief Manually triggers the event loop for all booked snapshots.
     * @details Call this in your macro after booking snapshots to ensure they write to disk immediately.
     */
    void RDFInterface::TriggerSnapshots() {
        if (_triggerSnapshots.empty()) return;

        std::cout << "[RDFInterface] Triggering " << _triggerSnapshots.size() << " snapshot(s)..." << std::endl;
        
        // Execute all registered triggers
        for (auto& trigger : _triggerSnapshots) {
            if (trigger) trigger(); 
        }
        
        // Clear them so we don't trigger them again in the destructor
        _triggerSnapshots.clear();
        std::cout << "[RDFInterface] Snapshots complete." << std::endl;
    }
  
    // --- Metadata & Utilities ---

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
    
    inline string RDFInterface::ColObjTypeString(const string& name){ 
        return CurrFrame().GetColumnType(name); 
    }

} // namespace rad
