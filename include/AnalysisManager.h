/**
 * @file AnalysisManager.h
 * @brief High-level orchestration of the RAD analysis chain.
 * @details
 * The AnalysisManager is the central entry point for the user. It manages:
 * 1. Data Streams (Rec, Truth, Systematics).
 * 2. Configuration Recipes (Kinematics, Cuts, Histograms).
 * 3. Execution phases (Initialization, Snapshotting, Running).
 * * It uses a "Three-Pass Initialization" strategy to handle dependencies:
 * - Pass 1: Define Kinematics (Variables).
 * - Pass 2: Compile Selections (Cuts depending on Variables).
 * - Pass 3: Book Histograms (Plots depending on Variables and Cuts).
 */

#pragma once

#include "ConfigReaction.h"
#include "KinematicsProcessor.h"
#include "PhysicsSelection.h"
#include "Histogrammer.h"
#include "CommonDefines.h"
#include "Indicing.h"
#include "Diagnostics.h" // Added for DiagnosticsPrinter

#include <memory>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <iostream>
#include <filesystem>
#include <functional>

namespace rad {

  /**
   * @class AnalysisManager
   * @brief The main driver class for RAD analysis.
   * @tparam ReactionClass The concrete reaction class (e.g. ePICReaction).
   * @tparam ProcessorClass The concrete processor class (e.g. KinematicsProcessor).
   */
  template <typename ReactionClass, typename ProcessorClass>
  class AnalysisManager {

  public:
    // Recipe Signatures
    using KineRecipe  = std::function<void(ProcessorClass&)>;
    using SelRecipe   = std::function<void(PhysicsSelection&)>;
    using HistoRecipe = std::function<void(histo::Histogrammer&)>;

    // =====================================================================
    // Setup
    // =====================================================================
    
    /**
     * @brief Constructor.
     * @param name Name of the analysis (used for output filenames).
     * @param treeName Input TTree name.
     * @param fileGlob Input file path/pattern.
     */
    AnalysisManager(const std::string& name, const std::string& treeName, const std::string& fileGlob);

    /** @brief Sets and creates the output directory. */
    void SetOutputDir(const std::string& dir);

    /**
     * @brief Define which data types to process.
     * @details The **First** type provided is treated as the **Primary Stream**.
     * The Primary Stream determines the Event Selection Mask for the Snapshot.
     * @param types Variadic list of stream prefixes (e.g. Rec(), Truth()).
     */
    template<typename... Args>
    void SetTypes(Args... types);

    /** @return Reference to the underlying Reaction object. */
    ReactionClass& Reaction();

    // =====================================================================
    // Configuration
    // =====================================================================

    /** @brief Apply a Kinematics recipe to ALL active streams. */
    void ConfigureKinematics(KineRecipe recipe);
    /** @brief Apply a Kinematics recipe to a SPECIFIC stream only. */
    void ConfigureKinematics(const std::string& streamName, KineRecipe recipe);

    /** @brief Apply a Selection recipe to ALL active streams. */
    void ConfigureSelection(SelRecipe recipe);
    /** @brief Apply a Selection recipe to a SPECIFIC stream only. */
    void ConfigureSelection(const std::string& streamName, SelRecipe recipe);

    /** @brief Apply a Histogram recipe to ALL active streams. */
    void ConfigureHistograms(HistoRecipe recipe);
    /** @brief Apply a Histogram recipe to a SPECIFIC stream only. */
    void ConfigureHistograms(const std::string& streamName, HistoRecipe recipe);

    // =====================================================================
    // Execution
    // =====================================================================

    /**
     * @brief Initializes the entire analysis chain.
     * @details Uses a 3-Pass system to respect dependencies:
     * 1. Variables (Kinematics) - Creates columns.
     * 2. Cuts (Selection) - Creates masks using variables.
     * 3. Histograms - Books actions using variables and masks.
     */
    void Init();

    /**
     * @brief Snapshots data to a SINGLE flat TTree.
     * @details 
     * - Merges columns from all streams (Rec + Truth).
     * - Uses the Primary Stream's mask to filter events.
     * - Uses the Primary Stream's combination index for alignment.
     * - Filename defaults to: {OutputDir}/{AnalysisName}_Tree.root
     * @param addCols Additional columns to save manually.
     * @param filename Optional override for the output filename.
     */
    void Snapshot(const ROOT::RDF::ColumnNames_t& addCols={}, const std::string& filename = "");

    /**
     * @brief Runs the analysis and writes histograms.
     * @details 
     * - Creates ONE file PER stream to avoid key clashes.
     * - Filename: {OutputDir}/{AnalysisName}_{Prefix}_{Suffix}
     * - Example:  output/Y4260_rec_Hist.root
     * @param suffix Suffix for the histogram file (default "Hist.root").
     */
    void Run(const std::string& suffix = "Hist.root");

    /**
     * @brief Print comprehensive diagnostics for the entire analysis setup.
     * @details Useful for debugging mapping, combinatorial, and observable issues.
     */
    void PrintDiagnostics() const;

  private:
    /**
     * @struct AnalysisStream
     * @brief Holds the components for a single data stream (e.g., "rec_" or "tru_").
     */
    struct AnalysisStream {
        std::string prefix;
        std::unique_ptr<ProcessorClass> kine;
        std::unique_ptr<PhysicsSelection> sel;
        std::unique_ptr<histo::Histogrammer> hist;

        AnalysisStream(ReactionClass* reaction, const std::string& p) 
            : prefix(p) 
        {
            kine = std::make_unique<ProcessorClass>(reaction, prefix);
            sel  = std::make_unique<PhysicsSelection>(*kine);
            hist = std::make_unique<histo::Histogrammer>(*kine, sel.get());
        }
    };

    ReactionClass _reaction;
    std::string _name;
    std::string _outputDir;
    bool _initialized = false;
    std::string _primaryStream;
    std::map<std::string, AnalysisStream> _streams;

    void AddStream(const std::string& prefix);
    bool CheckStream(const std::string& name);
    std::string MakePath(const std::string& filename);
    ROOT::RDF::ColumnNames_t CollectStreamColumns(const ProcessorClass& kine);
  };

  // ===========================================================================
  // IMPLEMENTATION
  // ===========================================================================

  template <typename R, typename P>
  inline AnalysisManager<R,P>::AnalysisManager(const std::string& name, const std::string& treeName, const std::string& fileGlob) 
      : _reaction(treeName, fileGlob), _name(name) {}

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::SetOutputDir(const std::string& dir) {
      _outputDir = dir;
      if (!_outputDir.empty() && !std::filesystem::exists(_outputDir)) {
          std::filesystem::create_directories(_outputDir);
      }
  }

  template <typename R, typename P>
  template<typename... Args>
  inline void AnalysisManager<R,P>::SetTypes(Args... types) {
      (AddStream(types), ...);
  }

  template <typename R, typename P>
  inline R& AnalysisManager<R,P>::Reaction() { return _reaction; }

  // --- Configuration ---

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureKinematics(KineRecipe recipe) {
      for(auto& [key, stream] : _streams) { if(stream.kine) recipe(*stream.kine); }
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureKinematics(const std::string& streamName, KineRecipe recipe) {
      if(CheckStream(streamName)) recipe(*_streams.at(streamName).kine);
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureSelection(SelRecipe recipe) {
      for(auto& [key, stream] : _streams) { if(stream.sel) recipe(*stream.sel); }
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureSelection(const std::string& streamName, SelRecipe recipe) {
      if(CheckStream(streamName)) recipe(*_streams.at(streamName).sel);
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureHistograms(HistoRecipe recipe) {
      for(auto& [key, stream] : _streams) { if(stream.hist) recipe(*stream.hist); }
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureHistograms(const std::string& streamName, HistoRecipe recipe) {
      if(CheckStream(streamName)) recipe(*_streams.at(streamName).hist);
  }

  // --- Execution ---

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::Init() {
      if(_initialized) return;
      if(_primaryStream.empty()) throw std::runtime_error("[AnalysisManager] No streams defined! Call SetTypes().");

      // PASS 1: Initialize Kinematics (Create Variables)
      for(auto& [key, stream] : _streams) {
          stream.kine->Init();
      }

      // PASS 2: Compile Selections (Create Masks)
      for(auto& [key, stream] : _streams) {
          if(stream.sel) stream.sel->Init();
      }

      // PASS 3: Initialize Histograms (Book Actions)
      for(auto& [key, stream] : _streams) {
          if(stream.hist) stream.hist->Init();
      }

      _initialized = true;
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::Snapshot(const ROOT::RDF::ColumnNames_t& addCols, const std::string& filename) {
      Init();
      
      std::string finalName = filename.empty() ? _name + "_Tree.root" : filename;
      std::string finalPath = MakePath(finalName);

      auto cols = addCols;
      
      // 1. Collect Primary Stream (Rec)
      if(_streams.count(_primaryStream)) {
          auto c = CollectStreamColumns(*_streams.at(_primaryStream).kine);
          cols.insert(cols.end(), c.begin(), c.end());
      }

      // 2. Collect Secondary Streams (Truth)
      for(auto& [key, stream] : _streams) {
          if(key == _primaryStream) continue; 
          auto c = CollectStreamColumns(*stream.kine);
          cols.insert(cols.end(), c.begin(), c.end());
      }

      // 3. Determine Mask (Primary Driven)
      std::string mask = "";
      auto& prime = _streams.at(_primaryStream);
      mask = prime.sel->GetMaskColumn();
      
      // Fallback: If no mask is defined (no cuts), create a "Pass All" mask
      if(mask.empty()) {
           auto pNames = prime.kine->Creator().GetParticleNames();
           if(!pNames.empty()) {
               mask = prime.prefix + "Analysis_AllIndices";
               // We need a reference column to know how many combinations exist (size)
               std::string ref = prime.prefix + pNames[0] + "_px";
               _reaction.Define(mask, [](const ROOT::RVecD& r){ 
                  return rad::util::EnumerateIndicesFrom(0, r.size()); 
               }, {ref});
           }
      }

      std::cout << "[AnalysisManager] Snapshotting combined tree to " << finalPath << std::endl;
      _reaction.BookSnapshotCombi(finalPath, "tree", cols, mask);
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::Run(const std::string& suffix) {
      Init();
      // If Snapshot() wasn't called, this triggers the event loop now.
      _reaction.TriggerSnapshots();
      
      if(!_streams.empty()) {
          for(auto& [key, stream] : _streams) {
              std::string fname = _name + "_" + stream.prefix + suffix;
              std::string finalPath = MakePath(fname);
              
              stream.hist->File(finalPath, "RECREATE");
              std::cout << "[AnalysisManager] Wrote " << stream.prefix << " histograms to " << finalPath << std::endl;
          }
      }
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::PrintDiagnostics() const {
    diag::DiagnosticsPrinter::PrintSectionHeader("ANALYSIS MANAGER DIAGNOSTICS", '=', 90);

    // Call PrintReactionDiagnostics on the Reaction object
    // Note: This assumes ReactionClass has this method.
    // If _reaction is const in this context, the method must be const.
    // We cast away constness if needed, but ideally PrintReactionDiagnostics is const.
    const_cast<R&>(_reaction).PrintReactionDiagnostics();

    diag::DiagnosticsPrinter::PrintBlank();
    std::cout << "Registered Streams: " << _streams.size() << std::endl;

    for (const auto& entry : _streams) {
      std::cout << "\nStream: " << entry.first << std::endl;
      
      // FIXED: Use .kine instead of .processor
      if(entry.second.kine) {
         entry.second.kine->PrintProcessorDiagnostics();
      }
    }
  }

  // --- Private Helpers ---

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::AddStream(const std::string& prefix) {
      if(_streams.find(prefix) == _streams.end()) {
          _streams.try_emplace(prefix, &_reaction, prefix);
          if(_primaryStream.empty()) _primaryStream = prefix; 
      }
  }

  template <typename R, typename P>
  inline bool AnalysisManager<R,P>::CheckStream(const std::string& name) {
      if(_streams.find(name) == _streams.end()) {
          std::cerr << "AnalysisManager Warning: Stream '" << name << "' not found. Recipe ignored." << std::endl;
          return false;
      }
      return true;
  }

  template <typename R, typename P>
  inline std::string AnalysisManager<R,P>::MakePath(const std::string& filename) {
      if(!_outputDir.empty()) {
          return (std::filesystem::path(_outputDir) / filename).string();
      }
      return filename;
  }

  template <typename R, typename P>
  inline ROOT::RDF::ColumnNames_t AnalysisManager<R,P>::CollectStreamColumns(const P& kine) {
      ROOT::RDF::ColumnNames_t cols;
      
      for(const auto& var : kine.GetDefinedNames()) {
          cols.push_back(kine.GetPrefix() + var + kine.GetSuffix());
      }

      // Add Signal Flag if it exists
      std::string sigCol = kine.GetPrefix() + rad::consts::TruthMatchedCombi();
      if(_reaction.ColumnExists(sigCol)) {
          cols.push_back(sigCol);
      }
      
      return cols;
  }

} // namespace rad
