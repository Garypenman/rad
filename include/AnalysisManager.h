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
#include "Diagnostics.h"

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

    // =====================================================================
    // Stream Management
    // =====================================================================

    /**
     * @brief Adds a specific analysis stream.
     * @details 
     * - Example: AddStream("rec", "loose") -> Creates stream "rec_loose".
     * - Variables will be named "rec_Mass_loose".
     * - Input data will be read from "rec_" columns.
     * @param dataSource The data prefix to read (e.g. "rec", "tru").
     * @param label Analysis label (e.g. "loose", "tight", "sysUp"). Leave empty for default.
     */
    void AddStream(const std::string& dataSource, const std::string& label = "");

    /**
     * @brief Shortcut to add standard default streams (label="").
     * @details SetTypes("rec", "tru") creates streams named "rec" and "tru".
     * @param types Variadic list of stream prefixes (e.g. Rec(), Truth()).
     */
    //  template<typename... Args>
    //  void SetTypes(const std::string& name,Args... types);

    /** @return Reference to the underlying Reaction object. */
    ReactionClass& Reaction();

    // =====================================================================
    // Configuration (Pattern Matching Enabled)
    // =====================================================================

    /** @brief Apply a Kinematics recipe to ALL active streams. */
    void ConfigureKinematics(KineRecipe recipe);
    /** * @brief Apply a Kinematics recipe to streams matching the pattern. 
     * @param pattern Can match the Full Name ("rec_loose"), the Source ("rec"), or the Label ("loose").
     * Exact matches only (e.g. "loose" will NOT match "loose2").
     */
    void ConfigureKinematics(const std::string& pattern, KineRecipe recipe);

    /** @brief Apply a Selection recipe to ALL active streams. */
    void ConfigureSelection(SelRecipe recipe);
    /** @brief Apply a Selection recipe to streams matching the pattern. */
    void ConfigureSelection(const std::string& pattern, SelRecipe recipe);

    /** @brief Apply a Histogram recipe to ALL active streams. */
    void ConfigureHistograms(HistoRecipe recipe);
    /** @brief Apply a Histogram recipe to streams matching the pattern. */
    void ConfigureHistograms(const std::string& pattern, HistoRecipe recipe);

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
     * @brief Snapshots data to separate flat TTrees (One per Stream).
     * @details 
     * - Iterates over ALL active streams.
     * - Creates file: {OutputDir}/{AnalysisName}_{StreamName}_Tree.root
     * - Uses that specific stream's Mask and Variables.
     * - This ensures 'loose' and 'tight' streams are stored safely with correct event counts.
     * @param addCols Additional columns to save (e.g. Truth Scalars) in ALL trees.
     * @param filenameBase Optional override for base filename (e.g. "MyOutput.root" -> "MyOutput_rec_0.root").
     */

    void Snapshot(const ROOT::RDF::ColumnNames_t& addCols={}, const std::string& filename = "");
  
    /**
     * @brief Runs the analysis and writes histograms.
     * @details 
     * - Creates ONE file PER stream to avoid key clashes.
     * - Filename: {OutputDir}/{AnalysisName}_{Prefix}_{Suffix}
     * - Example:  output/Y4260_rec_loose_Hist.root
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
     * @brief Holds the components for a single analysis pipeline.
     */
    struct AnalysisStream {
        std::string fullName;   // Unique ID (e.g. "rec_loose")
        std::string source;     // Data Source (e.g. "rec")
        std::string label;      // Variation Label (e.g. "loose")
        
        std::unique_ptr<ProcessorClass> kine;
        std::unique_ptr<PhysicsSelection> sel;
        std::unique_ptr<histo::Histogrammer> hist;

        AnalysisStream(ReactionClass* reaction, const std::string& src, const std::string& lbl) 
            : source(src), label(lbl) 
        {
            // Logic: 
            // 1. FullName = "rec" + "_" + "loose" (if label exists)
            // 2. InputPrefix = "rec_" (Must have underscore to find data)
            // 3. OutputSuffix = "_loose" (If label exists)
            
            fullName = src;
            std::string outputSuffix = "";
            if(!lbl.empty()) {
                fullName +=  lbl;
                outputSuffix = "_" + lbl;
            }
            
            // Ensure input prefix ends in "_" (e.g. "rec" -> "rec_")
            std::string inputPrefix = src;
            if(inputPrefix.back() != '_') inputPrefix += "_";

            kine = std::make_unique<ProcessorClass>(reaction, inputPrefix, outputSuffix);
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

    // Helper to check if a stream matches a pattern (Exact matching on Source, Label, or FullName)
    bool StreamMatches(const AnalysisStream& stream, const std::string& pattern);
    
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

  // --- Stream Management ---

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::AddStream(const std::string& dataSource, const std::string& label) {
      // Logic handled in struct constructor, just construct key here
      std::string key = dataSource;
      if(!label.empty()) key += "_" + label;

      if(_streams.find(key) != _streams.end()) {
          std::cerr << "AnalysisManager Warning: Stream '" << key << "' already exists." << std::endl;
          return;
      }

      _streams.try_emplace(key, &_reaction, dataSource, label);

      // First added stream becomes default Primary
      if(_primaryStream.empty()) _primaryStream = key;
  }

  // template <typename R, typename P>
  // template<typename... Args>
  // inline void AnalysisManager<R,P>::SetTypes(const std::string& name,Args... types) {
  //     // Just add basic streams with empty labels
  //     (AddStream(types, name), ...);
  // }

  template <typename R, typename P>
  inline R& AnalysisManager<R,P>::Reaction() { return _reaction; }

  // --- Configuration ---

  template <typename R, typename P>
  inline bool AnalysisManager<R,P>::StreamMatches(const AnalysisStream& stream, const std::string& pattern) {
      // Exact Match Logic prevents "loose" from matching "loose2"
      if(pattern == stream.fullName) return true; // Matches "rec_loose"
      if(pattern == stream.source)   return true; // Matches "rec"
      if(pattern == stream.label)    return true; // Matches "loose"
      return false;
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureKinematics(KineRecipe recipe) {
      for(auto& [key, stream] : _streams) { if(stream.kine) recipe(*stream.kine); }
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureKinematics(const std::string& pattern, KineRecipe recipe) {
      bool found = false;
      for(auto& [key, stream] : _streams) {
          if(StreamMatches(stream, pattern) && stream.kine) {
              recipe(*stream.kine);
              found = true;
          }
      }
      if(!found) std::cerr << "AnalysisManager Warning: No streams matched pattern '" << pattern << "'" << std::endl;
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureSelection(SelRecipe recipe) {
      for(auto& [key, stream] : _streams) { if(stream.sel) recipe(*stream.sel); }
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureSelection(const std::string& pattern, SelRecipe recipe) {
      bool found = false;
      for(auto& [key, stream] : _streams) {
          if(StreamMatches(stream, pattern) && stream.sel) {
              recipe(*stream.sel);
              found = true;
          }
      }
      if(!found) std::cerr << "AnalysisManager Warning: No streams matched pattern '" << pattern << "'" << std::endl;
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureHistograms(HistoRecipe recipe) {
      for(auto& [key, stream] : _streams) { if(stream.hist) recipe(*stream.hist); }
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureHistograms(const std::string& pattern, HistoRecipe recipe) {
      bool found = false;
      for(auto& [key, stream] : _streams) {
          if(StreamMatches(stream, pattern) && stream.hist) {
              recipe(*stream.hist);
              found = true;
          }
      }
      if(!found) std::cerr << "AnalysisManager Warning: No streams matched pattern '" << pattern << "'" << std::endl;
  }

  // --- Execution ---

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::Init() {
      if(_initialized) return;
      if(_primaryStream.empty()) throw std::runtime_error("[AnalysisManager] No streams defined! Call SetTypes() or AddStream().");

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
  void AnalysisManager<R,P>::Snapshot(const ROOT::RDF::ColumnNames_t& addCols, const std::string& filenameBase) {
    
    using namespace rad::consts::data_type;//for Rec() Truth()
    
      Init();
      
      // 1. Pre-collect Truth Columns (if any)
      // We want to add these to ALL "Rec" trees (e.g. Q2_true, x_true)
      ROOT::RDF::ColumnNames_t globalTruthCols;
      for(auto& [key, stream] : _streams) {
          // Heuristic: If source starts with "tru", it's a truth stream
	if(stream.source.find(Truth()) == 0) {
              auto tCols = CollectStreamColumns(*stream.kine);
              globalTruthCols.insert(globalTruthCols.end(), tCols.begin(), tCols.end());
          }
      }

      // 2. Loop over ALL streams and book a separate tree for each
      for(auto& [key, stream] : _streams) {
          
          // Construct Unique Filename
          std::string specificFile;
          if(filenameBase.empty()) {
              specificFile = _name + "_" + stream.fullName + "_Tree.root";
          } else {
               auto dotPos = filenameBase.find_last_of(".");
               if(dotPos != std::string::npos) {
                   std::string base = filenameBase.substr(0, dotPos);
                   std::string ext  = filenameBase.substr(dotPos);
                   specificFile = base + "_" + stream.fullName + ext;
               } else {
                   specificFile = filenameBase + "_" + stream.fullName + ".root";
               }
          }
          std::string finalPath = MakePath(specificFile);

          // Collect Columns for THIS stream
          auto cols = addCols;

	  // Always add the Event Counter
          cols.push_back("rdfentry_");

          auto streamCols = CollectStreamColumns(*stream.kine);
          cols.insert(cols.end(), streamCols.begin(), streamCols.end());

          // [FIX] If this is a Reconstruction stream, include the Truth columns too!
          // This ensures Rec trees have the "repeated truth values" you requested.
          if(stream.source.find(Rec()) == 0) {
              cols.insert(cols.end(), globalTruthCols.begin(), globalTruthCols.end());
          }

          // Determine Mask for THIS stream
          std::string mask = stream.sel->GetMaskColumn();
          
          // Fallback: If no cuts defined, create a "Pass All" mask
          if(mask.empty()) {
               auto pNames = stream.kine->Creator().GetParticleNames();
               if(!pNames.empty()) {
                   mask = stream.kine->GetPrefix() + "Analysis_AllIndices" + stream.kine->GetSuffix();
                   std::string ref = stream.kine->FullName(pNames[0] + "_"+ consts::NamePx());
                   if(!_reaction.ColumnExists(mask)) {
                       _reaction.Define(mask, [](const ROOT::RVecD& r){ 
                           return rad::util::EnumerateIndicesFrom(0, r.size()); 
                       }, {ref});
                   }
               }
          }

          std::cout << "[AnalysisManager] Snapshotting stream '" << stream.fullName 
                    << "' to " << finalPath << " (with Mask: " << mask << ")" << std::endl;
          
          _reaction.BookSnapshotCombi(finalPath, "tree", cols, mask);
      }
  }
  
  template <typename R, typename P>
  inline void AnalysisManager<R,P>::Run(const std::string& suffix) {
      Init();
      // If Snapshot() wasn't called, this triggers the event loop now.
      _reaction.TriggerSnapshots();
      
      if(!_streams.empty()) {
          for(auto& [key, stream] : _streams) {
              std::string fname = _name + "_" + stream.fullName + suffix;
              std::string finalPath = MakePath(fname);
              
              stream.hist->File(finalPath, "RECREATE");
              std::cout << "[AnalysisManager] Wrote " << stream.fullName << " histograms to " << finalPath << std::endl;
          }
      }
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::PrintDiagnostics() const {
    diag::DiagnosticsPrinter::PrintSectionHeader("ANALYSIS MANAGER DIAGNOSTICS", '=', 90);
    const_cast<R&>(_reaction).PrintReactionDiagnostics();
    diag::DiagnosticsPrinter::PrintBlank();
    std::cout << "Registered Streams: " << _streams.size() << std::endl;
    for (const auto& entry : _streams) {
      std::cout << "\nStream: " << entry.first << " (Source: " << entry.second.source << ", Label: " << entry.second.label << ")" << std::endl;
      if(entry.second.kine) entry.second.kine->PrintReactionMap(); // Use PrintReactionMap or PrintProcessorDiagnostics if available
    }
  }

  // --- Private Helpers ---

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
          cols.push_back(kine.FullName(var));
      }
      // Signal Flag uses Input Prefix (rec_) not stream suffix
      std::string sigCol = kine.GetPrefix() + rad::consts::TruthMatchedCombi();
      if(_reaction.ColumnExists(sigCol)) {
          cols.push_back(sigCol);
      }
      return cols;
  }

} // namespace rad
