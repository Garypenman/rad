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
     * @tparam ReactionType The concrete reaction class (e.g. ePICReaction).
     */
  template <typename ReactionClass,typename ProcessorClass>
    class AnalysisManager {
 /**
     * @struct AnalysisStream
     * @brief Holds the components for a single data stream (e.g., "rec_" or "tru_").
     * @details
     * A stream is a self-contained processing unit consisting of:
     * - ProcessorClass: Creates particles and variables.
     * - PhysicsSelection: Applies cuts.
     * - Histogrammer: Creates plots.
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
            // Selection depends on Kine
            sel  = std::make_unique<PhysicsSelection>(*kine);
            // Histograms depend on Kine and Selection
            hist = std::make_unique<histo::Histogrammer>(*kine, sel.get());
        }
    };
    
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
        AnalysisManager(const std::string& name, const std::string& treeName, const std::string& fileGlob) 
            : _reaction(treeName, fileGlob), _name(name)
        {
        }

        /** @brief Sets and creates the output directory. */
        void SetOutputDir(const std::string& dir) {
            _outputDir = dir;
            if (!_outputDir.empty() && !std::filesystem::exists(_outputDir)) {
                std::filesystem::create_directories(_outputDir);
            }
        }

        /**
         * @brief Define which data types to process.
         * @details The **First** type provided is treated as the **Primary Stream**.
         * The Primary Stream determines the Event Selection Mask for the Snapshot.
         * @param types Variadic list of stream prefixes (e.g. Rec(), Truth()).
         */
        template<typename... Args>
        void SetTypes(Args... types) {
            (AddStream(types), ...);
        }

        /** @return Reference to the underlying Reaction object. */
        ReactionClass& Reaction() { return _reaction; }

        // =====================================================================
        // Configuration
        // =====================================================================

        /** @brief Apply a Kinematics recipe to ALL active streams. */
        void ConfigureKinematics(KineRecipe recipe) {
            for(auto& [key, stream] : _streams) { if(stream.kine) recipe(*stream.kine); }
        }
        /** @brief Apply a Kinematics recipe to a SPECIFIC stream only. */
        void ConfigureKinematics(const std::string& streamName, KineRecipe recipe) {
            if(CheckStream(streamName)) recipe(*_streams.at(streamName).kine);
        }

        /** @brief Apply a Selection recipe to ALL active streams. */
        void ConfigureSelection(SelRecipe recipe) {
            for(auto& [key, stream] : _streams) { if(stream.sel) recipe(*stream.sel); }
        }
        /** @brief Apply a Selection recipe to a SPECIFIC stream only. */
        void ConfigureSelection(const std::string& streamName, SelRecipe recipe) {
            if(CheckStream(streamName)) recipe(*_streams.at(streamName).sel);
        }

        /** @brief Apply a Histogram recipe to ALL active streams. */
        void ConfigureHistograms(HistoRecipe recipe) {
            for(auto& [key, stream] : _streams) { if(stream.hist) recipe(*stream.hist); }
        }
        /** @brief Apply a Histogram recipe to a SPECIFIC stream only. */
        void ConfigureHistograms(const std::string& streamName, HistoRecipe recipe) {
            if(CheckStream(streamName)) recipe(*_streams.at(streamName).hist);
        }

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
        void Init() {
	  cout<< "AnalysisMAnager::Init "<< endl;
            if(_initialized) return;
            if(_primaryStream.empty()) throw std::runtime_error("[AnalysisManager] No streams defined! Call SetTypes().");

            // PASS 1: Initialize Kinematics (Create Variables)
            // This ensures all columns (e.g. "rec_Mass") exist in the dataframe.
            for(auto& [key, stream] : _streams) {
	      cout<< "AnalysisMAnager::Init Kine "<< endl;
                stream.kine->Init();
                stream.kine->PrintReactionMap();
                // Define Signal Flags (Truth Matching)
                // Depends on Kinematics being ready (to know topology)
		//    _reaction.DefineSignalFlag(stream.prefix); 
            }

            // PASS 2: Compile Selections (Create Masks)
            // This uses the variables created in Pass 1.
            for(auto& [key, stream] : _streams) {
 	  cout<< "AnalysisMAnager::Init Selections"<< endl;
                if(stream.sel) stream.sel->Init();
            }

            // PASS 3: Initialize Histograms (Book Actions)
            // This uses Variables (Pass 1) and Masks (Pass 2).
            for(auto& [key, stream] : _streams) {
	      cout<< "AnalysisMAnager::Init histograms "<< endl;
                if(stream.hist) stream.hist->Init();
            }

            _initialized = true;
        }

        /**
         * @brief Snapshots data to a SINGLE flat TTree.
         * @details 
         * - Merges columns from all streams (Rec + Truth).
         * - Uses the Primary Stream's mask to filter events.
         * - Uses the Primary Stream's combination index for alignment.
         * - Filename defaults to: {OutputDir}/{AnalysisName}_Tree.root
         * @param filename Optional override for the output filename.
         */
        void Snapshot( const ROOT::RDF::ColumnNames_t& addCols={},const std::string& filename = "") {
	  
            Init();
            
            // Default Name: "Y4260_Tree.root"
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
            // BookSnapshotCombi is an Immediate Action (it runs the event loop)
            _reaction.BookSnapshotCombi(finalPath, "tree", cols, mask);
        }

        /**
         * @brief Runs the analysis and writes histograms.
         * @details 
         * - Creates ONE file PER stream to avoid key clashes.
         * - Filename: {OutputDir}/{AnalysisName}_{Prefix}_{Suffix}
         * - Example:  output/Y4260_rec_Hist.root
         * @param suffix Suffix for the histogram file (default "Hist.root").
         */
        void Run(const std::string& suffix = "Hist.root") {
            Init();
            // If Snapshot() wasn't called, this triggers the event loop now.
            // If Snapshot() WAS called, this is a no-op (loop already ran).
            _reaction.TriggerSnapshots();
            
            if(!_streams.empty()) {
                for(auto& [key, stream] : _streams) {
                    // Construct unique filename per stream
                    std::string fname = _name + "_" + stream.prefix + suffix;
                    std::string finalPath = MakePath(fname);
                    
                    // Safe to RECREATE because each stream has its own file.
                    stream.hist->File(finalPath, "RECREATE");
                    
                    std::cout << "[AnalysisManager] Wrote " << stream.prefix << " histograms to " << finalPath << std::endl;
                }
            }
        }
      
    private:
        ReactionClass _reaction;
        std::string _name;
        std::string _outputDir;
        bool _initialized = false;
        std::string _primaryStream;
        std::map<std::string, AnalysisStream> _streams;

        void AddStream(const std::string& prefix) {
            if(_streams.find(prefix) == _streams.end()) {
                _streams.try_emplace(prefix, &_reaction, prefix);
                // First added stream becomes the Primary (Driver)
                if(_primaryStream.empty()) _primaryStream = prefix; 
            }
        }

        bool CheckStream(const std::string& name) {
            if(_streams.find(name) == _streams.end()) {
                std::cerr << "AnalysisManager Warning: Stream '" << name << "' not found. Recipe ignored." << std::endl;
                return false;
            }
            return true;
        }

        std::string MakePath(const std::string& filename) {
            if(!_outputDir.empty()) {
                return (std::filesystem::path(_outputDir) / filename).string();
            }
            return filename;
        }

        ROOT::RDF::ColumnNames_t CollectStreamColumns(const ProcessorClass& kine) {
            ROOT::RDF::ColumnNames_t cols;
            
            // 1. Collect Variables
            // GetDefinedNames() returns everything registered in the Processor:
            // - P4 components (ele_px...) added by DefineNewComponentVecs
            // - Calculations (MassJ...) added by RegisterCalc
            for(const auto& var : kine.GetDefinedNames()) {
                cols.push_back(kine.GetPrefix() + var + kine.GetSuffix());
            }

            // 2. Add Signal Flag if it exists
            // This flag is defined by ConfigReaction via DefineSignalFlag called in Init().
            // Uses standard constant "isTruth" (or "is_signal_combi")
            std::string sigCol = kine.GetPrefix() + rad::consts::TruthMatchedCombi();
            if(_reaction.ColumnExists(sigCol)) {
                cols.push_back(sigCol);
            }
            
            return cols;
        }
    };

} // namespace rad

/* /\** */
/*  * @file AnalysisManager.h */
/*  * @brief High-level orchestration of the RAD analysis chain. */
/*  *\/ */

/* #pragma once */

/* #include "ConfigReaction.h" */
/* #include "ProcessorClass.h" */
/* #include "PhysicsSelection.h" */
/* #include "Histogrammer.h" */

/* #include <memory> */
/* #include <vector> */
/* #include <map> */
/* #include <string> */
/* #include <stdexcept> */
/* #include <iostream> */
/* #include <filesystem> */
/* #include <functional> */

/* namespace rad { */

/*     /\** */
/*      * @struct AnalysisStream */
/*      * @brief Holds the components for a single data stream (e.g., "rec_" or "tru_"). */
/*      *\/ */
/*     struct AnalysisStream { */
/*         std::string prefix; */
/*         std::unique_ptr<ProcessorClass> kine; */
/*         std::unique_ptr<PhysicsSelection> sel; */
/*         std::unique_ptr<histo::Histogrammer> hist; */

/*         AnalysisStream(ConfigReaction* reaction, const std::string& p)  */
/*             : prefix(p)  */
/*         { */
/*             kine = std::make_unique<ProcessorClass>(reaction, prefix); */
/*             sel  = std::make_unique<PhysicsSelection>(*kine); */
/*             hist = std::make_unique<histo::Histogrammer>(*kine, sel.get()); */
/*         } */
/*     }; */

/*     template <typename ReactionClass> */
/*     class AnalysisManager { */
/*     public: */
/*         using KineRecipe  = std::function<void(ProcessorClass&)>; */
/*         using SelRecipe   = std::function<void(PhysicsSelection&)>; */
/*         using HistoRecipe = std::function<void(histo::Histogrammer&)>; */

/*         // ===================================================================== */
/*         // Setup */
/*         // ===================================================================== */
        
/*         AnalysisManager(const std::string& name, const std::string& treeName, const std::string& fileGlob)  */
/*             : _reaction(treeName, fileGlob), _name(name) */
/*         { */
/*         } */

/*         void SetOutputDir(const std::string& dir) { */
/*             _outputDir = dir; */
/*             if (!_outputDir.empty() && !std::filesystem::exists(_outputDir)) { */
/*                 std::filesystem::create_directories(_outputDir); */
/*             } */
/*         } */

/*         template<typename... Args> */
/*         void SetTypes(Args... types) { */
/*             (AddStream(types), ...); */
/*         } */

/*         ReactionClass& Reaction() { return _reaction; } */

/*         // ===================================================================== */
/*         // Configuration */
/*         // ===================================================================== */

/*         void ConfigureKinematics(KineRecipe recipe) { */
/*             for(auto& [key, stream] : _streams) { if(stream.kine) recipe(*stream.kine); } */
/*         } */
/*         void ConfigureKinematics(const std::string& streamName, KineRecipe recipe) { */
/*             if(CheckStream(streamName)) recipe(*_streams.at(streamName).kine); */
/*         } */

/*         void ConfigureSelection(SelRecipe recipe) { */
/*             for(auto& [key, stream] : _streams) { if(stream.sel) recipe(*stream.sel); } */
/*         } */
/*         void ConfigureSelection(const std::string& streamName, SelRecipe recipe) { */
/*             if(CheckStream(streamName)) recipe(*_streams.at(streamName).sel); */
/*         } */

/*         void ConfigureHistograms(HistoRecipe recipe) { */
/*             for(auto& [key, stream] : _streams) { if(stream.hist) recipe(*stream.hist); } */
/*         } */
/*         void ConfigureHistograms(const std::string& streamName, HistoRecipe recipe) { */
/*             if(CheckStream(streamName)) recipe(*_streams.at(streamName).hist); */
/*         } */

/*         // ===================================================================== */
/*         // Execution */
/*         // ===================================================================== */

/*         void Init() { */
/*             if(_initialized) return; */
/*             if(_primaryStream.empty()) throw std::runtime_error("[AnalysisManager] No streams defined!"); */

/*             for(auto& [key, stream] : _streams) { */
/*                 stream.kine->Init(); */
/* 		//  _reaction.DefineSignalFlag(stream.prefix);  */
/*                 if(stream.sel) stream.sel->Compile(); */
/*             } */
/*             _initialized = true; */
/*         } */

/*         /\** */
/*          * @brief Snapshots data to a SINGLE flat TTree. */
/*          * @details Merges columns from all streams.  */
/*          * Filename defaults to: {OutputDir}/{AnalysisName}_Tree.root */
/*          *\/ */
/*         void Snapshot(const std::string& filename = "") { */
/*             Init(); */
            
/*             // Default Name: "Y4260_Tree.root" */
/*             std::string finalName = filename.empty() ? _name + "_Tree.root" : filename; */
/*             std::string finalPath = MakePath(finalName); */

/*             ROOT::RVec<std::string> cols; */
            
/*             // 1. Collect Primary Stream (Rec) */
/*             if(_streams.count(_primaryStream)) { */
/*                 auto c = CollectStreamColumns(*_streams.at(_primaryStream).kine); */
/*                 cols.insert(cols.end(), c.begin(), c.end()); */
/*             } */

/*             // 2. Collect Secondary Streams (Truth) */
/*             for(auto& [key, stream] : _streams) { */
/*                 if(key == _primaryStream) continue;  */
/*                 auto c = CollectStreamColumns(*stream.kine); */
/*                 cols.insert(cols.end(), c.begin(), c.end()); */
/*             } */

/*             // 3. Determine Mask (Primary Driven) */
/*             std::string mask = ""; */
/*             auto& prime = _streams.at(_primaryStream); */
/*             mask = prime.sel->GetMaskColumn(); */
            
/*             if(mask.empty()) { */
/*                  auto pNames = prime.kine->Creator().GetParticleNames(); */
/*                  if(!pNames.empty()) { */
/*                      mask = prime.prefix + "Analysis_AllIndices"; */
/*                      std::string ref = prime.prefix + pNames[0] + "_px"; */
/*                      _reaction.Define(mask, [](const ROOT::RVecD& r){  */
/*                         return rad::util::EnumerateIndicesFrom(0, r.size());  */
/*                      }, {ref}); */
/*                  } */
/*             } */

/*             std::cout << "[AnalysisManager] Snapshotting combined tree to " << finalPath << std::endl; */
/*             _reaction.BookSnapshotCombi(finalPath, "tree", cols, mask); */
/*         } */

/*         /\** */
/*          * @brief Runs the analysis and writes histograms. */
/*          * @details Creates ONE file PER stream to avoid clashes. */
/*          * Filename: {OutputDir}/{AnalysisName}_{Prefix}_{Suffix} */
/*          * Example:  output/Y4260_rec_Hist.root */
/*          *\/ */
/*         void Run(const std::string& suffix = "Hist.root") { */
/*             Init(); */
/*             _reaction.TriggerSnapshots(); */
            
/*             if(!_streams.empty()) { */
/*                 for(auto& [key, stream] : _streams) { */
/*                     // Construct unique filename per stream */
/*                     std::string fname = _name + "_" + stream.prefix + suffix; */
/*                     std::string finalPath = MakePath(fname); */
                    
/*                     // Safe to RECREATE (Unique file per stream) */
/*                     stream.hist->File(finalPath, "RECREATE"); */
                    
/*                     std::cout << "[AnalysisManager] Wrote " << stream.prefix << " histograms to " << finalPath << std::endl; */
/*                 } */
/*             } */
/*         } */

/*     private: */
/*         ReactionClass _reaction; */
/*         std::string _name; */
/*         std::string _outputDir; */
/*         bool _initialized = false; */
/*         std::string _primaryStream; */
/*         std::map<std::string, AnalysisStream> _streams; */

/*         void AddStream(const std::string& prefix) { */
/*             if(_streams.find(prefix) == _streams.end()) { */
/*                 _streams.try_emplace(prefix, &_reaction, prefix); */
/*                 if(_primaryStream.empty()) _primaryStream = prefix;  */
/*             } */
/*         } */

/*         bool CheckStream(const std::string& name) { */
/*             if(_streams.find(name) == _streams.end()) { */
/*                 std::cerr << "AnalysisManager Warning: Stream '" << name << "' not found." << std::endl; */
/*                 return false; */
/*             } */
/*             return true; */
/*         } */

/*         std::string MakePath(const std::string& filename) { */
/*             if(!_outputDir.empty()) { */
/*                 return (std::filesystem::path(_outputDir) / filename).string(); */
/*             } */
/*             return filename; */
/*         } */

/*         ROOT::RVec<std::string> CollectStreamColumns(const ProcessorClass& kine) { */
/*             ROOT::RVec<std::string> cols; */
         
/*             for(const auto& var : kine.GetDefinedNames()) { */
/*                 cols.push_back(kine.GetPrefix() + var + kine.GetSuffix()); */
/*             } */
/*             if(_reaction.ColumnExists(kine.GetPrefix() + consts::TruthMatchedCombi())) cols.push_back(kine.GetPrefix() + consts::TruthMatchedCombi()); */
/*             return cols; */
/*         } */
/*     }; */

/* } // namespace rad */
