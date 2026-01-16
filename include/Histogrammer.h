/**
 * @file Histogrammer.h
 * @brief Manages N-Dimensional Histogramming with automatic "Unfolding" (Projection).
 * @details 
 * This class coordinates with KinematicsProcessor and PhysicsSelection to:
 * 1. Book N-Dimensional THnSparse histograms (Variable + Splits).
 * 2. Automatically apply "Good Candidate" masks.
 * 3. Write results to TFile, recursively projecting N-D histograms back into 
 * user-friendly 1D TH1D plots (e.g. "Mass_Sector1_Pt2").
 */

#pragma once

#include "ConfigReaction.h"
#include "KinematicsProcessor.h"
#include "PhysicsSelection.h"
#include "THnCombi.h" // Requires the THnCombi implementation provided previously

#include <THnSparse.h>
#include <TFile.h>
#include <vector>
#include <string>
#include <memory>
#include <iostream>

namespace rad {
namespace histo {

    /**
     * @struct SplitDef
     * @brief Configuration for a single splitting axis (category).
     */
    struct SplitDef {
        std::string name;    ///< Human-readable name (e.g., "Sector")
        std::string baseCol; ///< Column name in Processor (e.g., "loc_sector")
        int nbins;
        double min;
        double max;
    };

    /**
     * @class Histogrammer
     * @brief Manager for creating and persisting multidimensional kinematic histograms.
     */
    class Histogrammer {
    public:
        /**
         * @brief Constructor binding to a specific Processor stream.
         * @param proc The processor (Rec, Sim, etc.) used to resolve variable names.
         * @param sel Optional pointer to PhysicsSelection. If provided, histograms will 
         * automatically filter using the selection's compiled mask.
         */
        Histogrammer(KinematicsProcessor& proc, PhysicsSelection* sel = nullptr);

        /**
         * @brief Defines a split axis (category) for subsequent histograms.
         * @details All histograms created after calling this will have this extra dimension.
         * @param name Axis title (e.g. "Sector").
         * @param baseCol Variable name in RDF (e.g. "loc_sector"). Resolved via Processor.
         * @param nbins Number of bins.
         * @param min Axis min.
         * @param max Axis max.
         */
        void AddSplit(const std::string& name, const std::string& baseCol, int nbins, double min, double max);

        /**
         * @brief Books a multidimensional histogram (THnSparse).
         * @details
         * Creates a THnSparse with dimensions: [Var, Split1, Split2, ...].
         * Automatically attaches the PhysicsSelection mask if available.
         * * @param name Output histogram base name.
         * @param title Histogram title.
         * @param nbins Number of bins for the main variable.
         * @param min Minimum value for the main variable.
         * @param max Maximum value for the main variable.
         * @param varBaseName The RDF column name (e.g. "mm2"). Resolved via Processor.
         */
        void Create(const std::string& name, const std::string& title, 
                    int nbins, double min, double max, 
                    const std::string& varBaseName);

        /**
         * @brief Writes all histograms to a ROOT file.
         * @details
         * Performs "Unfolding": Recursively projects the N-D histograms into 1D TH1Ds 
         * for every bin combination of the split axes.
         * Creates TDirectories to organize output.
         * * @param filename Output filename.
         * * @param option file constructor option
         */
        void File(const std::string& filename, const std::string& option = "RECREATE");

    private:
        KinematicsProcessor& _proc;
        ConfigReaction& _rad; // Underlying RDF Manager
        std::string _maskCol;
        
        std::vector<SplitDef> _splits;
        
        // Store RResultPtrs to ensure RDF executes the actions
        std::vector<ROOT::RDF::RResultPtr<THnSparseD>> _results;

        /**
         * @brief Recursive helper to project N-D histogram to 1D slices.
         */
        void UnfoldAndWrite(THnSparseD* hn, TDirectory* dir, std::string current_suffix = "", int axis_depth = 1);
    };

    // =================================================================================
    // IMPLEMENTATION: Histogrammer
    // =================================================================================

    inline Histogrammer::Histogrammer(KinematicsProcessor& proc, PhysicsSelection* sel) 
        : _proc(proc), _rad(*proc.Reaction()), _maskCol("") 
    {
        if (sel) {
            _maskCol = sel->GetMaskColumn();
            if (_maskCol.empty()) {
                std::cerr << "Histogrammer Warning: Attached PhysicsSelection has no compiled mask. "
                          << "Histograms will use all combinations." << std::endl;
            }
        }
    }

    inline void Histogrammer::AddSplit(const std::string& name, const std::string& baseCol, int nbins, double min, double max) {
        _splits.push_back({name, baseCol, nbins, min, max});
    }

    inline void Histogrammer::Create(const std::string& name, const std::string& title, 
                int nbins, double min, double max, 
                const std::string& varBaseName) 
    {
        // 1. Resolve Column Names using Processor
        std::string fullVarName = _proc.FullName(varBaseName);
        std::vector<std::string> cols;
        cols.push_back(fullVarName);

        // 2. Setup Dimensions (1 Main Var + N Splits)
        int n_dims = 1 + _splits.size();
        std::vector<int> bins_vec(n_dims);
        std::vector<double> xmin(n_dims), xmax(n_dims);
        
        // Axis 0: The Main Variable
        bins_vec[0] = nbins; xmin[0] = min; xmax[0] = max;
        
        // Axis 1..N: The Splits
        for(size_t i=0; i<_splits.size(); ++i) {
            std::string fullSplitCol = _proc.FullName(_splits[i].baseCol);
            cols.push_back(fullSplitCol);
            
            bins_vec[i+1] = _splits[i].nbins;
            xmin[i+1] = _splits[i].min;
            xmax[i+1] = _splits[i].max;
        }

        // 3. Construct the THnSparse Prototype
        std::string hNameFull = _proc.GetPrefix() + name + _proc.GetSuffix();
        auto hist = std::make_shared<THnSparseD>(hNameFull.c_str(), title.c_str(), 
                                                 n_dims, bins_vec.data(), xmin.data(), xmax.data());
        
        // Label Axes
        hist->GetAxis(0)->SetName(name.c_str());
        for(size_t i=0; i<_splits.size(); ++i) {
            hist->GetAxis(i+1)->SetName(_splits[i].name.c_str());
        }

        // 4. Book the Action (Masked vs Unmasked)
        // Note: We assume all columns (Var & Splits) are of type ResultType_t (RVec<double>).
        // If splits are different types, generic lambdas or explicit casting in RDF is required.
	if(!_maskCol.empty()) {
        cols.push_back(_maskCol); 
        
        // Use <true> for Masked execution. No data types needed in template.
        THnCombi<true> action(hist);
        _results.push_back(_rad.CurrFrame().Book(std::move(action), cols));
	} 
	else {
	  // Use <false> for Standard execution.
	  THnCombi<false> action(hist);
	  _results.push_back(_rad.CurrFrame().Book(std::move(action), cols));
	}
	
    }

  inline void  Histogrammer::File(const std::string& filename, const std::string& option) {
            // Use the user-provided option (RECREATE wipes the corrupt file)
            auto file = std::unique_ptr<TFile>(TFile::Open(filename.c_str(), option.c_str())); 
            
            if (!file || file->IsZombie()) {
                std::cerr << "Histogrammer Error: Could not open file " << filename << std::endl;
                return;
            }

            // Trigger the Event Loop here if it hasn't run yet!
            for(auto& resultPtr : _results) {
                // Accessing the result (via Clone or *) triggers the RDataFrame loop
                auto hist_raw = resultPtr->Clone();
		auto hn = dynamic_cast<THnSparseD*>(hist_raw);
            
            if (hn) {
                // Create a directory for this variable (e.g. "rec_Mass_miss/")
                TDirectory* subdir = file->mkdir(hn->GetName());
                if (!subdir) subdir = file->GetDirectory(hn->GetName());
                subdir->cd();
                
                UnfoldAndWrite(hn, subdir);
            }
            delete hn; // Cleanup clone
        }
        file->Close();
    }

    inline void Histogrammer::UnfoldAndWrite(THnSparseD* hn, TDirectory* dir, std::string current_suffix, int axis_depth) {
        // Base Case: We have constrained all split axes. Project the Main Variable (Axis 0).
        if (axis_depth > (int)_splits.size()) {
            dir->cd();
            
            // Project Axis 0 (The Variable)
            // "E" option ensures errors are computed correctly
            TH1D* h1 = hn->Projection(0, "E"); 
            
            std::string baseName = hn->GetName(); 
            std::string full_name = baseName + current_suffix;
            
            h1->SetName(full_name.c_str());
            h1->SetTitle((std::string(hn->GetTitle()) + current_suffix).c_str());
            
            h1->Write();
            delete h1;
            return;
        }

        // Recursive Step: Iterate over all bins in the current Split Axis
        auto axis = hn->GetAxis(axis_depth);
        int nbins = axis->GetNbins();
        
        // Save range to restore later (crucial for recursion)
        int saved_min = axis->GetFirst();
        int saved_max = axis->GetLast();

        for(int i=1; i<=nbins; ++i) {
            // Restrict this axis to bin 'i'
            axis->SetRange(i, i);
            
            // Construct suffix: e.g. "_Sector1"
            // Note: Could be enhanced to use GetBinLabel() if labels are set
            std::string suffix = current_suffix + "_" + axis->GetName() + std::to_string(i);
            
            UnfoldAndWrite(hn, dir, suffix, axis_depth + 1);
        }
        
        // Restore Range
        axis->SetRange(saved_min, saved_max);
    }

} // namespace histo
} // namespace rad
