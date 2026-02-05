/**
 * @file SnapshotCombi.h
 * @brief Custom RDataFrame Action to write "Flat" TTrees from combinatorial data.
 * @details 
 * This action allows RDataFrame to save combinatorial data (Rec) and scalar data (Truth)
 * into a single flat TTree.
 * - Handles multi-threading via ROOT::TBufferMerger.
 * - Supports automatic broadcasting of scalar values (size 1) to match vector sizes.
 * - Performs safety checks on vector sizes to prevent index out-of-bounds errors.
 */

#pragma once

#include "RDFUtils.h" 

#include <ROOT/RDF/RActionImpl.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/TBufferMerger.hxx>
#include <TTree.h>
#include <TDirectory.h> 
#include <TFile.h>
#include <ROOT/TSeq.hxx>

#include <vector>
#include <string>
#include <tuple>
#include <memory>
#include <type_traits>
#include <variant>
#include <numeric>
#include <iostream>

namespace rad {
namespace io {

    template <typename T> struct IsRVec : std::false_type {};
    template <typename T> struct IsRVec<ROOT::VecOps::RVec<T>> : std::true_type {};
    template <typename T> inline constexpr bool IsRVec_v = IsRVec<T>::value;

    /**
     * @class SnapshotCombi
     * @brief RDataFrame Action that flattens combinatorial vectors into a TTree.
     */
    class SnapshotCombi : public ROOT::Detail::RDF::RActionImpl<SnapshotCombi> {
    public:
        using Result_t = unsigned long;

        /**
         * @brief Constructor.
         * @param filename Output ROOT file name.
         * @param treename Output TTree name.
         * @param colNames List of branch names to create.
         * @param colTypes List of types corresponding to the branches.
         */
        SnapshotCombi(const std::string& filename, const std::string& treename, 
                      const ROOT::RVec<std::string>& colNames,
                      const ROOT::RVec<ColType>& colTypes);

        SnapshotCombi(SnapshotCombi&&) = default;
        SnapshotCombi(const SnapshotCombi&) = delete;

        std::string GetActionName();
        
        /** @brief Resets counters and opens the TBufferMerger. */
        void Initialize();
        
        /** @brief Called at the start of each thread's processing slot. */
        void InitTask(TTreeReader*, unsigned int slot);
        
        /** @brief Called at the end of each thread's processing slot. Writes thread-local buffers. */
        void FinalizeTask(unsigned int slot); 
        
        /** @brief Called at the end of the event loop. Merges files and writes the final tree. */
        void Finalize();

        /** @brief Returns pointer to the total entry count. */
        std::shared_ptr<Result_t> GetResultPtr() const;

        /**
         * @brief Per-event execution function.
         * @details Iterates over the "Mask" column (last argument) and fills the tree 
         * for each valid combination index.
         * @param slot Thread slot ID.
         * @param args Variadic list of column values + Mask.
         */
        template <typename... Args>
        void Exec(unsigned int slot, const Args&... args);

    private:
        std::string _fileName;
        std::string _treeName;
        ROOT::RVec<std::string> _colNames;
        ROOT::RVec<ColType> _colTypes;

        std::shared_ptr<ROOT::TBufferMerger> _merger;
        std::shared_ptr<unsigned long> _totalCount; 
        
        ROOT::RVec<unsigned long> _slotCounts;      

        using Value_t = std::variant<double, float, int, unsigned int, short, bool, long long>;
        using Buffer_t = ROOT::RVec<Value_t>;

        ROOT::RVec<std::shared_ptr<ROOT::TBufferMergerFile>> _files; 
        ROOT::RVec<std::shared_ptr<TTree>> _trees;
        ROOT::RVec<std::unique_ptr<Buffer_t>> _buffers; 

        // --- Helpers ---
        template <typename Tuple, size_t... Is>
        void fill_buffers(unsigned int slot, int idx, const Tuple& t, std::index_sequence<Is...>);

        template <size_t I, typename T>
        void set_value(unsigned int slot, int idx, const T& input);

        /** @brief Helper for scalar inputs (returns value directly). */
        template <typename T>
        auto get_flat_val(const T& input, int idx);

        /** * @brief Helper for vector inputs (returns input[idx] or broadcasts input[0]). 
         * @throws std::runtime_error if index exceeds vector size (and size != 1).
         */
        template <typename T>
        auto get_flat_val(const ROOT::VecOps::RVec<T>& input, int idx);
    };


    // =========================================================================
    // IMPLEMENTATION
    // =========================================================================

    inline SnapshotCombi::SnapshotCombi(const std::string& filename, const std::string& treename, 
                                        const ROOT::RVec<std::string>& colNames,
                                        const ROOT::RVec<ColType>& colTypes)
        : _fileName(filename), _treeName(treename), 
          _colNames(colNames), _colTypes(colTypes),
          _totalCount(std::make_shared<unsigned long>(0))
    {}

    inline std::string SnapshotCombi::GetActionName() { return "SnapshotCombi"; }
    
    inline std::shared_ptr<unsigned long> SnapshotCombi::GetResultPtr() const { 
        return _totalCount; 
    }

    inline void SnapshotCombi::Initialize() {
        *_totalCount = 0;
        _merger = std::make_shared<ROOT::TBufferMerger>(_fileName.c_str(), "RECREATE");
        
        const auto nSlots = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;
        
        _files.resize(nSlots);
        _trees.resize(nSlots);
        _buffers.resize(nSlots);
        _slotCounts.assign(nSlots, 0);

        std::cout << "[SnapshotCombi] Initialized with " << nSlots << " slots." << std::endl;
    }

    inline void SnapshotCombi::FinalizeTask(unsigned int slot) {
        if (slot < _trees.size() && _trees[slot] && _trees[slot]->GetEntries() > 0) {
            _files[slot]->Write();
        }
        if(slot < _trees.size()) _trees[slot].reset();
        if(slot < _files.size()) _files[slot].reset();
    }

    inline void SnapshotCombi::Finalize() {
        std::cout << "[SnapshotCombi] Finalizing..." << std::endl;
        
        for (auto &file : _files) {
            if (file) {
                file->Write();
                file->Close();
            }
        }
        
        *_totalCount = std::accumulate(_slotCounts.begin(), _slotCounts.end(), 0ul);
        std::cout << "[SnapshotCombi] Total Entries: " << *_totalCount << std::endl;

        _trees.clear();  
        _files.clear();  
        _buffers.clear(); 
        _merger.reset();

        if (*_totalCount == 0) {
            TFile f(_fileName.c_str(), "UPDATE");
            if (!f.Get(_treeName.c_str())) {
                TTree* emptyTree = new TTree(_treeName.c_str(), "Flat Combinatorial Tree");
                emptyTree->Write();
                delete emptyTree;
            }
            f.Close();
        }
    }

    inline void SnapshotCombi::InitTask(TTreeReader*, unsigned int slot) {
        TDirectory::TContext c; 

        if (slot >= _files.size()) return;

        if (!_files[slot]) {
            _buffers[slot] = std::make_unique<Buffer_t>(_colNames.size());

            if(!_merger) throw std::runtime_error("[SnapshotCombi] InitTask called before Initialize");
            _files[slot] = _merger->GetFile(); 
            
            _trees[slot] = std::make_shared<TTree>(_treeName.c_str(), "Flat Combinatorial Tree");
            _trees[slot]->SetDirectory(_files[slot].get()); 
            _trees[slot]->SetImplicitMT(false); 
            _trees[slot]->SetAutoSave(0); 

            auto& buffer = *_buffers[slot]; 
            for(size_t i=0; i<_colNames.size(); ++i) {
                const char* name = _colNames[i].c_str();
                
                // [FIX] Force variant to hold the correct type BEFORE getting address
                switch(_colTypes[i]) {
                    case ColType::Double: 
                        buffer[i] = double(0); 
                        _trees[slot]->Branch(name, &std::get<double>(buffer[i])); 
                        break;
                    case ColType::Float:  
                        buffer[i] = float(0); 
                        _trees[slot]->Branch(name, &std::get<float>(buffer[i])); 
                        break;
                    case ColType::Int:    
                        buffer[i] = int(0); 
                        _trees[slot]->Branch(name, &std::get<int>(buffer[i])); 
                        break;
                    case ColType::UInt:   
                        buffer[i] = (unsigned int)(0); 
                        _trees[slot]->Branch(name, &std::get<unsigned int>(buffer[i])); 
                        break;
                    case ColType::Short:  
                        buffer[i] = (short)(0); 
                        _trees[slot]->Branch(name, &std::get<short>(buffer[i])); 
                        break;
                    case ColType::Bool:   
                        buffer[i] = false; 
                        _trees[slot]->Branch(name, &std::get<bool>(buffer[i])); 
                        break;
                    case ColType::Long:   
                        buffer[i] = (long long)(0); 
                        _trees[slot]->Branch(name, &std::get<long long>(buffer[i])); 
                        break;
                    default:              
                        buffer[i] = double(0);
                        _trees[slot]->Branch(name, &std::get<double>(buffer[i]));
                }
            }
        }
    }

    template <typename... Args>
    inline void SnapshotCombi::Exec(unsigned int slot, const Args&... args) {
        if constexpr (sizeof...(Args) > 0) {
            auto all_args = std::forward_as_tuple(args...);
            constexpr size_t NData = sizeof...(Args) - 1; 
            const auto& mask = std::get<NData>(all_args);

            if (!mask.empty()) {
                for (auto idx : mask) {
                    fill_buffers(slot, idx, all_args, std::make_index_sequence<NData>{});
                    _trees[slot]->Fill();
                    _slotCounts[slot]++; 
                }
            }
        }
    }

    // --- Private Helpers ---

    template <typename Tuple, size_t... Is>
    inline void SnapshotCombi::fill_buffers(unsigned int slot, int idx, const Tuple& t, std::index_sequence<Is...>) {
        ((set_value<Is>(slot, idx, std::get<Is>(t))), ...);
    }

    template <size_t I, typename T>
    inline void SnapshotCombi::set_value(unsigned int slot, int idx, const T& input) {
        auto val = get_flat_val(input, idx);
        auto& buffer = *_buffers[slot];

        switch(_colTypes[I]) {
            case ColType::Double: std::get<double>(buffer[I]) = val; break;
            case ColType::Float:  std::get<float>(buffer[I])  = val; break;
            case ColType::Int:    std::get<int>(buffer[I])    = val; break;
            case ColType::UInt:   std::get<unsigned int>(buffer[I]) = val; break;
            case ColType::Short:  std::get<short>(buffer[I])  = val; break;
            case ColType::Bool:   std::get<bool>(buffer[I])   = val; break;
            case ColType::Long:   std::get<long long>(buffer[I]) = val; break;
            default:              std::get<double>(buffer[I]) = val; break;
        }
    }

    template <typename T>
    inline auto SnapshotCombi::get_flat_val(const T& input, int idx) {
        if constexpr (!IsRVec_v<std::decay_t<T>>) return input; 
        else return input; 
    }
    
    template <typename T>
    inline auto SnapshotCombi::get_flat_val(const ROOT::VecOps::RVec<T>& input, int idx) {
        // CASE A: Broadcast (Scalar-in-Vector)
        // If the vector has size 1 (e.g. Truth variable), broadcast it to all Rec combinations.
        if (input.size() == 1) {
            return input[0];
        }

        // CASE B: Combinatorial Data
        // If the vector is larger, it must match the combinatorial indexing.
        // Safety Check:
        if (idx >= (int)input.size()) {
             throw std::runtime_error("[SnapshotCombi] Index Mismatch! "
                 "Combinatorial mask index " + std::to_string(idx) + 
                 " exceeds column size " + std::to_string(input.size()) + 
                 ". Check if a scalar column was not broadcast correctly.");
        }
        return input[idx];
    }

} // namespace io
} // namespace rad
