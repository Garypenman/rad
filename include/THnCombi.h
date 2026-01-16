/**
 * @file THnCombi.h
 * @brief RDataFrame Action for filling N-Dimensional Histograms (THnSparse) with combinatorial data.
 * @details
 * This class serves as the worker engine for the Histogrammer. It acts as a custom RDataFrame Action 
 * optimized for combinatorial analysis. Key features include:
 * * 1. **Multithreading:** Clones the target histogram for every worker thread slot. This allows lock-free filling 
 * during the event loop, followed by a merge step in `Finalize()`.
 * * 2. **Vectorization & Broadcasting:** Automatically handles mixed inputs of Scalars (per-event) and Vectors (per-combination). 
 * Scalars are "broadcast" (repeated) against the vector size.
 * * 3. **Masked Execution (Lazy Filtering):** Efficiently supports a "Good Combinations" mask (from `PhysicsSelection`). If enabled, 
 * it only loops over valid indices, preventing "Ghost Data" (NaNs/InvalidEntries) from 
 * polluting the histograms.
 */

#pragma once

#include <ROOT/RDF/RActionImpl.hxx>
#include <THnSparse.h>
#include <ROOT/RVec.hxx>
#include <ROOT/TSeq.hxx>

#include <vector>
#include <tuple>
#include <algorithm>
#include <memory>
#include <type_traits>

namespace rad {
  namespace histo {

    // =========================================================================
    // HELPER UTILITIES
    // =========================================================================
    // These generic helpers allow uniform access to both Scalar values (Event variables)
    // and Vector values (Particle/Combinatorial variables).

    /** @brief Returns size of a container (RVec). */
    template<typename T> 
    inline size_t get_size(const ROOT::RVec<T>& t) { return t.size(); }

    /** @brief Returns 1 for scalar types. */
    template<typename T> 
    inline size_t get_size(const T&) { return 1; }

    /** @brief Accesses value at index i for a container (RVec). */
    template<typename T> 
    inline double get_val(const ROOT::RVec<T>& t, size_t i) { return static_cast<double>(t[i]); }

    /** @brief Accesses value for a scalar (index i is ignored). */
    template<typename T> 
    inline double get_val(const T& t, size_t) { return static_cast<double>(t); }


    // =========================================================================
    // CLASS DEFINITION
    // =========================================================================

    /**
     * @class THnCombi
     * @brief RDataFrame Action for filling THnSparse histograms.
     * * @tparam HasMask Boolean flag controlling execution mode.
     * - If `true`: The **LAST** argument passed to `Exec` is treated as an `RVec<int>` mask of valid indices.
     * - If `false`: All entries in the input vectors are processed (0 to N).
     */
    template <bool HasMask>
    class THnCombi : public ROOT::Detail::RDF::RActionImpl<THnCombi<HasMask>> {
    
    public:
        /// Result type required by RDataFrame (The generic N-D histogram)
        using Result_t = THnSparseD;
        /// Shared pointer to the histogram
        using Hist_ptr = std::shared_ptr<THnSparseD>;

        /**
         * @brief Constructor.
         * @details Clones the prototype histogram for every worker thread slot to ensure 
         * thread safety without the need for mutexes during the high-frequency fill loop.
         * @param proto_hist The template histogram (empty) to be filled.
         */
        THnCombi(const Hist_ptr& proto_hist);

        THnCombi(THnCombi&&) = default;
        THnCombi(const THnCombi&) = delete; // Actions are move-only in RDF

        /** @brief Returns the final merged result (Slot 0). */
        Hist_ptr GetResultPtr() const;
        
        /// Required by RActionImpl API: Per-slot initialization (No-op here)
        void Initialize();
        /// Required by RActionImpl API: Per-tree initialization (No-op here)
        void InitTask(TTreeReader*, unsigned int);
        /// Required by RActionImpl API: Returns action name for reporting
        std::string GetActionName();

        /**
         * @brief Merges histograms from all thread slots into the main slot (Index 0).
         * @details Called automatically by RDataFrame at the end of the event loop.
         */
        void Finalize();

        /**
         * @brief The Per-Event Execution Function.
         * @details 
         * This generic template handles the unpacking of RDataFrame columns.
         * It supports any number of input columns (Var + Splits + Optional Mask).
         * * @tparam Args Variadic column types (deduced automatically by RDF).
         * @param slot The thread slot index (0 to NThreads-1).
         * @param args The column values for the current event.
         */
        template <typename... Args>
        void Exec(unsigned int slot, const Args&... args);

    private:
        /// Vector of histograms, one per thread slot.
        std::vector<Hist_ptr> _hists;

        /**
         * @brief Helper: Unpacks tuple values into the coordinate vector 'x'.
         * @details Uses a fold expression to iterate over the tuple elements at compile time.
         * @param x Buffer to fill with values.
         * @param idx The current index in the combinatorial arrays.
         * @param t The tuple of input arguments.
         */
        template <typename Tuple, size_t... Is>
        void fill_row(std::vector<double>& x, int idx, const Tuple& t, std::index_sequence<Is...>);

        /**
         * @brief Helper: Determines the maximum size among the data columns.
         * @details Essential for "Unmasked" mode to know how many combinations exist in the event.
         */
        template <typename Tuple, size_t... Is>
        size_t get_max_size(const Tuple& t, std::index_sequence<Is...>);
    };

    // =========================================================================
    // IMPLEMENTATION
    // =========================================================================

    template <bool HasMask>
    inline THnCombi<HasMask>::THnCombi(const Hist_ptr& proto_hist) {
        // Determine number of slots based on MT configuration
        const auto nSlots = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;
        _hists.resize(nSlots);
        
        for (auto i : ROOT::TSeqU(nSlots)) {
            // Clone the histogram for this slot.
            // Note: THnSparse::Clone returns TObject*, must cast back to THnSparseD.
            if (proto_hist) {
                _hists[i] = std::shared_ptr<THnSparseD>(
                    dynamic_cast<THnSparseD*>(proto_hist->Clone())
                );
            }
            (void)i; // Suppress unused warning in generic loops
        }
    }

    template <bool HasMask>
    inline typename THnCombi<HasMask>::Hist_ptr THnCombi<HasMask>::GetResultPtr() const { 
        return _hists[0]; 
    }

    template <bool HasMask>
    inline void THnCombi<HasMask>::Initialize() {}

    template <bool HasMask>
    inline void THnCombi<HasMask>::InitTask(TTreeReader*, unsigned int) {}

    template <bool HasMask>
    inline std::string THnCombi<HasMask>::GetActionName() { 
        return "THnCombi"; 
    }

    template <bool HasMask>
    inline void THnCombi<HasMask>::Finalize() {
        // Merge results from all threads into the primary histogram (Slot 0)
        if (_hists.empty()) return;
        for (size_t i = 1; i < _hists.size(); ++i) {
            if (_hists[i]) {
                _hists[0]->Add(_hists[i].get());
            }
        }
    }

    // --- TEMPLATE EXEC ---
    template <bool HasMask>
    template <typename... Args>
    inline void THnCombi<HasMask>::Exec(unsigned int slot, const Args&... args) {
        // SAFETY GUARD: RDataFrame internal SFINAE checks may instantiate Exec with 0 args.
        // We must prevent compilation of the logic body in that case.
        if constexpr (sizeof...(Args) > 0) {
            
            // Pack variadic arguments into a tuple for indexed access
            auto all_args = std::forward_as_tuple(args...);
            
            // Calculate NData dynamically based on provided arguments.
            // If HasMask is true, the last argument is the mask, not data.
            constexpr size_t NArgs = sizeof...(Args);

            // Logic separation based on Masking mode
            if constexpr (HasMask) {
                // --- MASKED MODE ---
                // We expect at least 1 arg (the mask).
                static_assert(NArgs >= 1, "THnCombi<true>: Exec called with no mask column.");
                
                constexpr size_t NData = NArgs - 1;
                const auto& indices = std::get<NData>(all_args); // Get Mask (Last Arg)

                // Static buffer to avoid reallocation. 
                std::vector<double> x(NData);

                // indices can be RVec<int>, RVec<unsigned long>, etc.
                // We loop ONLY over the valid indices provided by the mask.
                for (auto idx : indices) {
                    fill_row(x, idx, all_args, std::make_index_sequence<NData>{});
                    _hists[slot]->Fill(x.data());
                }
            } 
            else {
                // --- UNMASKED MODE ---
                constexpr size_t NData = NArgs;
                
                // Process all entries. We must determine the loop size.
                // The size is the maximum size of any vector input (handling scalar broadcasting).
                size_t n_entries = get_max_size(all_args, std::make_index_sequence<NData>{});
                std::vector<double> x(NData);
                
                for (size_t i = 0; i < n_entries; ++i) {
                    fill_row(x, i, all_args, std::make_index_sequence<NData>{});
                    _hists[slot]->Fill(x.data());
                }
            }
        }
    }

    // --- Private Helper Implementations ---

    template <bool HasMask>
    template <typename Tuple, size_t... Is>
    inline void THnCombi<HasMask>::fill_row(std::vector<double>& x, int idx, const Tuple& t, std::index_sequence<Is...>) {
        int col = 0;
        // Fold expression: executes get_val for each column index I in Is, populating the vector x
        ((x[col++] = get_val(std::get<Is>(t), idx)), ...);
    }

    template <bool HasMask>
    template <typename Tuple, size_t... Is>
    inline size_t THnCombi<HasMask>::get_max_size(const Tuple& t, std::index_sequence<Is...>) {
        // Finds the maximum size among all data columns (broadcasts scalars which have size 1)
        return std::max({ get_size(std::get<Is>(t))... });
    }

  } // namespace histo
} // namespace rad
