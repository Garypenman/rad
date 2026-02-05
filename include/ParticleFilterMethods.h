/**
 * @file ParticleFilterMethods.h
 * @brief Modifier strategies that invalidate particles based on auxiliary conditions.
 * @details
 * These classes inherit from ModifierBase and fit directly into the ParticleModifier
 * execution chain. If a filter condition fails, the particle's 4-vector is set to 
 * InvalidEntry (NaN), which signals the KinematicsProcessor to skip the combination.
 */

#pragma once

#include "ParticleModifierMethods.h"
#include <cmath>

namespace rad {

    /**
     * @class FilterBase
     * @brief Shared logic for all filters.
     */
    class FilterBase : public ModifierBase {
    public:
        // Helper to kill a particle (set to NaN)
        void Invalidate(Indice_t idx, ROOT::RVecD& px, ROOT::RVecD& py, 
                        ROOT::RVecD& pz, ROOT::RVecD& m) const 
        {
            px[idx] = rad::consts::InvalidEntry<double>();
            py[idx] = rad::consts::InvalidEntry<double>();
            pz[idx] = rad::consts::InvalidEntry<double>();
            m[idx]  = rad::consts::InvalidEntry<double>();
        }
        
        // Filters are clonable just like any Modifier
        virtual std::unique_ptr<ModifierBase> Clone() const = 0;
    };

    // =========================================================
    // Generic Filter Strategies
    // =========================================================

    /**
     * @class FilterRange
     * @brief Keeps particle if Aux variable is within [min, max].
     * @note Works on Double columns.
     */
    class FilterRange : public FilterBase {
    public:
        /** @brief Factory Ctor: params[0]=min, params[1]=max. */
        FilterRange(const Indices_t& dRows, const Indices_t&, const ROOT::RVec<double>& params);
        // Direct Ctor
        FilterRange(Indice_t row, double min, double max);

        std::unique_ptr<ModifierBase> Clone() const override;

        void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
                        ROOT::RVecD& pz, ROOT::RVecD& m,
                        const AuxCacheD& aux_d, const AuxCacheI&) const override;
    private:
        Indice_t _row;
        double _min;
        double _max;
    };

    /**
     * @class FilterExact
     * @brief Keeps particle if Aux variable equals target value.
     * @note Works on Int/Long columns (e.g., Detector ID).
     */
    class FilterExact : public FilterBase {
    public:
        /** @brief Factory Ctor: params[0]=target_val. */
        FilterExact(const Indices_t&, const Indices_t& iRows, const ROOT::RVec<double>& params);
        // Direct Ctor
        FilterExact(Indice_t row, long long target);

        std::unique_ptr<ModifierBase> Clone() const override;

        void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
                        ROOT::RVecD& pz, ROOT::RVecD& m,
                        const AuxCacheD&, const AuxCacheI& aux_i) const override;
    private:
        Indice_t _row;
        long long _target;
    };

    /**
     * @class FilterBitCheck
     * @brief Keeps particle if (val & mask) == result.
     * @note Works on Int/Long columns (e.g., Status Bits).
     */
    class FilterBitCheck : public FilterBase {
    public:
        /** @brief Factory Ctor: params[0]=mask, params[1]=result(default 0). */
        FilterBitCheck(const Indices_t&, const Indices_t& iRows, const ROOT::RVec<double>& params);
        // Direct Ctor
        FilterBitCheck(Indice_t row, long long mask, long long result);

        std::unique_ptr<ModifierBase> Clone() const override;

        void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
                        ROOT::RVecD& pz, ROOT::RVecD& m,
                        const AuxCacheD&, const AuxCacheI& aux_i) const override;
    private:
        Indice_t _row;
        long long _mask;
        long long _result;
    };

    // =========================================================
    // IMPLEMENTATION
    // =========================================================

    // --- FilterRange ---
    inline FilterRange::FilterRange(const Indices_t& dRows, const Indices_t&, const ROOT::RVec<double>& params) 
        : _row(dRows.empty() ? -1 : dRows[0]), 
          _min(params.size() > 0 ? params[0] : -1e99),
          _max(params.size() > 1 ? params[1] : 1e99) {}

    inline FilterRange::FilterRange(Indice_t row, double min, double max) 
        : _row(row), _min(min), _max(max) {}

    inline std::unique_ptr<ModifierBase> FilterRange::Clone() const {
        return std::make_unique<FilterRange>(_row, _min, _max);
    }

    inline void FilterRange::operator()(ROOT::RVecD& px, ROOT::RVecD& py, ROOT::RVecD& pz, ROOT::RVecD& m,
                                        const AuxCacheD& aux_d, const AuxCacheI&) const 
    {
        if (_target_idx < 0 || _row < 0) return;
        if (_target_idx >= static_cast<Indice_t>(aux_d[_row].size())) return; // Safety

        double val = aux_d[_row][_target_idx];
        
        // Check condition: if fail, invalidate
        if (rad::consts::IsInvalidEntry(val) || val < _min || val > _max) {
            Invalidate(_target_idx, px, py, pz, m);
        }
    }

    // --- FilterExact ---
    inline FilterExact::FilterExact(const Indices_t&, const Indices_t& iRows, const ROOT::RVec<double>& params)
        : _row(iRows.empty() ? -1 : iRows[0]), 
          _target(params.empty() ? 0 : static_cast<long long>(params[0])) {}

    inline FilterExact::FilterExact(Indice_t row, long long target) 
        : _row(row), _target(target) {}

    inline std::unique_ptr<ModifierBase> FilterExact::Clone() const {
        return std::make_unique<FilterExact>(_row, _target);
    }

    inline void FilterExact::operator()(ROOT::RVecD& px, ROOT::RVecD& py, ROOT::RVecD& pz, ROOT::RVecD& m,
                                        const AuxCacheD&, const AuxCacheI& aux_i) const 
    {
        if (_target_idx < 0 || _row < 0) return;
        if (_target_idx >= static_cast<Indice_t>(aux_i[_row].size())) return;

        long long val = aux_i[_row][_target_idx];
        
        if (val != _target) {
            Invalidate(_target_idx, px, py, pz, m);
        }
    }

    // --- FilterBitCheck ---
    inline FilterBitCheck::FilterBitCheck(const Indices_t&, const Indices_t& iRows, const ROOT::RVec<double>& params)
        : _row(iRows.empty() ? -1 : iRows[0]), 
          _mask(params.size() > 0 ? static_cast<long long>(params[0]) : 0),
          _result(params.size() > 1 ? static_cast<long long>(params[1]) : 0) {} // Default checks if masked bit is 0

    inline FilterBitCheck::FilterBitCheck(Indice_t row, long long mask, long long result) 
        : _row(row), _mask(mask), _result(result) {}

    inline std::unique_ptr<ModifierBase> FilterBitCheck::Clone() const {
        return std::make_unique<FilterBitCheck>(_row, _mask, _result);
    }

    inline void FilterBitCheck::operator()(ROOT::RVecD& px, ROOT::RVecD& py, ROOT::RVecD& pz, ROOT::RVecD& m,
                                           const AuxCacheD&, const AuxCacheI& aux_i) const 
    {
        if (_target_idx < 0 || _row < 0) return;
        if (_target_idx >= static_cast<Indice_t>(aux_i[_row].size())) return;

        long long val = aux_i[_row][_target_idx];

        // Example: (Status & 0x4) == 0x4
        if ((val & _mask) != _result) {
            Invalidate(_target_idx, px, py, pz, m);
        }
    }

} // namespace rad
