#pragma once
#include "CommonDefines.h"
#include <vector>
#include <memory>
#include <cmath>

namespace rad {

    using AuxCacheD = std::vector<ROOT::RVecD>; 
    using AuxCacheI = std::vector<ROOT::RVec<Long64_t>>;

    /**
     * @class ModifierBase
     * @brief Abstract Base Class for Momentum Modifiers.
     */
    class ModifierBase {
    public:
        virtual ~ModifierBase() = default;
        virtual std::unique_ptr<ModifierBase> Clone() const = 0;

        void SetIndex(Indice_t idx) { _target_idx = idx; }
        Indice_t GetIndex() const { return _target_idx; }

        virtual void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
                                ROOT::RVecD& pz, ROOT::RVecD& m,
                                const AuxCacheD& aux_d, 
                                const AuxCacheI& aux_i) const = 0;

    protected:
        Indice_t _target_idx = -1; 
    };

    // =========================================================
    // Standard Strategies
    // =========================================================

    /**
     * @brief Scales the 3-momentum (px, py, pz) by a fixed factor.
     * Useful for simple ad-hoc energy calibration.
     */
    class ModScaleMomentum : public ModifierBase {
    public:
        explicit ModScaleMomentum(double scale);
        ModScaleMomentum(const Indices_t&, const Indices_t&, double scale);

        std::unique_ptr<ModifierBase> Clone() const override;

        void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
                        ROOT::RVecD& pz, ROOT::RVecD& m,
                        const AuxCacheD&, const AuxCacheI&) const override;
    private:
        double _scale;
    };

    /**
     * @brief Constraints the mass of a particle to a fixed value.
     * Useful when the detector resolution suggests a specific particle (e.g. Proton)
     * but the mass spectrum is broad.
     */
    class ModFixMass : public ModifierBase {
    public:
        explicit ModFixMass(double mass);
        ModFixMass(const Indices_t&, const Indices_t&, double mass);

        std::unique_ptr<ModifierBase> Clone() const override;

        void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
                        ROOT::RVecD& pz, ROOT::RVecD& m,
                        const AuxCacheD&, const AuxCacheI&) const override;
    private:
        double _mass;
    };

    /**
     * @brief Overwrites momentum magnitude using a value from an auxiliary column.
     * Example: Replace track momentum with calorimeter energy for electrons.
     */
    class ModSetMomFromAux : public ModifierBase {
    public:
        explicit ModSetMomFromAux(Indice_t aux_row);
        ModSetMomFromAux(const Indices_t& dRows, const Indices_t&, int=0);

        std::unique_ptr<ModifierBase> Clone() const override;

        void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
                        ROOT::RVecD& pz, ROOT::RVecD& m,
                        const AuxCacheD& aux_d, const AuxCacheI&) const override;
    private:
        Indice_t _row;
    };

    // =========================================================
    // IMPLEMENTATION
    // =========================================================

    // --- ModScaleMomentum ---
    inline ModScaleMomentum::ModScaleMomentum(double scale) : _scale(scale) {}
    inline ModScaleMomentum::ModScaleMomentum(const Indices_t&, const Indices_t&, double scale) : _scale(scale) {}
    
    inline std::unique_ptr<ModifierBase> ModScaleMomentum::Clone() const {
        return std::make_unique<ModScaleMomentum>(_scale);
    }

    inline void ModScaleMomentum::operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
                    ROOT::RVecD& pz, ROOT::RVecD& m,
                    const AuxCacheD&, const AuxCacheI&) const 
    {
        px[_target_idx] *= _scale;
        py[_target_idx] *= _scale;
        pz[_target_idx] *= _scale;
    }

    // --- ModFixMass ---
    inline ModFixMass::ModFixMass(double mass) : _mass(mass) {}
    inline ModFixMass::ModFixMass(const Indices_t&, const Indices_t&, double mass) : _mass(mass) {}

    inline std::unique_ptr<ModifierBase> ModFixMass::Clone() const {
        return std::make_unique<ModFixMass>(_mass);
    }

    inline void ModFixMass::operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
                    ROOT::RVecD& pz, ROOT::RVecD& m,
                    const AuxCacheD&, const AuxCacheI&) const 
    {
        m[_target_idx] = _mass;
    }

    // --- ModSetMomFromAux ---
    inline ModSetMomFromAux::ModSetMomFromAux(Indice_t aux_row) : _row(aux_row) {}
    inline ModSetMomFromAux::ModSetMomFromAux(const Indices_t& dRows, const Indices_t&, int) : _row(dRows.at(0)) {}

    inline std::unique_ptr<ModifierBase> ModSetMomFromAux::Clone() const {
        return std::make_unique<ModSetMomFromAux>(_row);
    }

    inline void ModSetMomFromAux::operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
                    ROOT::RVecD& pz, ROOT::RVecD& m,
                    const AuxCacheD& aux_d, const AuxCacheI&) const 
    {
        double p2 = px[_target_idx]*px[_target_idx] + py[_target_idx]*py[_target_idx] + pz[_target_idx]*pz[_target_idx];
        if(p2 == 0) return;
        double old_mag = std::sqrt(p2);

        // Fetch new magnitude from the Aux Cache
        double new_mag = aux_d[_row][_target_idx];

        double scale = new_mag / old_mag;
        px[_target_idx] *= scale;
        py[_target_idx] *= scale;
        pz[_target_idx] *= scale;
    }

} // namespace rad


// #pragma once
// #include "CommonDefines.h"
// #include <vector>
// #include <memory>
// #include <cmath>

// namespace rad {

//     // Aliases for the Local Auxiliary Cache
//     // Structure: [Variable_Index][Particle_Index]
//     using AuxCacheD = std::vector<ROOT::RVecD>; 
//     using AuxCacheI = std::vector<ROOT::RVec<Long64_t>>;

//     /**
//      * @brief Abstract Base Class for Momentum Modifiers.
//      */
//     class ModifierBase {
//     public:
//         virtual ~ModifierBase() = default;

//         virtual std::unique_ptr<ModifierBase> Clone() const = 0;

//         /**
//          * @brief Sets the fixed array index of the particle to modify.
//          */
//         void SetIndex(Indice_t idx) { _target_idx = idx; }
//         Indice_t GetIndex() const { return _target_idx; }

//         /**
//          * @brief Execute the modification on the kinematic cache.
//          */
//         virtual void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
//                                 ROOT::RVecD& pz, ROOT::RVecD& m,
//                                 const AuxCacheD& aux_d, 
//                                 const AuxCacheI& aux_i) const = 0;

//     protected:
//         Indice_t _target_idx = -1; // The fixed index in the arrays to modify
//     };

//     // =========================================================
//     // Standard Strategies
//     // =========================================================

//     class ModScaleMomentum : public ModifierBase {
//     public:
//         explicit ModScaleMomentum(double scale) : _scale(scale) {}

//         // Factory Constructor: Accepts lists of indices (ignored here)
//         ModScaleMomentum(const Indices_t&, const Indices_t&, double scale) 
//             : _scale(scale) {}

//         std::unique_ptr<ModifierBase> Clone() const override {
//             return std::make_unique<ModScaleMomentum>(_scale);
//         }

//         void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
//                         ROOT::RVecD& pz, ROOT::RVecD& m,
//                         const AuxCacheD&, const AuxCacheI&) const override 
//         {
//             px[_target_idx] *= _scale;
//             py[_target_idx] *= _scale;
//             pz[_target_idx] *= _scale;
//         }
//     private:
//         double _scale;
//     };

//     class ModFixMass : public ModifierBase {
//     public:
//         explicit ModFixMass(double mass) : _mass(mass) {}

//         ModFixMass(const Indices_t&, const Indices_t&, double mass) 
//             : _mass(mass) {}

//         std::unique_ptr<ModifierBase> Clone() const override {
//             return std::make_unique<ModFixMass>(_mass);
//         }

//         void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
//                         ROOT::RVecD& pz, ROOT::RVecD& m,
//                         const AuxCacheD&, const AuxCacheI&) const override 
//         {
//             m[_target_idx] = _mass;
//         }
//     private:
//         double _mass;
//     };

//     class ModSetMomFromAux : public ModifierBase {
//     public:
//         explicit ModSetMomFromAux(Indice_t aux_row) : _row(aux_row) {}

//         // Factory Constructor: Expects 1 Double Column index in the list
//         ModSetMomFromAux(const Indices_t& dRows, const Indices_t&, int=0) 
//             : _row(dRows.at(0)) {}

//         std::unique_ptr<ModifierBase> Clone() const override {
//             return std::make_unique<ModSetMomFromAux>(_row);
//         }

//         void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
//                         ROOT::RVecD& pz, ROOT::RVecD& m,
//                         const AuxCacheD& aux_d, const AuxCacheI&) const override 
//         {
//             double p2 = px[_target_idx]*px[_target_idx] + py[_target_idx]*py[_target_idx] + pz[_target_idx]*pz[_target_idx];
//             if(p2 == 0) return;
//             double old_mag = std::sqrt(p2);

//             // Fetch new magnitude from the Aux Cache
//             double new_mag = aux_d[_row][_target_idx];

//             double scale = new_mag / old_mag;
//             px[_target_idx] *= scale;
//             py[_target_idx] *= scale;
//             pz[_target_idx] *= scale;
//         }
//     private:
//         Indice_t _row;
//     };

// }
