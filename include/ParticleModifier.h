#pragma once

#include "ParticleModifierMethods.h"
#include "ParticleCreator.h"
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <algorithm>

namespace rad {

    /**
     * @brief Manager class for applying kinematic corrections inside the event loop.
     * * The `ParticleModifier` acts as a container and executor for a sequence of modification strategies
     * (e.g., momentum scaling, mass fixing, energy calibration).
     */
    class ParticleModifier {
    public:
        ParticleModifier() = default;
        
        /**
         * @brief Deep Copy Constructor.
         * Required for the Twin Topology pattern (e.g., creating a Linked Processor).
         */
        ParticleModifier(const ParticleModifier& other);

        // --- Configuration ---

        /**
         * @brief Scales the momentum of a particle by a fixed factor.
         */
        void ScaleMomentum(const std::string& name, double scale);

        /**
         * @brief Forces the mass of a particle to a fixed value (modifies energy).
         */
        void FixMass(const std::string& name, double mass);

        /**
         * @brief Sets the particle's momentum magnitude from an auxiliary column.
         * Useful for correcting tracking momentum with calorimeter energy.
         */
        void SetMomentumFrom(const std::string& name, const std::string& colName);

        // --- Initialization & Execution ---

        /**
         * @brief Resolves particle names to fixed array indices.
         * Must be called after the ParticleCreator has finalized the ReactionMap.
         */
        void Init(const ParticleCreator& creator);

        /**
         * @brief Apply all registered modifications to the kinematic arrays.
         * This is called inside the high-performance event loop.
         */
        void Apply(ROOT::RVecD& px, ROOT::RVecD& py, 
                   ROOT::RVecD& pz, ROOT::RVecD& m, 
                   const AuxCacheD& aux_d, const AuxCacheI& aux_i) const;

        // --- Accessors ---

        const std::vector<std::string>& GetAuxDoubleCols() const;
        const std::vector<std::string>& GetAuxIntCols() const;
        Indice_t RegisterAuxDouble(const std::string& name);
        Indice_t RegisterAuxInt(const std::string& name);

    private:
        struct ModiConfig {
            std::string target_name;
            std::shared_ptr<ModifierBase> modifier;
        };

        std::vector<ModiConfig> _modi_configs; 
        std::vector<std::shared_ptr<ModifierBase>> _active_modifiers; 
        
        std::vector<std::string> _aux_double_cols; 
        std::vector<std::string> _aux_int_cols;    
    };

    // =================================================================================
    // IMPLEMENTATION
    // =================================================================================

    inline ParticleModifier::ParticleModifier(const ParticleModifier& other) 
        : _modi_configs(), // Deep copy manual step below
          _active_modifiers(), // Active modifiers are rebuilt via Init()
          _aux_double_cols(other._aux_double_cols),
          _aux_int_cols(other._aux_int_cols) 
    {
        // Deep copy the pending configurations
        for(const auto& conf : other._modi_configs) {
            _modi_configs.push_back({conf.target_name, conf.modifier->Clone()});
        }
    }

    inline void ParticleModifier::ScaleMomentum(const std::string& name, double scale) {
        auto mod = std::make_unique<ModScaleMomentum>(scale);
        _modi_configs.push_back({name, std::move(mod)});
    }

    inline void ParticleModifier::FixMass(const std::string& name, double mass) {
        auto mod = std::make_unique<ModFixMass>(mass);
        _modi_configs.push_back({name, std::move(mod)});
    }

    inline void ParticleModifier::SetMomentumFrom(const std::string& name, const std::string& colName) {
        Indice_t row = RegisterAuxDouble(colName);
        auto mod = std::make_unique<ModSetMomFromAux>(row);
        _modi_configs.push_back({name, std::move(mod)});
    }

    inline void ParticleModifier::Init(const ParticleCreator& creator) {
        _active_modifiers.clear();
        for(const auto& conf : _modi_configs) {
            // Clone the modifier to move it to the active list
            auto mod = conf.modifier->Clone();
            
            // Resolve the name to an index
            if(!creator.HasParticle(conf.target_name)) {
                std::cerr << "ParticleModifier Warning: Target '" << conf.target_name 
                          << "' not found in Creator. Skipping modification." << std::endl;
                continue;
            }
            mod->SetIndex(creator.GetReactionIndex(conf.target_name));
            _active_modifiers.push_back(std::move(mod));
        }
    }

    inline void ParticleModifier::Apply(ROOT::RVecD& px, ROOT::RVecD& py, 
                   ROOT::RVecD& pz, ROOT::RVecD& m, 
                   const AuxCacheD& aux_d, const AuxCacheI& aux_i) const 
    {
        for(const auto& mod : _active_modifiers) {
            (*mod)(px, py, pz, m, aux_d, aux_i);
        }
    }

    inline const std::vector<std::string>& ParticleModifier::GetAuxDoubleCols() const { return _aux_double_cols; }
    inline const std::vector<std::string>& ParticleModifier::GetAuxIntCols() const { return _aux_int_cols; }

    inline Indice_t ParticleModifier::RegisterAuxDouble(const std::string& name) {
        auto it = std::find(_aux_double_cols.begin(), _aux_double_cols.end(), name);
        if(it != _aux_double_cols.end()) return std::distance(_aux_double_cols.begin(), it);
        _aux_double_cols.push_back(name);
        return _aux_double_cols.size() - 1;
    }

    inline Indice_t ParticleModifier::RegisterAuxInt(const std::string& name) {
        auto it = std::find(_aux_int_cols.begin(), _aux_int_cols.end(), name);
        if(it != _aux_int_cols.end()) return std::distance(_aux_int_cols.begin(), it);
        _aux_int_cols.push_back(name);
        return _aux_int_cols.size() - 1;
    }

} // namespace rad
