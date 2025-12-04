#pragma once

#include "ParticleModifierMethods.h"
#include "ParticleCreator.h"
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <algorithm>

namespace rad {
namespace config {

    class ParticleModifier {
    public:
        ParticleModifier() = default;
        
        // Copy Constructor
        ParticleModifier(const ParticleModifier& other) 
            : _aux_double_cols(other._aux_double_cols),
              _aux_int_cols(other._aux_int_cols) 
        {
            for(const auto& cfg : other._modi_configs) {
                _modi_configs.push_back({cfg.target_name, cfg.modifier->Clone()});
            }
        }

        void Init(const ParticleCreator& creator) {
            _active_modifiers.clear();
            for(auto& cfg : _modi_configs) {
                if(!creator.HasParticle(cfg.target_name)) {
                    continue; 
                }
                
                Indice_t idx = creator.GetReactionIndex(cfg.target_name);
                cfg.modifier->SetIndex(idx);
                
                _active_modifiers.push_back(std::move(cfg.modifier));
            }
            // _modi_configs are preserved for potential re-Init/Copy
        }

        void Apply(ROOT::RVecD& px, ROOT::RVecD& py, ROOT::RVecD& pz, ROOT::RVecD& m,
                   const AuxCacheD& aux_d, const AuxCacheI& aux_i) const 
        {
            for(const auto& mod : _active_modifiers) {
                (*mod)(px, py, pz, m, aux_d, aux_i);
            }
        }

        // =====================================================================
        // Auxiliary Column Registration
        // =====================================================================

        Indice_t RegisterAuxDouble(const std::string& colName) {
            auto it = std::find(_aux_double_cols.begin(), _aux_double_cols.end(), colName);
            if(it != _aux_double_cols.end()) return (Indice_t)std::distance(_aux_double_cols.begin(), it);
            
            _aux_double_cols.push_back(colName);
            return (Indice_t)(_aux_double_cols.size() - 1);
        }

        Indice_t RegisterAuxInt(const std::string& colName) {
            auto it = std::find(_aux_int_cols.begin(), _aux_int_cols.end(), colName);
            if(it != _aux_int_cols.end()) return (Indice_t)std::distance(_aux_int_cols.begin(), it);
            
            _aux_int_cols.push_back(colName);
            return (Indice_t)(_aux_int_cols.size() - 1);
        }

        // Helper: Batch register a list of columns
        // Returns Indices_t (RVecI) instead of std::vector<int>
        Indices_t RegisterList(const std::vector<std::string>& cols, bool is_double) {
            Indices_t indices;
            indices.reserve(cols.size());
            for(const auto& c : cols) {
                if(is_double) indices.push_back(RegisterAuxDouble(c));
                else          indices.push_back(RegisterAuxInt(c));
            }
            return indices;
        }

        const std::vector<std::string>& GetAuxDoubleCols() const { return _aux_double_cols; }
        const std::vector<std::string>& GetAuxIntCols() const { return _aux_int_cols; }

        // =====================================================================
        // Generic Custom Factory
        // =====================================================================
        
        /**
         * @brief Registers any Custom Modifier class.
         * The class T must have a constructor: T(const Indices_t& dRows, const Indices_t& iRows, Args...).
         */
        template <typename T, typename... Args>
        void AddCustom(const std::string& name, 
                       const std::vector<std::string>& auxDoubleCols, 
                       const std::vector<std::string>& auxIntCols, 
                       Args&&... args) 
        {
            // 1. Register Columns & Get Row Indices (as Indices_t)
            auto dRows = RegisterList(auxDoubleCols, true);
            auto iRows = RegisterList(auxIntCols, false);

            // 2. Construct the Strategy
            auto mod = std::make_unique<T>(dRows, iRows, std::forward<Args>(args)...);
            
            // 3. Store
            _modi_configs.push_back({name, std::move(mod)});
        }

        // =====================================================================
        // Convenience Wrappers
        // =====================================================================

        void ScaleMomentum(const std::string& name, double scale) {
            auto mod = std::make_unique<ModScaleMomentum>(scale);
            _modi_configs.push_back({name, std::move(mod)});
        }

        void FixMass(const std::string& name, double mass) {
            auto mod = std::make_unique<ModFixMass>(mass);
            _modi_configs.push_back({name, std::move(mod)});
        }

        void SetMomentumFrom(const std::string& name, const std::string& colName) {
            Indice_t row = RegisterAuxDouble(colName);
            auto mod = std::make_unique<ModSetMomFromAux>(row);
            _modi_configs.push_back({name, std::move(mod)});
        }

    private:
        struct ModiConfig {
            std::string target_name;
            std::unique_ptr<ModifierBase> modifier;
        };

        std::vector<ModiConfig> _modi_configs; 
        std::vector<std::unique_ptr<ModifierBase>> _active_modifiers; 
        
        std::vector<std::string> _aux_double_cols;
        std::vector<std::string> _aux_int_cols;
    };

}} // namespace
