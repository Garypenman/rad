#pragma once
#include "CommonDefines.h"
#include "Random.h"
#include <vector>
#include <memory>
#include <cmath>

namespace rad {

  using AuxCacheD = ROOT::RVec<ROOT::RVecD>; 
  using AuxCacheI = ROOT::RVec<ROOT::RVec<Long64_t>>;

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
			    ROOT::RVecD& pz, ROOT::RVecD& e,
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
		    ROOT::RVecD& pz, ROOT::RVecD& e,
		    const AuxCacheD&, const AuxCacheI&) const override;
  private:
    double _scale;
  };

  /**
   * @brief Scales the 3-momentum (px, py, pz) by a random fluctuation.
   * scale factor is centered on 1 with given width
   */
  class ModSmearMomentum : public ModifierBase {
  public:
    explicit ModSmearMomentum(double width);
    ModSmearMomentum(const Indices_t&, const Indices_t&, double width);

    std::unique_ptr<ModifierBase> Clone() const override;

    void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
		    ROOT::RVecD& pz, ROOT::RVecD& e,
		    const AuxCacheD&, const AuxCacheI&) const override;
  private:
    double _width;
  };
  /**
   * @brief Scales the 3-momentum (px, py, pz) by a random fluctuation.
   * assuming calorimeter energy measurement
   * scale factor is centered on 1 with given width
   */
  class ModCalSmearMomentum : public ModifierBase {
  public:
    explicit ModCalSmearMomentum(double width);
    ModCalSmearMomentum(const Indices_t&, const Indices_t&, double width);

    std::unique_ptr<ModifierBase> Clone() const override;

    void operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
		    ROOT::RVecD& pz, ROOT::RVecD& e,
		    const AuxCacheD&, const AuxCacheI&) const override;
  private:
    double _width;
  };

  /**
   * @brief Constraints the mass of a particle to a fixed value, preserving the measured momentum
   * Useful when the detector resolution suggests a specific particle (e.g. Proton)
   * but the mass spectrum is broad.
   */
  class ModFixMassWithP : public ModifierBase {
  public:
    explicit ModFixMassWithP(double mass) : _mass(mass) {}
    std::unique_ptr<ModifierBase> Clone() const override {
      return std::make_unique<ModFixMassWithP>(_mass);
    }

    void operator()(ROOT::RVecD& px, ROOT::RVecD& py,
                    ROOT::RVecD& pz, ROOT::RVecD& e,
                    const AuxCacheD&, const AuxCacheI&) const override
    {
      if (_target_idx < 0 || _target_idx >= static_cast<Indice_t>(px.size())) return;
      
      double px_i = px[_target_idx];
      double py_i = py[_target_idx];
      double pz_i = pz[_target_idx];
      double p2 = px_i*px_i + py_i*py_i + pz_i*pz_i;
      e[_target_idx] = std::sqrt(p2 + _mass*_mass);

    }

  private:
    double _mass;
  };


  /**
   * @brief Constraints the mass of a particle to a fixed value, preserving the measured energy
   * Useful when the detector resolution suggests a specific particle (e.g. Proton)
   * but the mass spectrum is broad.
   */
  class ModFixMassWithE : public ModifierBase {
  public:
    explicit ModFixMassWithE(double mass) : _mass(mass) {}
    std::unique_ptr<ModifierBase> Clone() const override {
      return std::make_unique<ModFixMassWithE>(_mass);
    }

    void operator()(ROOT::RVecD& px, ROOT::RVecD& py,
                    ROOT::RVecD& pz, ROOT::RVecD& e,
                    const AuxCacheD&, const AuxCacheI&) const override
    {
      if (_target_idx < 0 || _target_idx >= static_cast<Indice_t>(px.size())) return;

      double px_i = px[_target_idx];
      double py_i = py[_target_idx];
      double pz_i = pz[_target_idx];
      double e_i = e[_target_idx];

      double p2 = px_i*px_i + py_i*py_i + pz_i*pz_i;
      // Compute new momentum magnitude
      double new_p = std::sqrt(std::max(0.0, e_i*e_i - _mass*_mass));
      
      if (p2 < 1e-12) return; // avoid division by zero

      double scale = new_p / std::sqrt(p2);

      px[_target_idx] *= scale;
      py[_target_idx] *= scale;
      pz[_target_idx] *= scale;
    }
    
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
		    ROOT::RVecD& pz, ROOT::RVecD& e,
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
					   ROOT::RVecD& pz, ROOT::RVecD& e,
					   const AuxCacheD&, const AuxCacheI&) const 
  {
    if (_target_idx < 0 || _target_idx >= static_cast<Indice_t>(px.size())) return;
    px[_target_idx] *= _scale;
    py[_target_idx] *= _scale;
    pz[_target_idx] *= _scale;
    e[_target_idx] *= _scale;
    
  }
  // --- ModSmearMomentum ---
  inline ModSmearMomentum::ModSmearMomentum(double width) : _width(width) {}
  inline ModSmearMomentum::ModSmearMomentum(const Indices_t&, const Indices_t&, double width) : _width(width) {}
    
  inline std::unique_ptr<ModifierBase> ModSmearMomentum::Clone() const {
    return std::make_unique<ModSmearMomentum>(_width);
  }

  inline void ModSmearMomentum::operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
					   ROOT::RVecD& pz, ROOT::RVecD& e,
					   const AuxCacheD&, const AuxCacheI&) const 
  {
    if (_target_idx < 0 || _target_idx >= static_cast<Indice_t>(px.size())) return;
    auto scale = random::Generator().Gaus(1,_width);
    px[_target_idx] *= scale;
    py[_target_idx] *= scale;
    pz[_target_idx] *= scale;
    e[_target_idx] *= scale;
    
  }
  // --- ModCalSmearMomentum ---
  inline ModCalSmearMomentum::ModCalSmearMomentum(double width) : _width(width) {}
  inline ModCalSmearMomentum::ModCalSmearMomentum(const Indices_t&, const Indices_t&, double width) : _width(width) {}
    
  inline std::unique_ptr<ModifierBase> ModCalSmearMomentum::Clone() const {
    return std::make_unique<ModCalSmearMomentum>(_width);
  }

  inline void ModCalSmearMomentum::operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
					      ROOT::RVecD& pz, ROOT::RVecD& e,
					      const AuxCacheD&, const AuxCacheI&) const 
  {
    if (_target_idx < 0 || _target_idx >= static_cast<Indice_t>(px.size())) return;

    
    double E_old = e[_target_idx];
    double E_new = random::Generator().Gaus(E_old, _width * std::sqrt(E_old));
    
    if (std::isnan(E_old) || E_old <= 1e-9) return; 
    
    double scale = E_new / E_old;
    
    // auto kin_energy = TMath::Sqrt(px[_target_idx]*px[_target_idx]+py[_target_idx]*py[_target_idx]+pz[_target_idx]*pz[_target_idx] +m[_target_idx]*m[_target_idx])-m[_target_idx];
	
    // auto scale_K = random::Generator().Gaus(1,_width/TMath::Sqrt(kin_energy));
    // kin_energy *=scale_K;

    // auto tot_energy = kin_energy + m[_target_idx];
    // auto momentum = TMath::Sqrt(tot_energy*tot_energy - m[_target_idx]*m[_target_idx]);
    // auto scale = momentum/TMath::Sqrt(px[_target_idx]*px[_target_idx]+py[_target_idx]*py[_target_idx]+pz[_target_idx]*pz[_target_idx]);
	
    px[_target_idx] *= scale;
    py[_target_idx] *= scale;
    pz[_target_idx] *= scale;
    e[_target_idx] = E_new;
    
  }

  // // --- ModFixMass ---
  // inline ModFixMass::ModFixMassWithP(double mass) : _mass(mass) {}
  // inline ModFixMass::ModFixMassWithP(const Indices_t&, const Indices_t&, double mass) : _mass(mass) {}

  // inline std::unique_ptr<ModifierBase> ModFixMassWithP::Clone() const {
  //   return std::make_unique<ModFixMass>(_mass);
  // }

  // inline void ModFixMass::operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
  // 				     ROOT::RVecD& pz, ROOT::RVecD& e,
  // 				     const AuxCacheD&, const AuxCacheI&) const 
  // {
  //   if (_target_idx < 0 || _target_idx >= static_cast<Indice_t>(m.size())) return;
  //   m[_target_idx] = _mass;
  // }

  // --- ModSetMomFromAux ---
  inline ModSetMomFromAux::ModSetMomFromAux(Indice_t aux_row) : _row(aux_row) {}
    
  // Hardened: Handle empty input vector safely
  inline ModSetMomFromAux::ModSetMomFromAux(const Indices_t& dRows, const Indices_t&, int) 
    : _row(dRows.empty() ? -1 : dRows[0]) {}

  inline std::unique_ptr<ModifierBase> ModSetMomFromAux::Clone() const {
    return std::make_unique<ModSetMomFromAux>(_row);
  }

  inline void ModSetMomFromAux::operator()(ROOT::RVecD& px, ROOT::RVecD& py, 
					   ROOT::RVecD& pz, ROOT::RVecD& e,
					   const AuxCacheD& aux_d, const AuxCacheI&) const 
  {
    // 1. Valid Index Checks (Bounds Safety)
    if (_target_idx < 0 || _target_idx >= static_cast<Indice_t>(px.size())) return;
    if (_row < 0 || _row >= static_cast<Indice_t>(aux_d.size())) return; 
    if (_target_idx >= static_cast<Indice_t>(aux_d[_row].size())) return;

    // 2. Fetch new magnitude
    double new_mag = aux_d[_row][_target_idx];

    // 3. Safety Check: Skip if NaN or <= 0 (e.g. missing cluster)
    if (std::isnan(new_mag) || new_mag <= 1e-9) return; 

    // 4. Calculate Scale Factor
    double p2 = px[_target_idx]*px[_target_idx] + py[_target_idx]*py[_target_idx] + pz[_target_idx]*pz[_target_idx];
    if (std::isnan(p2) || p2 <= 1e-9) return; 

    double old_mag = std::sqrt(p2);
    double scale = new_mag / old_mag;

    // 5. Apply
    px[_target_idx] *= scale;
    py[_target_idx] *= scale;
    pz[_target_idx] *= scale;
    e[_target_idx] *= scale;
    
  }
}
