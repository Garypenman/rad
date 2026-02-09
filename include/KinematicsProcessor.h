/**
 * @file KinematicsProcessor.h
 * @brief The Computational Engine for combinatorial kinematic analysis.
 * @details
 * This class orchestrates the combinatorial logic. It:
 * 1. Consumes the combinatorial indices from the `ParticleCreator`.
 * 2. Runs the "Event Loop" (via `operator()`) to generate 4-vectors for every valid combination.
 * 3. Provides an API to register calculations (Mass, Pt, etc.) on these combinations.
 * 4. Prepares flattened data columns for `SnapshotCombi`.
 */

#pragma once

#include "ConfigReaction.h"
#include "StringUtilities.h"
#include "RVecHelpers.h"
#include "ParticleCreator.h"
#include "KinematicsDispatch.h" // Contains DefineKinematicsProcessor template
#include "KineCalculation.h" 
#include "BasicKinematics.h"
#include "ParticleModifier.h" 

#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <iomanip>

namespace rad {

  using ROOT::RVec;
  using RVecRVecD = ROOT::RVec<ROOT::RVecD>;
  using RVecRVecI = ROOT::RVec<ROOT::RVecI>;

  /**
   * @class KinematicsProcessor
   * @brief Manages the kinematic calculations and RDF graph construction.
   * @details 
   * This class serves as the bridge between the physics topology (defined in `ParticleCreator`)
   * and the RDataFrame execution graph. It manages a registry of calculations 
   * so that the `AnalysisManager` knows exactly which columns to write to disk.
   */
  class KinematicsProcessor {

  public:
    /// Type alias for the complex nested vector structure [Combination][Component][Particle]
    using CombiOutputVec_t = RVec<RVec<RVecResultType>>;

    // =================================================================================
    // Lifecycle & Initialization
    // =================================================================================

    /**
     * @brief Constructor.
     * @param cr Pointer to the Reaction configuration.
     * @param prefix The data prefix (e.g., "rec_", "tru_").
     * @param suffix Optional suffix for systematics (e.g., "_sysUp").
     */
    KinematicsProcessor(ConfigReaction* cr, const std::string& prefix, const std::string& suffix = "");
    
    virtual ~KinematicsProcessor() = default;

    /** * @brief Initializes the processor and defines columns in the RDF graph. 
     * @details This function is idempotent (can be called multiple times safely).
     * It triggers the `ParticleCreator` setup, defines the main execution kernel via
     * `DefineKinematicsProcessor`, and runs definitions for all registered calculations.
     */
    void Init();
    
    // =================================================================================
    // Core Execution (Functor)
    // =================================================================================

    /**
     * @brief The main computational kernel.
     * @details Generates Px, Py, Pz, M for every combination defined in `indices`.
     * Applies Pre/Post modifiers (Momentum corrections, Energy loss, etc.).
     * * @tparam Tp Type of Momentum input (RVecF or RVecD).
     * @tparam Tm Type of Mass input (RVecF or RVecD).
     * @param indices The combinatorial indices [Particle][Combination].
     * @param px Input Px vector.
     * @param py Input Py vector.
     * @param pz Input Pz vector.
     * @param m Input Mass vector.
     * @param aux_pre_d Auxiliary Doubles (Pre-creation).
     * @param aux_pre_i Auxiliary Ints (Pre-creation).
     * @param aux_post_d Auxiliary Doubles (Post-creation).
     * @param aux_post_i Auxiliary Ints (Post-creation).
     * @return The nested structure containing kinematics for all combinations.
     */
    template<typename Tp, typename Tm> 
    CombiOutputVec_t operator()(const RVecIndices& indices, 
                                const Tp& px, const Tp& py, const Tp& pz, const Tm& m,
                                const RVecRVecD& aux_pre_d, const RVecRVecI& aux_pre_i,
                                const RVecRVecD& aux_post_d, const RVecRVecI& aux_post_i) const;

    /** * @brief Defines flattened columns for Px, Py, Pz, M for every particle. 
     * @details 
     * Essential for `SnapshotCombi`. This creates scalar RVecs (one entry per combination)
     * effectively flattening the nested `CombiOutputVec_t` structure.
     * It registers these variables (e.g., "ele_px") so `AnalysisManager` knows to save them.
     */
    void DefineNewComponentVecs();

    // =================================================================================
    // Configuration Helpers
    // =================================================================================

    /** @brief Defines the Auxiliary data pack columns in the RDF graph. */
    void DefineAux();

    /** @brief Define a generic particle group override. */
    void SetGroup(const std::string& groupName, const ROOT::RVec<std::string>& particles);

    /** @brief Override the standard "Mesons" group. */
    void SetMesonParticles(const ROOT::RVec<std::string>& particles);

    /** @brief Override the standard "Baryons" group. */
    void SetBaryonParticles(const ROOT::RVec<std::string>& particles);

    /** @brief Re-target this processor to a new Reaction (used internally). */
    void SetReaction(ConfigReaction* cr) { 
      _reaction = cr; 
      _creator.SetReaction(cr); 
    }
    
    // =================================================================================
    // Accessors
    // =================================================================================

    ParticleCreator& Creator();
    const ParticleCreator& Creator() const;
    
    ParticleModifier& PreModifier();
    ParticleModifier& PostModifier();

    ConfigReaction* Reaction() const;
    std::string GetSuffix() const;
    std::string GetPrefix() const;
    
    /** @brief Construct a full column name: prefix + base + suffix. */
    std::string FullName(const std::string& baseName) const;

    /** * @brief Get list of all variables defined via RegisterCalc/Mass/Pt etc.
     * @details This is used by `AnalysisManager::Snapshot` to auto-detect columns.
     */
    const ROOT::RVec<std::string>& GetDefinedNames() const { 
      return _registered_vars; 
    }

    // =================================================================================
    // Definition API & Calculation Registration
    // =================================================================================
    
    /** * @brief Registers a calculation based on the RVecIndexMap (e.g. Mass, Missing Mass).
     * @param name The output variable name (base name only).
     * @param func The kernel function to execute.
     */
    void RegisterCalc(const std::string& name, KineCalculation::MapKernel func);

    /** * @brief Registers a calculation based on explicit particle indices. 
     * @param name The output variable name.
     * @param func The kernel function.
     * @param particles List of particle names needed by the kernel.
     */
    void RegisterCalc(const std::string& name, KineCalculation::IndexKernel func, ROOT::RVec<ParticleNames_t> particles);

    // --- Legacy String-Based Definition ---
    
    void Define(const std::string& name, const std::string& func);
    void Define(const std::string& name, const std::string& func, const ROOT::RVec<ParticleNames_t>& particles);

    /** @brief Helper to define an arbitrary lambda kernel in RDF. */
    template <typename Lambda>
    void DefineKernel(const std::string& name, Lambda&& func);

    // =================================================================================
    // Physics Shortcuts
    // =================================================================================

    void Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});
    void Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});
    void Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

    void ParticleTheta(const ParticleNames_t& particles);
    void ParticlePhi(const ParticleNames_t& particles);
    void ParticleP(const ParticleNames_t& particles);
    void ParticleEta(const ParticleNames_t& particles);
    
    void PrintReactionMap() const;


    //////////////////////
    ///Diagnostics
    ///\brief Print registered calculations.
    void PrintCalculations() const;

    ///\brief Print registered output variables.
    void PrintRegisteredVariables() const;

    ///\brief Print group overrides.
    void PrintGroupOverrides() const;

    ///\brief Print comprehensive processor diagnostics.
    void PrintProcessorDiagnostics() const;

  private:
    ConfigReaction* _reaction = nullptr;
    std::string _prefix; 
    std::string _suffix; 
    
    bool _isInitialized = false; 
    
    ParticleCreator _creator;
    ParticleModifier _preModifier;
    ParticleModifier _postModifier;

    ROOT::RVec<KineCalculation> _calculations;
    ROOT::RVec<std::string> _registered_vars; // Registry of variable names for Snapshot
    
    struct GroupOverride {
        std::string name;
        ROOT::RVec<std::string> particles;
    };
    ROOT::RVec<GroupOverride> _groupOverrides;

       void ApplyGroupOverrides();
  };

  // =================================================================================
  // IMPLEMENTATION: KinematicsProcessor
  // =================================================================================

  inline KinematicsProcessor::KinematicsProcessor(ConfigReaction* cr, const std::string& prefix, const std::string& suffix) 
    : _reaction{cr}, _prefix{prefix}, _suffix{suffix}, _creator{cr, prefix, suffix} 
  {}

  inline void KinematicsProcessor::Init() {
    if (_isInitialized) return; 
    _isInitialized = true;

    _reaction->ValidateType(_prefix);
    // Clear the registry to ensure no stale/duplicate entries from config phase
    _registered_vars.clear();

    ApplyGroupOverrides();
    Creator().InitMap(); 
    _preModifier.Init(Creator());
    _postModifier.Init(Creator());
    DefineAux();

    // Calls KinematicsDispatch::DefineKinematicsProcessor
    // IMPORTANT: This internally calls this->DefineNewComponentVecs() !
    DefineKinematicsProcessor(*_reaction, *this, _prefix);

    for(auto& calc : _calculations) {
      calc.Define(this); 
    }
  }

  inline void KinematicsProcessor::ApplyGroupOverrides() {
      for(const auto& group : _groupOverrides) {
          Creator().OverrideGroup(group.name, group.particles); 
      }
  }

  inline void KinematicsProcessor::SetGroup(const std::string& groupName, const ROOT::RVec<std::string>& particles) {
      _groupOverrides.push_back({groupName, particles});
  }

  inline void KinematicsProcessor::SetMesonParticles(const ROOT::RVec<std::string>& particles) {
      SetGroup(rad::consts::Mesons(), particles);
  }

  inline void KinematicsProcessor::SetBaryonParticles(const ROOT::RVec<std::string>& particles) {
      SetGroup(rad::consts::Baryons(), particles);
  }

  inline void KinematicsProcessor::RegisterCalc(const std::string& name, KineCalculation::MapKernel func) {
      _calculations.emplace_back(name, func);
  }

  inline void KinematicsProcessor::RegisterCalc(const std::string& name, KineCalculation::IndexKernel func, ROOT::RVec<ParticleNames_t> particles) {
    
      _calculations.emplace_back(name, func, particles);
  }

  inline void KinematicsProcessor::DefineAux() {
      auto define_pack = [&](const std::string& name, const ROOT::RVec<std::string>& cols, bool is_int) {
          if(cols.empty()) {
              if(is_int) _reaction->Define(name, [](){ return RVecRVecI{}; }, {});
              else       _reaction->Define(name, [](){ return RVecRVecD{}; }, {});
              return;
          }
          _reaction->Define(name, rad::util::createPackVectorString(cols));
      };
      define_pack(_prefix + "aux_pre_d" + _suffix + DoNotWriteTag(), _preModifier.GetAuxDoubleCols(), false);
      define_pack(_prefix + "aux_pre_i" + _suffix + DoNotWriteTag(), _preModifier.GetAuxIntCols(), true);
      define_pack(_prefix + "aux_post_d" + _suffix + DoNotWriteTag(), _postModifier.GetAuxDoubleCols(), false);
      define_pack(_prefix + "aux_post_i" + _suffix + DoNotWriteTag(), _postModifier.GetAuxIntCols(), true);
  }
  
  inline ParticleCreator& KinematicsProcessor::Creator() { return _creator; }
  inline const ParticleCreator& KinematicsProcessor::Creator() const { return _creator; }
  
  inline ParticleModifier& KinematicsProcessor::PreModifier() { return _preModifier; }
  inline ParticleModifier& KinematicsProcessor::PostModifier() { return _postModifier; }

  inline ConfigReaction* KinematicsProcessor::Reaction() const { return _reaction; }
  inline std::string KinematicsProcessor::GetSuffix() const { return _suffix; }
  inline std::string KinematicsProcessor::GetPrefix() const { return _prefix; }
  
  inline std::string KinematicsProcessor::FullName(const std::string& baseName) const { 
      return _prefix + baseName + _suffix; 
  }

  // --- Core Operator ---
  template<typename Tp, typename Tm> 
  inline KinematicsProcessor::CombiOutputVec_t KinematicsProcessor::operator()(
        const RVecIndices& indices, const Tp& px, const Tp& py, const Tp& pz, const Tm& m,
        const RVecRVecD& aux_pre_d, const RVecRVecI& aux_pre_i,
        const RVecRVecD& aux_post_d, const RVecRVecI& aux_post_i) const 
  {
    const auto Ncomponents = 4; // x, y, z, m
    const auto Nparticles0 = indices.size(); // Number of input particles
    const auto Nparticles = Nparticles0 + _creator.GetNCreated(); 
    
    if (Nparticles == 0) return CombiOutputVec_t(Ncomponents); 
          
    const auto Ncombis = indices[0].size(); 
    CombiOutputVec_t result(Ncombis, RVec<RVecResultType>(Ncomponents, RVecResultType(Nparticles)));

    ROOT::RVecD temp_px(Nparticles, consts::InvalidEntry<double>());
    ROOT::RVecD temp_py(Nparticles, consts::InvalidEntry<double>());
    ROOT::RVecD temp_pz(Nparticles, consts::InvalidEntry<double>());
    ROOT::RVecD temp_m(Nparticles, consts::InvalidEntry<double>());

    AuxCacheD cache_pre_d(aux_pre_d.size(), ROOT::RVecD(Nparticles));
    AuxCacheI cache_pre_i(aux_pre_i.size(), ROOT::RVecI(Nparticles));
    AuxCacheD cache_post_d(aux_post_d.size(), ROOT::RVecD(Nparticles));
    AuxCacheI cache_post_i(aux_post_i.size(), ROOT::RVecI(Nparticles));
    
    for (size_t icombi = 0; icombi < Ncombis; ++icombi) {
      
      for (size_t ip = 0; ip < Nparticles0; ++ip) {
        size_t iparti = _creator.GetReactionIndex(ip);                
        const int original_index = indices[ip][icombi];     

        temp_px[iparti] = px[original_index];
        temp_py[iparti] = py[original_index];
        temp_pz[iparti] = pz[original_index];
        temp_m[iparti]  = m[original_index];

        for(size_t v=0; v<aux_pre_d.size(); ++v) cache_pre_d[v][iparti] = aux_pre_d[v][original_index];
        for(size_t v=0; v<aux_pre_i.size(); ++v) cache_pre_i[v][iparti] = aux_pre_i[v][original_index];
      }

      _preModifier.Apply(temp_px, temp_py, temp_pz, temp_m, cache_pre_d, cache_pre_i);
      
      _creator.ApplyCreation(temp_px, temp_py, temp_pz, temp_m);
      
      _postModifier.Apply(temp_px, temp_py, temp_pz, temp_m, cache_post_d, cache_post_i);
 
      result[icombi][OrderX()] = temp_px;
      result[icombi][OrderY()] = temp_py;
      result[icombi][OrderZ()] = temp_pz;
      result[icombi][OrderM()] = temp_m;
    }

    return result;
  }

  // --- Snapshot Support ---
  inline void KinematicsProcessor::DefineNewComponentVecs() {

    // 1. Get particle names
       auto particle_names = _creator.GetParticleNames();
       
       // 2. Define Components mapping
       const ROOT::RVec<std::pair<std::string, int>> components = {
         {"_px", 0}, {"_py", 1}, {"_pz", 2}, {"_m", 3}
       };
       
       std::string resultColName = _prefix + consts::KineComponents() + _suffix;

       for (const auto& pName : particle_names) {
           size_t idx = _creator.GetReactionIndex(pName);
           std::string full_pName = FullName(pName);

           for (const auto& comp : components) {
               std::string suffix = comp.first;
               auto compIdx = comp.second;
        
               // IMPORTANT: Register these so AnalysisManager::CollectStreamColumns 
               // knows they exist and adds them to the Snapshot list.
               _registered_vars.push_back(pName+suffix); 

               // Direct Component Copy
               _reaction->Define(full_pName + suffix, 
                 [idx,compIdx](const CombiOutputVec_t& res) {
                   ROOT::RVec<double> out(res.size());
                   for(size_t i=0; i<res.size(); ++i) {
                     // res[i] = List of Components (RVec<RVec<double>>)
                     // res[i][compIdx] = List of Particles (RVec<double>)
                     // res[i][compIdx][idx] = Value (double)
                     out[i] = res[i][compIdx][idx];
                   }
                   return out;
                 }, 
                 {resultColName}
                 );
           }
       }
  }
  
  // --- Definitions ---

  inline void KinematicsProcessor::Define(const std::string& name, const std::string& func) {
      std::string colName = FullName(name); 
      _reaction->Define(colName, util::createFunctionCallStringFromVec("rad::util::ApplyKinematics", 
                                       {func, Creator().GetMapName(), (_prefix + consts::KineComponents() + _suffix)}));
      _registered_vars.push_back(name);
  }
  
  template <typename Lambda>
  inline void KinematicsProcessor::DefineKernel(const std::string& name, Lambda&& func) {
    
      ROOT::RDF::ColumnNames_t cols = { Creator().GetMapName(), _prefix + consts::KineComponents() + _suffix };
      auto apply_func = [func](const RVecIndexMap& map, const ROOT::RVec<ROOT::RVec<RVecResultType>>& comps){
        return rad::util::ApplyKinematics(func, map, comps);
      };
      _reaction->Define(FullName(name), apply_func, cols);
      
      // calculated variables should be registered!
      _registered_vars.push_back(name);
  }

 
  // --- Shortcuts ---
  inline void KinematicsProcessor::Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    RegisterCalc(name, rad::FourVectorMassCalc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
  }
  inline void KinematicsProcessor::Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    RegisterCalc(name, rad::FourVectorMass2Calc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
  }
  inline void KinematicsProcessor::Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    RegisterCalc(name, rad::FourVectorPtCalc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
  }
  inline void KinematicsProcessor::ParticleTheta(const ParticleNames_t& particles) {
    //here we actually perform a loop over combies
    //so we need to call a dunction which returns
    //a single entry each time.
    //rad::ThreeVectorTheta has an overwite which only uses the
    //first entry of the first entry in the given RVecIndices list
    for(const auto& p: particles){
      RegisterCalc(p+"_theta", rad::ThreeVectorTheta, {{p}});
    }
  }
 inline void KinematicsProcessor::ParticlePhi(const ParticleNames_t& particles) {
    for(const auto& p: particles){
      RegisterCalc(p+"_phi", rad::ThreeVectorPhi, {{p}});
    }
  }
 inline void KinematicsProcessor::ParticleP(const ParticleNames_t& particles) {
    for(const auto& p: particles){
      RegisterCalc(p+"_pmag", rad::ThreeVectorMag, {{p}});
    }
  }
  inline void KinematicsProcessor::ParticleEta(const ParticleNames_t& particles) {
    for(const auto& p: particles){
      RegisterCalc(p+"_eta", rad::ThreeVectorPhi, {{p}});
    }
  }

  inline void KinematicsProcessor::PrintReactionMap() const {
    std::cout << "\n=== KinematicsProcessor [" << _prefix << "] " << _suffix << " Reaction Map ===" << std::endl;
    std::cout << std::left << std::setw(20) << "Particle Name" << "Index" << std::endl;
    std::cout << std::string(30, '-') << std::endl;
    for(const auto& [name, idx] : _creator.GetIndexMap()) {
      std::cout << std::left << std::setw(20) << name << idx << std::endl;
    }
    std::cout << "========================================\n" << std::endl;
  }

  // =================================================================================
  // IMPLEMENTATION: KineCalculation::Define
  // =================================================================================
  
  inline void KineCalculation::Define(KinematicsProcessor* processor) {
    if (_kern_type == KernelType::Map) {
          processor->DefineKernel(_name, _mapFunc);
      } 
      else if (_kern_type == KernelType::Index) {
          RVecIndices resolved_indices;
          for(const auto& group : _particles) {
              Indices_t idxs;
              for(const auto& pname : group) {
                  idxs.push_back(processor->Creator().GetReactionIndex(pname));
              }
              resolved_indices.push_back(idxs);
          }

          auto func_ptr = _indexFunc; 
          auto adapter = [resolved_indices, func_ptr](const RVecIndexMap&, 
                                                      const RVecResultType& px, const RVecResultType& py, 
                                                      const RVecResultType& pz, const RVecResultType& m) 
          {
              return func_ptr(resolved_indices, px, py, pz, m);
          };

          processor->DefineKernel(_name, adapter);
      }
  }

} // namespace rad
#include "Diagnostics_KinematicProcessor.hxx"
