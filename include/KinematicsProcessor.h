/**
 * @file KinematicsProcessor.h
 * @brief The Computational Engine for combinatorial kinematic analysis.
 */

#pragma once

#include "ConfigReaction.h"
#include "StringUtilities.h"
#include "RVecHelpers.h"
#include "ParticleCreator.h"
#include "KinematicsDispatch.h"
#include "KineCalculation.h" 
#include "BasicKinematics.h"
#include "ParticleModifier.h" 

#include <vector>
#include <string>
#include <memory>

namespace rad {

  using ROOT::RVec;
  using RVecRVecD = ROOT::RVec<ROOT::RVecD>;
  using RVecRVecI = ROOT::RVec<ROOT::RVecI>;

  /**
   * @class KinematicsProcessor
   * @brief Manages the "Event Loop" for kinematic calculations and RDF graph construction.
   * @details 
   * This class coordinates the setup and execution of kinematic analysis. It acts as the 
   * central manager for:
   * 1. **Topology Definition:** Uses `ParticleCreator` to define inputs and composite particles.
   * 2. **Calculation Management:** Registers and executes physics kernels (e.g., Mass, Missing Mass).
   * 3. **Multi-Stream Handling:** Manages "Clones" to apply the same analysis logic to different 
   * data streams (Rec vs. Truth) or hypotheses (Detected vs. Missing) automatically.
   * 4. **Execution:** Provides the functor `operator()` used by `RDataFrame::Define` to process events.
   */
  class KinematicsProcessor {

  public:
    using CombiOutputVec_t = RVec<RVec<RVecResultType>>;

    // =================================================================================
    // Lifecycle & Initialization
    // =================================================================================

    /**
     * @brief Constructor for the Master processor.
     * @param cr Pointer to the ConfigReaction object managing the RDF nodes.
     * @param prefix Data type prefix (e.g., "rec_" for reconstructed data).
     * @param suffix Optional suffix to distinguish this processor's outputs (e.g., "_miss").
     */
    KinematicsProcessor(ConfigReaction* cr, const std::string& prefix, const std::string& suffix = "");
    
    /** * @brief Copy Constructor used internally for cloning. 
     * @param other The source processor to copy configuration from.
     * @param new_suffix The new suffix for the cloned instance.
     */
    KinematicsProcessor(const KinematicsProcessor& other, const std::string& new_suffix);
    
    virtual ~KinematicsProcessor() = default;

    /** * @brief Creates a Clone for a different data type (e.g., Rec -> Truth).
     * @details Copies the configuration but swaps the prefix (e.g., "rec_" -> "mc_").
     * The clone is stored internally so that calling `Init()` on the master automatically 
     * initializes this clone.
     * @param new_prefix The prefix for the new stream (e.g., "mc_").
     * @param copyModifiers Whether to copy smearers/modifiers (default false).
     * @return Shared pointer to the new processor.
     */
    std::shared_ptr<KinematicsProcessor> CloneForType(const std::string& new_prefix, bool copyModifiers = false);

    /** * @brief Creates a Clone for a different physics hypothesis on the same data type.
     * @details Uses a unique suffix (e.g., "_miss") to create parallel branches in the RDF graph.
     * Useful for comparing "Detected Recoil" vs "Missing Recoil" topologies side-by-side.
     * @param suffix Unique suffix for the new stream.
     * @return Shared pointer to the new processor.
     */
    std::shared_ptr<KinematicsProcessor> CloneLinked(const std::string& suffix);

    /** * @brief Initializes this processor and cascades initialization to all registered clones.
     * @details 
     * 1. Initializes self (Creator, Modifiers, RDF columns).
     * 2. Calls `Init()` on all Type Clones (e.g., Truth).
     * 3. Calls `InitLinked()` on all Linked Clones (e.g., Missing Recoil).
     */
    void Init();
    
    /** * @brief Initializes a Linked Clone by adopting the topology from its Master.
     * @details This is called automatically by the Master's `Init()`. It copies the index maps
     * and group definitions from the Master to ensure consistency.
     * @param master The master processor to link against.
     */
    void InitLinked(const KinematicsProcessor& master);

    // =================================================================================
    // Core Execution (Functor)
    // =================================================================================

    /** * @brief The main execution functor called by RDataFrame for every event.
     * @details 
     * Iterates over all particle combinations:
     * 1. Loads input 4-vectors.
     * 2. Applies Pre-Modifiers (e.g., Energy Loss corrections).
     * 3. Creates intermediate particles (Sum/Diff).
     * 4. Applies Post-Modifiers (e.g., Smearing).
     * 5. Fills result vectors.
     * * @param indices Combinatorial indices for input particles.
     * @param px Input Px vector (RVec).
     * @param py Input Py vector (RVec).
     * @param pz Input Pz vector (RVec).
     * @param m  Input Mass vector (RVec).
     * @param aux_pre_d Auxiliary double columns for Pre-Modifiers.
     * @param aux_pre_i Auxiliary int columns for Pre-Modifiers.
     * @param aux_post_d Auxiliary double columns for Post-Modifiers.
     * @param aux_post_i Auxiliary int columns for Post-Modifiers.
     * @return A vector of 4-vectors (Px, Py, Pz, M) for all combinations.
     */
    template<typename Tp, typename Tm> 
    CombiOutputVec_t operator()(const RVecIndices& indices, 
                                const Tp& px, const Tp& py, const Tp& pz, const Tm& m,
                                const RVecRVecD& aux_pre_d, const RVecRVecI& aux_pre_i,
                                const RVecRVecD& aux_post_d, const RVecRVecI& aux_post_i) const;

    /** @brief Hook for defining custom component vectors if needed by derived classes. */
    void DefineNewComponentVecs();

    // =================================================================================
    // Configuration Helpers
    // =================================================================================

    /** @brief Defines auxiliary columns in RDF required for Modifiers (smearing, etc.). */
    void DefineAux();

    /** * @brief Overrides a group definition LOCALLY for this processor. 
     * @param groupName The abstract group name (e.g., "Baryons").
     * @param particles List of particle names to assign to this group.
     */
    void SetGroup(const std::string& groupName, const std::vector<std::string>& particles);

    /** @brief Helper to define the Meson system (for Mandelstam t/u/s and Angles). */
    void SetMesonParticles(const std::vector<std::string>& particles);

    /** @brief Helper to define the Baryon system (for Mandelstam t/u/s and Angles). */
    void SetBaryonParticles(const std::vector<std::string>& particles);

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
    
    /** @brief Helper to construct fully qualified column names (prefix + name + suffix). */
    std::string FullName(const std::string& baseName) const;

    // =================================================================================
    // Definition API & Calculation Registration
    // =================================================================================
    
    /** * @brief Register a calculation kernel using the standard Map signature.
     * @param name The output branch name (e.g., "Q2").
     * @param func The function to execute.
     */
    void RegisterCalc(const std::string& name, KineCalculation::MapKernel func);
    
    /** * @brief Register a calculation kernel using a list of specific particle names.
     * @param name The output branch name.
     * @param func The function to execute.
     * @param particles List of particle names required by the function.
     */
    void RegisterCalc(const std::string& name, KineCalculation::IndexKernel func, std::vector<ParticleNames_t> particles);

    // --- Legacy String-Based Definition ---
    
    void Define(const std::string& name, const std::string& func);
    void Define(const std::string& name, const std::string& func, const std::vector<ParticleNames_t>& particles);

    template <typename Lambda>
    void DefineKernel(const std::string& name, Lambda&& func);

    // =================================================================================
    // Physics Shortcuts
    // =================================================================================

    void Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});
    void Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});
    void Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});
    
     
    /** @brief Prints the resolved reaction map indices to stdout for debugging. */
    void PrintReactionMap() const;
 
  private:
    ConfigReaction* _reaction = nullptr;
    std::string _prefix; 
    std::string _suffix; 
    
    // Prevents double-initialization when cascading
    bool _isInitialized = false; 
    
    ParticleCreator _creator;
    ParticleModifier _preModifier;
    ParticleModifier _postModifier;

    std::vector<KineCalculation> _calculations;
    
    // Store clones to allow cascading initialization
    std::vector<std::shared_ptr<KinematicsProcessor>> _type_clones;
    std::vector<std::shared_ptr<KinematicsProcessor>> _linked_clones;

    struct GroupOverride {
        std::string name;
        std::vector<std::string> particles;
    };
    std::vector<GroupOverride> _groupOverrides;

    void ApplyGroupOverrides();
  };

  // =================================================================================
  // IMPLEMENTATION: KinematicsProcessor
  // =================================================================================

  inline KinematicsProcessor::KinematicsProcessor(ConfigReaction* cr, const std::string& prefix, const std::string& suffix) 
    : _reaction{cr}, _prefix{prefix}, _suffix{suffix}, _creator{cr, prefix, suffix} 
  {}

  inline KinematicsProcessor::KinematicsProcessor(const KinematicsProcessor& other, const std::string& new_suffix)
    : _reaction(other._reaction),
      _prefix(other._prefix),
      _suffix(new_suffix),
      _creator(other._creator, new_suffix),
      _preModifier(other._preModifier),   
      _postModifier(other._postModifier),
      _calculations(other._calculations), 
      _groupOverrides(other._groupOverrides) 
  {}

  inline std::shared_ptr<KinematicsProcessor> KinematicsProcessor::CloneForType(const std::string& new_prefix, bool copyModifiers) {
      auto clone = std::make_shared<KinematicsProcessor>(_reaction, new_prefix, _suffix);
      
      // 1. Copy Configuration
      clone->Creator().CopyConfigurationFrom(_creator); 
      
      // 2. Port Groups (Critical): Converts "rec_scat_ele_group" -> "mc_scat_ele_group"
      clone->Creator().PortGroupsFrom(_creator); 

      clone->_calculations = _calculations;
      if(copyModifiers) {
          clone->_preModifier = _preModifier;
          clone->_postModifier = _postModifier;
      }
      clone->_groupOverrides = _groupOverrides;
      
      // 3. Store in list for auto-initialization
      _type_clones.push_back(clone); 
      return clone;
  }

  inline std::shared_ptr<KinematicsProcessor> KinematicsProcessor::CloneLinked(const std::string& suffix) {
      auto clone = std::make_shared<KinematicsProcessor>(*this, suffix);
      
      // Store in separate list. Linked clones need InitLinked(), not just Init().
      _linked_clones.push_back(clone); 
      return clone;
  }

  inline void KinematicsProcessor::Init() {
    // Idempotency Check: Prevent double-init if called manually or via cascade
    if (_isInitialized) return; 
    _isInitialized = true;

    // 1. Initialize Self
    ApplyGroupOverrides();
    Creator().InitMap(); 
    _preModifier.Init(Creator());
    _postModifier.Init(Creator());
    DefineAux();
    DefineKinematicsProcessor(*_reaction, *this, _prefix);

    // Register all calculations into RDF
    for(auto& calc : _calculations) {
        calc.Define(this); 
    }

    // 2. Cascade: Initialize Type Clones (e.g., Truth)
    // These are independent processors with just a different prefix
    for(auto& clone : _type_clones) {
        clone->Init();
    }

    // 3. Cascade: Initialize Linked Clones (e.g., Missing Topology)
    // These depend on *this* processor's map, so we call InitLinked
    for(auto& clone : _linked_clones) {
        clone->InitLinked(*this);
    }
  }

  inline void KinematicsProcessor::InitLinked(const KinematicsProcessor& master) {
    // Idempotency Check
    if (_isInitialized) return; 
    
    // Safety: Linked processors must distinguish themselves via Suffix or Prefix
    if(_suffix == master.GetSuffix() && _prefix == master.GetPrefix()) 
        throw std::runtime_error("InitLinked Error: Linked processor must have different signature than Master.");
    
    _isInitialized = true;

    // 1. Adopt Topology from Master
    ApplyGroupOverrides(); 
    _creator.AdoptIndices(master.Creator());
    _creator.RebuildReactionMap(); // Rebuild map for the new suffix
    
    _preModifier.Init(Creator());
    _postModifier.Init(Creator());
    DefineAux();
    DefineKinematicsProcessor(*_reaction, *this, _prefix);
    
    for(auto& calc : _calculations) {
        calc.Define(this); 
    }

    // 2. Cascade: Initialize Type Clones OF THIS Linked Clone
    // e.g. If this is "Rec Missing", initialize "Truth Missing"
    for(auto& clone : _type_clones) {
        clone->Init();
    }
  }

  inline void KinematicsProcessor::ApplyGroupOverrides() {
      for(const auto& group : _groupOverrides) {
          Creator().OverrideGroup(group.name, group.particles); 
      }
  }

  inline void KinematicsProcessor::SetGroup(const std::string& groupName, const std::vector<std::string>& particles) {
      _groupOverrides.push_back({groupName, particles});
  }

  inline void KinematicsProcessor::SetMesonParticles(const std::vector<std::string>& particles) {
      SetGroup(rad::consts::Mesons(), particles);
  }

  inline void KinematicsProcessor::SetBaryonParticles(const std::vector<std::string>& particles) {
      SetGroup(rad::consts::Baryons(), particles);
  }

  inline void KinematicsProcessor::RegisterCalc(const std::string& name, KineCalculation::MapKernel func) {
      _calculations.emplace_back(name, func);
  }

  inline void KinematicsProcessor::RegisterCalc(const std::string& name, KineCalculation::IndexKernel func, std::vector<ParticleNames_t> particles) {
      _calculations.emplace_back(name, func, particles);
  }

  inline void KinematicsProcessor::DefineAux() {
      auto define_pack = [&](const std::string& name, const std::vector<std::string>& cols, bool is_int) {
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

  inline void KinematicsProcessor::DefineNewComponentVecs() {}

  // --- Core Operator ---
  template<typename Tp, typename Tm> 
  inline KinematicsProcessor::CombiOutputVec_t KinematicsProcessor::operator()(
        const RVecIndices& indices, const Tp& px, const Tp& py, const Tp& pz, const Tm& m,
        const RVecRVecD& aux_pre_d, const RVecRVecI& aux_pre_i,
        const RVecRVecD& aux_post_d, const RVecRVecI& aux_post_i) const 
  {
    const auto Ncomponents = 4; // x, y, z, m
    const auto Nparticles0 = indices.size(); // Number of input particles
    const auto Nparticles = Nparticles0 + _creator.GetNCreated(); // Total particles (Input + Created)
    
    // Safety check for empty events
    if (Nparticles == 0) return CombiOutputVec_t(Ncomponents); 
          
    const auto Ncombis = indices[0].size(); // Number of combinations in this event
    CombiOutputVec_t result(Ncombis, RVec<RVecResultType>(Ncomponents, RVecResultType(Nparticles)));

    // Temporary storage for 4-vectors of a single combination
    ROOT::RVecD temp_px(Nparticles, consts::InvalidEntry<double>());
    ROOT::RVecD temp_py(Nparticles, consts::InvalidEntry<double>());
    ROOT::RVecD temp_pz(Nparticles, consts::InvalidEntry<double>());
    ROOT::RVecD temp_m(Nparticles, consts::InvalidEntry<double>());

    // Cache buffers for Modifiers
    AuxCacheD cache_pre_d(aux_pre_d.size(), ROOT::RVecD(Nparticles));
    AuxCacheI cache_pre_i(aux_pre_i.size(), ROOT::RVecI(Nparticles));
    AuxCacheD cache_post_d(aux_post_d.size(), ROOT::RVecD(Nparticles));
    AuxCacheI cache_post_i(aux_post_i.size(), ROOT::RVecI(Nparticles));
    
    // --- Loop over Combinations ---
    for (size_t icombi = 0; icombi < Ncombis; ++icombi) {
      
      // 1. Load Input Particles
      for (size_t ip = 0; ip < Nparticles0; ++ip) {
        size_t iparti = _creator.GetReactionIndex(ip);                
        const int original_index = indices[ip][icombi];    

        temp_px[iparti] = px[original_index];
        temp_py[iparti] = py[original_index];
        temp_pz[iparti] = pz[original_index];
        temp_m[iparti]  = m[original_index];

        // Load Aux Data
        for(size_t v=0; v<aux_pre_d.size(); ++v) cache_pre_d[v][iparti] = aux_pre_d[v][original_index];
        for(size_t v=0; v<aux_pre_i.size(); ++v) cache_pre_i[v][iparti] = aux_pre_i[v][original_index];
      }

      // 2. Apply Pre-Modifiers (e.g. Energy Loss)
      _preModifier.Apply(temp_px, temp_py, temp_pz, temp_m, cache_pre_d, cache_pre_i);
      
      // 3. Create Intermediate Particles (e.g. P4 = P1 + P2)
      _creator.ApplyCreation(temp_px, temp_py, temp_pz, temp_m);
      
      // 4. Apply Post-Modifiers (e.g. Smearing)
      _postModifier.Apply(temp_px, temp_py, temp_pz, temp_m, cache_post_d, cache_post_i);
  
      // 5. Store Results
      result[icombi][OrderX()] = temp_px;
      result[icombi][OrderY()] = temp_py;
      result[icombi][OrderZ()] = temp_pz;
      result[icombi][OrderM()] = temp_m;
    }

    return result;
  }

  // --- Definitions ---

  inline void KinematicsProcessor::Define(const std::string& name, const std::string& func) {
      std::string colName = FullName(name); 
      _reaction->Define(colName, util::createFunctionCallStringFromVec("rad::util::ApplyKinematics", 
                                      {func, Creator().GetMapName(), (_prefix + consts::KineComponents() + _suffix)}));
  }
  
  template <typename Lambda>
  inline void KinematicsProcessor::DefineKernel(const std::string& name, Lambda&& func) {
      ROOT::RDF::ColumnNames_t cols = { Creator().GetMapName(), _prefix + consts::KineComponents() + _suffix };
      auto apply_func = [func](const RVecIndexMap& map, const ROOT::RVec<ROOT::RVec<RVecResultType>>& comps){
        return rad::util::ApplyKinematics(func, map, comps);
      };

      _reaction->Define(FullName(name), apply_func, cols);
  }

  inline void KinematicsProcessor::Define(const std::string& name, const std::string& func, const std::vector<ParticleNames_t>& particles) {
    std::string kine_parts = "{";
    for(const auto& pnames : particles){
      ParticleNames_t suffixed_names;
      for(const auto& p : pnames) suffixed_names.push_back(p + _suffix);
      kine_parts += util::combineVectorToString(util::prependToAll(suffixed_names, consts::data_type::Kine()));
      kine_parts += ",";
    }
    kine_parts.pop_back(); kine_parts += "}";
    
    _reaction->Define(FullName(name), util::createFunctionCallStringFromVec("rad::util::ApplyKinematics", 
                                      {func, kine_parts, (_prefix + consts::KineComponents() + _suffix)}));
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
      if (_type == Type::Map) {
          processor->DefineKernel(_name, _mapFunc);
      } 
      else if (_type == Type::Index) {
          // 1. Resolve Strings -> Indices
          RVecIndices resolved_indices;
          for(const auto& group : _particles) {
              Indices_t idxs;
              for(const auto& pname : group) {
                  idxs.push_back(processor->Creator().GetReactionIndex(pname));
              }
              resolved_indices.push_back(idxs);
          }

          // 2. Adapt to MapKernel Signature
          auto func_ptr = _indexFunc; 
          auto adapter = [resolved_indices, func_ptr](const RVecIndexMap&, 
                                                      const RVecResultType& px, const RVecResultType& py, 
                                                      const RVecResultType& pz, const RVecResultType& m) 
          {
              return func_ptr(resolved_indices, px, py, pz, m);
          };

          // 3. Register
          processor->DefineKernel(_name, adapter);
      }
  }

} // namespace rad

// namespace rad {

//   using ROOT::RVec;
  
//   // Aliases for Packed Auxiliary Columns passed from RDataFrame
//   using RVecRVecD = ROOT::RVec<ROOT::RVecD>;
//   using RVecRVecI = ROOT::RVec<ROOT::RVecI>;

//   /**
//    * @class KinematicsProcessor
//    * @brief The Computational Engine for combinatorial kinematic analysis.
//    * * @details
//    * This class acts as the bridge between the physics definition (Topology) and the 
//    * computational execution (RDataFrame). It manages the "Event Loop" for a specific 
//    * data tier (e.g., "rec_" or "tru_").
//    * * **Key Architectural Patterns:**
//    * - **Twin Topology (Master/Linked):** Supports creating a "Master" processor that runs the 
//    * combinatorial loop, and a "Linked" processor that reuses the Master's data for 
//    * alternative hypothesis testing (e.g., Missing Neutron vs. Detected Neutron) at zero extra cost.
//    * - **Structure of Arrays (SoA):** Transposes input data into fixed-index arrays for 
//    * efficient, cache-friendly SIMD processing.
//    */
//   class KinematicsProcessor {

//   public:
//     /**
//      * @brief Output Data Structure.
//      * Dimensions: `[Combination_Index][Component_Index][Particle_Index]`
//      * - Component_Index 0: Px
//      * - Component_Index 1: Py
//      * - Component_Index 2: Pz
//      * - Component_Index 3: Mass
//      */
//     using CombiOutputVec_t = RVec<RVec<RVecResultType>>;

//     // =================================================================================
//     // Lifecycle & Initialization
//     // =================================================================================

//     /**
//      * @brief Master Constructor.
//      * @param cr Pointer to the ConfigReaction interface.
//      * @param prefix The data tier prefix (e.g., "rec_", "tru_"). 
//      * @param suffix Optional suffix for output columns (e.g., "_miss"). Use this to avoid name collisions.
//      */
//     KinematicsProcessor(ConfigReaction* cr, const std::string& prefix, const std::string& suffix = "");
    
//     /**
//      * @brief Fork Constructor (Copy & Rename).
//      * Creates a deep copy of the configuration but assigns a new suffix. 
//      * Useful for systematic variations (e.g., "rec_" -> "rec_sys_miss").
//      */
//     KinematicsProcessor(const KinematicsProcessor& other, const std::string& new_suffix);
    
//     virtual ~KinematicsProcessor() = default;

//     /**
//      * @brief Fully initializes the processor in "Master" mode.
//      * * @details
//      * 1. **Topology:** Triggers `ParticleCreator::InitMap` to assign fixed indices to all particles.
//      * 2. **Modifiers:** Initializes `PreModifier` and `PostModifier` (resolves names to indices).
//      * 3. **Data:** Registers packed auxiliary columns (e.g., `rec_aux_cal_energy`).
//      * 4. **Execution:** Registers the `DefineKinematicsProcessor` loop with RDataFrame.
//      */
//     void Init();
    
//     /**
//      * @brief Initializes the processor in "Linked" mode (Optimization).
//      * * @details
//      * Adopts the index structure from a Master processor to ensure topological alignment.
//      * This allows the Linked processor to analyze the *same* combinations as the Master
//      * but with different definitions (e.g., swapping a detected neutron for a missing one).
//      * * @param master The fully initialized Master processor to adopt topology from.
//      */
//     void InitLinked(const KinematicsProcessor& master);

//     // =================================================================================
//     // Core Execution
//     // =================================================================================

//     /**
//      * @brief The Core Event Loop (Functor).
//      * * @details
//      * Called by RDataFrame for every event. It performs the transformation:
//      * `Raw Tracks -> Combinations -> Physics Candidates`.
//      * * Pipeline:
//      * 1. **Extract:** Maps random input tracks to the fixed Reaction Cache.
//      * 2. **Pre-Mod:** Applies calibration/corrections to inputs.
//      * 3. **Create:** Builds intermediate particles (Sum/Diff).
//      * 4. **Post-Mod:** Applies constraints to created particles.
//      */
//     template<typename Tp, typename Tm> 
//     CombiOutputVec_t operator()(const RVecIndices& indices, 
//                                 const Tp& px, const Tp& py, const Tp& pz, const Tm& m,
//                                 const RVecRVecD& aux_pre_d, const RVecRVecI& aux_pre_i,
//                                 const RVecRVecD& aux_post_d, const RVecRVecI& aux_post_i) const;

//     /**
//      * @brief Reserved hook for explicit component unpacking.
//      * Required by KinematicsDispatch.
//      */
//     void DefineNewComponentVecs();

//     // =================================================================================
//     // Configuration Helpers
//     // =================================================================================

//     /** @brief Internal: Defines the packed auxiliary data columns for RDataFrame. */
//     void DefineAux();

//     // =================================================================================
//     // Accessors
//     // =================================================================================

//     /** @brief Access the Topology Architect. */
//     ParticleCreator& Creator();
//     const ParticleCreator& Creator() const;
    
//     /** @brief Access the Input Modifier (Calibration) Manager. */
//     ParticleModifier& PreModifier();

//     /** @brief Access the Output Modifier (Constraint) Manager. */
//     ParticleModifier& PostModifier();

//     ConfigReaction* Reaction() const;
//     std::string GetSuffix() const;
//     std::string GetPrefix() const;
    
//     /** @brief Returns `{prefix}{baseName}{suffix}` (e.g., "rec_Q2_miss"). */
//     std::string FullName(const std::string& baseName) const;

//     // =================================================================================
//     // Definition API (User-Facing)
//     // =================================================================================
    
//     /**
//      * @brief Registers a kinematic calculation using a JIT-string Kernel.
//      * @param name Name of the output column (e.g., "Q2").
//      * @param func Name of the global C++ function to call (e.g., "rad::CalculateQ2").
//      */
//     void Define(const std::string& name, const std::string& func);

//     /**
//      * @brief Registers a kinematic calculation using a JIT-string Kernel with explicit particles.
//      * @param name Name of the output column.
//      * @param func Name of the global C++ function.
//      * @param particles List of particle names to pass to the function.
//      */
//     void Define(const std::string& name, const std::string& func, const std::vector<ParticleNames_t>& particles);

//     /**
//      * @brief Registers a C++ Lambda/Functor as a kinematic kernel.
//      * Preferred for performance and compile-time safety.
//      */
//     template <typename Lambda>
//     void DefineKernel(const std::string& name, Lambda&& func);

//     // =================================================================================
//     // Physics Shortcuts
//     // =================================================================================

//     /** @brief Calculates Invariant Mass: `Sqrt( (Sum P_pos - Sum P_neg)^2 )`. */
//     void Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

//     /** @brief Calculates Missing Mass Squared: `(Sum P_pos - Sum P_neg)^2`. */
//     void Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});
    
//     /** @brief Calculates Transverse Momentum (Pt). */
//     void Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

//     /** @brief Prints the Particle->Index mapping to stdout. */
//     void PrintReactionMap() const;
 
//   private:
//     ConfigReaction* _reaction = nullptr;
//     std::string _prefix; 
//     std::string _suffix; 
    
//     ParticleCreator _creator;
//     ParticleModifier _preModifier;
//     ParticleModifier _postModifier;
//   };

//   // =================================================================================
//   // IMPLEMENTATION
//   // =================================================================================

//   inline KinematicsProcessor::KinematicsProcessor(ConfigReaction* cr, const std::string& prefix, const std::string& suffix) 
//     : _reaction{cr}, _prefix{prefix}, _suffix{suffix}, _creator{cr, prefix, suffix} 
//   {}

//   inline KinematicsProcessor::KinematicsProcessor(const KinematicsProcessor& other, const std::string& new_suffix)
//     : _reaction(other._reaction),
//       _prefix(other._prefix),
//       _suffix(new_suffix),
//       _creator(other._creator, new_suffix),
//       _preModifier(other._preModifier),   
//       _postModifier(other._postModifier)
//   {}

//   inline void KinematicsProcessor::Init() {
//     // 1. Build Map (assign fixed indices) for THIS prefix
//     // No loops over types here. We only care about our own prefix.
//     Creator().InitMap(); 
    
//     // 2. Initialize Modifiers
//     _preModifier.Init(Creator());
//     _postModifier.Init(Creator());

//     // 3. Define Aux Columns (e.g. rec_aux_pre_d)
//     DefineAux();

//     // 4. Register Processor execution loop for THIS prefix
//     DefineKinematicsProcessor(*_reaction, *this, _prefix);
//   }

//   inline void KinematicsProcessor::DefineAux() {
//       auto define_pack = [&](const std::string& name, const std::vector<std::string>& cols, bool is_int) {
//           if(cols.empty()) {
//               if(is_int) _reaction->Define(name, [](){ return RVecRVecI{}; }, {});
//               else       _reaction->Define(name, [](){ return RVecRVecD{}; }, {});
//               return;
//           }
//           _reaction->Define(name, rad::util::createPackVectorString(cols));
//       };

//       // Define packs using _prefix
//       define_pack(_prefix + "aux_pre_d" + _suffix + DoNotWriteTag(), _preModifier.GetAuxDoubleCols(), false);
//       define_pack(_prefix + "aux_pre_i" + _suffix + DoNotWriteTag(), _preModifier.GetAuxIntCols(), true);
      
//       define_pack(_prefix + "aux_post_d" + _suffix + DoNotWriteTag(), _postModifier.GetAuxDoubleCols(), false);
//       define_pack(_prefix + "aux_post_i" + _suffix + DoNotWriteTag(), _postModifier.GetAuxIntCols(), true);
//   }

//   inline void KinematicsProcessor::InitLinked(const KinematicsProcessor& master) {
//     if(_suffix == master.GetSuffix() && _prefix == master.GetPrefix()) 
//         throw std::runtime_error("InitLinked Error: Linked processor must have different signature than Master.");

//     // 1. Adopt Indices: Enforce 1:1 Topology with Master
//     _creator.AdoptIndices(master.Creator());
    
//     // 2. Rebuild Map: Build the specific map for THIS prefix
//     _creator.RebuildReactionMap();

//     // 3. Initialize Modifiers 
//     _preModifier.Init(Creator());
//     _postModifier.Init(Creator());

//     // 4. Define Aux and Execution Loop
//     DefineAux();
//     DefineKinematicsProcessor(*_reaction, *this, _prefix);
//   }
  
//   inline ParticleCreator& KinematicsProcessor::Creator() { return _creator; }
//   inline const ParticleCreator& KinematicsProcessor::Creator() const { return _creator; }
  
//   inline ParticleModifier& KinematicsProcessor::PreModifier() { return _preModifier; }
//   inline ParticleModifier& KinematicsProcessor::PostModifier() { return _postModifier; }

//   inline ConfigReaction* KinematicsProcessor::Reaction() const { return _reaction; }
//   inline std::string KinematicsProcessor::GetSuffix() const { return _suffix; }
//   inline std::string KinematicsProcessor::GetPrefix() const { return _prefix; }
  
//   inline std::string KinematicsProcessor::FullName(const std::string& baseName) const { 
//       // e.g. "Q2" -> "rec_Q2_miss"
//       return _prefix + baseName + _suffix; 
//   }

//   inline void KinematicsProcessor::DefineNewComponentVecs() {}

//   // --- Core Operator ---
//   template<typename Tp, typename Tm> 
//   inline KinematicsProcessor::CombiOutputVec_t KinematicsProcessor::operator()(
//         const RVecIndices& indices, const Tp& px, const Tp& py, const Tp& pz, const Tm& m,
//         const RVecRVecD& aux_pre_d, const RVecRVecI& aux_pre_i,
//         const RVecRVecD& aux_post_d, const RVecRVecI& aux_post_i) const 
//   {
//     const auto Ncomponents = 4; 
//     const auto Nparticles0 = indices.size(); 
//     const auto Nparticles = Nparticles0 + _creator.GetNCreated(); 
//     if (Nparticles == 0) return CombiOutputVec_t(Ncomponents); 
          
//     const auto Ncombis = indices[0].size(); 
//     CombiOutputVec_t result(Ncombis, RVec<RVecResultType>(Ncomponents, RVecResultType(Nparticles)));

//     ROOT::RVecD temp_px(Nparticles, consts::InvalidEntry<double>());
//     ROOT::RVecD temp_py(Nparticles, consts::InvalidEntry<double>());
//     ROOT::RVecD temp_pz(Nparticles, consts::InvalidEntry<double>());
//     ROOT::RVecD temp_m(Nparticles, consts::InvalidEntry<double>());

//     AuxCacheD cache_pre_d(aux_pre_d.size(), ROOT::RVecD(Nparticles));
//     AuxCacheI cache_pre_i(aux_pre_i.size(), ROOT::RVecI(Nparticles));
//     AuxCacheD cache_post_d(aux_post_d.size(), ROOT::RVecD(Nparticles));
//     AuxCacheI cache_post_i(aux_post_i.size(), ROOT::RVecI(Nparticles));
    
//     for (size_t icombi = 0; icombi < Ncombis; ++icombi) {
      
//       for (size_t ip = 0; ip < Nparticles0; ++ip) {
// 	size_t iparti = _creator.GetReactionIndex(ip);                
// 	const int original_index = indices[ip][icombi];   

//         temp_px[iparti] = px[original_index];
//         temp_py[iparti] = py[original_index];
//         temp_pz[iparti] = pz[original_index];
//         temp_m[iparti]  = m[original_index];

//         for(size_t v=0; v<aux_pre_d.size(); ++v) cache_pre_d[v][iparti] = aux_pre_d[v][original_index];
//         for(size_t v=0; v<aux_pre_i.size(); ++v) cache_pre_i[v][iparti] = aux_pre_i[v][original_index];
//       }

//       //Apply corrections to particle momenta, before creating other (non-detected) particles
//       _preModifier.Apply(temp_px, temp_py, temp_pz, temp_m, cache_pre_d, cache_pre_i);

//       //Create other particles
//       _creator.ApplyCreation(temp_px, temp_py, temp_pz, temp_m);

//       //Apply corrections to particle momenta after creation of other particles
//       _postModifier.Apply(temp_px, temp_py, temp_pz, temp_m, cache_post_d, cache_post_i);
  
//       result[icombi][OrderX()] = temp_px;
//       result[icombi][OrderY()] = temp_py;
//       result[icombi][OrderZ()] = temp_pz;
//       result[icombi][OrderM()] = temp_m;
//     }

//     return result;
//   }

//   // =================================================================================
//   // Definitions & Shortcuts
//   // =================================================================================

//   inline void KinematicsProcessor::Define(const std::string& name, const std::string& func) {
//       // Define specific column: "rec_Q2" using "rec_ReactionMap" and "rec_KineComponents"
//       std::string colName = FullName(name); // adds prefix and suffix
//       _reaction->Define(colName, util::createFunctionCallStringFromVec("rad::util::ApplyKinematics", 
//                                       {func, Creator().GetMapName(), (_prefix + consts::KineComponents() + _suffix)}));
//   }
  
//   template <typename Lambda>
//   inline void KinematicsProcessor::DefineKernel(const std::string& name, Lambda&& func) {
//       ROOT::RDF::ColumnNames_t cols = { Creator().GetMapName(), _prefix + consts::KineComponents() + _suffix };
//       auto apply_func = [func](const RVecIndexMap& map, const ROOT::RVec<ROOT::RVec<RVecResultType>>& comps){
//         return rad::util::ApplyKinematics(func, map, comps);
//       };

//       _reaction->Define(FullName(name), apply_func, cols);
//   }

//   inline void KinematicsProcessor::Define(const std::string& name, const std::string& func, const std::vector<ParticleNames_t>& particles) {
//     std::string kine_parts = "{";
//     for(const auto& pnames : particles){
//       ParticleNames_t suffixed_names;
//       // Prepend Kine tag and append suffix to look up index columns
//       // Note: We use the logical names here, ApplyKinematics expects index aliases
//       for(const auto& p : pnames) suffixed_names.push_back(p + _suffix);
//       kine_parts += util::combineVectorToString(util::prependToAll(suffixed_names, consts::data_type::Kine()));
//       kine_parts += ",";
//     }
//     kine_parts.pop_back(); kine_parts += "}";
    
//     _reaction->Define(FullName(name), util::createFunctionCallStringFromVec("rad::util::ApplyKinematics", 
//                                       {func, kine_parts, (_prefix + consts::KineComponents() + _suffix)}));
//   }

//   // --- Shortcuts ---
//   inline void KinematicsProcessor::Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
//     KineCalculation calc(*this, name, rad::FourVectorMassCalc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
//   }
//   inline void KinematicsProcessor::Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
//     KineCalculation calc(*this, name, rad::FourVectorMass2Calc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
//   }
//   inline void KinematicsProcessor::Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
//     KineCalculation calc(*this, name, rad::FourVectorPtCalc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
//   }
  
//   inline void KinematicsProcessor::PrintReactionMap() const {
//     std::cout << "\n=== KinematicsProcessor [" << _prefix << "] " << _suffix << " Reaction Map ===" << std::endl;
//     std::cout << std::left << std::setw(20) << "Particle Name" << "Index" << std::endl;
//     std::cout << std::string(30, '-') << std::endl;
//     for(const auto& [name, idx] : _creator.GetIndexMap()) {
//       std::cout << std::left << std::setw(20) << name << idx << std::endl;
//     }
//     std::cout << "========================================\n" << std::endl;
//   }
// } // namespace rad

