#pragma once

#include "ConfigReaction.h"
#include "StringUtilities.h"
#include "RVecHelpers.h"
#include "ParticleCreator.h"
#include "KinematicsDispatch.h"
#include "KineCalculation.h"
#include "BasicKinematics.h"
#include "ParticleModifier.h" 

namespace rad {

  using ROOT::RVec;
  
  // Aliases for Packed Auxiliary Columns passed from RDataFrame
  using RVecRVecD = ROOT::RVec<ROOT::RVecD>;
  using RVecRVecI = ROOT::RVec<ROOT::RVecI>;

  /**
   * @class KinematicsProcessor
   * @brief The Computational Engine for combinatorial kinematic analysis.
   * * @details
   * This class acts as the bridge between the physics definition (Topology) and the 
   * computational execution (RDataFrame). It manages the "Event Loop" for a specific 
   * data tier (e.g., "rec_" or "tru_").
   * * **Key Architectural Patterns:**
   * - **Twin Topology (Master/Linked):** Supports creating a "Master" processor that runs the 
   * combinatorial loop, and a "Linked" processor that reuses the Master's data for 
   * alternative hypothesis testing (e.g., Missing Neutron vs. Detected Neutron) at zero extra cost.
   * - **Structure of Arrays (SoA):** Transposes input data into fixed-index arrays for 
   * efficient, cache-friendly SIMD processing.
   */
  class KinematicsProcessor {

  public:
    /**
     * @brief Output Data Structure.
     * Dimensions: `[Combination_Index][Component_Index][Particle_Index]`
     * - Component_Index 0: Px
     * - Component_Index 1: Py
     * - Component_Index 2: Pz
     * - Component_Index 3: Mass
     */
    using CombiOutputVec_t = RVec<RVec<RVecResultType>>;

    // =================================================================================
    // Lifecycle & Initialization
    // =================================================================================

    /**
     * @brief Master Constructor.
     * @param cr Pointer to the ConfigReaction interface.
     * @param prefix The data tier prefix (e.g., "rec_", "tru_"). 
     * @param suffix Optional suffix for output columns (e.g., "_miss"). Use this to avoid name collisions.
     */
    KinematicsProcessor(ConfigReaction* cr, const std::string& prefix, const std::string& suffix = "");
    
    /**
     * @brief Fork Constructor (Copy & Rename).
     * Creates a deep copy of the configuration but assigns a new suffix. 
     * Useful for systematic variations (e.g., "rec_" -> "rec_sys_miss").
     */
    KinematicsProcessor(const KinematicsProcessor& other, const std::string& new_suffix);
    
    virtual ~KinematicsProcessor() = default;

    /**
     * @brief Fully initializes the processor in "Master" mode.
     * * @details
     * 1. **Topology:** Triggers `ParticleCreator::InitMap` to assign fixed indices to all particles.
     * 2. **Modifiers:** Initializes `PreModifier` and `PostModifier` (resolves names to indices).
     * 3. **Data:** Registers packed auxiliary columns (e.g., `rec_aux_cal_energy`).
     * 4. **Execution:** Registers the `DefineKinematicsProcessor` loop with RDataFrame.
     */
    void Init();
    
    /**
     * @brief Initializes the processor in "Linked" mode (Optimization).
     * * @details
     * Adopts the index structure from a Master processor to ensure topological alignment.
     * This allows the Linked processor to analyze the *same* combinations as the Master
     * but with different definitions (e.g., swapping a detected neutron for a missing one).
     * * @param master The fully initialized Master processor to adopt topology from.
     */
    void InitLinked(const KinematicsProcessor& master);

    // =================================================================================
    // Core Execution
    // =================================================================================

    /**
     * @brief The Core Event Loop (Functor).
     * * @details
     * Called by RDataFrame for every event. It performs the transformation:
     * `Raw Tracks -> Combinations -> Physics Candidates`.
     * * Pipeline:
     * 1. **Extract:** Maps random input tracks to the fixed Reaction Cache.
     * 2. **Pre-Mod:** Applies calibration/corrections to inputs.
     * 3. **Create:** Builds intermediate particles (Sum/Diff).
     * 4. **Post-Mod:** Applies constraints to created particles.
     */
    template<typename Tp, typename Tm> 
    CombiOutputVec_t operator()(const RVecIndices& indices, 
                                const Tp& px, const Tp& py, const Tp& pz, const Tm& m,
                                const RVecRVecD& aux_pre_d, const RVecRVecI& aux_pre_i,
                                const RVecRVecD& aux_post_d, const RVecRVecI& aux_post_i) const;

    /**
     * @brief Reserved hook for explicit component unpacking.
     * Required by KinematicsDispatch.
     */
    void DefineNewComponentVecs();

    // =================================================================================
    // Configuration Helpers
    // =================================================================================

    /** @brief Internal: Defines the packed auxiliary data columns for RDataFrame. */
    void DefineAux();

    // =================================================================================
    // Accessors
    // =================================================================================

    /** @brief Access the Topology Architect. */
    ParticleCreator& Creator();
    const ParticleCreator& Creator() const;
    
    /** @brief Access the Input Modifier (Calibration) Manager. */
    ParticleModifier& PreModifier();

    /** @brief Access the Output Modifier (Constraint) Manager. */
    ParticleModifier& PostModifier();

    ConfigReaction* Reaction() const;
    std::string GetSuffix() const;
    std::string GetPrefix() const;
    
    /** @brief Returns `{prefix}{baseName}{suffix}` (e.g., "rec_Q2_miss"). */
    std::string FullName(const std::string& baseName) const;

    // =================================================================================
    // Definition API (User-Facing)
    // =================================================================================
    
    /**
     * @brief Registers a kinematic calculation using a JIT-string Kernel.
     * @param name Name of the output column (e.g., "Q2").
     * @param func Name of the global C++ function to call (e.g., "rad::CalculateQ2").
     */
    void Define(const std::string& name, const std::string& func);

    /**
     * @brief Registers a kinematic calculation using a JIT-string Kernel with explicit particles.
     * @param name Name of the output column.
     * @param func Name of the global C++ function.
     * @param particles List of particle names to pass to the function.
     */
    void Define(const std::string& name, const std::string& func, const std::vector<ParticleNames_t>& particles);

    /**
     * @brief Registers a C++ Lambda/Functor as a kinematic kernel.
     * Preferred for performance and compile-time safety.
     */
    template <typename Lambda>
    void DefineKernel(const std::string& name, Lambda&& func);

    // =================================================================================
    // Physics Shortcuts
    // =================================================================================

    /** @brief Calculates Invariant Mass: `Sqrt( (Sum P_pos - Sum P_neg)^2 )`. */
    void Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

    /** @brief Calculates Missing Mass Squared: `(Sum P_pos - Sum P_neg)^2`. */
    void Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});
    
    /** @brief Calculates Transverse Momentum (Pt). */
    void Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

    /** @brief Prints the Particle->Index mapping to stdout. */
    void PrintReactionMap() const;
 
  private:
    ConfigReaction* _reaction = nullptr;
    std::string _prefix; 
    std::string _suffix; 
    
    ParticleCreator _creator;
    ParticleModifier _preModifier;
    ParticleModifier _postModifier;
  };

  // =================================================================================
  // IMPLEMENTATION
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
      _postModifier(other._postModifier)
  {}

  inline void KinematicsProcessor::Init() {
    // 1. Build Map (assign fixed indices) for THIS prefix
    // No loops over types here. We only care about our own prefix.
    Creator().InitMap(); 
    
    // 2. Initialize Modifiers
    _preModifier.Init(Creator());
    _postModifier.Init(Creator());

    // 3. Define Aux Columns (e.g. rec_aux_pre_d)
    DefineAux();

    // 4. Register Processor execution loop for THIS prefix
    DefineKinematicsProcessor(*_reaction, *this, _prefix);
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

      // Define packs using _prefix
      define_pack(_prefix + "aux_pre_d" + _suffix + DoNotWriteTag(), _preModifier.GetAuxDoubleCols(), false);
      define_pack(_prefix + "aux_pre_i" + _suffix + DoNotWriteTag(), _preModifier.GetAuxIntCols(), true);
      
      define_pack(_prefix + "aux_post_d" + _suffix + DoNotWriteTag(), _postModifier.GetAuxDoubleCols(), false);
      define_pack(_prefix + "aux_post_i" + _suffix + DoNotWriteTag(), _postModifier.GetAuxIntCols(), true);
  }

  inline void KinematicsProcessor::InitLinked(const KinematicsProcessor& master) {
    if(_suffix == master.GetSuffix() && _prefix == master.GetPrefix()) 
        throw std::runtime_error("InitLinked Error: Linked processor must have different signature than Master.");

    // 1. Adopt Indices: Enforce 1:1 Topology with Master
    _creator.AdoptIndices(master.Creator());
    
    // 2. Rebuild Map: Build the specific map for THIS prefix
    _creator.RebuildReactionMap();

    // 3. Initialize Modifiers 
    _preModifier.Init(Creator());
    _postModifier.Init(Creator());

    // 4. Define Aux and Execution Loop
    DefineAux();
    DefineKinematicsProcessor(*_reaction, *this, _prefix);
  }
  
  inline ParticleCreator& KinematicsProcessor::Creator() { return _creator; }
  inline const ParticleCreator& KinematicsProcessor::Creator() const { return _creator; }
  
  inline ParticleModifier& KinematicsProcessor::PreModifier() { return _preModifier; }
  inline ParticleModifier& KinematicsProcessor::PostModifier() { return _postModifier; }

  inline ConfigReaction* KinematicsProcessor::Reaction() const { return _reaction; }
  inline std::string KinematicsProcessor::GetSuffix() const { return _suffix; }
  inline std::string KinematicsProcessor::GetPrefix() const { return _prefix; }
  
  inline std::string KinematicsProcessor::FullName(const std::string& baseName) const { 
      // e.g. "Q2" -> "rec_Q2_miss"
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
    const auto Ncomponents = 4; 
    const auto Nparticles0 = indices.size(); 
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

      //Apply corrections to particle momenta, before creating other (non-detected) particles
      _preModifier.Apply(temp_px, temp_py, temp_pz, temp_m, cache_pre_d, cache_pre_i);

      //Create other particles
      _creator.ApplyCreation(temp_px, temp_py, temp_pz, temp_m);

      //Apply corrections to particle momenta after creation of other particles
      _postModifier.Apply(temp_px, temp_py, temp_pz, temp_m, cache_post_d, cache_post_i);
  
      result[icombi][OrderX()] = temp_px;
      result[icombi][OrderY()] = temp_py;
      result[icombi][OrderZ()] = temp_pz;
      result[icombi][OrderM()] = temp_m;
    }

    return result;
  }

  // =================================================================================
  // Definitions & Shortcuts
  // =================================================================================

  inline void KinematicsProcessor::Define(const std::string& name, const std::string& func) {
      // Define specific column: "rec_Q2" using "rec_ReactionMap" and "rec_KineComponents"
      std::string colName = FullName(name); // adds prefix and suffix
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
      // Prepend Kine tag and append suffix to look up index columns
      // Note: We use the logical names here, ApplyKinematics expects index aliases
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
    KineCalculation calc(*this, name, rad::FourVectorMassCalc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
  }
  inline void KinematicsProcessor::Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    KineCalculation calc(*this, name, rad::FourVectorMass2Calc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
  }
  inline void KinematicsProcessor::Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    KineCalculation calc(*this, name, rad::FourVectorPtCalc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
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
} // namespace rad

