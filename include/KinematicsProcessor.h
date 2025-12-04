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
  // These represent a vector of columns, where each column is an RVec
  using RVecRVecD = ROOT::RVec<ROOT::RVecD>;
  using RVecRVecI = ROOT::RVec<ROOT::RVecI>;

  /**
   * @brief Master functor for executing combinatorial kinematic analysis on a single event.
   * * This class acts as the **Computational Engine** of the framework. It is responsible for:
   * 1. **Orchestration**: Managing the `ParticleCreator` to define the static reaction map.
   * 2. **Execution**: Running the thread-safe event loop that transforms input indices into kinematic data.
   * 3. **Modification**: Applying Pre-creation (input scaling) and Post-creation (mass constraint) modifiers.
   * 4. **Output**: Returning transposed, structure-of-arrays data for efficient RDataFrame analysis.
   * * @note This is a **Base Class**. It is topology-agnostic. Specific physics scenarios 
   * (e.g., Electro-production vs. Photoproduction) should inherit from this class to configure 
   * the specific particle groups and creation rules.
   */
  class KinematicsProcessor {

  public:
    /**
     * @brief Output Type Alias.
     * Structure: `[Combination_Index][Component_Index][Particle_Index]`
     * Component Index maps to: 0=Px, 1=Py, 2=Pz, 3=Mass.
     */
    using CombiOutputVec_t = RVec<RVec<RVecResultType>>;

    // =================================================================================
    // Lifecycle & Initialization
    // =================================================================================

    /**
     * @brief Standard Constructor.
     * Initializes the engine with a connection to the RDataFrame configuration.
     * @param cr Pointer to the ConfigReaction instance.
     * @param suffix Optional suffix string (e.g., "_miss") appended to all output columns.
     * Use this when running multiple hypothesis processors side-by-side to avoid name collisions.
     */
    KinematicsProcessor(config::ConfigReaction* cr, const std::string& suffix = "");
    
    /**
     * @brief Copy Constructor (The "Fork").
     * Creates a new processor that clones the configuration of an existing one but assigns a new suffix.
     * * **Use Case:** "Twin Topology" Analysis.
     * 1. Define a standard topology (Master).
     * 2. Fork it using this constructor (Linked).
     * 3. Override specific groups (e.g. Baryons) in the Linked processor.
     * 4. Use `InitLinked()` to share data between them.
     */
    KinematicsProcessor(const KinematicsProcessor& other, const std::string& new_suffix);
    
    virtual ~KinematicsProcessor() = default;

    /**
     * @brief Initialize the Processor, calculate the ReactionMap, and register with RDataFrame.
     * * This method:
     * 1. Triggers `ParticleCreator::InitMap` to assign fixed indices to all particles.
     * 2. Initializes `ParticleModifier` to resolve string names to those indices.
     * 3. Registers the `DefineAux` columns (packed auxiliary data).
     * 4. Registers this functor via `DefineKinematicsProcessor` for execution.
     */
    void Init();
    
    /**
     * @brief Initialize as a "Linked" topology (Optimization).
     * * Instead of running the expensive combinatorial loop again, this method:
     * 1. Adopts the integer indices from the Master processor.
     * 2. Rebuilds the ReactionMap based on *local* group overrides (e.g. "Baryon" = "n_miss").
     * 3. **Aliases** the component data columns from the Master to this processor's suffix.
     * * @warning Modifiers defined on this processor will **NOT** run, because the event loop 
     * is skipped. All data creation/modification must happen in the Master.
     * @param master The fully initialized Master processor to borrow data from.
     */
    void InitLinked(const KinematicsProcessor& master);

    // =================================================================================
    // Core Execution
    // =================================================================================

    /**
     * @brief Functor Operator: The Core Event Loop.
     * * This function is called by RDataFrame for every event. It performs the heavy lifting:
     * 1. **Transposition**: Maps input tracks (random order) to the Reaction Cache (fixed order).
     * 2. **Pre-Modification**: Applies scaling/correction to input tracks.
     * 3. **Creation**: Calculates intermediate particles (e.g. Virtual Photon, J/psi).
     * 4. **Post-Modification**: Applies constraints to created particles.
     * * @tparam Tp Input momentum type (e.g. RVecF or RVecD).
     * @tparam Tm Input mass type.
     * @param indices Combinatorial indices `[Role][Combination]`.
     * @param px, py, pz, m Raw input component vectors.
     * @param aux_pre_d Packed double columns for Pre-Modifiers.
     * @param aux_pre_i Packed int columns for Pre-Modifiers.
     * @param aux_post_d Packed double columns for Post-Modifiers.
     * @param aux_post_i Packed int columns for Post-Modifiers.
     * @return CombiOutputVec_t The fully calculated, transposed kinematic data.
     */
    template<typename Tp, typename Tm> 
    CombiOutputVec_t operator()(const RVecIndices& indices, 
                                const Tp& px, const Tp& py, const Tp& pz, const Tm& m,
                                const RVecRVecD& aux_pre_d, const RVecRVecI& aux_pre_i,
                                const RVecRVecD& aux_post_d, const RVecRVecI& aux_post_i) const;

    /**
     * @brief Hook for unpacking components if needed.
     * Currently unused, but reserved for defining explicit `rec_px`, `rec_py` columns 
     * from the monolithic output if required.
     */
    void DefineNewComponentVecs();

    // =================================================================================
    // Configuration Helpers
    // =================================================================================

    /**
     * @brief Defines the packed auxiliary columns in RDataFrame.
     * Inspects the registered modifiers to see which columns (e.g. "rec_cal_energy") 
     * are needed, and packs them into `RVec<RVec>` for efficient passing to the operator.
     */
    void DefineAux(const std::string& type);

    // =================================================================================
    // Accessors
    // =================================================================================

    /** @brief Access the ParticleCreator (Architect). */
    config::ParticleCreator& Creator();
    const config::ParticleCreator& Creator() const;
    
    /** @brief Access the Pre-Creation Modifier Manager (Inputs). */
    config::ParticleModifier& PreModifier();
    
    /** @brief Access the Post-Creation Modifier Manager (Intermediates). */
    config::ParticleModifier& PostModifier();

    config::ConfigReaction* Reaction() const;
    std::string GetSuffix() const;
    
    /** @brief Format a base name with the processor's suffix (e.g. "Q2" -> "Q2_miss"). */
    std::string FullName(const std::string& baseName) const;

    // =================================================================================
    // Definition API (Adding Kinematics)
    // =================================================================================
    
    /**
     * @brief Define a kinematic variable using a JIT-string Kernel.
     * @param name Name of the output column.
     * @param func String name of the function (e.g., "rad::ElS_Q2").
     */
    void Define(const std::string& name, const std::string& func);

    /**
     * @brief Define a kinematic variable using a C++ Functor/Lambda.
     * Compiles the logic directly into the analysis loop for maximum performance.
     * @param name Name of the output column.
     * @param func A callable object or lambda matching the kinematic signature.
     */
    template <typename Lambda>
    void DefineKernel(const std::string& name, Lambda&& func);

    /**
     * @brief Define using JIT-string kernel with explicit input particles.
     * Automatically handles resolving particle names to the underlying suffixed index columns.
     */
    void Define(const std::string& name, const std::string& func, const std::vector<ParticleNames_t>& particles);

    // =================================================================================
    // Generic Shortcuts
    // =================================================================================

    /** @brief Calculate Invariant Mass of the specified group. */
    void Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});
    
    /** @brief Calculate Invariant Mass Squared (M2). */
    void Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});
    
    /** @brief Calculate Transverse Momentum (Pt). */
    void Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

  private:
    config::ConfigReaction* _reaction = nullptr;
    config::ParticleCreator _creator;
    
    // The Modifier Managers
    config::ParticleModifier _preModifier;
    config::ParticleModifier _postModifier;

    std::string _suffix;
  };

  // =================================================================================
  // IMPLEMENTATION
  // =================================================================================

  inline KinematicsProcessor::KinematicsProcessor(config::ConfigReaction* cr, const std::string& suffix) 
    : _reaction{cr}, _creator{cr, suffix}, _suffix{suffix} 
  {}

  inline KinematicsProcessor::KinematicsProcessor(const KinematicsProcessor& other, const std::string& new_suffix)
    : _reaction(other._reaction),
      _creator(other._creator, new_suffix),
      _preModifier(other._preModifier),   // Copy Config
      _postModifier(other._postModifier), // Copy Config
      _suffix(new_suffix)
  {}

  inline void KinematicsProcessor::Init() {
    auto types = _reaction->GetTypes();
    
    // 1. Build Map (assign fixed indices)
    Creator().InitMap(types);
    
    // 2. Initialize Modifiers (resolve names to those indices)
    _preModifier.Init(Creator());
    _postModifier.Init(Creator());

    for(const auto& type : types){
      // 3. Define Packed Aux Columns (must exist before processor is defined)
      DefineAux(type);
      // 4. Register Processor execution loop
      DefineKinematicsProcessor(*_reaction, *this, type);
    }
  }

  inline void KinematicsProcessor::DefineAux(const std::string& type) {
      
      // Helper to define RVec<RVec> column from a list of names
      auto define_pack = [&](const std::string& name, const std::vector<std::string>& cols, bool is_int) {
          if(cols.empty()) {
              // Define dummy empty column if no aux data needed
              if(is_int) _reaction->Define(name, [](){ return RVecRVecI{}; }, {});
              else       _reaction->Define(name, [](){ return RVecRVecD{}; }, {});
              return;
          }
          _reaction->Define(name, rad::utils::createPackVectorString(cols));
      };

      // Define 4 Packs (Pre/Post x Double/Int) with unique suffix names
      define_pack(type + "aux_pre_d" + _suffix+config::DoNotWriteTag(), _preModifier.GetAuxDoubleCols(), false);
      define_pack(type + "aux_pre_i" + _suffix+config::DoNotWriteTag(), _preModifier.GetAuxIntCols(), true);
      
      define_pack(type + "aux_post_d" + _suffix+config::DoNotWriteTag(), _postModifier.GetAuxDoubleCols(), false);
      define_pack(type + "aux_post_i" + _suffix+config::DoNotWriteTag(), _postModifier.GetAuxIntCols(), true);
  }

  inline void KinematicsProcessor::InitLinked(const KinematicsProcessor& master) {
    if(_suffix == master.GetSuffix()) 
        throw std::runtime_error("InitLinked Error: Linked processor must have a different suffix ('" 
                                 + _suffix + "') than Master.");

    // 1. Adopt Indices: Use Master's integer mappings for consistency
    _creator.AdoptIndices(master.Creator());
    
    // 2. Rebuild Map: Use Local groups (e.g. n_miss vs n) to define the new topology map
    _creator.RebuildReactionMap();

    // 3. Alias Data Columns: Point our names (e.g. rec_kine_comps_miss) to Master's data
    auto types = _reaction->GetTypes();
    for(const auto& type : types) {
        std::string src_comps = type + names::KineComponents() + master.GetSuffix();
        std::string dest_comps = type + names::KineComponents() + _suffix;
        _reaction->Define(dest_comps, src_comps);
        
        std::string src_inds = type + names::KineIndices() + master.GetSuffix();
        std::string dest_inds = type + names::KineIndices() + _suffix;
        _reaction->Define(dest_inds, src_inds);

        // Alias helper columns (e.g. kine_Jpsi_miss -> kine_Jpsi)
        for(auto const& [name, idx] : _creator.GetIndexMap()) {
             std::string h_src = names::data_type::Kine() + name + master.GetSuffix();
             std::string h_dest = names::data_type::Kine() + name + _suffix;
             _reaction->Define(h_dest, h_src);
        }
    }
  }

  inline config::ParticleCreator& KinematicsProcessor::Creator() { return _creator; }
  inline const config::ParticleCreator& KinematicsProcessor::Creator() const { return _creator; }
  
  inline config::ParticleModifier& KinematicsProcessor::PreModifier() { return _preModifier; }
  inline config::ParticleModifier& KinematicsProcessor::PostModifier() { return _postModifier; }

  inline config::ConfigReaction* KinematicsProcessor::Reaction() const { return _reaction; }
  inline std::string KinematicsProcessor::GetSuffix() const { return _suffix; }
  inline std::string KinematicsProcessor::FullName(const std::string& baseName) const { return baseName + _suffix; }

  inline void KinematicsProcessor::DefineNewComponentVecs() {}

  // --- Core Operator (Template) ---
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

    // 1. Momentum Cache (Stack Allocated for Speed)
    // Used for in-place modification before writing to result
    ROOT::RVecD temp_px(Nparticles, constant::InvalidEntry<double>());
    ROOT::RVecD temp_py(Nparticles, constant::InvalidEntry<double>());
    ROOT::RVecD temp_pz(Nparticles, constant::InvalidEntry<double>());
    ROOT::RVecD temp_m(Nparticles, constant::InvalidEntry<double>());

    // 2. Aux Cache (Transposed Local Buffer)
    // Structure: [Variable_Index][Fixed_Particle_Index]
    // We explicitly transpose these to ensure modifiers have O(1) access to correct values.
    config::AuxCacheD cache_pre_d(aux_pre_d.size(), ROOT::RVecD(Nparticles));
    config::AuxCacheI cache_pre_i(aux_pre_i.size(), ROOT::RVecI(Nparticles));
    config::AuxCacheD cache_post_d(aux_post_d.size(), ROOT::RVecD(Nparticles));
    config::AuxCacheI cache_post_i(aux_post_i.size(), ROOT::RVecI(Nparticles));
    
    for (size_t icombi = 0; icombi < Ncombis; ++icombi) {
      
      // --- Load Inputs & Transpose ---
      for (size_t ip = 0; ip < Nparticles0; ++ip) {
        size_t iparti = _creator.GetReactionIndex(ip);    // Fixed Target Index                
        const int original_index = indices[ip][icombi];   // Source Input Index

        // Copy Momentum
        temp_px[iparti] = px[original_index];
        temp_py[iparti] = py[original_index];
        temp_pz[iparti] = pz[original_index];
        temp_m[iparti]  = m[original_index];

        // Transpose Pre-Aux Data
        // Maps unordered detector data (original_index) to ordered reaction data (iparti)
        for(size_t v=0; v<aux_pre_d.size(); ++v) cache_pre_d[v][iparti] = aux_pre_d[v][original_index];
        for(size_t v=0; v<aux_pre_i.size(); ++v) cache_pre_i[v][iparti] = aux_pre_i[v][original_index];
      }

      // --- Execution Pipeline ---

      // A. Pre-Modification: Modify Inputs (Scale, Calibrate)
      _preModifier.Apply(temp_px, temp_py, temp_pz, temp_m, cache_pre_d, cache_pre_i);

      // B. Creation: Build Intermediate Particles (using modified inputs)
      _creator.ApplyCreation(temp_px, temp_py, temp_pz, temp_m);

      // C. Post-Modification: Modify Created Particles (Mass Constraints)
      _postModifier.Apply(temp_px, temp_py, temp_pz, temp_m, cache_post_d, cache_post_i);

      // --- Finalize: Write to Output ---
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
    auto types = _reaction->GetTypes();
    for(const auto& type : types){
      // Use FullName(name) to append suffix
      _reaction->Define(type + FullName(name), utils::createFunctionCallStringFromVec("rad::util::ApplyKinematics", 
                                   {func, Creator().GetMapName(), (type+names::KineComponents() + _suffix)}));
    }
  }
  
  template <typename Lambda>
  inline void KinematicsProcessor::DefineKernel(const std::string& name, Lambda&& func) {
    auto types = _reaction->GetTypes();
    for(const auto& type : types){
      ROOT::RDF::ColumnNames_t cols = { Creator().GetMapName(), type + names::KineComponents() + _suffix };
      auto apply_func = [func](const RVecIndexMap& map, const ROOT::RVec<ROOT::RVec<RVecResultType>>& comps){
        return rad::util::ApplyKinematics(func, map, comps);
      };
      _reaction->Define(type + FullName(name), apply_func, cols);
    }
  }

  inline void KinematicsProcessor::Define(const std::string& name, const std::string& func, const std::vector<ParticleNames_t>& particles) {
    auto types = _reaction->GetTypes();
    std::string kine_parts = "{";
    for(const auto& pnames : particles){
      ParticleNames_t suffixed_names;
      for(const auto& p : pnames) suffixed_names.push_back(p + _suffix);
      kine_parts += utils::combineVectorToString(utils::prependToAll(suffixed_names, names::data_type::Kine()));
      kine_parts += ",";
    }
    kine_parts.pop_back(); kine_parts += "}";
    
    for(const auto& type : types){
      _reaction->Define(type + FullName(name), utils::createFunctionCallStringFromVec("rad::util::ApplyKinematics", 
                                   {func, kine_parts, (type + names::KineComponents() + _suffix)}));
    }
  }

  // --- Shortcuts ---
  inline void KinematicsProcessor::Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    // KineCalculation constructor handles name resolution and kernel compilation
    KineCalculation calc(*this, name, rad::FourVectorMassCalc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
  }
  inline void KinematicsProcessor::Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    KineCalculation calc(*this, name, rad::FourVectorMass2Calc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
  }
  inline void KinematicsProcessor::Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    KineCalculation calc(*this, name, rad::FourVectorPtCalc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
  }

} // namespace rad



// namespace rad {

//   using ROOT::RVec;
 
//   /**
//    * @brief Master functor for executing combinatorial kinematic analysis on a single event.
//    * * This class acts as the "Engine" of the combinatorial system. It:
//    * 1. Orchestrates particle index grouping.
//    * 2. Defines the static ReactionMap (fixed output order).
//    * 3. Executes the thread-safe creation loop (modifying momenta, creating intermediate particles).
//    * 4. Returns transposed momentum components for downstream RDataFrame processing.
//    * * @note This is the BASE CLASS. It contains no specific physics topology (like Electro/Photo production).
//    * Specific reaction topologies should inherit from this class.
//    */
//   class KinematicsProcessor {

//   public:
        
//     // Aliases for cleaner code within the class
//     // Output format: [Combination][Component(x,y,z,m)][ParticleIndex]
//     using CombiOutputVec_t = RVec<RVec<RVecResultType>>;

//     // =================================================================================
//     // Lifecycle & Initialization
//     // =================================================================================

//     /**
//      * @brief Constructs the processor.
//      * @param cr Pointer to the ConfigReaction instance for RDataFrame access.
//      */
//     KinematicsProcessor(config::ConfigReaction* cr);
    
//     virtual ~KinematicsProcessor() = default;

//     /**
//      * @brief Initialise Processor Creator and Register with RDataFrame.
//      * * This method:
//      * 1. Triggers the ParticleCreator to calculate the static ReactionMap.
//      * 2. Registers this functor instance via `DefineKinematicsProcessor` for each 
//      * data type (reconstructed, truth, etc.) found in ConfigReaction.
//      */
//     void Init();

//     // =================================================================================
//     // Core Execution (Functor)
//     // =================================================================================

//     /**
//      * @brief Functor operator: The Core Event Loop.
//      * * Executes once per event to process all combinations.
//      * This function performs the core calculation: transposition, modification, and creation 
//      * using temporary stack-allocated vectors for cache-friendliness.
//      * * @tparam Tp Type of momentum component RVecs (e.g., RVec<float> or RVec<double>).
//      * @tparam Tm Type of mass component RVec.
//      * @param indices RVecIndices[role][combo]: The index for each particle/combo pair.
//      * @param px, py, pz, m The primary component vectors read from the input tree.
//      * @return CombiOutputVec_t The transposed and expanded momentum data.
//      */
//     template<typename Tp, typename Tm> 
//     CombiOutputVec_t operator()(const RVecIndices& indices, const Tp& px, const Tp& py, const Tp& pz, const Tm& m) const;

//     /**
//      * @brief Placeholder for unpacking components if needed.
//      * Can be used to define individual columns (e.g. "rec_px") from the monolithic output.
//      */
//     void DefineNewComponentVecs();

//     // =================================================================================
//     // Accessors
//     // =================================================================================

//     config::ParticleCreator& Creator();
//     const config::ParticleCreator& Creator() const;
//     config::ConfigReaction* Reaction() const;

//     // =================================================================================
//     // Definition Interfaces (API for adding Kinematics)
//     // =================================================================================
    
//     /**
//      * @brief Define a kinematic variable using a string-based Kernel name.
//      * * Example: Define("Q2", "rad::ElS_Q2<...>")
//      * This uses JIT compilation via RDataFrame's string parsing.
//      * * @param name The name of the output column (e.g., "Q2").
//      * @param func The string name of the function to apply.
//      */
//     void Define(const std::string& name, const std::string& func);

//     /**
//      * @brief Define a kinematic variable using a C++ Functor/Lambda (Compiled).
//      * * Use this for maximum performance or inline logic. 
//      * The functor is captured by value into the RDataFrame loop.
//      * * @tparam Lambda The type of the functor (deduced automatically).
//      * @param name The name of the output column.
//      * @param func A callable with signature: 
//      * ResultType f(const RVecIndexMap&, const RVecD& px, py, pz, m)
//      */
//     template <typename Lambda>
//     void DefineKernel(const std::string& name, Lambda&& func);

//     /**
//      * @brief Define using string kernel with specific input particles.
//      * * Helper that constructs the argument string for ApplyKinematics using 
//      * the provided particle lists.
//      */
//     void Define(const std::string& name, const std::string& func, const std::vector<ParticleNames_t>& particles);

//     // =================================================================================
//     // Generic Shortcuts (Topology Agnostic)
//     // =================================================================================

//     /**
//      * @brief Calculates the Invariant Mass of the specified particle group.
//      * * Uses `DefineKernel` to compile the specific indices into the loop.
//      * @param name Name of the output column (e.g., "JpsiMass").
//      * @param particles_pos List of particles to ADD (Sum P4).
//      * @param particles_neg List of particles to SUBTRACT (Diff P4).
//      */
//     void Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

//     /**
//      * @brief Calculates the Invariant Mass Squared (M2) of the specified group.
//      */
//     void Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

//     /**
//      * @brief Calculates the Transverse Momentum (Pt) of the specified group.
//      */
//     void Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

//   private:
//     config::ConfigReaction* _reaction = nullptr;
//     config::ParticleCreator _creator;
//   };

//   // =================================================================================
//   // IMPLEMENTATION
//   // =================================================================================

//   inline KinematicsProcessor::KinematicsProcessor(config::ConfigReaction* cr) 
//     : _reaction{cr}, _creator{cr} 
//   {
//     // Constructor purely initializes the Engine.
//     // Physics topology (e.g. VirtGamma creation) must be done in derived classes.
//   }

//   inline void KinematicsProcessor::Init() {
//     // 1. Init the indices and reaction map based on Creator's definitions
//     auto types = _reaction->GetTypes();
//     Creator().InitMap(types);

//     // 2. Define processor for each type (rec_, tru_, etc.)
//     for(const auto& type : types){
//       DefineKinematicsProcessor(*_reaction, *this, type);
//     }
//   }

//   inline config::ParticleCreator& KinematicsProcessor::Creator() { return _creator; }
//   inline const config::ParticleCreator& KinematicsProcessor::Creator() const { return _creator; }
//   inline config::ConfigReaction* KinematicsProcessor::Reaction() const { return _reaction; }

//   inline void KinematicsProcessor::DefineNewComponentVecs() {
//     // Logic to unpack back to individual arrays could go here if needed.
//     // e.g. _reaction->Define("rec_px", ...);
//   }

//   // --- Core Operator (Template) ---
//   template<typename Tp, typename Tm> 
//   inline KinematicsProcessor::CombiOutputVec_t KinematicsProcessor::operator()(const RVecIndices& indices, const Tp& px, const Tp& py, const Tp& pz, const Tm& m) const {
          
//     const auto Ncomponents = 4; // px, py, pz, m
//     const auto Nparticles0 = indices.size(); // Number of Input particle roles
//     const auto Nparticles = Nparticles0 + _creator.GetNCreated(); // Total including created
    
//     if (Nparticles == 0) {
//       // Return structured array with zero size if no particles/combos are found.
//       return CombiOutputVec_t(Ncomponents); 
//     }
          
//     // Assuming indices[i_role] is consistent in size across all roles.
//     const auto Ncombis = indices[0].size(); 

//     // Initialize transposed result: [combination][components][particles]
//     // This pre-allocation is critical for performance.
//     CombiOutputVec_t result(Ncombis, RVec<RVecResultType>(Ncomponents, RVecResultType(Nparticles)));

//     // Temporary component vectors
//     // Allocated on stack/local scope for cache friendliness during modification loop.
//     // We use InvalidEntry to detect uninitialized access if necessary.
//     ROOT::RVecD temp_px(Nparticles, constant::InvalidEntry<double>());
//     ROOT::RVecD temp_py(Nparticles, constant::InvalidEntry<double>());
//     ROOT::RVecD temp_pz(Nparticles, constant::InvalidEntry<double>());
//     ROOT::RVecD temp_m(Nparticles, constant::InvalidEntry<double>());
    
//     // --- Loop over Combinations ---
//     for (size_t icombi = 0; icombi < Ncombis; ++icombi) {
      
//       // 1. Transpose and Copy Original Components
//       // For each combi fill px[iparti]...
//       for (size_t ip = 0; ip < Nparticles0; ++ip) {
//         size_t iparti = _creator.GetReactionIndex(ip); // Get fixed target index                   
//         const int original_index = indices[ip][icombi]; // Get source index for this combo

//         // Copy primary components from input arrays (px, py, pz, m)
//         temp_px[iparti] = px[original_index];
//         temp_py[iparti] = py[original_index];
//         temp_pz[iparti] = pz[original_index];
//         temp_m[iparti]  = m[original_index];
//       }

//       // 2. Apply Particle Creator (Modifications/Intermediate Particles)
//       // This calculates things like J/psi or VirtGamma and updates the temp vectors in place.
//       _creator.ApplyCreation(temp_px, temp_py, temp_pz, temp_m);

//       // 3. Finalize: Record components in correct order
//       // Copy the local cache into the complex output structure.
//       result[icombi][OrderX()] = temp_px;
//       result[icombi][OrderY()] = temp_py;
//       result[icombi][OrderZ()] = temp_pz;
//       result[icombi][OrderM()] = temp_m;
//     }

//     // Return the transposed and expanded momentum data: [combination][components][particles]
//     return result;
//   }

//   // =================================================================================
//   // Definitions (String & Kernel)
//   // =================================================================================

//   inline void KinematicsProcessor::Define(const std::string& name, const std::string& func) {
//     auto types = _reaction->GetTypes();
//     for(const auto& type : types){
//       // Use string manipulation to build the C++ call: ApplyKinematics(func, map, comps)
//       _reaction->Define(type+name, utils::createFunctionCallStringFromVec("rad::util::ApplyKinematics", 
//                                    {func, names::ReactionMap(), (type+names::KineComponents())}));
//     }
//   }
  
//   template <typename Lambda>
//   inline void KinematicsProcessor::DefineKernel(const std::string& name, Lambda&& func) {
//     auto types = _reaction->GetTypes();
//     for(const auto& type : types){
      
//       // Define the input columns required by the wrapper
//       ROOT::RDF::ColumnNames_t cols = {
//         names::ReactionMap(), 
//         type + names::KineComponents()
//       };
      
//       // Create the Wrapper
//       // Capture 'func' by VALUE [func] for thread safety.
//       auto apply_func = [func](const RVecIndexMap& map, const ROOT::RVec<ROOT::RVec<RVecResultType>>& comps){
//         return rad::util::ApplyKinematics(func, map, comps);
//       };
      
//       // Define via ConfigReaction
//       _reaction->Define(type+name, apply_func, cols);
//     }
//   }

//   inline void KinematicsProcessor::Define(const std::string& name, const std::string& func, const std::vector<ParticleNames_t>& particles) {
//     auto types = _reaction->GetTypes();
    
//     // Convert the vector of particle names into a string representation for the JIT compiler
//     std::string kine_parts = "{";
//     for(const auto& pnames : particles){
//       kine_parts += utils::combineVectorToString(utils::prependToAll(pnames, names::data_type::Kine()));
//       kine_parts += ",";
//     }
//     kine_parts.pop_back();
//     kine_parts += "}";
    
//     for(const auto& type : types){
//       _reaction->Define(type+name, utils::createFunctionCallStringFromVec("rad::util::ApplyKinematics", 
//                                    {func, kine_parts, (type+names::KineComponents())}));
//     }
//   }

//   // =================================================================================
//   // Generic Shortcuts
//   // =================================================================================

//  inline void KinematicsProcessor::Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
//     // Explicitly instantiate the template so it matches IndexKernel signature
//     KineCalculation calc(*this, name, 
//                          rad::FourVectorMassCalc<rad::RVecResultType, rad::RVecResultType>, 
//                          {particles_pos, particles_neg});
//   }

//   inline void KinematicsProcessor::Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
//     KineCalculation calc(*this, name, 
//                          rad::FourVectorMass2Calc<rad::RVecResultType, rad::RVecResultType>, 
//                          {particles_pos, particles_neg});
//   }

//   inline void KinematicsProcessor::Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
//     KineCalculation calc(*this, name, 
//                          rad::FourVectorPtCalc<rad::RVecResultType, rad::RVecResultType>, 
//                          {particles_pos, particles_neg});
//   }
// } // namespace rad





