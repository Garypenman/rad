#pragma once

#include "ConfigReaction.h"
#include "StringUtilities.h"
#include "RVecHelpers.h"
#include "ParticleCreator.h"
#include "KinematicsDispatch.h"
#include "KineCalculation.h"
#include "BasicKinematics.h"

namespace rad {

  using ROOT::RVec;
 
  /**
   * @brief Master functor for executing combinatorial kinematic analysis on a single event.
   * * This class acts as the "Engine" of the combinatorial system. It:
   * 1. Orchestrates particle index grouping.
   * 2. Defines the static ReactionMap (fixed output order).
   * 3. Executes the thread-safe creation loop (modifying momenta, creating intermediate particles).
   * 4. Returns transposed momentum components for downstream RDataFrame processing.
   * * @note This is the BASE CLASS. It contains no specific physics topology (like Electro/Photo production).
   * Specific reaction topologies should inherit from this class.
   */
  class KinematicsProcessor {

  public:
        
    // Aliases for cleaner code within the class
    // Output format: [Combination][Component(x,y,z,m)][ParticleIndex]
    using CombiOutputVec_t = RVec<RVec<RVecResultType>>;

    // =================================================================================
    // Lifecycle & Initialization
    // =================================================================================

    /**
     * @brief Constructs the processor.
     * @param cr Pointer to the ConfigReaction instance for RDataFrame access.
     */
    KinematicsProcessor(config::ConfigReaction* cr);
    
    virtual ~KinematicsProcessor() = default;

    /**
     * @brief Initialise Processor Creator and Register with RDataFrame.
     * * This method:
     * 1. Triggers the ParticleCreator to calculate the static ReactionMap.
     * 2. Registers this functor instance via `DefineKinematicsProcessor` for each 
     * data type (reconstructed, truth, etc.) found in ConfigReaction.
     */
    void Init();

    // =================================================================================
    // Core Execution (Functor)
    // =================================================================================

    /**
     * @brief Functor operator: The Core Event Loop.
     * * Executes once per event to process all combinations.
     * This function performs the core calculation: transposition, modification, and creation 
     * using temporary stack-allocated vectors for cache-friendliness.
     * * @tparam Tp Type of momentum component RVecs (e.g., RVec<float> or RVec<double>).
     * @tparam Tm Type of mass component RVec.
     * @param indices RVecIndices[role][combo]: The index for each particle/combo pair.
     * @param px, py, pz, m The primary component vectors read from the input tree.
     * @return CombiOutputVec_t The transposed and expanded momentum data.
     */
    template<typename Tp, typename Tm> 
    CombiOutputVec_t operator()(const RVecIndices& indices, const Tp& px, const Tp& py, const Tp& pz, const Tm& m) const;

    /**
     * @brief Placeholder for unpacking components if needed.
     * Can be used to define individual columns (e.g. "rec_px") from the monolithic output.
     */
    void DefineNewComponentVecs();

    // =================================================================================
    // Accessors
    // =================================================================================

    config::ParticleCreator& Creator();
    const config::ParticleCreator& Creator() const;
    config::ConfigReaction* Reaction() const;

    // =================================================================================
    // Definition Interfaces (API for adding Kinematics)
    // =================================================================================
    
    /**
     * @brief Define a kinematic variable using a string-based Kernel name.
     * * Example: Define("Q2", "rad::ElS_Q2<...>")
     * This uses JIT compilation via RDataFrame's string parsing.
     * * @param name The name of the output column (e.g., "Q2").
     * @param func The string name of the function to apply.
     */
    void Define(const std::string& name, const std::string& func);

    /**
     * @brief Define a kinematic variable using a C++ Functor/Lambda (Compiled).
     * * Use this for maximum performance or inline logic. 
     * The functor is captured by value into the RDataFrame loop.
     * * @tparam Lambda The type of the functor (deduced automatically).
     * @param name The name of the output column.
     * @param func A callable with signature: 
     * ResultType f(const RVecIndexMap&, const RVecD& px, py, pz, m)
     */
    template <typename Lambda>
    void DefineKernel(const std::string& name, Lambda&& func);

    /**
     * @brief Define using string kernel with specific input particles.
     * * Helper that constructs the argument string for ApplyKinematics using 
     * the provided particle lists.
     */
    void Define(const std::string& name, const std::string& func, const std::vector<ParticleNames_t>& particles);

    // =================================================================================
    // Generic Shortcuts (Topology Agnostic)
    // =================================================================================

    /**
     * @brief Calculates the Invariant Mass of the specified particle group.
     * * Uses `DefineKernel` to compile the specific indices into the loop.
     * @param name Name of the output column (e.g., "JpsiMass").
     * @param particles_pos List of particles to ADD (Sum P4).
     * @param particles_neg List of particles to SUBTRACT (Diff P4).
     */
    void Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

    /**
     * @brief Calculates the Invariant Mass Squared (M2) of the specified group.
     */
    void Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

    /**
     * @brief Calculates the Transverse Momentum (Pt) of the specified group.
     */
    void Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

  private:
    config::ConfigReaction* _reaction = nullptr;
    config::ParticleCreator _creator;
  };

  // =================================================================================
  // IMPLEMENTATION
  // =================================================================================

  inline KinematicsProcessor::KinematicsProcessor(config::ConfigReaction* cr) 
    : _reaction{cr}, _creator{cr} 
  {
    // Constructor purely initializes the Engine.
    // Physics topology (e.g. VirtGamma creation) must be done in derived classes.
  }

  inline void KinematicsProcessor::Init() {
    // 1. Init the indices and reaction map based on Creator's definitions
    auto types = _reaction->GetTypes();
    Creator().InitMap(types);

    // 2. Define processor for each type (rec_, tru_, etc.)
    for(const auto& type : types){
      DefineKinematicsProcessor(*_reaction, *this, type);
    }
  }

  inline config::ParticleCreator& KinematicsProcessor::Creator() { return _creator; }
  inline const config::ParticleCreator& KinematicsProcessor::Creator() const { return _creator; }
  inline config::ConfigReaction* KinematicsProcessor::Reaction() const { return _reaction; }

  inline void KinematicsProcessor::DefineNewComponentVecs() {
    // Logic to unpack back to individual arrays could go here if needed.
    // e.g. _reaction->Define("rec_px", ...);
  }

  // --- Core Operator (Template) ---
  template<typename Tp, typename Tm> 
  inline KinematicsProcessor::CombiOutputVec_t KinematicsProcessor::operator()(const RVecIndices& indices, const Tp& px, const Tp& py, const Tp& pz, const Tm& m) const {
          
    const auto Ncomponents = 4; // px, py, pz, m
    const auto Nparticles0 = indices.size(); // Number of Input particle roles
    const auto Nparticles = Nparticles0 + _creator.GetNCreated(); // Total including created
    
    if (Nparticles == 0) {
      // Return structured array with zero size if no particles/combos are found.
      return CombiOutputVec_t(Ncomponents); 
    }
          
    // Assuming indices[i_role] is consistent in size across all roles.
    const auto Ncombis = indices[0].size(); 

    // Initialize transposed result: [combination][components][particles]
    // This pre-allocation is critical for performance.
    CombiOutputVec_t result(Ncombis, RVec<RVecResultType>(Ncomponents, RVecResultType(Nparticles)));

    // Temporary component vectors
    // Allocated on stack/local scope for cache friendliness during modification loop.
    // We use InvalidEntry to detect uninitialized access if necessary.
    ROOT::RVecD temp_px(Nparticles, constant::InvalidEntry<double>());
    ROOT::RVecD temp_py(Nparticles, constant::InvalidEntry<double>());
    ROOT::RVecD temp_pz(Nparticles, constant::InvalidEntry<double>());
    ROOT::RVecD temp_m(Nparticles, constant::InvalidEntry<double>());
    
    // --- Loop over Combinations ---
    for (size_t icombi = 0; icombi < Ncombis; ++icombi) {
      
      // 1. Transpose and Copy Original Components
      // For each combi fill px[iparti]...
      for (size_t ip = 0; ip < Nparticles0; ++ip) {
        size_t iparti = _creator.GetReactionIndex(ip); // Get fixed target index                   
        const int original_index = indices[ip][icombi]; // Get source index for this combo

        // Copy primary components from input arrays (px, py, pz, m)
        temp_px[iparti] = px[original_index];
        temp_py[iparti] = py[original_index];
        temp_pz[iparti] = pz[original_index];
        temp_m[iparti]  = m[original_index];
      }

      // 2. Apply Particle Creator (Modifications/Intermediate Particles)
      // This calculates things like J/psi or VirtGamma and updates the temp vectors in place.
      _creator.ApplyCreation(temp_px, temp_py, temp_pz, temp_m);

      // 3. Finalize: Record components in correct order
      // Copy the local cache into the complex output structure.
      result[icombi][OrderX()] = temp_px;
      result[icombi][OrderY()] = temp_py;
      result[icombi][OrderZ()] = temp_pz;
      result[icombi][OrderM()] = temp_m;
    }

    // Return the transposed and expanded momentum data: [combination][components][particles]
    return result;
  }

  // =================================================================================
  // Definitions (String & Kernel)
  // =================================================================================

  inline void KinematicsProcessor::Define(const std::string& name, const std::string& func) {
    auto types = _reaction->GetTypes();
    for(const auto& type : types){
      // Use string manipulation to build the C++ call: ApplyKinematics(func, map, comps)
      _reaction->Define(type+name, utils::createFunctionCallStringFromVec("rad::util::ApplyKinematics", 
                                   {func, names::ReactionMap(), (type+names::KineComponents())}));
    }
  }
  
  template <typename Lambda>
  inline void KinematicsProcessor::DefineKernel(const std::string& name, Lambda&& func) {
    auto types = _reaction->GetTypes();
    for(const auto& type : types){
      
      // Define the input columns required by the wrapper
      ROOT::RDF::ColumnNames_t cols = {
        names::ReactionMap(), 
        type + names::KineComponents()
      };
      
      // Create the Wrapper
      // Capture 'func' by VALUE [func] for thread safety.
      auto apply_func = [func](const RVecIndexMap& map, const ROOT::RVec<ROOT::RVec<RVecResultType>>& comps){
        return rad::util::ApplyKinematics(func, map, comps);
      };
      
      // Define via ConfigReaction
      _reaction->Define(type+name, apply_func, cols);
    }
  }

  inline void KinematicsProcessor::Define(const std::string& name, const std::string& func, const std::vector<ParticleNames_t>& particles) {
    auto types = _reaction->GetTypes();
    
    // Convert the vector of particle names into a string representation for the JIT compiler
    std::string kine_parts = "{";
    for(const auto& pnames : particles){
      kine_parts += utils::combineVectorToString(utils::prependToAll(pnames, names::data_type::Kine()));
      kine_parts += ",";
    }
    kine_parts.pop_back();
    kine_parts += "}";
    
    for(const auto& type : types){
      _reaction->Define(type+name, utils::createFunctionCallStringFromVec("rad::util::ApplyKinematics", 
                                   {func, kine_parts, (type+names::KineComponents())}));
    }
  }

  // =================================================================================
  // Generic Shortcuts
  // =================================================================================

 inline void KinematicsProcessor::Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    // Explicitly instantiate the template so it matches IndexKernel signature
    KineCalculation calc(*this, name, 
                         rad::FourVectorMassCalc<rad::RVecResultType, rad::RVecResultType>, 
                         {particles_pos, particles_neg});
  }

  inline void KinematicsProcessor::Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    KineCalculation calc(*this, name, 
                         rad::FourVectorMass2Calc<rad::RVecResultType, rad::RVecResultType>, 
                         {particles_pos, particles_neg});
  }

  inline void KinematicsProcessor::Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    KineCalculation calc(*this, name, 
                         rad::FourVectorPtCalc<rad::RVecResultType, rad::RVecResultType>, 
                         {particles_pos, particles_neg});
  }
} // namespace rad








// #pragma once
// #include "ConfigReaction.h"
// #include "StringUtilities.h"
// #include "RVecHelpers.h"
// #include "ParticleCreator.h"
// #include "KinematicsDispatch.h"
// #include "KineCalculation.h"
// #include "ElectroIonReaction.h"
// #include "BasicKinematics.h"
// #include "ElectronScatterKinematics.h"

// namespace rad {

//   using ROOT::RVec;
 
//   /**
//    * @brief Master functor for executing combinatorial kinematic analysis on a single event.
//    * * This class orchestrates particle index grouping, static map definition, 
//    * momentum modification, and intermediate particle creation, returning the transposed 
//    * momentum components needed for downstream RDataFrame processing.
//    */
//   class KinematicsProcessor {

//   public:
        
//     // Aliases for cleaner code within the class
//     using CombiOutputVec_t = RVec<RVec<RVecResultType>>;

//     /**
//      * @brief Constructs the processor and defines the fixed (static) ReactionMap.
//      * * This defines the final order of particles in the output momentum arrays, 
//      * ensuring fixed lookup positions regardless of combinatorial variations.
//      * @param cr Pointer to the ConfigReaction instance for RDataFrame access.
//      */
//     KinematicsProcessor(config::ConfigReaction* cr) : _reaction{cr},_creator{cr} {
//       //First created particle is virtual photon
//       Creator().Diff(names::VirtGamma(),{{names::BeamEle()},{names::ScatEle()}});
 
   
//     }

//     /**
//      * Initialise Processor Creator
//      * Pass functor instance to RDF::Define
//      */
//     void Init(){

//       //Init the indices and reaction map
//       auto types = _reaction->GetTypes();
//       Creator().InitMap(types);

//       //Define processor for each type
//       //
//       for(const auto& type:types){
// 	DefineKinematicsProcessor(*_reaction,*this,type);
//       }
    
//     }
//     /**
//          * @brief Functor operator: Executes once per event to process all combinations.
//          * * This function performs the core calculation: transposition, modification, and creation 
//          * using temporary vectors for cache-friendliness.
//          * * @tparam Tp Type of momentum component RVecs (e.g., RVec<double>).
//          * @tparam Tm Type of mass component RVec.
//          * @param indices RVecIndices[role][combo]: The index for each particle/combo pair.
//          * @param px, py, pz, m The primary component vectors read from the input tree.
//          * @return CombiOutputVec_t The transposed and expanded momentum data: [combination][components][particles].
//          */
//     template<typename Tp, typename Tm> 
//     CombiOutputVec_t operator()(const RVecIndices& indices, const Tp& px, const Tp& py, const Tp& pz, const Tm& m) const {
            
//       const auto Ncomponents = 4; // px, py, pz, m
//       const auto Nparticles0 = indices.size();
//       const auto Nparticles = Nparticles0 + _creator.GetNCreated();
//       if (Nparticles == 0) {
// 	// Return structured array with zero size if no particles/combos are found.
// 	return CombiOutputVec_t(Ncomponents); 
//       }
            
//       // Assuming indices[i_role] is consistent in size across all roles.
//       const auto Ncombis = indices[0].size(); 

//       // Initialize transposed result: [combination][components][particles]
//       CombiOutputVec_t result(Ncombis, RVec<RVecResultType>(Ncomponents, RVecResultType(Nparticles)));

//       //temporary component vectors
//       //for cache freindliness for ParticleModifier etc 
//       ROOT::RVecD temp_px(Nparticles,constant::InvalidEntry<double>());
//       ROOT::RVecD temp_py(Nparticles,constant::InvalidEntry<double>());
//       ROOT::RVecD temp_pz(Nparticles,constant::InvalidEntry<double>());
//       ROOT::RVecD temp_m(Nparticles,constant::InvalidEntry<double>());
      
//       // --- 1. Transpose and Copy Original Components ---
//       for (size_t icombi = 0; icombi < Ncombis; ++icombi) {
// 	//For each combi fill px[iparti],...
	
// 	for (size_t ip = 0; ip < Nparticles0; ++ip) {
// 	  size_t iparti = _creator.GetReactionIndex(ip);                    
// 	  const int original_index = indices[ip][icombi];
// 	  // cout<<"kine process "<<ip<<" "<<iparti<<" "<<original_index<<endl;
// 	 // Copy primary components from input arrays (px, py, pz, m) using the combo index.
// 	  temp_px[iparti]=px[original_index];
// 	  temp_py[iparti]=py[original_index];
// 	  temp_pz[iparti]=pz[original_index];
// 	  temp_m[iparti]=m[original_index];
// 	  // --- 2. Apply Particle Modifications (TO BE IMPLEMENTED) ---
// 	  // e.g., result[0][iparti][icombi] = _modifier.apply(result[0][iparti][icombi]);
	  
// 	}
// 	// --- 3. Apply Particle Creator---
// 	_creator.ApplyCreation(temp_px,temp_py,temp_pz,temp_m);

// 	//--- 4. This combi components are now finalised-> record them in correct order
// 	result[icombi][OrderX()] = temp_px;
// 	result[icombi][OrderY()] = temp_py;
// 	result[icombi][OrderZ()] = temp_pz;
// 	result[icombi][OrderM()] = temp_m;

	
//       }
//       // CombiOutputVec_t The transposed and expanded momentum data: [combination][components][particles]
//       return result;
//     }

//     /**
//      * @brief Defines the final component columns by unpacking the functor's multi-dimensional output.
//      * * This is the canonical way to convert the Functor's nested RVec return into 
//      * separate, high-performance RDataFrame columns (combi_px, combi_py, etc.).
//      */
//     void DefineNewComponentVecs() {
//       //_reaction->Define("kine_px", [](const CombiOutputVec_t& comps){ return comps[0]; }, {names::KineComponents()});
//       // _reaction->Define("kine_py", [](const CombiOutputVec_t& comps){ return comps[1]; }, {names::KineComponents()});
//       // _reaction->Define("kine_pz", [](const CombiOutputVec_t& comps){ return comps[2]; }, {names::KineComponents()});
//       // _reaction->Define("kine_m", [](const CombiOutputVec_t& comps){ return comps[3]; }, {names::KineComponents()});
//     }

//     config::ParticleCreator& Creator(){return _creator;}

//     /**
//      * Simplify Kinematic Calculation Define calls
//      * see examples, Q2 etc below
//      * Asssume indices from ReactionMap
//      */
//     void Define(const string& name,const string& func){
//       auto types = _reaction->GetTypes();
//       for(const auto& type:types){
// 	_reaction->Define(type+name,utils::createFunctionCallStringFromVec("rad::util::ApplyKinematics",{func,names::ReactionMap(),(type+names::KineComponents())}));
//       }
//     }
    
//     /**
//      * @brief Defines a kinematic variable using a C++ Functor/Lambda as the kernel.
//      * Use this when you want to write logic inline or pass a complex functor
//      * rather than a string representation of a function.
//      * * @tparam Lambda The type of the functor (deduced automatically).
//      * @param name The name of the output column.
//      * @param func The functor to apply. 
//      * Signature: ResultType func(const RVecIndexMap& map, const RVecD& px, const RVecD& py, const RVecD& pz, const RVecD& m)
//      */

//     template <typename Lambda>
//     void DefineLambda(const string& name, Lambda&& func){
      
//       auto types = _reaction->GetTypes();
//       for(const auto& type:types){
// 	// Define the input columns required by the wrapper
// 	// We need the static map and the specific component vector for this type (rec_, tru_, etc.)
// 	ROOT::RDF::ColumnNames_t cols = {
// 	  names::ReactionMap(), 
// 	  type + names::KineComponents()
// 	};
// 	//create ApplyKinematics wrapper for func
// 	auto apply_func = [func](const RVecIndexMap& map, const ROOT::RVec<ROOT::RVec<RVecResultType>>& comps){
// 	  return rad::util::ApplyKinematics(func, map, comps);
// 	};
// 	//Define via Lambda 
// 	_reaction->Define(type+name,apply_func,cols);
//       }
//     }
    
//     /**
//      * Simplify Kinematic Calculation Define calls
//      * see examples, Q2 etc below
//      * Supply relevent indices in particles
//      */
//     // void Define(const string& name,const string& func, const ParticleNames_t& particles){
//     void Define(const string& name,const string& func, const std::vector<ParticleNames_t>& particles){
//       auto types = _reaction->GetTypes();
//       cout<<"ParticleCreator::Define " <<name<<endl;
//       //change {{"p1"},{"p2"}} to "{{kine_p1}, {kine_p2}}"
//       string kine_parts ="{";
//       for(const auto& pnames:particles){
// 	kine_parts+=utils::combineVectorToString(utils::prependToAll(pnames,names::data_type::Kine()));
// 	kine_parts+=",";
//       }
//       kine_parts.pop_back();
//       kine_parts+="}";
//       cout<<"ParticleCreator::Define "<<name<<" "<<kine_parts<<endl;
//       for(const auto& type:types){
// 	_reaction->Define(type+name,utils::createFunctionCallStringFromVec("rad::util::ApplyKinematics",{func,kine_parts,(type+names::KineComponents())}));
//       }
//     }

//     //User Shortcuts
//     //Note template types are all RVecResultType
//     //for ReactionMap indices
//     void Q2(){
//       KineCalculation Q2(*this,"Q2",rad::ElS_Q2);
//     }
    
//     void CosThetaCM(){
//       Define("CosThetaCM","rad::ElS_CosThetaCM< rad::RVecResultType, rad::RVecResultType>");
//     }
//     void PhiCM(){
//       Define("PhiCM","rad::ElS_PhiCM< rad::RVecResultType, rad::RVecResultType>");
//     }

//     //For given set of indices
//     void Mass(const string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={}){
//       //Note can take 2 sets of indices, 1 for +ve 1 for -ve 4-momenta
//       Define(name,"rad::FourVectorMassCalc<ROOT::RVecD,ROOT::RVecD>",{particles_pos,particles_neg});
//     }
//     void Mass2(const string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={}){
//       //Note can take 2 sets of indices, 1 for +ve 1 for -ve 4-momenta
//       Define(name,"rad::FourVectorMass2Calc<ROOT::RVecD,ROOT::RVecD>",{particles_pos,particles_neg});
//     }
//     void Pt(const string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={}){
//       //Note can take 2 sets of indices, 1 for +ve 1 for -ve 4-momenta
//       Define(name,"rad::FourVectorPtCalc<ROOT::RVecD,ROOT::RVecD>",{particles_pos,particles_neg});
//     }

//   private :
//     RVecIndexMap _mapIndices;
//     std::vector<string> _reactionNames;
//     config::ConfigReaction* _reaction = nullptr;
//     config::ParticleCreator _creator;
    
//   };
  
// }
