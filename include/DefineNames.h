#pragma once


//!  Helper functions for names of reaction parts

namespace rad{
  namespace names{
    ///Comments not using string_view as need to return a string
    // returning static const seems to be slower, so stick with
    // returning string which should actually be slower 
    //constexpr const std::string_view  ScatEle() {return "scat_ele";}
    //Lilely due to some small string optimisations
    const std::string  ScatEle() {return "scat_ele";}
    const std::string  VirtGamma() {return "virt_gam";}
    const std::string  Dependencies() {return "dependencies";}
 
    const std::string  BeamEle() {return "beam_ele";}
    const std::string  BeamIon() {return "beam_ion";}
    const std::string    TargetIon() {return "tar_ion";}
    const std::string   BeamGamma() {return "beam_gamma";}
    const std::string  P4BeamEle() {return "p4beam_ele";}
    const std::string  P4BeamIon() {return "p4beam_ion";}
    const std::string  P4TargetIon() {return "p4tar_ion";}
    const std::string  P4BeamGamma() {return "p4beam_gamma";}
    
    /**
     * Names used to identify reaction components
     */
    const std::string    ReactionCombos()  {return "reaction_combos__dnwtag"; }//note DoNotWriteTag()
    const std::string  KineIndices()  {return "kine_indices__dnwtag"; }//note DoNotWriteTag()
    const std::string  KineComponents()  {return "kine_comps__dnwtag"; }//note DoNotWriteTag()
    const std::string    ReactionMap()  {return "reaction_map"; }
    const std::string  Mesons()  {return "meson__dnwtag"; }
    const std::string    Baryons() {return "baryon__dnwtag";}
    const std::string    Beams() {return "beams__dnwtag";}
    const std::string    Createds() {return "create_parts__dnwtag";}
    const std::string    Depends() {return "dependency_parts__dnwtag";}
    const std::string    ScatGroup() { return "scat_ele_group__dnwtag"; } 
    /**
     * @brief Defines the fixed order of particle roles in the consolidated momentum arrays.
     * * This order is essential for correctly building the static ReactionMap and accessing 
     * components in the KinematicsProcessor functor.
     */
 
    enum class ParticleGroupOrder{Beams, Baryons, Mesons, ScatEle, Deps, Createds};
  
    constexpr int OrderBeams() {return static_cast<int> (ParticleGroupOrder::Beams); }
    constexpr int OrderBeamIon() {return 0; }//First in beams
    constexpr int OrderBeamEle() {return 1; }//Second in beams
    
    constexpr int OrderBaryons() {return static_cast<int> (ParticleGroupOrder::Baryons); }
    constexpr int OrderMesons()  {return static_cast<int> (ParticleGroupOrder::Mesons); }
    constexpr int OrderScatEle() {return static_cast<int> (ParticleGroupOrder::ScatEle); }
    constexpr int OrderDependencies() {return static_cast<int> (ParticleGroupOrder::Deps); }
    constexpr int OrderCreated() {return static_cast<int> (ParticleGroupOrder::Createds); }
    constexpr int OrderVirtGamma() {return 0; }//First in Createds

    /**
     * @brief Returns the standard string name corresponding to a fixed particle index.
     * @param idx The integer value of the ParticleGroupOrder enum (e.g., OrderBeamIon()).
     * @return const std::string The human-readable particle name.
     * @throws std::runtime_error if the provided index is unknown.
     */
    std::string ParticleGroupFromIdx(const int idx) {
    
      // Cast the integer index back to the enum class for clear switch handling.
      switch (static_cast<ParticleGroupOrder>(idx)) {
        
      case ParticleGroupOrder::Beams:
	return Beams();
                       
      case ParticleGroupOrder::Baryons:
	return Baryons();
            
      case ParticleGroupOrder::Mesons:
	return Mesons();
            
      case ParticleGroupOrder::ScatEle:
	return ScatEle();
            
      case ParticleGroupOrder::Deps:
	return Depends();
            
      case ParticleGroupOrder::Createds:
	return Createds();
            
      default:
	// Ensures safety against indices outside the defined enum range.
	throw std::runtime_error("ParticleNameFromIdx: Unknown ParticleGroupOrder index: " + std::to_string(idx));
      }
    }
    
    /**
     * Names for common groups of columns 
     */
    const std::string  P4Components() {return "components_p4";}
    const std::string  P3Components() {return "components_p3";}
   

    /**
     * Links to reaction component links for users c++ functions
     */
    // enum class InitGroup{Bot,Top}; //ordering below must match this
    // enum class ElectroGroup{ BeamIon,BeamEle,ScatEle=4,VirtGam}; //4=>order in particleMap (after baryons and mesons)
    // enum class PhotoGroup{ TarIon,BeamGam}; //ordering below must match this
    // enum class FinalGroup{ Baryons=2,Mesons}; //ordering below must match this

    
    // constexpr uint  InitialTopIdx() {return static_cast<uint>(InitGroup::Top);}
    // constexpr uint  InitialBotIdx() {return static_cast<uint>(InitGroup::Bot);}

    // constexpr uint  ElectroIonIdx() {return InitialBotIdx();}
    // constexpr uint  ElectroEleIdx() {return InitialTopIdx();}
    
    // constexpr uint  ScatEleIdx() {return static_cast<uint>(ElectroGroup::ScatEle);}
    // constexpr uint  VirtGammaIdx() {return static_cast<uint>(ElectroGroup::VirtGam);}
    
    // constexpr uint  BaryonsIdx() {return static_cast<uint>(FinalGroup::Baryons);}
    // constexpr uint  MesonsIdx() {return static_cast<uint>(FinalGroup::Mesons);}
    
    // constexpr uint  PhotoGammaIdx() {return InitialTopIdx();}
    // constexpr uint  PhotoIonIdx() {return InitialBotIdx();}

    /**
     * Types of data
     */
    namespace data_type{
      const std::string  Rec() {return "rec_";}
      const std::string  Truth() {return "tru_";}
      const std::string  MC() {return "mc_";}
      const std::string  Kine() {return "kine_";}
    }
    
  }//names
}//rad
