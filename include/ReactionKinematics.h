#pragma once

#include "Beams.h"
#include "BasicKinematics.h"
#include "ConfigReaction.h"
#include "DefineNames.h" 
#include "CommonDefines.h"

namespace rad {

  // Forward declaration so TTop can compile before ElectronScatterKinematics.h is fully read
  // inline PxPyPzMVector PhotoFourVector(const RVecIndexMap& react, 
  //                                      const RVecResultType& px, const RVecResultType& py, 
  //                                      const RVecResultType& pz, const RVecResultType& m);

  // --- Helper: Initial State ---

  /**
   * @brief Reconstructs the total initial state 4-vector (Beam + Target).
   * Uses the specific beam indices provided.
   */
  inline PxPyPzMVector InitialStateFourVector(const RVecIndexMap& react, 
                                              const RVecResultType& px, const RVecResultType& py, 
                                              const RVecResultType& pz, const RVecResultType& m) 
  {
    // Sums all particles in the "Beams" group
    return FourVector(react[names::OrderBeams()], px, py, pz, m);
  }

  // --- Missing Kinematics (Initial - Subtracted) ---

  /**
   * @brief Calculates Missing Transverse Momentum (Pt).
   */
  inline double FourVectorMissPtCalc(const RVecIndexMap& react, const Indices_t& ineg,
                                     const RVecResultType& px, const RVecResultType& py, 
                                     const RVecResultType& pz, const RVecResultType& m) 
  {
    auto psum = InitialStateFourVector(react, px, py, pz, m);
    SubtractFourVector(psum, ineg, px, py, pz, m);
    return psum.Pt();
  }

  inline double FourVectorMissPzCalc(const RVecIndexMap& react, const Indices_t& ineg,
                                     const RVecResultType& px, const RVecResultType& py, 
                                     const RVecResultType& pz, const RVecResultType& m) 
  {
    auto psum = InitialStateFourVector(react, px, py, pz, m);
    SubtractFourVector(psum, ineg, px, py, pz, m);
    return psum.Pz();
  }

  inline double FourVectorMissECalc(const RVecIndexMap& react, const Indices_t& ineg,
                                    const RVecResultType& px, const RVecResultType& py, 
                                    const RVecResultType& pz, const RVecResultType& m) 
  {
    auto psum = InitialStateFourVector(react, px, py, pz, m);
    SubtractFourVector(psum, ineg, px, py, pz, m);
    return psum.E();
  }

  inline double FourVectorMissMass2Calc(const RVecIndexMap& react, const Indices_t& ineg,
                                        const RVecResultType& px, const RVecResultType& py, 
                                        const RVecResultType& pz, const RVecResultType& m) 
  {
    auto psum = InitialStateFourVector(react, px, py, pz, m);
    SubtractFourVector(psum, ineg, px, py, pz, m);
    return psum.M2();
  }

  inline double FourVectorMissMassCalc(const RVecIndexMap& react, const Indices_t& ineg,
                                       const RVecResultType& px, const RVecResultType& py, 
                                       const RVecResultType& pz, const RVecResultType& m) 
  {
    auto psum = InitialStateFourVector(react, px, py, pz, m);
    SubtractFourVector(psum, ineg, px, py, pz, m);
    return psum.M();
  }


  // --- Final State Reconstruction ---

  inline PxPyPzMVector CMVectorFinal(const RVecIndexMap& react, 
                                     const RVecResultType& px, const RVecResultType& py, 
                                     const RVecResultType& pz, const RVecResultType& m) 
  {
    Indices_t icm = react[names::OrderBaryons()];
    const auto& mesons = react[names::OrderMesons()];
    icm.insert(icm.end(), mesons.begin(), mesons.end());
    return FourVector(icm, px, py, pz, m);
  }

  inline PxPyPzMVector BaryonFourVector(const RVecIndexMap& react, 
                                        const RVecResultType& px, const RVecResultType& py, 
                                        const RVecResultType& pz, const RVecResultType& m) 
  {
    return FourVector(react[names::OrderBaryons()], px, py, pz, m);
  }

  inline PxPyPzMVector MesonFourVector(const RVecIndexMap& react, 
                                       const RVecResultType& px, const RVecResultType& py, 
                                       const RVecResultType& pz, const RVecResultType& m) 
  {
    return FourVector(react[names::OrderMesons()], px, py, pz, m);
  }


  // --- Mandelstam Variables ---

  /**
   * @brief Calculates t from the Bottom Vertex (Target - Baryon).
   */
  inline double TBot(const RVecIndexMap& react, 
                     const RVecResultType& px, const RVecResultType& py, 
                     const RVecResultType& pz, const RVecResultType& m) 
  {
    // Uses specific BeamIon index as requested
    auto p_target = FourVector(react[names::OrderBeams()][names::OrderBeamIon()], px, py, pz, m);
    
    // Subtract Baryon(s)
    SubtractFourVector(p_target, react[names::OrderBaryons()], px, py, pz, m);
    
    // t = M^2 of the transfer vector
    return -p_target.M2(); 
  }

  /**
   * @brief Calculates t from the Top Vertex (Photon - Meson).
   */
  inline double TTop(const RVecIndexMap& react, 
                     const RVecResultType& px, const RVecResultType& py, 
                     const RVecResultType& pz, const RVecResultType& m) 
  {
    // Get photon 4-vector (Forward declared)
    auto phot = PhotoFourVector(react, px, py, pz, m);
    
    // Get meson 4-vector
    auto meso = MesonFourVector(react, px, py, pz, m);
    
    // t = (gamma - meson)^2
    auto psum = phot - meso;
    return -psum.M2();
  }

  /**
   * @brief Calculates t0 (Minimum t) in the Center of Mass frame.
   */
  inline double T0(const RVecIndexMap& react, 
                   const RVecResultType& px, const RVecResultType& py, 
                   const RVecResultType& pz, const RVecResultType& m) 
  {
    // Target (Bottom)
    auto tar = FourVector(react[names::OrderBeams()][names::OrderBeamIon()], px, py, pz, m);
    // Baryon (Bottom Out)
    auto bar = FourVector(react[names::OrderBaryons()], px, py, pz, m);
    
    // Generate Final State CM (Mesons + Baryons)
    auto cm = CMVectorFinal(react, px, py, pz, m);
    auto cmBoost = cm.BoostToCM();
    
    auto CMTar = boost(tar, cmBoost);
    auto CMBar = boost(bar, cmBoost);
    
    // t0 formula assuming collinear emission
    double t0 = CMBar.M2() + CMTar.M2() - 2 * (CMBar.E() * CMTar.E() - CMBar.P() * CMTar.P());
    return t0;
  }

  /**
   * @brief Calculates t' = t - t0 (Bottom Vertex).
   */
  inline double TPrimeBot(const RVecIndexMap& react, 
                          const RVecResultType& px, const RVecResultType& py, 
                          const RVecResultType& pz, const RVecResultType& m) 
  {
    return TBot(react, px, py, pz, m) - T0(react, px, py, pz, m); 
  }

  /**
   * @brief Calculates t' = t - t0 (Top Vertex).
   */
  inline double TPrimeTop(const RVecIndexMap& react, 
                          const RVecResultType& px, const RVecResultType& py, 
                          const RVecResultType& pz, const RVecResultType& m) 
  {
    return TTop(react, px, py, pz, m) - T0(react, px, py, pz, m); 
  }

} // namespace rad







// #pragma once
// #include "Beams.h"
// #include "BasicKinematics.h"
// #include "ConfigReaction.h"

// void ReactionKinematics(){}

// namespace rad{

    
 

// //     ///\brief missing transverse momentum of reaction = top+bot -  neg[i]
// //     //const config::RVecIndexMap react must be copied for thread safety.
// //    template<typename Tp, typename Tm>
// //     double FourVectorMissPtCalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
// //     {
// //       auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
// //       psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
// //       SubtractFourVector(psum,ineg,px,py,pz,m);
// //       return psum.Pt();
// //     }

// //    ///\brief missing transverse momentum of reaction = top+bot -  neg[i]
// //     //const config::RVecIndexMap react must be copied for thread safety.
// //    template<typename Tp, typename Tm>
// //     double FourVectorMissPzCalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
// //     {
// //       auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
// //       psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
// //       SubtractFourVector(psum,ineg,px,py,pz,m);
// //       return psum.Pz();
// //     }

// //   ///\brief missing transverse momentum of reaction = top+bot -  neg[i]
// //     //const config::RVecIndexMap react must be copied for thread safety.
// //     template<typename Tp, typename Tm>
// //     double FourVectorMissThetaCalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
// //     {
// //       auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
// //       psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
// //       SubtractFourVector(psum,ineg,px,py,pz,m);
// //       return psum.Theta();
// //     }
    
// //     ///\brief missing transverse momentum of reaction = top+bot -  neg[i]
// //     //const config::RVecIndexMap react must be copied for thread safety.
// //    template<typename Tp, typename Tm>
// //     double FourVectorMissPhiCalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
// //     {
// //       auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
// //       psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
// //       SubtractFourVector(psum,ineg,px,py,pz,m);
// //       return psum.Phi();
// //     }


// //     ///\brief missing transverse momentum of reaction = top+bot -  neg[i]
// //     //const config::RVecIndexMap react must be copied for thread safety.
// //    template<typename Tp, typename Tm>
// //     double FourVectorMissPCalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
// //     {
// //       auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
// //       psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
// //       SubtractFourVector(psum,ineg,px,py,pz,m);
// //       return psum.P();
// //     }

// //     ///\brief missing transverse momentum of reaction = top+bot -  neg[i]
// //     //const config::RVecIndexMap react must be copied for thread safety.
// //    template<typename Tp, typename Tm>
// //     double FourVectorMissECalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
// //     {
// //       auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
// //       psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
// //       SubtractFourVector(psum,ineg,px,py,pz,m);
// //       return psum.E();
// //     }

//   ///\brief functions to compute standard reaction kinematics
//   template<typename Tp, typename Tm>
//   PxPyPzMVector CMVectorFinal(const RVecIndexMap& react, const Tp &px, const Tp &py, const Tp &pz, const Tm &m){
    
//     Indices_t icm=react[names::OrderBaryons()];
//     icm.insert(icm.end(),react[names::OrderMesons()].begin(),react[names::OrderMesons()].end());
    
//     return FourVector(icm,px,py,pz,m);
//   }

 

// //     /**
// //      * reaction baryon 4-vector
// //      */
// //     template<typename Tp, typename Tm>
// //     PxPyPzMVector BaryonFourVector(const config::RVecIndexMap& react, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
// //       return FourVector(react[names::BaryonsIdx()],px,py,pz,m);
// //     }
// //     /**
// //      * reaction meson 4-vector
// //      */
// //     template<typename Tp, typename Tm>
// //     PxPyPzMVector MesonFourVector(const config::RVecIndexMap& react, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
// //       return FourVector(react[names::MesonsIdx()],px,py,pz,m);
// //     }

   
// //     template<typename Tp, typename Tm>
// //     double T0(const config::RVecIndexMap& react,
// // 	  const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){

// //       auto tar = beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
// //       auto bar = FourVector(react[names::BaryonsIdx()],px,py,pz,m);
// //       //generate CM from sum of final state meson and baryon particles
// //       auto cm = CMVector(react,px,py,pz,m);
// //       auto cmBoost = cm.BoostToCM();
// //       PxPyPzMVector  CMTar=boost(tar,cmBoost);
// //       PxPyPzMVector  CMBar=boost(bar,cmBoost);
    
// //       //return  M1*M1 + M3*M3  - 2 * ( E1*E3 -p1*p3*costh );
// //       Double_t t0 = CMBar.M2() + CMTar.M2() - 2*(CMBar.E()*CMTar.E() - CMBar.P()*CMTar.P() ) ;
// //       /* cout << "Inside T0 Func" << endl; */
// //       /* cout << "CMBar: " << CMBar << endl; */
// //       /* cout << "CMTar: " << CMTar << endl; */
// //       return t0;
// //     }
  
   
// //     ///\brief return 4 momentum transfer squared of "in particles" - "out particles" on bottom vertex
// //     template<typename Tp, typename Tm>
// //     double TBot(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
// //       auto psum = beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
// //       /* cout << "Inside TBot Func" << endl; */
// //       /* cout << "pbeam: " << psum << endl; */
// //       SubtractFourVector(psum,react[names::BaryonsIdx()],px,py,pz,m); 
// //       /* cout << "pbeam_sub_pprime: " << psum << endl; */

// //       return - (psum.M2());
// //     }
// //     ///\brief return 4 momentum transfer squared of "in particles" - "out particles" on top vertex
// //     template<typename Tp, typename Tm>
// //     double TTop(const config::RVecIndexMap& react,
// // 	    const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
// //       //Get photon 4-vector
// //       auto phot=PhotoFourVector(react,px,py,pz,m);
// //       //Get meson (Top) 4-vector
// //       auto meso=FourVector(react[names::MesonsIdx()],px,py,pz,m);
// //       //subtract
// //       /* cout << "Inside TTop Func" << endl; */
// //       /* cout << "phot: " << phot << endl; */
// //       /* cout << "meso: " << meso << endl; */
// //       auto psum = phot-meso;
// //       //return t
// //       return - (psum.M2());
// //     }

// //       //For some reason this calculation only works when boosted into
// //       //proton beam rest frame......
// //       //auto pbeam = beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
// //       //auto boov=pbeam.BoostToCM();
// //       // auto phot=boost(PhotoFourVector(react,px,py,pz,m),boov);
// //       //auto meso=boost(FourVector(react[names::MesonsIdx()],px,py,pz,m),boov);

  
// //     ///\brief return 4 momentum transfer squared, t, minus t0 (or tmin) of "in particles" - "out particles"
// //     template<typename Tp, typename Tm>
// //     double TPrimeBot(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
// //       auto tbot = TBot(react,px,py,pz,m);
// //       auto t0 = T0(react,px,py,pz,m);
// //       /* cout << "TBot: " << tbot << endl; */
// //       /* cout << "T0: " << t0 << endl; */
// //       /* cout << "TpBot; " << tbot+t0 << endl; */
// //       /* cout << endl; */
// //       return tbot + t0;
// //     //return TBot(react,px,py,pz,m) + T0(react,px,py,pz,m);
// //     }
// //     ///\brief return 4 momentum transfer squared, t, minus t0 (or tmin) of "in particles" - "out particles"
// //     template<typename Tp, typename Tm>
// //     double TPrimeTop(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
// //       auto ttop = TTop(react,px,py,pz,m);
// //       auto t0 = T0(react,px,py,pz,m);
// //       /* cout << "TTop: " << ttop << endl; */
// //       /* cout << "T0: " << t0 << endl; */
// //       /* cout << "TpTot; " << ttop+t0 << endl; */
// //       /* cout << endl; */
// //       return ttop + t0;
// //       //return TTop(react,px,py,pz,m) + T0(react,px,py,pz,m);
// //     }
  


// //   }

// }
