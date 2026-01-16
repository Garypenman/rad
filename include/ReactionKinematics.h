#pragma once

#include "BasicKinematics.h"
#include "ConfigReaction.h"
#include "DefineNames.h" 
#include "CommonDefines.h"

namespace rad {
  namespace physics {

    /**
     * @brief Reconstructs the total initial state 4-vector (Beam + Target).
     * Uses the specific beam indices provided.
     */
    inline PxPyPzMVector InitialStateFourVector(const RVecIndexMap& react, 
						const RVecResultType& px, const RVecResultType& py, 
						const RVecResultType& pz, const RVecResultType& m) 
    {
      // Sums all particles in the "Beams" group
      return FourVector(react[consts::OrderBeams()], px, py, pz, m);
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
      Indices_t icm = react[consts::OrderBaryons()];
      const auto& mesons = react[consts::OrderMesons()];
      icm.insert(icm.end(), mesons.begin(), mesons.end());
      return FourVector(icm, px, py, pz, m);
    }

    inline PxPyPzMVector BaryonFourVector(const RVecIndexMap& react, 
					  const RVecResultType& px, const RVecResultType& py, 
					  const RVecResultType& pz, const RVecResultType& m) 
    {
      return FourVector(react[consts::OrderBaryons()], px, py, pz, m);
    }

    inline PxPyPzMVector MesonFourVector(const RVecIndexMap& react, 
					 const RVecResultType& px, const RVecResultType& py, 
					 const RVecResultType& pz, const RVecResultType& m) 
    {
      return FourVector(react[consts::OrderMesons()], px, py, pz, m);
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
      auto p_target = FourVector(react[consts::OrderBeams()][consts::OrderBeamIon()], px, py, pz, m);
    
      // Subtract Baryon(s)
      SubtractFourVector(p_target, react[consts::OrderBaryons()], px, py, pz, m);
    
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
      auto tar = FourVector(react[consts::OrderBeams()][consts::OrderBeamIon()], px, py, pz, m);
      // Baryon (Bottom Out)
      auto bar = FourVector(react[consts::OrderBaryons()], px, py, pz, m);
    
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

  } // namespace physics
} // namespace rad






