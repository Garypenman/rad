#pragma once

#include "ReactionKinematics.h"  // FourVector, boost, PxPyPzMVector, etc.
#include "BasicKinematics.h"     // If InitialFourVector / helpers live here
#include "ConfigReaction.h"      // RVecIndexMap, names:: indices
#include "CommonDefines.h"       // Ensures RVecResultType, ResultType_t, etc.
#include "TMath.h"

namespace rad {
  namespace gn2s0s0s12 {
    //namespace physics {

    // --- Structures ---

    /**
     * @brief Result structure for generic decay angles.
     */
    struct DecayAngles_t {
      double cosTheta = 0.; ///< Cosine of the polar decay angle.
      double theta = 0.;    ///< Polar decay angle.
      double phi = 0.;      ///< Azimuthal decay angle.
    };

    
    // --- Center-of-Mass (CM) Kinematics ---

    /**
     * @brief Calculates the initial CM four-momentum for photoproduction:
     *        P_CM = P_initial(bottom beam) + q (photon).
     *
     * @details
     * Assumes `PhotoFourVector(...)` is defined by the reaction configuration
     * (e.g. PhotoIonReaction.h or ElectroIonReaction.h).
     *
     * @param react The fixed reaction index map.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return PxPyPzMVector The CM four-momentum vector.
     */
    inline PxPyPzMVector PhotoCMVector(const RVecIndexMap& react,
                                       const RVecResultType& px, const RVecResultType& py,
                                       const RVecResultType& pz, const RVecResultType& m)
    {
      // Assumes OrderBeams() [0] = ion, [1] = electron
      // Assuming ion is at pos 0 in the beam group
      auto ion = FourVector(react[consts::OrderBeams()][consts::OrderBeamIon()], px, py, pz, m); 
      auto phot = PhotoFourVector(react, px, py, pz, m);
      
      return ion + phot;
    }

    /**
     * @brief Calculates decay angles in the overall CM frame of the initial state.
     *
     * @param react The fixed reaction index map.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return DecayAngles_t {cosTheta, phi} of the meson system in the CM frame.
     */
    inline DecayAngles_t PhotoCMDecay(const RVecIndexMap& react,
                                      const RVecResultType& px, const RVecResultType& py,
                                      const RVecResultType& pz, const RVecResultType& m)
    {
      const auto cm = PhotoCMVector(react, px, py, pz, m);
      const auto cmBoost = cm.BoostToCM();

      auto mes = FourVector(react[consts::OrderMesons()], px, py, pz, m); // Assumes single particle or combined meson vector
      const PxPyPzMVector cmMes = boost(mes, cmBoost);

      DecayAngles_t result;
      result.cosTheta = TMath::Cos(cmMes.Theta());
      result.phi = cmMes.Phi();
      return result;
    }

    // --- Helicity Frame Decay Angles ---

    /**
     * @brief Calculates helicity decay angles in the meson rest frame.
     *
     * @details
     * Helicity frame convention:
     *   - z-axis along **-baryon** direction in the meson rest frame.
     *   - y-axis along (baryon × photon).
     *   - x-axis completes the right-handed set: x = y × z.
     *
     * @param react The fixed reaction index map.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return DecayAngles_t {cosTheta, phi} for the first meson child (index 0).
     */
    inline DecayAngles_t PhotoHelicityDecay(const RVecIndexMap& react,
                                            const RVecResultType& px, const RVecResultType& py,
                                            const RVecResultType& pz, const RVecResultType& m)
    {
      const auto baryon = FourVector(react[consts::OrderBaryons()], px, py, pz, m);
      //const auto meson = MesonFourVector(react, px, py, pz, m);
      const auto meson = FourVector(react[consts::OrderMesons()], px, py, pz, m);
      const auto photon = PhotoFourVector(react, px, py, pz, m);
      
      const auto decBoost = meson.BoostToCM();

      // Four-vectors in the meson rest frame
      const auto decBar   = boost(baryon, decBoost);
      const auto decGamma = boost(photon, decBoost);

      // Helicity axes
      const XYZVector zV = -decBar.Vect().Unit();
      const XYZVector yV = decBar.Vect().Cross(decGamma.Vect()).Unit();
      const XYZVector xV = yV.Cross(zV).Unit();

      // First decay product of the meson
      const auto child1    = FourVector(react[consts::OrderMesons()][0], px, py, pz, m);
      const auto decChild1 = boost(child1, decBoost);

      // Projections and angles
      const XYZVector proj(decChild1.Vect().Dot(xV),
                           decChild1.Vect().Dot(yV),
                           decChild1.Vect().Dot(zV));

      DecayAngles_t result;
      result.cosTheta = TMath::Cos(proj.Theta());
      result.phi      = proj.Phi();
      return result;
    }

    /**
     * @brief Convenience: returns cos(theta) in the helicity frame.
     */
    inline ResultType_t CosThetaHel(const RVecIndexMap& react,
                                    const RVecResultType& px, const RVecResultType& py,
                                    const RVecResultType& pz, const RVecResultType& m)
    {
      return PhotoHelicityDecay(react, px, py, pz, m).cosTheta;
    }
    
    /**
     * @brief Convenience: returns theta in the helicity frame.
     */
    inline ResultType_t ThetaHel(const RVecIndexMap& react,
                               const RVecResultType& px, const RVecResultType& py,
                               const RVecResultType& pz, const RVecResultType& m)
    {
      return PhotoHelicityDecay(react, px, py, pz, m).theta;
    }

    /**
     * @brief Convenience: returns phi in the helicity frame.
     */
    inline ResultType_t PhiHel(const RVecIndexMap& react,
                               const RVecResultType& px, const RVecResultType& py,
                               const RVecResultType& pz, const RVecResultType& m)
    {
      return PhotoHelicityDecay(react, px, py, pz, m).phi;
    }
    
    // --- GottfriedJackson (GJ) Frame Decay Angles ---

    /**
     * @brief Calculates GottfriedJackson (GJ) decay angles in the meson rest frame.
     *
     * @details
     * GJ frame convention:
     *   - z-axis along **photon** direction in the meson rest frame.
     *   - y-axis along (baryon × photon).
     *   - x-axis completes the right-handed set: x = y × z.
     *
     * @param react The fixed reaction index map.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return DecayAngles_t {cosTheta, phi} for the first meson child (index 0).
     */
    inline DecayAngles_t PhotoGJDecay(const RVecIndexMap& react,
                                      const RVecResultType& px, const RVecResultType& py,
                                      const RVecResultType& pz, const RVecResultType& m)
    {
      const auto baryon = FourVector(react[consts::OrderBaryons()], px, py, pz, m);
      const auto meson = FourVector(react[consts::OrderMesons()], px, py, pz, m);
      
      const auto decBoost = meson.BoostToCM();

      // Four-vectors in the meson rest frame
      const auto decBar   = boost(baryon, decBoost);
      const auto decGamma = boost(PhotoFourVector(react, px, py, pz, m), decBoost);

      // GJ axes
      const XYZVector zV = decGamma.Vect().Unit();
      const XYZVector yV = decBar.Vect().Cross(decGamma.Vect()).Unit();
      const XYZVector xV = yV.Cross(zV).Unit();

      // First decay product of the meson
      const auto child1    = FourVector(react[consts::OrderMesons()][0], px, py, pz, m);
      const auto decChild1 = boost(child1, decBoost);

      // Projections and angles
      const XYZVector proj(decChild1.Vect().Dot(xV),
                           decChild1.Vect().Dot(yV),
                           decChild1.Vect().Dot(zV));

      DecayAngles_t result;
      result.cosTheta = TMath::Cos(proj.Theta());
      result.phi      = proj.Phi();
      return result;
    }

    /**
     * @brief Convenience: returns cos(theta) in the GJ frame.
     */
    inline ResultType_t CosThetaGJ(const RVecIndexMap& react,
                                   const RVecResultType& px, const RVecResultType& py,
                                   const RVecResultType& pz, const RVecResultType& m)
    {
      return PhotoGJDecay(react, px, py, pz, m).cosTheta;
    }

    /**
     * @brief Convenience: returns phi in the GJ frame.
     */
    inline ResultType_t PhiGJ(const RVecIndexMap& react,
                              const RVecResultType& px, const RVecResultType& py,
                              const RVecResultType& pz, const RVecResultType& m)
    {
      return PhotoGJDecay(react, px, py, pz, m).phi;
    }

  } // namespace gn2s0s0s12
} // namespace rad
