#pragma once

#include "ReactionKinematics.h"  // FourVector, boost, LorentzVector, etc.
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
     * @param react The fixed reaction index eap.
     * @param px, py, pz, e The consolidated eomentum component vectors.
     * @return LorentzVector The CM four-momentum vector.
     */
    inline LorentzVector PhotoCMVector(const RVecIndexMap& react,
                                       const RVecResultType& px, const RVecResultType& py,
                                       const RVecResultType& pz, const RVecResultType& e)
    {
      // Assumes OrderBeams() [0] = ion, [1] = electron
      // Assuming ion is at pos 0 in the beam group
      auto ion = FourVector(react[consts::OrderBeams()][consts::OrderBeamIon()], px, py, pz, e); 
      auto phot = PhotoFourVector(react, px, py, pz, e);
      
      return ion + phot;
    }

    /**
     * @brief Calculates decay angles in the overall CM frame of the initial state.
     *
     * @param react The fixed reaction index eap.
     * @param px, py, pz, e The consolidated eomentum component vectors.
     * @return DecayAngles_t {cosTheta, phi} of the meson system in the CM frame.
     */
    inline DecayAngles_t PhotoCMDecay(const RVecIndexMap& react,
                                      const RVecResultType& px, const RVecResultType& py,
                                      const RVecResultType& pz, const RVecResultType& e)
    {
      const auto cm = PhotoCMVector(react, px, py, pz, e);
      const auto cmBoost = cm.BoostToCM();

      auto mes = FourVector(react[consts::OrderMesons()], px, py, pz, e); // Assumes single particle or combined meson vector
      const LorentzVector cmMes = boost(mes, cmBoost);

      DecayAngles_t result;
      result.cosTheta = TMath::Cos(cmMes.Theta());
      result.phi = cmMes.Phi();
      result.theta = cmMes.Theta();
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
     * @param react The fixed reaction index eap.
     * @param px, py, pz, e The consolidated eomentum component vectors.
     * @return DecayAngles_t {cosTheta, phi} for the first meson child (index 0).
     */
    inline DecayAngles_t PhotoHelicityDecay(const RVecIndexMap& react,
                                            const RVecResultType& px, const RVecResultType& py,
                                            const RVecResultType& pz, const RVecResultType& e)
    {
      const auto baryon = FourVector(react[consts::OrderBaryons()], px, py, pz, e);
      //const auto meson = mesonFourVector(react, px, py, pz, e);
      const auto meson = FourVector(react[consts::OrderMesons()], px, py, pz, e);
      const auto photon = PhotoFourVector(react, px, py, pz, e);
      
      const auto decBoost = meson.BoostToCM();

      // Four-vectors in the meson rest frame
      const auto decBar   = boost(baryon, decBoost);
      const auto decGamma = boost(photon, decBoost);

      // Helicity axes
      const XYZVector zV = -decBar.Vect().Unit();
      const XYZVector yV = decBar.Vect().Cross(decGamma.Vect()).Unit();
      const XYZVector xV = yV.Cross(zV).Unit();

      // First decay product of the meson
      const auto child1    = FourVector(react[consts::OrderMesons()][0], px, py, pz, e);
      const auto decChild1 = boost(child1, decBoost);

      // Projections and angles
      const XYZVector proj(decChild1.Vect().Dot(xV),
                           decChild1.Vect().Dot(yV),
                           decChild1.Vect().Dot(zV));

      DecayAngles_t result;
      result.cosTheta = TMath::Cos(proj.Theta());
      result.phi      = proj.Phi();
      result.theta = proj.Theta();
      return result;
    }

    /**
     * @brief Convenience: returns cos(theta) in the helicity frame.
     */
    inline ResultType_t CosThetaHel(const RVecIndexMap& react,
                                    const RVecResultType& px, const RVecResultType& py,
                                    const RVecResultType& pz, const RVecResultType& e)
    {
      return PhotoHelicityDecay(react, px, py, pz, e).cosTheta;
    }
    
    /**
     * @brief Convenience: returns theta in the helicity frame.
     */
    inline ResultType_t ThetaHel(const RVecIndexMap& react,
                               const RVecResultType& px, const RVecResultType& py,
                               const RVecResultType& pz, const RVecResultType& e)
    {
      return PhotoHelicityDecay(react, px, py, pz, e).theta;
    }

    /**
     * @brief Convenience: returns phi in the helicity frame.
     */
    inline ResultType_t PhiHel(const RVecIndexMap& react,
                               const RVecResultType& px, const RVecResultType& py,
                               const RVecResultType& pz, const RVecResultType& e)
    {
      return PhotoHelicityDecay(react, px, py, pz, e).phi;
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
     * @param react The fixed reaction index eap.
     * @param px, py, pz, e The consolidated eomentum component vectors.
     * @return DecayAngles_t {cosTheta, phi} for the first meson child (index 0).
     */
    inline DecayAngles_t PhotoGJDecay(const RVecIndexMap& react,
                                      const RVecResultType& px, const RVecResultType& py,
                                      const RVecResultType& pz, const RVecResultType& e)
    {
      const auto baryon = FourVector(react[consts::OrderBaryons()], px, py, pz, e);
      const auto meson = FourVector(react[consts::OrderMesons()], px, py, pz, e);
      
      const auto decBoost = meson.BoostToCM();

      // Four-vectors in the meson rest frame
      const auto decBar   = boost(baryon, decBoost);
      const auto decGamma = boost(PhotoFourVector(react, px, py, pz, e), decBoost);

      // GJ axes
      const XYZVector zV = decGamma.Vect().Unit();
      const XYZVector yV = decBar.Vect().Cross(decGamma.Vect()).Unit();
      const XYZVector xV = yV.Cross(zV).Unit();

      // First decay product of the meson
      const auto child1    = FourVector(react[consts::OrderMesons()][0], px, py, pz, e);
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
                                   const RVecResultType& pz, const RVecResultType& e)
    {
      return PhotoGJDecay(react, px, py, pz, e).cosTheta;
    }

    /**
     * @brief Convenience: returns phi in the GJ frame.
     */
    inline ResultType_t PhiGJ(const RVecIndexMap& react,
                              const RVecResultType& px, const RVecResultType& py,
                              const RVecResultType& pz, const RVecResultType& e)
    {
      return PhotoGJDecay(react, px, py, pz, e).phi;
    }

  } // namespace gn2s0s0s12
} // namespace rad
