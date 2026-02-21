#pragma once

#include "ReactionKinematics.h"
#include "CommonDefines.h" // Ensures RVecResultType is defined
#include "TMath.h"
#include <algorithm> // For std::clamp

namespace rad {
  namespace physics {

    // --- Structures ---

    /**
     * @brief Result structure for scattering kinematics in the proton/ion rest frame.
     */
    struct ScatterInProtonRestResult {
      double theta = 0;  ///< Scattering angle in the proton rest frame (radians).
      double energy = 0; ///< Virtual photon energy in the proton rest frame (GeV).
      double Q2 = 0;     ///< Momentum transfer squared (Q^2 > 0).
    };

    /**
     * @brief Result structure for decay angles.
     */
    struct DecayAngles_t {
      double cosTheta = 0.; ///< Cosine of the decay angle (polar angle).
      double phi = 0.;      ///< Azimuthal angle.
    };

 

    /**
     * @brief Calculates the squared four-momentum transfer, Q^2 = -q^2.
     * @param react The fixed ReactionMap.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return ResultType_t The Q^2 value (always positive).
     */
    inline ResultType_t ElS_Q2(const RVecIndexMap& react, 
			       const RVecResultType& px, const RVecResultType& py, 
			       const RVecResultType& pz, const RVecResultType& m) 
    {
      // Q^2 is defined as the negative invariant mass squared of the virtual photon.
      auto phot = PhotoFourVector(react, px, py, pz, m);
      return -phot.M2();
    }
    
    /**
     * @brief Calculates the struck quark momentum fraction, xbj = Q^2/(2p.q).
     * @param react The fixed ReactionMap.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return ResultType_t The xbj value (always positive).
     */
    inline ResultType_t ElS_xbj(const RVecIndexMap& react, 
			       const RVecResultType& px, const RVecResultType& py, 
			       const RVecResultType& pz, const RVecResultType& m) 
    {
      auto ion = FourVector(react[consts::OrderBeams()][consts::OrderBeamIon()], px, py, pz, m); 
      auto phot = PhotoFourVector(react, px, py, pz, m);
      return -phot.M2() / (2*ion.Dot(phot));
    }
    
    /**
     * @brief Calculates the inelasticity parameter, y = p.q/p.k
     * @param react The fixed ReactionMap.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return ResultType_t The y value (always positive).
     */
    inline ResultType_t ElS_y(const RVecIndexMap& react, 
			       const RVecResultType& px, const RVecResultType& py, 
			       const RVecResultType& pz, const RVecResultType& m) 
    {
      auto ion = FourVector(react[consts::OrderBeams()][consts::OrderBeamIon()], px, py, pz, m); 
      auto ebeam = FourVector(react[consts::OrderBeams()][consts::OrderBeamEle()], px, py, pz, m);
      auto phot = PhotoFourVector(react, px, py, pz, m);
      return ion.Dot(phot) / ion.Dot(ebeam);
    }
    
    /**
     * @brief Calculates the virtual photon energy (ion rest frame), nu = p.q / M
     * @param react The fixed ReactionMap.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return ResultType_t The nu value (always positive).
     */
    inline ResultType_t ElS_nu(const RVecIndexMap& react, 
			       const RVecResultType& px, const RVecResultType& py, 
			       const RVecResultType& pz, const RVecResultType& m) 
    {
      auto ion = FourVector(react[consts::OrderBeams()][consts::OrderBeamIon()], px, py, pz, m); 
      auto phot = PhotoFourVector(react, px, py, pz, m);
      return ion.Dot(phot) / ion.M();
    }
    
    /**
     * @brief Calculates the electroproduction variable, tau  = Q^2/(4M^2).
     * @param react The fixed ReactionMap.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return ResultType_t The tau value (always positive).
     */
    inline ResultType_t ElS_tau(const RVecIndexMap& react, 
			       const RVecResultType& px, const RVecResultType& py, 
			       const RVecResultType& pz, const RVecResultType& m) 
    {
      auto ion = FourVector(react[consts::OrderBeams()][consts::OrderBeamIon()], px, py, pz, m); 
      auto phot = PhotoFourVector(react, px, py, pz, m);
      //auto Q2 = -phot.M2(); //call ElS_Q2() here or just compute by hand?
      return -phot.M2() / (4*ion.M2());
    }
    
    /**
     * @brief Calculates the  variable, tau' = Q'^2/(2p.q).
     * @param react The fixed ReactionMap.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return ResultType_t The tau' value (always positive).
     */
    inline ResultType_t ElS_tauprime(const RVecIndexMap& react, 
			       const RVecResultType& px, const RVecResultType& py, 
			       const RVecResultType& pz, const RVecResultType& m) 
    {
      // Q^2 is defined as the negative invariant mass squared of the virtual photon.
      auto ion = FourVector(react[consts::OrderBeams()][consts::OrderBeamIon()], px, py, pz, m); 
      auto phot = PhotoFourVector(react, px, py, pz, m);
      auto mes = FourVector(react[consts::OrderMesons()], px, py, pz, m); // Assumes single particle or combined meson vector
      return mes.M2() / (2*ion.Dot(phot));
    }
    
    
    /**
     * @brief Calculates the four-momentum of the initial state Center-of-Mass (CM) system.
     * @param react The fixed ReactionMap.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return PxPyPzMVector The CM four-momentum vector (P_ion + P_virtual_photon).
     */
    inline PxPyPzMVector CMVectorInitial(const RVecIndexMap& react, 
					 const RVecResultType& px, const RVecResultType& py, 
					 const RVecResultType& pz, const RVecResultType& m) 
    {
      // Assumes OrderBeams() [0] = ion, [1] = electron
      // Assuming ion is at pos 0 in the beam group
      auto ion = FourVector(react[consts::OrderBeams()][consts::OrderBeamIon()], px, py, pz, m); 
      auto phot = PhotoFourVector(react, px, py, pz, m);

      return ion + phot;
    }

    // --- Proton Rest Frame Kinematics ---

    /**
     * @brief Calculates key scattering kinematics in the rest frame of the proton/ion target.
     * @param react The fixed ReactionMap.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return ScatterInProtonRestResult {theta, energy, Q2} in the rest frame.
     */
    inline ScatterInProtonRestResult ScatterInProtonRest(const RVecIndexMap& react,
							 const RVecResultType& px, const RVecResultType& py,
							 const RVecResultType& pz, const RVecResultType& m) 
    {
      // Assumes necessary indices are defined in react
      const auto pbeam = FourVector(react[consts::OrderBeams()][consts::OrderBeamIon()], px, py, pz, m);
      const auto ebeam = FourVector(react[consts::OrderBeams()][consts::OrderBeamEle()], px, py, pz, m);
      const auto scatele = FourVector(react[consts::OrderScatEle()][0], px, py, pz, m);

      // Boost frame: Target ion CM
      const auto pboost = pbeam.BoostToCM();

      // Boost momenta
      const auto prbeam = boost(ebeam, pboost);
      const auto prscat = boost(scatele, pboost);
      const auto prgamstar = prbeam - prscat; // Virtual photon in proton rest frame

      // Calculate scattering angle (theta)
      const double num = prbeam.Vect().Dot(prscat.Vect());
      const double denom = prbeam.Vect().R() * prscat.Vect().R();

      double theta = consts::InvalidEntry<double>();

      if (denom != 0.0) {
	double ratio = num / denom;
	ratio = std::clamp(ratio, -1.0, 1.0); // Ensure ratio is [-1, 1] for acos
	theta = TMath::ACos(ratio);
      }

      const double energy = prgamstar.E();
      const double Q2 = -prgamstar.M2();

      return {theta, energy, Q2};
    }
    
    /**
     * @brief Calculates the polarization parameter (epsilon) for the virtual photon.
     * @param react The fixed ReactionMap.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return ResultType_t The polarization parameter (epsilon).
     */
    inline ResultType_t ElS_PolVirtPhot(const RVecIndexMap& react, 
					const RVecResultType& px, const RVecResultType& py, 
					const RVecResultType& pz, const RVecResultType& m) 
    {
      // Calculate necessary variables in the proton rest frame
      const auto prvec = ScatterInProtonRest(react, px, py, pz, m);

      // Renaming for clarity in calculation
      const auto ElScatTh = prvec.theta;
      const auto GammaE = prvec.energy;
      const auto Q2 = prvec.Q2;

      // Formula for virtual photon polarization (epsilon)
      const auto pol = 1. / (1. + 2. * (1. + GammaE * GammaE / Q2) * TMath::Tan(ElScatTh / 2.) * TMath::Tan(ElScatTh / 2.));

      return pol;
    }
    
    /**
     * @brief Calculates the circular polarization parameter for the virtual photon sqrt(1-epsilon^2).
     * @param react The fixed ReactionMap.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return ResultType_t The circular polarization parameter.
     */
    inline ResultType_t ElS_CircPolVirtPhot(const RVecIndexMap& react, 
					const RVecResultType& px, const RVecResultType& py, 
					const RVecResultType& pz, const RVecResultType& m) 
    {
      auto eps = ElS_PolVirtPhot(react,px,py,pz,m);
      return sqrt(1-eps*eps);
    }

    // --- Decay Frame Kinematics (CM and Proton Rest) ---

    /**
     * @brief Calculates decay angles in the Center-of-Mass (CM) frame of the electron-virtual photon system.
     * This frame is defined using the electron beam and the virtual photon direction.
     * @param react The fixed ReactionMap.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return DecayAngles_t {CosTheta, Phi}.
     */
    inline DecayAngles_t ElectroCMDecay(const RVecIndexMap& react, 
					const RVecResultType& px, const RVecResultType& py, 
					const RVecResultType& pz, const RVecResultType& m) 
    {
      auto cm = CMVectorInitial(react, px, py, pz, m);
      auto cmBoost = cm.BoostToCM();

      // Get relevant particles
      auto beam = FourVector(react[consts::OrderBeams()][consts::OrderBeamEle()], px, py, pz, m);
      auto mes = FourVector(react[consts::OrderMesons()], px, py, pz, m); // Assumes single particle or combined meson vector
      auto photon = PhotoFourVector(react, px, py, pz, m);

      // Boost to CM frame
      PxPyPzMVector CMBeam = boost(beam, cmBoost);
      PxPyPzMVector CMMes = boost(mes, cmBoost);
      PxPyPzMVector CMGamma = boost(photon, cmBoost);

      // Define Helicity Coordinate System (z-axis along virtual photon)
      XYZVector zV = CMGamma.Vect().Unit();
      XYZVector yV = CMGamma.Vect().Cross(CMBeam.Vect()).Unit();
      XYZVector xV = yV.Cross(zV).Unit();

      // Project meson decay vector onto the coordinate axes
      XYZVector angles(CMMes.Vect().Dot(xV), CMMes.Vect().Dot(yV), CMMes.Vect().Dot(zV));

      DecayAngles_t result;
      // Cos(theta) is projection onto z-axis (magnitude of z-component)
      result.cosTheta = TMath::Cos(angles.Theta());
      result.phi = angles.Phi();

      return result;
    }

    /**
     * @brief Calculates the decay angle Cos(Theta) in the CM frame.
     * @param react, px, py, pz, m The inputs.
     * @return ResultType_t Cosine of the polar decay angle.
     */
    inline ResultType_t ElS_CosThetaCM(const RVecIndexMap& react, 
				       const RVecResultType& px, const RVecResultType& py, 
				       const RVecResultType& pz, const RVecResultType& m) 
    {
      return ElectroCMDecay(react, px, py, pz, m).cosTheta;
    }

    /**
     * @brief Calculates the decay angle Phi in the CM frame.
     * @param react, px, py, pz, m The inputs.
     * @return ResultType_t Azimuthal decay angle (Phi).
     */
    inline ResultType_t ElS_PhiCM(const RVecIndexMap& react, 
				  const RVecResultType& px, const RVecResultType& py, 
				  const RVecResultType& pz, const RVecResultType& m) 
    {
      return ElectroCMDecay(react, px, py, pz, m).phi;
    }

    /**
     * @brief Calculates decay angles in the rest frame of the produced baryon (Proton Rest Frame).
     * This frame is defined using the initial proton beam direction.
     * @param react The fixed ReactionMap.
     * @param px, py, pz, m The consolidated momentum component vectors.
     * @return DecayAngles_t {CosTheta, Phi}.
     */
    inline DecayAngles_t ElectroProtonRestDecay(const RVecIndexMap& react, 
						const RVecResultType& px, const RVecResultType& py, 
						const RVecResultType& pz, const RVecResultType& m) 
    {
      // Boost frame: Initial proton rest frame
      auto pr = FourVector(react[consts::OrderBeams()][consts::OrderBeamIon()], px, py, pz, m);
      auto prBoost = pr.BoostToCM();

      // Get relevant particles
      auto beam = FourVector(react[consts::OrderBeams()][consts::OrderBeamEle()], px, py, pz, m);
      auto mes = FourVector(react[consts::OrderMesons()], px, py, pz, m);
      auto photon = PhotoFourVector(react, px, py, pz, m);

      // Boost to Proton Rest Frame
      PxPyPzMVector prBeam = boost(beam, prBoost);
      PxPyPzMVector prMes = boost(mes, prBoost);
      PxPyPzMVector prGamma = boost(photon, prBoost);

      // Define Helicity Coordinate System (z-axis opposite virtual photon, for convention)
      XYZVector zV = -prGamma.Vect().Unit(); // Often defined opposite q*
      XYZVector yV = prGamma.Vect().Cross(prBeam.Vect()).Unit();
      XYZVector xV = yV.Cross(zV).Unit();

      // Project meson decay vector onto the coordinate axes
      XYZVector angles(prMes.Vect().Dot(xV), prMes.Vect().Dot(yV), prMes.Vect().Dot(zV));

      DecayAngles_t result;
      result.cosTheta = TMath::Cos(angles.Theta());
      result.phi = angles.Phi();

      return result;
    }

    /**
     * @brief Calculates the decay angle Cos(Theta) in the Proton Rest frame.
     * @param react, px, py, pz, m The inputs.
     * @return ResultType_t Cosine of the polar decay angle.
     */
    inline ResultType_t ElS_CosThetaProtonRest(const RVecIndexMap& react, 
					       const RVecResultType& px, const RVecResultType& py, 
					       const RVecResultType& pz, const RVecResultType& m) 
    {
      return ElectroProtonRestDecay(react, px, py, pz, m).cosTheta;
    }

    /**
     * @brief Calculates the decay angle Phi in the Proton Rest frame.
     * @param react, px, py, pz, m The inputs.
     * @return ResultType_t Azimuthal decay angle (Phi).
     */
    inline ResultType_t ElS_PhiProtonRest(const RVecIndexMap& react, 
					  const RVecResultType& px, const RVecResultType& py, 
					  const RVecResultType& pz, const RVecResultType& m) 
    {
      return ElectroProtonRestDecay(react, px, py, pz, m).phi;
    }

  } // namespace physics
} // namespace rad

