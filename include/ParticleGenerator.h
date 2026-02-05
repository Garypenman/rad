#pragma once

#include "ParticleCreator.h"
#include "ConfigReaction.h"
#include "StringUtilities.h"
#include "Random.h"

/**
 * @file ParticleGenerator.h
 * @brief Contains functionality for generating new particles in RDataFrame operations based on 4-vector data.
 *
 * This file provides utilities to add new particles, intended to be used
 * within RDataFrame Define operations.
 * The current implementation offers methods for two-body decays.
 */

namespace rad{
  namespace generator{

    using ROOT::RVecD;

    /**
     * @brief Calculate the 4-vector of a two-body decay product in the lab frame.
     *
     * Given a parent particle index, its 4-momentum components, and decay product masses,
     * this function calculates the 4-vector of the first decay product after a two-body decay,
     * applying random decay kinematics in the parent rest frame and boosting back to the lab frame.
     *
     * @tparam Tp Type for momentum components (e.g., double, float).
     * @tparam Tm Type for mass components (e.g., double, float).
     * @param pidx Index of the parent particle in the input vectors.
     * @param masses Masses of the two decay products as a vector of length 2.
     * @param px X-components of all particle momenta.
     * @param py Y-components of all particle momenta.
     * @param pz Z-components of all particle momenta.
     * @param m Masses of all particles.
     * @return RVec<double> The (x, y, z, m) 4-vector of the first decay product in the lab frame.
     */
    template<typename Tp, typename Tm>
    RVec<double> CalculateTwoBody(const int &pidx, const RVec<Tm> &masses, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m) {

      //get boost
      auto parent = FourVector(pidx,px,py,pz,m);
      auto decBoost = parent.BoostToCM();

      //pmag of 2body decay
      auto Mpar = parent.M();

      //Calculate 2-body breakup rest-frame momentum
      auto sum_masses = masses[0] + masses[1];
      auto dif_masses = masses[0] - masses[1];
      auto term1 = Mpar*Mpar - sum_masses*sum_masses;
      auto term2 = Mpar*Mpar - dif_masses*dif_masses;
      auto p = TMath::Sqrt(term1*term2)/(2*Mpar);

      //random thetaphi in helicity frame
      double costheta = rad::random::Generator().Uniform(-1.,1.);
      auto sintheta = TMath::Sqrt(1 - costheta*costheta);
      auto phi =  rad::random::Generator().Uniform(0, 2.0*TMath::Pi());

      //define momentum components
      auto dpx = p * sintheta * cos(phi);
      auto dpy = p * sintheta * sin(phi);
      auto dpz = p * costheta;

      //define first decay product
      PxPyPzMVector cmpart1 = {dpx,dpy,dpz,masses[0]};

      //boost back to lab frame
      auto part1 = boost(cmpart1,-decBoost);

      RVec<double> result={part1.X(),part1.Y(),part1.Z(),part1.M()};
      return result;
    }

    /**
     * @brief Appends a two-body decay product to momentum and mass vectors.
     *
     * Utilizes CalculateTwoBody to compute the decay product, then appends its components to the relevant vectors.
     *
     * @tparam Tp Type for momentum components (e.g., double, float).
     * @tparam Tm Type for mass components (e.g., double, float).
     * @param pidx Index of the parent particle in the input vectors.
     * @param masses Masses of the two decay products as a vector of length 2.
     * @param px X-components of all particle momenta (will be appended to).
     * @param py Y-components of all particle momenta (will be appended to).
     * @param pz Z-components of all particle momenta (will be appended to).
     * @param m Masses of all particles (will be appended to).
     * @param iafter Indices at which to insert the new particle, required to make sure prior particles are calculated first
     * @return int Index of the newly created particle in the vectors.
     */
    template<typename Tp, typename Tm>
    int ParticleCreateTwoBody(const int &pidx, const RVec<Tm>& masses, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m, const RVecI &iafter) {

      auto result = CalculateTwoBody(pidx,masses,px,py,pz,m);

      int idx = px.size();

      px.push_back(result[0]);
      py.push_back(result[1]);
      pz.push_back(result[2]);
      m.push_back(result[3]);

      return idx;
    }

    /**
     * @class ParticleGenerator
     * @brief Helper class to generate new particles for RDataFrame Define operations.
     *
     * Inherits from rad::config::ParticleCreator and provides methods to add new particles,
     * currently supporting only two-body decays.
     */
    class ParticleGenerator : public rad::config::ParticleCreator{

    public:

      /**
       * @brief Default constructor.
       */
      ParticleGenerator() = default;

      /**
       * @brief Construct a ParticleGenerator with a ConfigReaction and seed.
       * @param cr Reference to the configuration for the reaction.
       * @param seed Seed for the random number generator.
       */
      ParticleGenerator(rad::config::ConfigReaction &cr, size_t seed=1):
        rad::config::ParticleCreator{cr}
      {
        Reaction()->InitRandom(seed);
      };

      /**
       * @brief Generate a two-body decay and define the resulting particles.
       *
       * This method is designed to be used in RDataFrame Define operations. It creates the first decay product
       * using a user-defined function, and the second decay particle by difference from the parent.
       *
       * @param names Names of the daughter particles (length 2).
       * @param masses Masses of the daughter particles.
       * @param parent Name of the parent particle.
       */
      void GenerateTwoBody(const ROOT::RVec<std::string> &names,
                           const ROOT::RVecD& masses, const string &parent){

        //combine masses vector to string
        std::string smasses=rad::utils::combineAnyVectorToString(masses);

        //create Define function string
        auto expr = Form( "rad::generator::ParticleCreateTwoBody(%s,%s",parent.data(), smasses.data());

        //create first decay particle
        DefineParticle(names[0],ROOT::RVec<string>(),expr);

        //create second decay particle
        Diff(names[1],{parent},{names[0]});
      }

    };// end Class ParticleGenerator

  }//end namespace generator
}//end namespace rad
