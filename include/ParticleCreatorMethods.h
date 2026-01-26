#pragma once
#include "CommonDefines.h"
#include "BasicKinematics.h"

namespace rad{

  /**
   * @brief Creates a particle by Difference (P_new = P_pos - P_neg).
   * @param position The index in the arrays where the new particle will be stored.
   * @param indices [0] = Indices to Add, [1] = Indices to Subtract.
   */
  inline void ParticleCreateByDiff(const Indice_t position,const RVecIndices &indices, ROOT::RVecD &px, ROOT::RVecD &py, ROOT::RVecD &pz, ROOT::RVecD &m) {
      const auto& ipos = indices[0];
      const auto& ineg = indices[1];
      
      auto p4 = FourVector(ipos,px,py,pz,m);
      SubtractFourVector(p4,ineg,px,py,pz,m);
      
      px[position] = p4.X();
      py[position] = p4.Y();
      pz[position] = p4.Z();
      m[position]  = p4.M();
  }

  /**
   * @brief Creates a particle by Summing (P_new = Sum P_i).
   * @param position The index in the arrays where the new particle will be stored.
   * @param indices [0] = Indices to Add.
   */
  inline void ParticleCreateBySum(const Indice_t position,const RVecIndices &indices, ROOT::RVecD &px, ROOT::RVecD &py, ROOT::RVecD &pz, ROOT::RVecD &m) {
      const auto& ipos = indices[0];
      auto p4 = FourVector(ipos,px,py,pz,m);

      px[position] = p4.X();
      py[position] = p4.Y();
      pz[position] = p4.Z();
      m[position]  = p4.M();
  }
}
