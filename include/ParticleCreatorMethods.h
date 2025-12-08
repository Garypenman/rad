#pragma once
#include "CommonDefines.h"

namespace rad{

  ///Create a particle as the diffence between ipos particles and ineg
    inline void ParticleCreateByDiff(const Indice_t position,const RVecIndices &indices, ROOT::RVecD &px, ROOT::RVecD &py, ROOT::RVecD &pz, ROOT::RVecD &m) {//,const RVecI& iafter){
      const auto& ipos = indices[0];
      const auto& ineg = indices[1];
      
      //sum the postive 4-vectors
      auto p4 = FourVector(ipos,px,py,pz,m);
      //subtract the negative 4-vectors
      SubtractFourVector(p4,ineg,px,py,pz,m);
      px[position] =p4.X();
      py[position] =p4.Y();
      pz[position] =p4.Z();
      m[position] =p4.M();
    }

  ///Create a particle as the sum of ipos particles
    inline void ParticleCreateBySum(const Indice_t position,const RVecIndices &indices, ROOT::RVecD &px, ROOT::RVecD &py, ROOT::RVecD &pz, ROOT::RVecD &m) {//,const RVecI& iafter){
      const auto& ipos = indices[0];
      
      //sum the postive 4-vectors
      auto p4 = FourVector(ipos,px,py,pz,m);

      px[position] =p4.X();
      py[position] =p4.Y();
      pz[position] =p4.Z();
      m[position] =p4.M();
    }




}
