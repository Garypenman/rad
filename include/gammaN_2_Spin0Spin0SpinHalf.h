#pragma once

//#include "Beams.h"
#include "BasicKinematics.h"
#include "ConfigReaction.h"


namespace rad{
  namespace gn2s0s0s12{
    
   struct Angles_t{
     double Theta=0;
     double CosTheta=0.;
     double Phi=0.;
    };

   ///\brief functions to compute standard reaction kinematics

    /**
    * calculate CM kinematics from beam
    */
   template<typename Tp, typename Tm>
    PxPyPzMVector PhotoCMVector(const config::RVecIndexMap& react, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){

     //Note PhotoFourVector must be defined in the reaction config file
     //e.g. in PhotoIonReaction.h or ElectroIonReaction.h
     return beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m) +
       PhotoFourVector(react,px,py,pz,m);
 
     }

   /**
    * calculate CM decay angles
    */
    template<typename Tp, typename Tm>
    Angles_t PhotoCMDecay(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){

      auto cm = PhotoCMVector(react,px,py,pz,m);
      auto cmBoost = cm.BoostToCM();
      auto mes = FourVector(react[names::MesonsIdx()],px,py,pz,m);
      PxPyPzMVector cm_mes=boost(mes,cmBoost);

      Angles_t result;
      result.CosTheta=TMath::Cos(cm_mes.Theta());
      result.Theta=cm_mes.Theta();
      result.Phi=cm_mes.Phi();
      return result;
     }
    /**
    * calculate Helicity decay angles
    * z-axis along -baryon in meson rest frame
    */
    template<typename Tp, typename Tm>
    Angles_t PhotoHelicityDecay(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){

      auto baryon = reactkine::BaryonFourVector(react,px,py,pz,m);
      auto meson = reactkine::MesonFourVector(react,px,py,pz,m);

      auto decBoost = meson.BoostToCM();
      //vectors in rest/decay frame of meson
      auto decBar=boost(baryon,decBoost);
      auto decGamma=boost(PhotoFourVector(react,px,py,pz,m),decBoost);
      
      XYZVector  zV=-decBar.Vect().Unit();
      XYZVector  yV=decBar.Vect().Cross(decGamma.Vect()).Unit();
      XYZVector  xV=yV.Cross(zV).Unit();

      //four vector of first [0] decay product
      auto child1 = FourVector(react[names::MesonsIdx()][0],px,py,pz,m);
      auto decChild1=boost(child1,decBoost);

      //calculate decay angles
      XYZVector angles(decChild1.Vect().Dot(xV),decChild1.Vect().Dot(yV),decChild1.Vect().Dot(zV));
      //store in angles struct
      Angles_t result;
      result.CosTheta=TMath::Cos(angles.Theta());
      result.Theta=angles.Theta();
      result.Phi=angles.Phi();
      return result;
    }
    template<typename Tp, typename Tm>
    double CosThetaHel(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = PhotoHelicityDecay(react,px,py,pz,m);
      return angles.CosTheta;
    }
    template<typename Tp, typename Tm>
      double ThetaHel(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = PhotoHelicityDecay(react,px,py,pz,m);
      return angles.Theta;
    }
    template<typename Tp, typename Tm>
      double PhiHel(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = PhotoHelicityDecay(react,px,py,pz,m);
      return angles.Phi;
    }
    /**
    * calculate GJ decay angles
    * z-axis along gamma direction in meson rest frame
    */
    template<typename Tp, typename Tm>
    Angles_t PhotoGJDecay(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){

      auto baryon = reactkine::BaryonFourVector(react,px,py,pz,m);
      auto meson = reactkine::MesonFourVector(react,px,py,pz,m);

      auto decBoost = meson.BoostToCM();
      //vectors in rest/decay frame of meson
      auto decBar=boost(baryon,decBoost);
      auto decGamma=boost(PhotoFourVector(react,px,py,pz,m),decBoost);
      
      XYZVector  zV=decGamma.Vect().Unit();
      XYZVector  yV=decBar.Vect().Cross(decGamma.Vect()).Unit();
      XYZVector  xV=yV.Cross(zV).Unit();

      //four vector of first [0] decay product
      auto child1 = FourVector(react[names::MesonsIdx()][0],px,py,pz,m);
      auto decChild1=boost(child1,decBoost);

      //calculate decay angles
      XYZVector angles(decChild1.Vect().Dot(xV),decChild1.Vect().Dot(yV),decChild1.Vect().Dot(zV));
      //store in angles struct
      Angles_t result;
      result.CosTheta=TMath::Cos(angles.Theta());
      result.Theta=angles.Theta();
      result.Phi=angles.Phi();
      return result;
    }
    template<typename Tp, typename Tm>
    double CosThetaGJ(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = PhotoGJDecay(react,px,py,pz,m);
      return angles.CosTheta;
    }
     template<typename Tp, typename Tm>
    double PhiGJ(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = PhotoGJDecay(react,px,py,pz,m);
      return angles.Phi;
    }
     
     template<typename Tp, typename Tm>
       double TrentoPhi(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      
       auto ebeam = beams::InitialFourVector(react[names::ElectroEleIdx()][0],px,py,pz,m);
      auto escat = FourVector(react[names::ScatEleIdx()],px,py,pz,m);
      
      auto pbeam = beams::InitialFourVector(react[names::ElectroIonIdx()][0],px,py,pz,m);
      auto baryon = reactkine::BaryonFourVector(react,px,py,pz,m);
      
      // Virtual photon q = k - k'    
      auto q = ebeam - escat;    
      // Boost everything to gamma*?p CM (recommended for numerical stability)
      auto cm   = q + pbeam;
      auto Bcm  = cm.BoostToCM();    
      auto k_cm   = boost(ebeam,  Bcm);    
      auto kp_cm  = boost(escat,  Bcm);    
      auto q_cm   = boost(q,      Bcm);    
      auto Pp_cm  = boost(baryon, Bcm);    
      // Extract 3-vectors    
      auto k    = k_cm.Vect();    
      auto kp   = kp_cm.Vect();    
      auto qv   = q_cm.Vect();    
      auto Ppv  = Pp_cm.Vect();    
      // Plane normals    
      auto n_lep = (k.Cross(kp)).Unit();  
      // lepton-scattering plane  
      auto n_had = (qv.Cross(Ppv)).Unit();  
      // hadron plane    
      auto qhat  = qv.Unit();    
      // Trento azimuthal angle    
      double cosphi = n_lep.Dot(n_had);  
      double sinphi = qhat.Dot(n_lep.Cross(n_had));
      return std::atan2(sinphi, cosphi);
     }
     
  }
}
