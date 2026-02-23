#pragma once
#include "gammaN_2_Spin0Spin0SpinHalf.h"


namespace rad{
  namespace rdf{
     namespace gn2s0s0s12{
 
       void PhotoCMAngles(config::ConfigReaction& cr,const string& name){
	 cr.DefineForAllTypes(name, Form("rad::gn2s0s0s12::PhotoCMDecay(%s,components_p4)",names::ReactionMap().data()));
       }
       
       void CosThetaHel(config::ConfigReaction& cr,const string& name,const string& convention){
	 cr.DefineForAllTypes(name, Form("rad::gn2s0s0s12::CosTheta%s(%s,components_p4)",convention.data(),names::ReactionMap().data()));
       }
       void ThetaHel(config::ConfigReaction& cr,const string& name,const string& convention){
	 cr.DefineForAllTypes(name, Form("rad::gn2s0s0s12::Theta%s(%s,components_p4)",convention.data(),names::ReactionMap().data()));
       }
       void PhiHel(config::ConfigReaction& cr,const string& name,const string& convention){
	 cr.DefineForAllTypes(name, Form("rad::gn2s0s0s12::Phi%s(%s,components_p4)",convention.data(),names::ReactionMap().data()));
       }
    
       void HelicityAngles(config::ConfigReaction& cr,const string& name){
	 CosThetaHel(cr,name+"_CosTheta","Hel");
	 ThetaHel(cr,name+"_Theta","Hel");
	 PhiHel(cr,name+"_Phi","Hel");
       }
       void GJAngles(config::ConfigReaction& cr,const string& name){
	 CosThetaHel(cr,name+"_CosTheta","GJ");
	 PhiHel(cr,name+"_Phi","GJ");
      }
       
       void PhiTrento(config::ConfigReaction& cr){
	 cr.DefineForAllTypes("TrentoPhi", Form("rad::gn2s0s0s12::TrentoPhi(%s,components_p4)",names::ReactionMap().data()));
       }
    
       
    
     }
  }
}
