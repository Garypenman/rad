#pragma once

//!  Derived class to configure HepMC root files for electroproduction

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for HepMC files with fixed particle order
*/
#include "ElectroIonReaction.h"


namespace rad{
  namespace config{
   using rad::names::data_type::MC;

    //! Class definition

    class HepMCElectro : public ElectroIonReaction {


    public:

    HepMCElectro(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={}) : ElectroIonReaction{treeName,fileNameGlob,columns} {
	
      }
    HepMCElectro(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ElectroIonReaction{treeName,filenames,columns} {

      }

     void AliasMomentumComponents(){
      AddType(MC());
      
      ParticleInjector injector(this);
      injector.DefineParticleInfo({"px", "py", "pz", "m", "pid"});
    
      injector.AddSource(MC(), 
			 {"particles.momentum.m_v1","particles.momentum.m_v2","particles.momentum.m_v3","particles.mass","particles.pid"}); 
      
      
      injector.CreateUnifiedVectors();
      
      // setBranchAlias("particles.momentum.m_v1",MC()+"px");
      // setBranchAlias("particles.momentum.m_v2",MC()+"py");
      // setBranchAlias("particles.momentum.m_v3",MC()+"pz");
      //setBranchAlias("particles.mass",MC()+"m");
      //setBranchAlias("particles.pid",MC()+"pid");

      DefineForAllTypes("phi", Form("rad::ThreeVectorPhi(components_p3)"));
      DefineForAllTypes("theta", Form("rad::ThreeVectorTheta(components_p3)"));
      DefineForAllTypes("eta", Form("rad::ThreeVectorEta(components_p3)"));
      DefineForAllTypes("pmag", Form("rad::ThreeVectorMag(components_p3)"));
      // DefineForAllTypes("pmag", Form("rad::ThreeVectorMag"),{},names::P3Components());

      // setMyBaseFrame();
     }

     void AliasStableMomentumComponents(){
       AliasMomentumComponents();
       setBranchAlias("particles.status",MC()+"genStat");

       Define(MC()+"final_pid",[](const ROOT::RVecI& pid,const ROOT::RVecI& stat){
	 auto n = pid.size();
	 auto final_pid = pid;
	 final_pid*=(stat==1);
	 // for(size_t i=0;i<n;++i){
	 //   final_id[i]=pid[i]*(stat[i]==1);
	 // }
	 return final_pid;

       },{MC()+"pid",MC()+"genStat"});//simID points from rec to tru

     }
     };
  }
}
