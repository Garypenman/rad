#pragma once

//!  Derived class to configure HepMC root files for electroproduction

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for HepMC files with fixed particle order
*/
#include "PhotoIonReaction.h"


namespace rad{
  namespace config{
    using rad::names::data_type::MC;

    //! Class definition

    class HepMCPhoto : public PhotoIonReaction {


    public:

    HepMCPhoto(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={}) : PhotoIonReaction{treeName,fileNameGlob,columns} {
	
      }
    HepMCPhoto(const std::string_view treeName, const ROOT::RVec<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} ) : PhotoIonReaction{treeName,filenames,columns} {

      }

     void AliasMomentumComponents(){
      AddType("mc");
      setBranchAlias("particles.momentum.m_v1",MC()+"px");
      setBranchAlias("particles.momentum.m_v2",MC()+"py");
      setBranchAlias("particles.momentum.m_v3",MC()+"pz");
      setBranchAlias("particles.mass",MC()+"m");
      setBranchAlias("particles.pid",MC()+"pid");

      DefineForAllTypes("phi", Form("rad::ThreeVectorPhi(components_p3)"));
      DefineForAllTypes("theta", Form("rad::ThreeVectorTheta(components_p3)"));
      DefineForAllTypes("eta", Form("rad::ThreeVectorEta(components_p3)"));
      DefineForAllTypes("pmag", Form("rad::ThreeVectorMag(components_p3)"));
    }

    };
  }
}
