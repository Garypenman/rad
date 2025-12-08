#pragma once

//!  Derived class to configure ePIC root files

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for Electron scattering reactions
  Derive from this if your experiment is electron scattering
*/
#include "ConfigReaction.h"
#include "RVecHelpers.h"
#include "BasicKinematics.h"
#include "Beams.h"
#include "CommonDefines.h"
#include "ParticleInjector.h" 

//#include "KinematicsProcessor.h"

namespace rad{

  namespace electroion{
    const std::string  BeamIndices() {return Form("%s,%s",names::BeamIon().data(),names::BeamEle().data()); }//"beam_ele,beam_ion";}
  }
}
using rad::electroion::BeamIndices;

//Following include needs previous line.
//#include "ParticleCreator.h"

namespace rad{
  namespace config{
    
 
    //! Class definition

    class ElectroIonReaction : public ConfigReaction {


    public:

      ElectroIonReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ConfigReaction{treeName,fileNameGlob,columns} {

      }
     ElectroIonReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ConfigReaction{treeName,filenames,columns} {

      }
      ElectroIonReaction(ROOT::RDataFrame rdf) : ConfigReaction{rdf} {
      }
      /**
       * Make map that links particle names to indices in user functions
       * in C++ functions you can use the RVecIndexMap object indexed by 
       * name of the reaction component you need
       */
      void makeParticleMap() override {
	//now we have beam and scattered electron
	DefineVirtualPhoton();

	
	//note, ordering in arguments, map and names must be maintained
  	Define(names::ReactionMap().data(),
	       [](const int& beamion, const int& beamel,
		  const RVecI& baryons,const RVecI& mesons,const int& scatel,const int& virtgamma){
		 return RVecIndexMap{{beamion},{beamel},baryons,mesons,{scatel},{virtgamma}};},
	       {names::BeamIon().data(),names::BeamEle().data(),
		names::Baryons().data(),names::Mesons().data(),
		names::ScatEle().data(),names::BeamGamma().data()});

	ConfigReaction::makeParticleMap();
      }
      // /**
      //  * Make map that links particle names to indices in user functions
      //  * in C++ functions you can use the RVecIndexMap object indexed by 
      //  * name of the reaction component you need
      //  */
      // void makeBeamIndices() override {
      // 	//note, ordering in arguments, map and names must be maintained
      // 	Define(names::ReactionMap().data(),
      // 	       [](const int& beamion, const int& beamel){
      // 		 return ROOT::RVecI{{beamion},{beamel} } };
      // 	       );
      // }
      /** 
       * Set constant index for beam electron, scattered electron and beam ion
       * This assumes constant position in collection (e.g in some HepMC3 files)
       * and update the current frame to the aliased one
       */
      void setBeamElectronIndex(const int idx){
	setParticleIndex(names::BeamEle().data(),idx);
      }
      void setScatElectronIndex(const int idx){
	setParticleIndex(names::ScatEle().data(),idx,11);
      }
      void setScatElectronCandidates(const Indices_t& idx){
	setParticleCandidates(names::ScatEle().data(),idx);
      }
      template<typename Lambda>
      void setScatElectronIndex(Lambda&& func,const ROOT::RDF::ColumnNames_t & columns = {}){
	setParticleIndex(names::ScatEle().data(),func,columns,11);
      }
      void setBeamIonIndex(const int idx){
	setParticleIndex(names::BeamIon().data(),idx);
      }
      /**
       * Allow variable index for scattered electron
       */
       template<typename Lambda>
      void setScatElectron(Lambda&& func,const ROOT::RDF::ColumnNames_t & columns = {} ){
	 setParticleIndex(names::ScatEle().data(),func, columns, 11);
      }
      /**
       * Allow variable index for scattered electron
       */
       template<typename Lambda>
      void setScatElectronCandidates(Lambda&& func,const ROOT::RDF::ColumnNames_t & columns = {} ){
	 setParticleCandidates(names::ScatEle().data(),func, columns, 11);
       }
      /**
       * Allow variable index for beam electron
       */
      template<typename Lambda>
      void setBeamElectron(Lambda&& func,const ROOT::RDF::ColumnNames_t & columns = {} ){
	setParticleIndex(names::BeamEle().data(),func, columns);
      }
      /**
       * Allow variable index for beam ion
       */
      template<typename Lambda>
      void setBeamIon(Lambda&& func,const ROOT::RDF::ColumnNames_t & columns = {} ){
	setParticleIndex(names::BeamIon().data(),func, columns);
      }
      /**
       * set fixed P4 for electron beam
       */
      void setBeamElectron(double x,double y,double z){
	_p4el_beam = PxPyPzMVector{x,y,z,0.000510999};
      }
     /**
       * set fixed P4 for ion beam
       */
      void setBeamIon(double x,double y,double z,double m=0.938272){
	_p4ion_beam = PxPyPzMVector{x,y,z,m};
     }
      void DefineBeamElectron(){
	if( !_useBeamsFromMC ) return;
	//add to particles lists, i.e. components of _p4el_beam to rec_px etc
	//note copying p4 so return will never change
	auto p4=_p4el_beam;
	Define(rad::names::P4BeamEle(),[p4](){return p4;},{});
	//	Particles().Beam(rad::names::BeamEle().data(),rad::names::P4BeamEle());
      }
      void DefineBeamIon(){
	if( !_useBeamsFromMC ) return;
	//add to particles lists, i.e. components of _p4ion_beam to rec_px etc
	auto p4=_p4ion_beam;
	//note copying p4 so return will never change
	Define(rad::names::P4BeamIon(),[p4](){return p4;},{});
	// 	Particles().Beam(rad::names::BeamIon().data(),rad::names::P4BeamIon());

      }
      void DefineVirtualPhoton(){
	//add to particles lists, i.e. components of _p4el_beam to rec_px etc
	//note copying p4 so return will never change
	//Particles().Diff(rad::names::BeamGamma().data(),{rad::names::BeamEle()},{rad::names::ScatEle()});
      }
      void FixBeamElectronMomentum(double x,double y,double z){
	setBeamElectron(x,y,z);
	_useBeamsFromMC=true;
	DefineBeamElectron();
      }
      void FixBeamIonMomentum(double x,double y,double z,double m=0.938272){
	setBeamIon(x,y,z,m);
	_useBeamsFromMC=true;
	DefineBeamIon();
      }
       /**
       * Get the particle creator project to add intermediate
       * beam or missing particles
       * Set myself as reaction to be sure 
       * and avoid having to define copy constuctor 
       */
      //const ParticleCreator& Particles() {_particles.SetReaction(this);return _particles;};

      PxPyPzMVector P4BeamIon()const {return _p4ion_beam;}
      PxPyPzMVector P4BeamEle()const {return _p4el_beam;}

    
    protected:
      PxPyPzMVector _p4el_beam;
      PxPyPzMVector _p4ion_beam;
    private:
      /**
       *
       */
      //ParticleCreator _particles{*this};

      
    };//ElectroIonReaction

  }//config
  namespace electroion{
    /**
     * virtual photon 4-vector
     */
    template<typename Tp, typename Tm>
    PxPyPzMVector PhotoFourVector(const RVecIndexMap& react, const Tp &px, const Tp &py, const Tp &pz, const Tm &m){
      return FourVector(react[names::OrderCreated()][names::OrderVirtGamma()],px,py,pz,m);
      
    }
 
  }
  
}//rad

//Declare we are using this PhotoFourVector in kinematics
using rad::electroion::PhotoFourVector;
