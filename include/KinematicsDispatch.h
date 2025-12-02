#pragma once

//#include "KinematicsProcessor.h"


namespace rad{
  // Helper to define the KinematicsProcessor call with concrete types (Tp=double, Tm=double)
  // This is required to pass a functor to RDF::Define
  using rad::config::DeduceColumnVectorType;
  using rad::config::ColType;

  //Make KinemticsProcessor a template so can #include this
  //But this could also become a generic Discpatch anyway
  template <typename T>
  void DefineKinematicsProcessor(config::ConfigReaction& cr,
				 T& processor, //T->KinematicsProcessor
				 //				 const std::string& name, 
				 const std::string& type) 
  {
    //construct required column names
    //Note KineIndices should be independent ot type
    //Particle order is defined by the ParticleCreator
    auto indicesName = type + names::KineIndices();
    ROOT::RDF::ColumnNames_t columns = {indicesName ,type+"px",type+"py",type+"pz",type+"m"};
    auto name = type + names::KineComponents();

    //Get the column types for deducing template arguments
    auto vecTypeP = DeduceColumnVectorType(&cr, type + "px");
    auto vecTypeM = DeduceColumnVectorType(&cr, type + "m");
    
    if (vecTypeP == ColType::Float &&
	vecTypeM == ColType::Double) {

      cr.Define(name, 
		// This is the concrete, non-templated callable RDataFrame needs:
		//must capture processor by value (copy) for thread safety
		[processor](const RVecIndices& indices, 
			    const ROOT::RVecF& px, const ROOT::RVecF& py, 
			    const ROOT::RVecF& pz, const ROOT::RVecD& m) 
		{
		  // The C++ compiler resolves the template types here (double, double)
		  return processor.template operator()< ROOT::RVecF, ROOT::RVecD>(indices, px, py, pz, m);
		}, 
		columns
		);
      //done Define
    }
    else if (vecTypeP == ColType::Double &&
	     vecTypeM == ColType::Double) {

      cr.Define(name, 
		// This is the concrete, non-templated callable RDataFrame needs:
		//must capture processor by value (copy) for thread safety
		[processor](const RVecIndices& indices, 
			    const ROOT::RVecD& px, const ROOT::RVecD& py, 
			    const ROOT::RVecD& pz, const ROOT::RVecD& m) 
		{
		  // The C++ compiler resolves the template types here (double, double)
		  return processor.template operator()< ROOT::RVecD, ROOT::RVecD>(indices, px, py, pz, m);
		}, 
		columns
		);
      //done Define
    }
    else if (vecTypeP == ColType::Double &&
	     vecTypeM == ColType::Float) {

      cr.Define(name, 
		// This is the concrete, non-templated callable RDataFrame needs:
		//must capture processor by value (copy) for thread safety
		[processor](const RVecIndices& indices, 
			    const ROOT::RVecD& px, const ROOT::RVecD& py, 
			    const ROOT::RVecD& pz, const ROOT::RVecF& m) 
		{
		  // The C++ compiler resolves the template types here (double, double)
		  return processor.template operator()< ROOT::RVecD, ROOT::RVecF>(indices, px, py, pz, m);
		}, 
		columns
		);
      //done Define
    }
    else if (vecTypeP == ColType::Float &&
	     vecTypeM == ColType::Float) {

      cr.Define(name, 
		// This is the concrete, non-templated callable RDataFrame needs:
		//must capture processor by value (copy) for thread safety
		[processor](const RVecIndices& indices, 
			    const ROOT::RVecF& px, const ROOT::RVecF& py, 
			    const ROOT::RVecF& pz, const ROOT::RVecF& m) 
		{
		  // The C++ compiler resolves the template types here (double, double)
		  return processor.template operator()< ROOT::RVecF, ROOT::RVecF>(indices, px, py, pz, m);
		}, 
		columns
		);
      //done Define
    }
    else{
      throw std::runtime_error(std::string("KinematicsDispatch cannot deduce template Types for KinemticsProcessor, ")+type+" "+name);
    }

    //unpack the created RVec back into component arrays
    processor.DefineNewComponentVecs();
     
    return;
    
  }
}
