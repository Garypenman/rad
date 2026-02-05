#pragma once

#include <TMath.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <Math/VectorUtil.h>
#include <ROOT/RVec.hxx>


namespace rad{

  using ROOT::Math::PxPyPzMVector ;
  using ROOT::Math::XYZVector ;

  //kinematics functions etc should return common type
  using ResultType_t = Double_t;
  using RVecResultType = ROOT::RVec<ResultType_t>;
    
  using Indice_t = Int_t;
  using Indices_t = ROOT::RVecI;
  using RVecCombis = ROOT::RVec<Indices_t>; //indices : RVecCombis[particle][combo]
  using RVecIndices = ROOT::RVec<Indices_t>; //indices : RVecIndices[set][particle]
  using RVecIndexMap = ROOT::RVec<Indices_t>;

  enum class ComponentsOrder{X, Y, Z, M};
  constexpr int OrderX() {return static_cast<int> (ComponentsOrder::X); }
  constexpr int OrderY() {return static_cast<int> (ComponentsOrder::Y); }
  constexpr int OrderZ() {return static_cast<int> (ComponentsOrder::Z); }
  constexpr int OrderM() {return static_cast<int> (ComponentsOrder::M); }
 
  using ParticleNames_t = ROOT::RVec<std::string>;
 
}
//#pragma link C++ class ROOT::VecOps::RVec<ROOT::VecOps::RVec<int> >+;
#pragma link C++ class rad::RVecIndices+;
#pragma link C++ class ROOT::RVec<rad::RVecResultType>+;
//#pragma link C++ class ROOT::VecOps::RVec<ROOT::VecOps::RVec<ROOT::VecOps::RVec<rad::ResultType_t> > >+;
