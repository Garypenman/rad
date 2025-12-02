#pragma once
#include <string>
#include <vector>
#include <functional>
#include "CommonDefines.h" 

namespace rad {

  class KineCalculation {
  public:
    // Signature 1: Standard Kernel (Uses full Map)
    using MapKernel = ResultType_t(*)(const RVecIndexMap&, 
                                      const RVecResultType&, const RVecResultType&, 
                                      const RVecResultType&, const RVecResultType&);

    // Signature 2: Indexed Kernel (Uses specific resolved indices)
    // This matches: ResultType_t Func(const RVecIndices&, const RVecResultType&, ...)
    using IndexKernel = ResultType_t(*)(const RVecIndices&, 
                                        const RVecResultType&, const RVecResultType&, 
                                        const RVecResultType&, const RVecResultType&);


    // -------------------------------------------------------------------
    // Constructor 1: Standard Map Kernel
    // -------------------------------------------------------------------
    template <typename ProcessorT>
    KineCalculation(ProcessorT& kine, std::string name, MapKernel func) {
       kine.DefineKernel(name, func);
    }

    // -------------------------------------------------------------------
    // Constructor 2: Indexed Kernel Adapter
    // -------------------------------------------------------------------
    template <typename ProcessorT>
    KineCalculation(ProcessorT& kine, std::string name, IndexKernel func, 
                    std::vector<ParticleNames_t> particle_groups) 
    {
       // 1. Resolve String Names -> Integer Indices
       RVecIndices resolved_indices;
       for(const auto& group : particle_groups) {
           Indices_t idxs;
           for(const auto& pname : group) {
               idxs.push_back(kine.Creator().GetReactionIndex(pname));
           }
           resolved_indices.push_back(idxs);
       }

       // 2. Capture and Adapt
       auto adapter_lambda = [resolved_indices, func](const RVecIndexMap&, 
                                                      const RVecResultType& px, const RVecResultType& py, 
                                                      const RVecResultType& pz, const RVecResultType& m) 
       {
           return func(resolved_indices, px, py, pz, m);
       };

       // 3. Register
       kine.DefineKernel(name, adapter_lambda);
    }

  };

} // namespace rad

// #pragma once
// #include <string>
// #include <functional>
// #include "CommonDefines.h" 

// namespace rad {

//   class KineCalculation {
//   public:
//     // The exact signature of your non-templated physics kernels
//     using KernelSignature = ResultType_t(*)(const RVecIndexMap&, 
//                                             const RVecResultType&, const RVecResultType&, 
//                                             const RVecResultType&, const RVecResultType&);

//     // We template 'ProcessorT' to avoid circular include issues. 
//     // It works as long as 'kine' has a DefineKernel method.
//   template <typename ProcessorT>
//     KineCalculation( ProcessorT& kine,std::string name, KernelSignature func )
//       : _name(std::move(name)), _func(func) 
//     {
//        // Immediate registration
//        kine.DefineLambda(_name, *this);
//     }


//     // -------------------------------------------------------------------
//     // The Functor Operator (Must be const)
//     // -------------------------------------------------------------------
//     ResultType_t operator()(const RVecIndexMap& map, 
//                             const RVecResultType& px, const RVecResultType& py, 
//                             const RVecResultType& pz, const RVecResultType& m) const 
//     {
//         return _func(map, px, py, pz, m);
//     }

//   private:
//     std::string _name;
//     KernelSignature _func; 
//   };

// } // namespace rad
