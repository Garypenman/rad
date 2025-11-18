#pragma once

#include "ConfigReaction.h"
#include "Constants.h"
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

namespace rad{
  namespace config{


    class ParticleModifier {
    public:
      // Default constructor: Initializes _modifications to empty and _reaction to nullptr.
      // Use this if you plan to set the ConfigReaction object later using a setter.
      ParticleModifier() = default;

      // Constructor: Takes a non-const reference to ConfigReaction and stores its address.
      // 'explicit' prevents unintended implicit conversions.
      // Per your guarantee: 'cr' is always a valid object, and '_reaction' will not be null after this.
      explicit ParticleModifier(rad::config::ConfigReaction& cr) : _reaction{&cr} {
        // _reaction is now guaranteed to be a valid pointer to 'cr'.
      }

      /**
       * @brief Applies all queued momentum modifications to the RDataFrame chain
       * by constructing a single Filter expression and applying it.
       * @param filter_label An optional string label for this Filter node in the
       * RDataFrame computation graph. This is not used to specify
       * input columns; rather, it's for identification/debugging.
       *
       * This method leverages RDataFrame's lazy evaluation and JIT compilation.
       * By calling Filter with an expression that contains side-effecting function calls
       * (e.g., `rad::epic::ModifyMomentum`), the operations are forced to execute
       * at this specific point in the RDataFrame chain.
       * The filter condition is always made to evaluate to 'true' to keep all events.
       */
      void Apply(const std::string& filter_label) {
        // If no modifications have been queued, there's nothing to apply.
        if (_modifications.empty()) {
	  std::cout << "ParticleModifier::Apply: No modifications added. Skipping Apply." << std::endl;
	  return;
        }

        // _reaction is guaranteed to be valid per class design, so no nullptr check is needed here.

        // Start building the combined filter string expression.
        // The outer parentheses are crucial to ensure the arithmetic evaluation of
        // the concatenated boolean results (from functions like ModifyMomentum)
        // is performed correctly before being multiplied by 'false'.
        std::string filter_expression_str = "(";
        for (const auto& func_call_str : _modifications) { // Iterate through the stored string function calls
	  filter_expression_str += func_call_str;
	  filter_expression_str += "+"; // Concatenate results using '+' (booleans cast to int: 0 for false, 1 for true)
        }
        filter_expression_str.pop_back(); // Remove the last '+' character (e.g., converts "(func1+func2+)" to "(func1+func2)")

        // Append the "always true" component to the filter expression:
        // "*false" converts the integer sum of boolean results to 0.
        // "+ true" then adds 1 (as integer 1) to make the final expression always evaluate to 'true' (boolean 1).
        // This ensures the filter keeps all events while guaranteeing the execution of the side-effecting 'func_call_str'.
        filter_expression_str += ")*false + true";

        std::cout << "ParticleModifier::Apply: Applying filter with expression: '"
                  << filter_expression_str << "' and label: '" << filter_label << "'" << std::endl;

        // Apply all queued modifications using a single Filter operation.
        // Filter signature used: RInterface<Proxied, DS_t> Filter(std::string_view expression, std::string_view name = "")
        //   - 'filter_expression_str' (your 'apply_all'): is the actual condition string RDataFrame evaluates.
        //     Columns referenced within this string (e.g., "rec_px", "rec_py", "rec_pz" inside ModifyMomentum)
        //     will be automatically inferred as dependencies and materialized by RDataFrame.
        //   - 'filter_label' (your 'name'): is used purely as a label for this filter node in the RDataFrame graph.
	_reaction->Filter(filter_expression_str, filter_label);

	//clear modifications so can use again
	_modifications.clear();
      }

      /**
       * @brief Queues a modification to change the momentum of a specific particle.
       * This prepares a string representing a call to the 'rad::epic::ModifyMomentum' function.
       * @param particle A string that represents the particle's index or a way to access it
       * (e.g., "0" for the first particle, or "particle_index_column[some_idx]").
       * @param new_mag A string that represents the new desired magnitude for the momentum
       * (e.g., "5.0", or "target_magnitude_column").
       *
       * This method assumes that 'rec_px', 'rec_py', 'rec_pz' are existing RVec<float> columns
       * in the RDataFrame, and that 'rad::epic::ModifyMomentum' is a C++ function
       * designed to modify these RVecs in-place when called through RDataFrame's JIT.
       */
      void ChangeRecMomentumOfTo(const std::string& particle, const std::string& new_mag) {
        // Assemble the arguments for the 'rad::epic::ModifyMomentum' C++ function.
        // These are strings that RDataFrame's JIT compiler will resolve to actual values/RVecs.
        std::vector<std::string> vars = {particle, new_mag, Rec()+"px", Rec()+"py", Rec()+"pz"};

        // Use a utility function to
        // construct the full string representation of the function call.
        // Example output: "rad::epic::ModifyMomentum(0, 5.0, rec_px, rec_py, rec_pz)"
        auto func_call_str = rad::utils::createFunctionCallStringFromVec("rad::config::ModifyMomentumTo", vars);

        // Store the generated function call string. Multiple such strings will be combined in Apply().
        _modifications.push_back(func_call_str);
      }

      void FixMassTo(const std::string& particle, double mass){
	// Assemble the arguments for the 'rad::epic::ModifyMomentum' C++ function.
        // These are strings that RDataFrame's JIT compiler will resolve to actual values/RVecs.
	std::vector<std::string> vars = {particle, Form("%lf",mass), "rec_m"};

        // Use a utility function to
        // construct the full string representation of the function call.
        // Example output: "rad::config::Fix4VectorMass(0, 0.938, rec_px, rec_py, rec_pz, rec_m)"
        auto func_call_str = rad::utils::createFunctionCallStringFromVec("rad::config::Fix4VectorMass", vars);

        // Store the generated function call string. Multiple such strings will be combined in Apply().
        _modifications.push_back(func_call_str);
  
      }
      
      void OverwriteRecParticle(const string oldpart, const string newpart){
	std::vector<std::string> vars = {oldpart, newpart, Rec()+"px", Rec()+"py", Rec()+"pz", Rec()+"m", Rec()+"pid"};
	auto func_call_str = rad::utils::createFunctionCallStringFromVec("rad::config::OverwriteIfNoRecon", vars);
	
	_modifications.push_back(func_call_str);
	
      }
      
      void ModifyMomentum(const string oldpart, const string newpart){
	std::vector<std::string> vars = {oldpart, newpart, MC()+"px", MC()+"py", MC()+"pz", MC()+"m"};
	
	auto func_call_str = rad::utils::createFunctionCallStringFromVec("rad::config::ModifyMomentumParts", vars);
	
	_modifications.push_back(func_call_str);

      }
      
    protected:
      void addModification(const std::string& mod){
	_modifications.push_back(mod);
      }

      private:
      // Pointer to the ConfigReaction object, which encapsulates the RDataFrame chain.
      // It is guaranteed to be valid and non-null after ParticleModifier construction (if explicit constructor is used).
      // Not 'mutable' as 'Apply' is not a const method, and 'Filter' modifies the chain.
      rad::config::ConfigReaction* _reaction = nullptr;

      // A collection of string representations of function calls.
      // Each string corresponds to a single momentum modification operation to be applied.
      std::vector<std::string> _modifications;
    };
    
   
    /**
     * @brief Changes the magnitude of a 3-vector momentum at a specified index within RVecs.
     *
     * This function scales the x, y, and z components of a 3-vector such that its
     * new magnitude matches the provided 'new_mag'. It operates directly on the
     * RVecs, modifying them in place.
     *
     * @param index The index of the 3-vector to be modified within the RVecs.
     * It's crucial that this index is valid for all provided RVecs.
     * @param new_mag The desired new magnitude for the 3-vector at 'index'.
     * @param x A reference to the ROOT::RVec containing the x-components of the momenta.
     * This RVec will be modified.
     * @param y A reference to the ROOT::RVec containing the y-components of the momenta.
     * This RVec will be modified.
     * @param z A reference to the ROOT::RVec containing the z-components of the momenta.
     * This RVec will be modified.
     *
     * @note This function assumes that 'x', 'y', 'z', RVecs are all of the
     * same size and that 'index' is a valid, in-bounds index for all of them.
     * No bounds checking is performed within the function for performance.
     * Caller is responsible for ensuring valid indices.
     *
     *
     * @warning Modifies the input RVecs 'x', 'y', and 'z' in place.
     */
    template <typename Tn, typename To>
    bool ModifyMomentumTo(int index, Tn new_mag, ROOT::RVec<To> &x, ROOT::RVec<To> &y, ROOT::RVec<To> &z ){
      // Calculate the scaling factor needed to achieve the new magnitude.
      // This is done by dividing the desired new magnitude by the current magnitude.
      auto mag = TMath::Sqrt(x[index]*x[index] + y[index]*y[index] + z[index]*z[index]); 
      auto scale = static_cast<To>(new_mag/mag);
      // Apply the scaling factor to each component of the 3-vector at the specified index.
      // The *= operator modifies the RVec elements in place.
      x[index]*=scale;
      y[index]*=scale;
      z[index]*=scale;
      return true;//always return true 
    }
    
    template <typename Tm>
    bool Fix4VectorMass(int index, double mass,ROOT::RVec<Tm> &m ){
      m[index]=mass;
      return true;//always return true 
    }
    
    
    template <typename Tp, typename Tm, typename Tpdg>
      bool OverwriteIfNoRecon(const int oldpart, const int newpart, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m, RVec<Tpdg> &pid){
      
      /* std::cout << "Old: " << oldpart << " " << "New: " << newpart << " " << std::endl; */
      /* std::cout << "Old PDG: " << pid[oldpart] << "New PDG: " << pid[newpart] << " " << std::endl; */
      /* std:: cout << endl; */
      std::cout << "Old px: " << px[oldpart] << "New px: " << px[newpart] << " " << std::endl;
      
      //doesnt work while pid vector bugged / offset from momentum vectors
      /* if(pid[oldpart] == rad::constant::InvalidEntry<int>() || pid[oldpart] == -1){ */
      //use zero momentum as proxy for good pid in particle creator for now.
      if( px[oldpart] == 0 || std::isnan(px[oldpart]) ){
	px[oldpart] = px[newpart];
	py[oldpart] = py[newpart];
	pz[oldpart] = pz[newpart];
	m[oldpart] = m[newpart];
	pid[oldpart] = pid[newpart];
      }
      return true;
    }
    
    template <typename Tp, typename Tm>
      bool ModifyMomentumParts(const int oldpart, const int newpart, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m){
      
      px[oldpart] = px[newpart];
      py[oldpart] = py[newpart];
      pz[oldpart] = pz[newpart];
      m[oldpart] = m[newpart];
      return true;
    }
    
  }
}
