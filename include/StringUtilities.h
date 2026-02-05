#pragma once

#include <ROOT/RVec.hxx> // For ROOT::RVecS alias (which is RVec<std::string>)
#include <iostream>    // Required for input/output operations (e.g., std::cout)
#include <string>      // Required for std::string
#include <sstream>     // Required for std::stringstream (for efficient string building)
#include <utility>     // Required for std::forward (for perfect forwarding in templates)
#include <algorithm> // Required for std::transform and std::tolower
#include <cctype>    // Required for std::tolower


namespace rad{
  namespace util{
    
    using RVecS = ROOT::RVec<std::string>;

    /**
     * @brief Concatenates an arbitrary number of string containers into a single ROOT::RVec<std::string>.
     * * @tparam Containers A pack of containers (e.g., ROOT::RVec<std::string>, std::initializer_list<std::string>)
     * where the element type is convertible to std::string.
     * @param containers The containers to be concatenated.
     * @return ROOT::RVec<std::string> The consolidated vector of strings.
     */
    template <typename... Containers>
    ROOT::RVec<std::string> concatenateStringVectors(const Containers&... containers) {
        
        // --- 1. Calculate the total required size (optional but recommended for reserve) ---
        size_t total_size = (0 + ... + containers.size());
        
        ROOT::RVec<std::string> result;
        result.reserve(total_size);

        // --- 2. Concatenate using a C++17 Fold Expression ---
        // The fold expression iterates over the parameter pack and inserts elements.
        ( (result.insert(result.end(), containers.begin(), containers.end())), ... );
        
        return result;
    }
    /**
     * @brief Helper struct to append arguments with a comma separator.
     * This is used internally by the createFunctionCallString template.
     * @tparam T The type of the current argument.
     */
    template <typename T>
      struct ArgAppender {
	std::stringstream& ss; // Reference to the stringstream where the string is being built
	bool& firstArg;        // Reference to a boolean flag to track if it's the first argument

	/**
	 * @brief Constructor for ArgAppender.
	 * @param stream The stringstream to append to.
	 * @param isFirst Reference to the first argument flag.
	 */
      ArgAppender(std::stringstream& stream, bool& isFirst) : ss(stream), firstArg(isFirst) {}

	/**
	 * @brief Overloaded operator() to append an argument.
	 * It adds a comma before the argument if it's not the first one.
	 * @param arg The argument to append.
	 */
	void operator()(T&& arg) const {
	  // If it's not the first argument, add a comma and a space for readability.
	  if (!firstArg) {
            ss << ", ";
	  }
	  // Append the argument to the stringstream.
	  ss << std::forward<T>(arg);
	  // After appending, it's no longer the first argument.
	  firstArg = false;
	}
      };

    /**
     * @brief Creates a string representation of a function call.
     * This version handles functions with zero arguments.
     * @param funcName The name of the function.
     * @return A string representing the function call (e.g., "myFunc()").
     */
    std::string createFunctionCallString(const std::string& funcName) {
      return funcName + "()";
    }

    /**
     * @brief Creates a string representation of a function call.
     * This template version handles functions with one or more arguments.
     * @tparam Args The types of the arguments (deduced automatically).
     * @param funcName The name of the function.
     * @param args The arguments to be included in the function call string.
     * These can be any types that can be streamed to std::stringstream
     * (e.g., std::string, int, double, etc.).
     * @return A string representing the function call (e.g., "myFunc(arg1, arg2)").
     */
    template <typename... Args>
      std::string createFunctionCallString(const std::string& funcName, Args&&... args) {
      std::stringstream ss; // Create a stringstream to build the result string.
      ss << funcName << "("; // Start with the function name and opening parenthesis.

      bool firstArg = true; // Flag to manage comma separation for arguments.

      // Use a fold expression (C++17 and later) to iterate over arguments.
      // For each argument, an ArgAppender object is created and invoked,
      // which handles adding commas and appending the argument to the stringstream.
      (ArgAppender<Args>(ss, firstArg)(std::forward<Args>(args)), ...);

      ss << ")"; // End with the closing parenthesis.
      return ss.str(); // Return the built string.
    }

   /**
     * @brief Creates a string representation of a function call.
     * This template version handles a vector of string arguments.
     * @param funcName The name of the function.
     * @param args The arguments to be included in the function call string.
     * These can be any types that can be streamed to std::stringstream
     * (e.g., std::string, int, double, etc.).
     * @return A string representing the function call (e.g., "myFunc(arg1, arg2)").
     */
    std::string createFunctionCallStringFromVec(const std::string& funcName, const  ROOT::RVec<std::string>&  args) {
      std::stringstream ss; // Create a stringstream to build the result string.
      ss << funcName << "("<<args[0]; // Start with the function name and opening parenthesis.
      if(args.size()>1){
	std::for_each(args.begin()+1, args.end(),
		       [&ss](const std::string& s) {
			 ss << "," << s;
		       });
      }
      ss << ")"; // End with the closing parenthesis.
      return ss.str(); // Return the built string.
    }

    /**
     * @brief Replaces all occurrences of a specified substring within a string.
     *
     * This function iterates through the input string, finding all instances of
     * 'oldSubstr' and replacing them with 'newSubstr'.
     *
     * @param str The original string in which replacements will be made.
     * @param oldSubstr The substring to be replaced.
     * @param newSubstr The substring to replace 'oldSubstr' with.
     * @return A new string with all occurrences replaced.
     */
    std::string replaceAll(std::string& str, const std::string& oldSubstr, const std::string& newSubstr) {
      // Start searching from the beginning of the string.
      size_t pos = 0;

      // Loop until no more occurrences of oldSubstr are found.
      while ((pos = str.find(oldSubstr, pos)) != std::string::npos) {
        // Replace the found occurrence.
        // str.replace(position, length_of_old_substring, new_substring_content)
        str.replace(pos, oldSubstr.length(), newSubstr);

        // Advance the search position by the length of the new substring.
        // This is crucial to avoid infinite loops if newSubstr contains oldSubstr
        // (e.g., replacing "a" with "aa") and to continue searching
        // after the newly inserted text.
        pos += newSubstr.length();
      }
      // Return the modified string.
      return str;
    }

/**
 * @brief Combines a vector of strings into a single string, formatted as a
 * comma-separated list enclosed in curly braces.
 *
 * @param stringVector The vector of strings to combine.
 * @return A single string representing the combined vector (e.g., "{item1, item2, item3}").
 */
std::string combineVectorToString(const ROOT::RVec<std::string>& stringVector) {
    std::stringstream ss; // Create a stringstream for efficient string building.
    ss << "{";            // Prepend with an opening curly brace.

    // Use a loop to iterate through the vector elements.
    for (size_t i = 0; i < stringVector.size(); ++i) {
        ss << stringVector[i]; // Append the current string.

        // If it's not the last element, append a comma and a space.
        if (i < stringVector.size() - 1) {
            ss << ", ";
        }
    }

    ss << "}"; // Append with a closing curly brace.

    return ss.str(); // Return the final combined string.
}


/**
 * @brief Combines a vector of strings into a single string, formatted as a
 * comma-separated list enclosed in curly braces, with each item given quotes.
 *
 * @param stringVector The vector of strings to combine.
 * @return A single string representing the combined vector (e.g., "{"item1", "item2", "item3"}").
 */

    inline std::string combineVectorToQuotedString(const ROOT::RVec<std::string>& parts) {
   std::string result = "{";
   for (const auto& part : parts) {
     result += "\"" + part + "\",";
   }
   if (!parts.empty()) result.pop_back(); // remove trailing comma
   result += "}";
   return result;
 }
    
/**
 * @brief Combines a vector of strings into a single string, formatted as a
 * comma-separated list enclosed in curly braces.
 *
 * @param stringVector The vector of strings to combine.
 * @return A single string representing the combined vector (e.g., "{item1, item2, item3}").
 */
  template<typename T>
  inline std::string combineAnyVectorToString(const ROOT::RVec<T>& vec){
      
      std::string result = "{";
      
      //go through vec and make string of each element
      for (auto iter:vec){
	auto p = std::to_string(iter);
	result.append(p);
	result.append(",");
      }

      //remove last "," easily
      if(!result.empty()){
	result.pop_back();
      }
      //close the curly brackets {}
      result.append("}");
      return result;
    }
    
// Helper function to convert a string to lowercase
// This is useful for case-insensitive comparisons
std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

    /**
     * @brief Filters a vector of strings, returning only those that contain the specified substring.
     *
     * @param stringList A constant reference to the vector of strings to be filtered.
     * @param substring The substring to search for within each string.
     * @param caseSensitive If true, the search is case-sensitive. If false, the search is case-insensitive.
     * @return A new vector containing only the strings that contain the substring.
     */
    ROOT::RVec<std::string> filterStrings(
					   const ROOT::RVec<std::string>& stringList,
					   const std::string& substring,
					   bool caseSensitive = true)
    {
      ROOT::RVec<std::string> filteredList;

      if (substring.empty()) { // If substring is empty, all strings contain it (or none, depending on interpretation).
	// Here, we'll return all strings if the substring is empty.
        return stringList;
      }

      if (caseSensitive) {
        for (const std::string& s : stringList) {
	  // std::string::find returns std::string::npos if the substring is not found
	  if (s.find(substring) != std::string::npos) {
	    filteredList.push_back(s);
	  }
        }
      } else { // Case-insensitive search
        std::string lowerSubstring = toLower(substring);
        for (const std::string& s : stringList) {
	  std::string lowerS = toLower(s);
	  if (lowerS.find(lowerSubstring) != std::string::npos) {
	    filteredList.push_back(s);
	  }
        }
      }


      return filteredList;
    }
    /**
     * @brief Filters a ROOT::RVecS (RVec of strings), returning only those that contain the specified substring.
     * This function acts as a wrapper, converting RVecS to ROOT::RVec<std::string> for filtering,
     * and then converting the result back to RVecS.
     *
     * @param rvecS A constant reference to the ROOT::RVecS to be filtered.
     * @param substring The substring to search for within each string.
     * @param caseSensitive If true, the search is case-sensitive. If false, the search is case-insensitive.
     * @return A new ROOT::RVecS containing only the strings that contain the substring.
     */
    /* RVecS filterStrings(
			const RVecS& rvecS,
			const std::string& substring,
			bool caseSensitive = true)
    {
      // 1. Convert ROOT::RVecS to ROOT::RVec<std::string>
      // RVecs can be directly constructed from iterators, or implicitly converted
      // to ROOT::RVec in many contexts. Explicit construction is clear.
      ROOT::RVec<std::string> stdVec(rvecS.begin(), rvecS.end());

      // 2. Call the core filtering function
      auto filteredStdVec = filterStrings(stdVec, substring, caseSensitive );
      
      // 3. Convert the result back to ROOT::RVecS
      // RVec can be constructed from a ROOT::RVec.
      RVecS filteredRVecS(filteredStdVec);

      return filteredRVecS;
    }
    */
    /**
     * @brief Removes elements from the destination vector if they exist in the removal vector.
     * * The order of elements in the destination vector is preserved.
     * Uses std::unordered_set for O(1) average time lookups, resulting in faster overall execution.
     *
     * @param dest The vector to be filtered (modified in place).
     * @param keys_to_remove The vector containing the keys that should be removed from 'dest'.
     */
    void removeExistingStrings(
        ROOT::RVec<std::string>& dest, 
        const ROOT::RVec<std::string>& keys_to_remove) 
    {
        // 1. Create a hash set for fast O(1) average-time lookups.
        // This is much faster than repeatedly using std::find (which is O(N) per call).
        std::unordered_set<std::string> removal_set(
            keys_to_remove.begin(), 
            keys_to_remove.end());

        // 2. Use the Erase-Remove Idiom with a lambda to filter the destination vector.
        // std::remove_if moves elements to be kept to the front of the vector.
        dest.erase(
            std::remove_if(dest.begin(), dest.end(),
                // The lambda checks if the element 's' exists in our fast removal_set.
                [&removal_set](const std::string& s) {
                    return removal_set.count(s) > 0;
                }),
            dest.end()); // std::erase removes the elements moved to the end.
    }
  
    /**
     * @brief Returns a vector containing all common entries (intersection) of two input vectors.
     * * The input vectors are sorted internally to allow for efficient set intersection.
     *
     * @param vec1 The first vector<string>.
     * @param vec2 The second vector<string>.
     * @return ROOT::RVec<std::string> A new vector containing the unique common elements.
     */
    ROOT::RVec<std::string> getCommonStrings(
					      const ROOT::RVec<std::string>& vec1,
					      const ROOT::RVec<std::string>& vec2) 
    {
      // 1. Create mutable copies for sorting (required by std::set_intersection).
      ROOT::RVec<std::string> sorted_vec1 = vec1;
      ROOT::RVec<std::string> sorted_vec2 = vec2;

      // 2. Sort both input vectors. This is the necessary prerequisite for set_intersection.
      // Complexity: O(N log N)
      std::sort(sorted_vec1.begin(), sorted_vec1.end());
      std::sort(sorted_vec2.begin(), sorted_vec2.end());

      // 3. Determine the maximum possible size for the result vector and reserve space.
      // The intersection size cannot be larger than the smallest input vector.
      size_t max_size = std::min(sorted_vec1.size(), sorted_vec2.size());
      ROOT::RVec<std::string> result;
      result.reserve(max_size);

      // 4. Perform the set intersection.
      // Complexity: O(N + M) (linear time after sorting)
      std::set_intersection(
			    sorted_vec1.begin(), sorted_vec1.end(),
			    sorted_vec2.begin(), sorted_vec2.end(),
			    std::back_inserter(result)
			    );

      return result;
    }

    // using ColumnNames_t_Std = ROOT::RVec<std::string>;
    //using NestedColumnNames_t = ROOT::RVec<ROOT::RVec<std::string>>;
    using ColumnNames_t_Std = ROOT::RVec<std::string>;
    using NestedColumnNames_t = ROOT::RVec<ROOT::RVec<std::string>>;

    /**
     * @brief Flattens a 2D vector structure (vector<vector<string>>) into a single 1D vector.
     * * @param nested_names The NestedColumnNames_t structure (vector<vector<string>>) to flatten.
     * @return ColumnNames_t_Std A single ROOT::RVec<std::string> containing all entries.
     */
    ColumnNames_t_Std flattenColumnNames(const NestedColumnNames_t& nested_names) {
        
        ColumnNames_t_Std flat_names;
        size_t total_size = 0;

        // 1. Calculate total size (Optional, but good for reserving memory)
        for (const auto& inner_vec : nested_names) {
            total_size += inner_vec.size();
        }
        flat_names.reserve(total_size);

        // 2. Insert the contents of each inner vector into the output vector.
        for (const auto& inner_vec : nested_names) {
            // Use ROOT::RVec::insert for efficient appending of the range
            flat_names.insert(flat_names.end(), inner_vec.begin(), inner_vec.end());
        }

        return flat_names;
    }

    /**
     * @brief Appends a specified suffix string to every element in a vector of strings.
     * * @param input_vector The vector<string> whose elements will be modified.
     * @param suffix The string to append to the end of every element.
     * @return ROOT::RVec<std::string> A new vector containing the modified strings.
     */
    ROOT::RVec<std::string> appendSuffixToAll(
        const ROOT::RVec<std::string>& input_vector,
        const std::string& suffix) 
    {
        ROOT::RVec<std::string> result;
        result.reserve(input_vector.size()); // Optimize memory allocation
        
        // Use std::transform to apply the lambda function to every element.
        std::transform(input_vector.begin(), input_vector.end(), 
                       std::back_inserter(result),
                       [&suffix](const std::string& s) {
                           return s + suffix;
                       });
                       
        return result;
    }
  /**
     * @brief Appends a specified suffix string to every element in a vector of strings.
     * * @param input_vector The vector<string> whose elements will be modified.
     * @param suffix The string to append to the end of every element.
     * @return ROOT::RVec<std::string> A new vector containing the modified strings.
     */
    ROOT::RVec<std::string> prependToAll(
        const ROOT::RVec<std::string>& input_vector,
        const std::string& suffix) 
    {
        ROOT::RVec<std::string> result;
        result.reserve(input_vector.size()); // Optimize memory allocation
        
        // Use std::transform to apply the lambda function to every element.
        std::transform(input_vector.begin(), input_vector.end(), 
                       std::back_inserter(result),
                       [&suffix](const std::string& s) {
                           return suffix + s;
                       });
                       
        return result;
    }
  
  /**
     * @brief Generates a JIT string to pack multiple RDataFrame columns into a single RVec<RVec>.
     * * Example Output: "rad::util::PackColumns(rec_px, rec_py, rec_pz)"
     * * @param cols The list of column names to pack.
     * @return std::string The function call string.
     */
    inline std::string createPackVectorString(const ROOT::RVec<std::string>& cols) {
        if (cols.empty()) {
            return ""; // Should be handled by caller, but safe default
        }
        
        // We use a variadic C++ helper function 'rad::util::PackColumns'
        // which must be available to the Cling JIT.
        return createFunctionCallStringFromVec("rad::util::PackColumns", cols);
    }

    // Helper to replace invalid chars (., :, etc.) with underscores
    inline std::string MakeValidName(std::string s) {
      std::replace(s.begin(), s.end(), '.', '_');
      std::replace(s.begin(), s.end(), ':', '_');
      std::replace(s.begin(), s.end(), '/', '_');
      return s;
    };
    
    //////////////////////////////////////////////////////////////////
    std::string ColumnsToString(const ROOT::RDF::ColumnNames_t &cols) {
      if(cols.empty()==true) return "{}";
      
      string toString ="{";
      for(const auto& p:cols){
	toString=(toString + p + ",");
      }
      toString.pop_back(); //remove last ,
      toString+='}';
      return toString;
    }
    //////////////////////////////////////////////////////////////////
    /**
     * @brief Joins column names with commas but WITHOUT surrounding braces.
     * Essential for passing arguments to C++ variadic templates in JIT code.
     * e.g. returns "col1,col2,col3" instead of "{col1,col2,col3}"
     */
    inline std::string ColumnsToStringNoBraces(const ROOT::RDF::ColumnNames_t &cols) {
      if(cols.empty()) return "";
      
      std::string toString = "";
      for(const auto& p : cols){
         toString += (p + ",");
      }
      if (!toString.empty()) toString.pop_back(); // remove last comma
       
      return toString;
    }
  }
}
