#pragma once
#include "Constants.h"
#include <ROOT/RVec.hxx>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <limits>

namespace rad{
  namespace util{
    
    using ROOT::RVecI;
    
    // =========================================================================
    //  Index Finding Utilities
    // =========================================================================

    /**
     * @brief Finds the index of the next occurrence of a value in a vector.
     * @param vec The vector to search.
     * @param value The value to find.
     * @param vecstart Optional pointer to the start of the search window (default: beginning).
     * @return size_t Index of the value relative to the vector begin, or vec.size() if not found.
     */
    inline size_t findIndex(const ROOT::RVecI& vec, const int value, const int* vecstart=nullptr) {
      if(vecstart==nullptr) vecstart = vec.begin();
      return std::distance(vec.begin(), std::find(vecstart, vec.end(), value));
    }

     /**
     * @brief Find the Nth occurrence of a value in a vector.
     * @param n_occurance The 1-based occurrence number (1 finds the first instance).
     * @param value The value to search for.
     * @return int Index of the Nth occurrence, or -1 if not found/invalid.
     */
    inline int findNthIndex(const ROOT::RVecI& vec, int n_occurance, const int value){
      if(n_occurance <= 0) return -1;
      
      auto currIdx = findIndex(vec, value);
      auto vsize = vec.size();
      
      if(currIdx == vsize) return -1;

      // Iterate to find subsequent occurrences
      n_occurance--;
      while((n_occurance--) > 0){
        // Search starting from the element AFTER the current find
        currIdx = findIndex(vec, value, &vec[currIdx+1]);
        if(currIdx == vsize) return -1;
      }
      return currIdx;
    }
  
    // =========================================================================
    //  Vector Manipulation & Reordering
    // =========================================================================

    /**
     * @brief Creates a Reverse Index lookup (Argsort for integer sequences).
     * If the input is `[3, 0, 2, 1]`, the output is `[1, 3, 2, 0]`.
     * Useful for mapping "Which track is Particle X?" back to "What is Particle X's track index?".
     */
    template<typename T>
    ROOT::VecOps::RVec<T> ReverseIndex(const ROOT::VecOps::RVec<T>& vec){
      if(vec.empty()) return {};
      
      // Allocate result big enough to hold the largest index found
      auto max_val = ROOT::VecOps::Max(vec);
      ROOT::VecOps::RVec<T> result(max_val + 1, -1); // Initialize with -1
      
      T entry = 0;
      for(auto idx : vec){
        if(idx >= 0) result[idx] = entry;
        ++entry;
      }
      return result;
    }
   
   /**
     * @brief Rearranges `vec` according to the indices in `imatch`.
     * Equivalent to `ROOT::VecOps::Take`, but explicit name helps clarify intent.
     * @param vec The source data vector.
     * @param imatch The vector of indices to select.
     */
    template<typename T, typename Ti>
    ROOT::VecOps::RVec<T> Rearrange(const ROOT::VecOps::RVec<T>& vec, const ROOT::RVec<Ti>& imatch){
      return ROOT::VecOps::Take(vec, imatch);
    }

    /**
     * @brief Complex reordering: Moves entries from `iorder0` mapping to `iorder1` mapping.
     * Creates a new vector of size `n`. 
     * Missing elements are set to `InvalidEntry`.
     */
    template<typename T, typename T0, typename T1>
    ROOT::VecOps::RVec<T> Reorder(const ROOT::VecOps::RVec<T>& vec0, 
                                  const ROOT::RVec<T0>& iorder0, 
                                  const ROOT::RVec<T1>& iorder1, 
                                  const size_t n) {
      
      ROOT::VecOps::RVec<T> vec1(n, rad::consts::InvalidEntry<T>()); 
      size_t target_size = iorder1.size();
      
      for(size_t i=0; i<target_size; ++i){
        if(iorder1[i] < 0) continue;
        if(iorder0[i] < 0) continue;
        // Map data from Source Order -> Target Order
        vec1[iorder1[i]] = vec0[iorder0[i]];
      }
      return vec1;
    }
    
    /** @brief Returns [0, 1, 2, ... size-1]. */
    template <typename T>
    ROOT::RVec<T> Enumerate(size_t size) {
      ROOT::RVec<T> ret;
      ret.reserve(size);
      for (auto i = 0UL; i < size; ++i) ret.emplace_back(static_cast<T>(i));
      return ret;
    }

    /** @brief Returns [istart, istart+1, ... istart+size-1]. */
    inline ROOT::RVecI EnumerateIndicesFrom(int istart, size_t size) {
        ROOT::RVecI temp_vec(size); 
        std::iota(temp_vec.begin(), temp_vec.end(), istart);
        return temp_vec;
    }

    // =========================================================================
    //  Statistical & Inspection Helpers
    // =========================================================================
    
    template<typename T>
    bool Contains(const ROOT::VecOps::RVec<T>& vec, const T& val){   
       return std::find(vec.begin(), vec.end(), val) != vec.end();
    }
    
    template<typename T>
    size_t Count(const ROOT::VecOps::RVec<T>& vec, const T& val){
      return std::count(vec.begin(), vec.end(), val);
    }
    
    /** @brief Safe access to first element. Returns InvalidEntry if empty. */
    template<typename T>
    T First(const ROOT::RVec<T>& values){
      return values.empty() ? rad::consts::InvalidEntry<T>() : values[0];
    }
    
    template<typename T>
    T Mean(const ROOT::RVec<T>& values){
      return values.empty() ? rad::consts::InvalidEntry<T>() : ROOT::VecOps::Mean(values);
    }
    
    template<typename T>
    T Sum(const ROOT::RVec<T>& values){
      return values.empty() ? rad::consts::InvalidEntry<T>() : ROOT::VecOps::Sum(values);
    }
    
    /** @brief Returns the Max value in a vector, or InvalidIndex if empty. */
    inline int MaxIndex(const ROOT::RVecI& indices) {
        if (indices.empty()) return consts::InvalidIndex();
        return ROOT::VecOps::Max(indices);
    }

    /** @brief Returns the Max value in a nested vector of indices. */
    inline int MaxIndex(const ROOT::RVec<ROOT::RVecI>& nested_indices) {
      if (nested_indices.empty()) return consts::InvalidIndex();
      
      int global_max = std::numeric_limits<int>::min(); 
      bool found = false;

      for (const auto& indices : nested_indices) {
        if (!indices.empty()) {
          int current_max = ROOT::VecOps::Max(indices);
          if (current_max > global_max) global_max = current_max;
          found = true;
        }
      }
      return found ? global_max : -1;
    }

    // =========================================================================
    //  Grouping & Concatenation
    // =========================================================================

    /**
     * @brief Groups a set of individual values into a single RVec<T>.
     * * This variadic function is useful for creating vectors on the fly inside 
     * RDataFrame `Define` calls (e.g., grouping index columns or creating fixed 
     * lists of PIDs).
     * * Usage: `rad::util::Group<int>(idx1, idx2, idx3)`
     * * @tparam T The target type stored in the resulting RVec.
     * @tparam ColumnValues Variadic types of the input arguments (must be convertible to T).
     * @param values The parameter pack of values to be grouped.
     * @return ROOT::RVec<T> A vector containing the input values cast to T.
     */
    template<typename T, typename... ColumnValues>
    ROOT::RVec<T> Group(ColumnValues... values){
      ROOT::RVec<T> group{static_cast<T>(values)...};
      return group;
    }

    /**
     * @brief Variadic Grouping: Packs multiple RVecs into a single RVec<RVec>.
     * * This is essential for `ConfigReaction::setGroupParticles`. It takes N columns
     * (e.g., ele_indices, pos_indices) and bundles them into a single column 
     * (Meson_Candidates) for the ParticleCreator.
     * @tparam T The value type (e.g., Indices_t). All inputs must be compatible.
     * @param first The first vector.
     * @param args The rest of the vectors.
     * @return ROOT::RVec<T> A vector containing the inputs as elements.
     */
    template <typename T, typename... Args>
    ROOT::RVec<T> Group(const T& first, const Args&... args) {
        ROOT::RVec<T> result;
        result.reserve(1 + sizeof...(args));
        result.push_back(first);
        (result.push_back(args), ...); // Fold expression to push all args
        return result;
    }

    /**
     * @brief Alias for Group, often used in JIT string generation for kernels.
     */
    template <typename T, typename... Args>
    ROOT::RVec<ROOT::RVec<T>> PackColumns(const ROOT::RVec<T>& first, const Args&... args) {
        return Group<ROOT::RVec<T>>(first, args...);
    }

    // --- Concatenation Helpers ---

    template <typename T>
    ROOT::RVec<T> Concatenate(const ROOT::RVec<T> &v0, const ROOT::RVec<T> &v1) {
       ROOT::RVec<T> res;
       res.reserve(v0.size() + v1.size());
       res.insert(res.end(), v0.begin(), v0.end());
       res.insert(res.end(), v1.begin(), v1.end());
       return res;
    }

    // Variadic Recursion helpers for Concatenate
    template<typename ResT, typename HeadT, typename... TailT>
    void append_to(ResT& res, const HeadT& head, const TailT&... tail) {
        for(const auto& val : head) res.push_back(static_cast<typename ResT::value_type>(val));
        if constexpr (sizeof...(tail) > 0) append_to(res, tail...);
    }
    
    inline size_t get_total_size() { return 0; }
    template<typename T, typename... Args>
    size_t get_total_size(const T& v, const Args&... args) {
        return v.size() + get_total_size(args...);
    }

    /**
     * @brief Unified Variadic Concatenate.
     * Usage: Concatenate<double>(vec1, vec2, vec3).
     */
    template <typename ResT = double, typename... Args>
    ROOT::RVec<ResT> Concatenate(const Args&... args) {
       ROOT::RVec<ResT> res;
       res.reserve(get_total_size(args...));
       append_to(res, args...);
       return res;
    }
    
    /**
     *  @brief Truncate vec at n
     */
    template<typename T>
    ROOT::VecOps::RVec<T> Truncate(const ROOT::VecOps::RVec<T>& vec,const size_t size){
      // std::cout<<"Truncate "<<vec<<size<<std::endl;
      ROOT::RVec<T> ret;
      ret.reserve(size);
      for (auto i = 0UL; i < size; ++i) {
	ret.emplace_back(vec[i]);
      }
      return ret;
    }
  }//util
}//rad

// #pragma once
// #include "Constants.h"
// #include <ROOT/RVec.hxx>
// #include <algorithm>
// #include <iostream>

// namespace rad{
//   namespace util{
    
//      //! Code simplifications
  
//     using ROOT::RVecI;
    
//     /**
//      * Find the next occurance of value in vec after vecstart
//      */
//     inline size_t findIndex(const ROOT::RVecI& vec,const int value,const int* vecstart=nullptr) {
//       if(vecstart==nullptr) vecstart = vec.begin();
//       return distance(vec.begin(), find(vecstart, vec.end(), value ));
//     }

//      /**
//      * Find the nth occurance of value in vec. Note n_occurance should be 1 for first occurance
//      */
//     inline int findNthIndex(const ROOT::RVecI& vec, int n_occurance,const int value){
//       if(n_occurance == 0 ) return -1;
//       //find first occurance
//       auto currIdx = findIndex(vec,value);
//       auto vsize= vec.size();
//       //check if we reached the end of the vector
//       if(currIdx==vsize) return -1;

//       //Look for the n_occurance case
//       n_occurance--;
//       while( (n_occurance--) >0){
// 	currIdx=findIndex(vec,value,&vec[currIdx+1]);
// 	//check if we reached the end of the vector
// 	if(currIdx==vsize) return -1;
//       }
//       return currIdx;
//     }
  
//     /**
//      * Rearrange vec in order of its own elements
//      * e.g. [1,3,2,0] -> [3,0,2,1]
//      */
//     template<typename T>
//     ROOT::VecOps::RVec<T> ReverseIndex(const ROOT::VecOps::RVec<T>& vec){
//       ROOT::VecOps::RVec<T> result(vec.size());
//       T entry = 0;
//       for(auto idx:vec){
// 	result[idx] = entry;
// 	++entry;
//       }
//       return result;
//     }
   
//    /**
//      * Rearrange vec in order of elements in imatch
//      */
//     template<typename T,typename Ti>
//     ROOT::VecOps::RVec<T> Rearrange(const ROOT::VecOps::RVec<T>& vec,const ROOT::RVec<Ti>& imatch){
//       //std::cout<<"Rearrange "<<vec<<" "<<imatch<<" "<<ROOT::VecOps::Take(vec,imatch)<<std::endl;
//       return ROOT::VecOps::Take(vec,imatch);
//     }
//     /**
//      * Reorder vec0 moving entries at iorder0 to iorder1, 
//      * with missing elements set to InvalidEntry()
//      */
//     template<typename T,typename T0,typename T1>
//     //ROOT::VecOps::RVec<T> Reorder(const ROOT::VecOps::RVec<T>& vec0,const ROOT::RVecU& iorder0,const ROOT::RVecU& iorder1,const size_t n){
//     ROOT::VecOps::RVec<T> Reorder(const ROOT::VecOps::RVec<T>& vec0,const ROOT::RVec<T0>& iorder0,const ROOT::RVec<T1>& iorder1,const size_t n){

//       //create new vector size  n
//       ROOT::VecOps::RVec<T> vec1(n,rad::consts::InvalidEntry<T>()); //create new vector size n=iorder1.size
//       size_t target_size = iorder1.size();
//       //need to loop over order0
//       for(size_t i=0;i<target_size;++i){
// 	//add value of vec0 at iorder1
// 	if(iorder1[i]<0) continue;
// 	if(iorder0[i]<0) continue;
// 	vec1[iorder1[i]]=vec0[iorder0[i]]; //give vec1[iorder1] value of vec0[iorder0]
//       }
//       // std::cout<<"reorder done "<<vec1<<std::endl;
//       return vec1;
//     }
    
//     template <typename T>
//      ROOT::RVec<T> Enumerate(size_t size)
//      {
//       ROOT::RVec<T> ret;
//       ret.reserve(size);
//       for (auto i = 0UL; i < size; ++i) {
// 	ret.emplace_back(i);
//       }
//       return ret;
//     }
//     /**
//      * @brief Initializes an Indices_t (ROOT::RVecI) with a sequence of enumerated integers.
//      * * The vector is filled with values: [istart, istart+1, ..., istart + size - 1].
//      *
//      * @param istart The starting integer value of the sequence.
//      * @param size The number of elements to generate in the sequence.
//      * @return Indices_t The initialized RVecI containing the enumerated sequence.
//      */
//     Indices_t EnumerateIndicesFrom(int istart, size_t size) {
        
//         // 1. Create a temporary ROOT::RVec to hold the sequence.
//         // We use ROOT::RVec because std::iota is defined for it.
//         ROOT::RVec<int> temp_vec(size); 
        
//         // 2. Fill the temporary vector with the sequence.
//         // std::iota fills the range [temp_vec.begin(), temp_vec.end()) 
//         // with sequentially increasing values, starting at istart.
//         std::iota(temp_vec.begin(), temp_vec.end(), istart);

//         // 3. Construct and return the Indices_t (ROOT::RVecI).
//         // RVec has a constructor that efficiently copies/moves data from a ROOT::RVec.
//         return Indices_t(std::move(temp_vec));
//     }
//    /**
//      * Truncate vec at n
//      */
//     template<typename T>
//     ROOT::VecOps::RVec<T> Truncate(const ROOT::VecOps::RVec<T>& vec,const size_t size){
//       // std::cout<<"Truncate "<<vec<<size<<std::endl;
//       ROOT::RVec<T> ret;
//       ret.reserve(size);
//       for (auto i = 0UL; i < size; ++i) {
// 	ret.emplace_back(vec[i]);
//       }
//       return ret;
//     }
//     /**
//      * Check if vec contains element == val
//      */
//      template<typename T>
//      bool Contains(const ROOT::VecOps::RVec<T>& vec,const T& val){   
//        return std::find(vec.begin(), vec.end(), val) != vec.end() ? true :false;
//      }
//     /**
//      * count the instances of val in vec
//      */
//     template<typename T>
//     size_t Count(const ROOT::VecOps::RVec<T>& vec,const T& val){
//       return std::count(vec.begin(), vec.end(), val);
//     }
//     /**
//      * increment the values in vec
//      */
//     template<typename T>
//     void Increment(ROOT::VecOps::RVec<T>& vec, long off){
//       if(vec.size()==0) return ;
//       std::for_each(vec.begin(), vec.end(), [&off](T &val) { val+=off; });
     
//     }
//     /**
//      * increment the values in vec
//      */
//     template<typename T>
//     void Append(const T& val, ROOT::VecOps::RVec<T>& vec){
//       vec.push_back(val);
//     }
//     /**
//      * increment the values in vec
//      */
//     template<typename T>
//     ROOT::VecOps::RVec<T> AppendToCopy(const T& val, const ROOT::VecOps::RVec<T>& vec){
//       auto result = vec;
//       result.push_back(val);
//       return result;
//     }
//     /**
//      * increment the values in vec and return a copy
//      */
//     template<typename T>
//     ROOT::VecOps::RVec<T> IncrementCopy(ROOT::VecOps::RVec<T> vec, long off){
//       if(vec.size()==0) return vec;
//       std::for_each(vec.begin(), vec.end(), [&off](T &val) { val+=off; });
//       return vec;
//     }

//    template <typename T>
//    ROOT::RVec<T> Concatenate(const ROOT::RVec<T> &v0, const ROOT::RVec<T> &v1)
//    {
  
//      ROOT::RVec<T> res;
//      res.reserve(v0.size() + v1.size());
//      std::copy(v0.begin(), v0.end(), std::back_inserter(res));
//      std::copy(v1.begin(), v1.end(), std::back_inserter(res));
//      return res;
//    }
//    template <typename T>
//    ROOT::RVec<T> Concatenate(const ROOT::RVec<T> &v0, const ROOT::RVec<T> &v1,const ROOT::RVec<T> &v2)
//    {
//      ROOT::RVec<T> res;
//      res.reserve(v0.size() + v1.size() + v2.size() );
//      std::copy(v0.begin(), v0.end(), std::back_inserter(res));
//      std::copy(v1.begin(), v1.end(), std::back_inserter(res));
//      std::copy(v2.begin(), v2.end(), std::back_inserter(res));
//      return res;
//    }
//     template <typename T>
//    ROOT::RVec<T> Concatenate(const ROOT::RVec<T> &v0, const ROOT::RVec<T> &v1,const ROOT::RVec<T> &v2, const ROOT::RVec<T> &v3)
//    {

//      ROOT::RVec<T> res;
//      res.reserve(v0.size() + v1.size() + v2.size() + v3.size());
//      std::copy(v0.begin(), v0.end(), std::back_inserter(res));
//      std::copy(v1.begin(), v1.end(), std::back_inserter(res));
//      std::copy(v2.begin(), v2.end(), std::back_inserter(res));
//      std::copy(v3.begin(), v3.end(), std::back_inserter(res));
//      return res;
//    }
//     template <typename T>
//     ROOT::VecOps::RVec<ROOT::RVec<T>> Concatenate(const ROOT::VecOps::RVec<ROOT::RVec<T>> vec_all)
//     {
//      const size_t nelements = vec_all.size();
//      size_t length=1;
//      for(auto i =0UL ; i<nelements;++i){
//        length+=vec_all[i].size();
//      }
 
//      ROOT::RVec<T> res(length);

//      size_t entry=0;
//      for(auto i =0UL ; i<nelements;++i){
//        size_t curr_size=vec_all[i].size();
//        for(auto j =0UL ; j<curr_size;++j){
// 	 res[entry++]=vec_all[i][j];
//        }
//      }

//      return res;
//     }

//     // Helper to calculate total size (C++11/14 compatible)
//     // (If you are on C++17 you can use fold expressions, but this is safer for ROOT JIT)
//   // Base case for recursion
//     inline size_t get_total_size() { return 0; }

//     // template<typename T>
//     // size_t get_size(const T& v) { return v.size(); }

//     template<typename T, typename... Args>
//     size_t get_total_size(const T& v, const Args&... args) {
//         return v.size() + get_total_size(args...);
//     }
  

//     // Helper to append vectors recursively
//     template<typename ResT>
//     void append_to(ResT&) {} // Base case

//     template<typename ResT, typename HeadT, typename... TailT>
//     void append_to(ResT& res, const HeadT& head, const TailT&... tail) {
//         // static_cast handles the float -> double conversion efficiently
//         for(const auto& val : head) res.push_back(static_cast<typename ResT::value_type>(val));
//         append_to(res, tail...);
//     }

//     /**
//      * @brief Unified Concatenate function.
//      * * Usage:
//      * - Concatenate(a, b)        -> Returns RVec<double> (Default)
//      * - Concatenate<int>(a, b)   -> Returns RVec<int>
//      * - Concatenate<float>(a, b) -> Returns RVec<float>
//      */
//     template <typename ResT = double, typename... Args>
//     ROOT::RVec<ResT> Concatenate(const Args&... args) {
//        ROOT::RVec<ResT> res;
//        res.reserve(get_total_size(args...));
       
//        // append_to handles the static_cast automatically
//        append_to(res, args...);
       
//        return res;
//     }
   
    
//     /** combine v0 and v1 into 1 vector with elements
//      * which are vectors of length 2.
//      */
//     template <typename T>
//   ROOT::RVec<T> CombineElements(const ROOT::RVec<T> &v0, const ROOT::RVec<T> &v1)
//    {
//      const size_t length = v0.size();
//      ROOT::VecOps::RVec<ROOT::RVec<T>> res(length, ROOT::RVec<T>(2));
//      for(auto i =0UL ; i<length;++i){
//        res[i][0]=v0[i];
//        res[i][1]=v1[i];
//      }
//      return res;
//    }
//     /** combine the n vectors in vec_all into 1 vector with elements
//      * which are vectors of length n.
//      */
//     template <typename T>
//     ROOT::VecOps::RVec<ROOT::RVec<T>> CombineElements(const ROOT::VecOps::RVec<ROOT::RVec<T>> vec_all)
//    {
//      const size_t nelements = vec_all.size();
//      const size_t length = vec_all[0].size();
//      ROOT::VecOps::RVec<ROOT::RVec<T>> res(length, ROOT::RVec<T>(nelements));
//      for(auto i =0UL ; i<length;++i){
//        for(auto j =0UL ; j<nelements;++j){
// 	 res[i][j]=vec_all[j][i];
//        }
//      }
//      return res;
//    }
  
 
//     /**
//      * Return the positions where c is true
//      */
//     ROOT::RVecU PositionsWhere(const ROOT::RVecB& c)
//     {
//       const size_t length = c.size();
//       ROOT::RVecU r(Count(c,true));//new vector sized to non-zero elements
//       uint entry=0;
//       for (auto i=0UL; i<length; ++i) {
// 	if(c[i]==true)
// 	  r[entry++]=i;
//       }
//       return r;
//     }

//     /**
//     * Group a set of types into an RVec
//     */
//     template<typename T, typename... ColumnValues>
//     ROOT::RVec<T> Group(ColumnValues... values){
//       ROOT::RVec<T> group{static_cast<T>(values)...};
//       return group;
//     }
    
//     /**
//     * Return first entry in RVec
//     */
//     template<typename T>
//     T First(const ROOT::RVec<T>& values){
//       return values.empty() ? rad::consts::InvalidEntry<T>() : values[0];
//     }
    
//     /**
//     * Return mean value of RVec
//     */
//     template<typename T>
//       T Mean(const ROOT::RVec<T>& values){
//       return values.empty() ? rad::consts::InvalidEntry<T>() : ROOT::VecOps::Mean(values);
//     }
    
//     /**
//     * Return sum value of RVec
//     */
//     template<typename T>
//       T Sum(const ROOT::RVec<T>& values){
//       return values.empty() ? rad::consts::InvalidEntry<T>() : ROOT::VecOps::Sum(values);
//     }
    
//     /**
//      * @brief Returns the highest integer value within a single Indices_t (RVecI).
//      * @param indices The Indices_t vector to search.
//      * @return int The maximum index value. Returns consts::InvalidIndex() if the vector is empty.
//      */
//     inline int MaxIndex(const Indices_t& indices) {
//         if (indices.empty()) {
//             // Assuming consts::InvalidIndex() is defined and indicates an invalid/empty state.
// 	  return consts::InvalidIndex();
//         }
//         // Use ROOT's optimized vectorized maximum calculation.
//         return ROOT::VecOps::Max(indices);
//     }
//     /**
//      * @brief Returns the highest integer value across all vectors within an RVecIndices structure.
//      * @param nested_indices The RVecIndices structure to search (RVec<RVecI>).
//      * @return int The overall maximum index value. Returns consts::InvalidIndex() if structure is empty.
//      */
//     inline int MaxIndex(const RVecIndices& nested_indices) {
//       if (nested_indices.empty()) {
// 	return consts::InvalidIndex();
//       }

//       int global_max = std::numeric_limits<int>::min(); // Start with the smallest possible integer
//       bool found_any_element = false;

//       for (const auto& indices : nested_indices) {
// 	if (!indices.empty()) {
// 	  // Find the max within the current inner vector
// 	  int current_max = ROOT::VecOps::Max(indices);
                
// 	  // Update the global maximum
// 	  if (current_max > global_max) {
// 	    global_max = current_max;
// 	  }
// 	  found_any_element = true;
// 	}
//       }

//       if (!found_any_element) {
// 	return -1; // Or consts::InvalidIndex()
//       }
//       return global_max;
//     }
//      /**
//      * @brief Variadic helper to pack multiple RVecs into a single RVec<RVec>.
//      * Called by RDataFrame via JIT string.
//      * * @tparam T The value type (e.g., double). All inputs must be RVec<T>.
//      * @tparam Args Variadic pack of RVec<T>.
//      */
//     template <typename T, typename... Args>
//     ROOT::RVec<ROOT::RVec<T>> PackColumns(const ROOT::RVec<T>& first, const Args&... args) {
//         // 1. Create result vector
//         // Size = 1 (first) + sizeof...(args)
//         ROOT::RVec<ROOT::RVec<T>> result;
//         result.reserve(1 + sizeof...(args));

//         // 2. Add first element
//         result.push_back(first);

//         // 3. Fold expression to push back the rest
//         (result.push_back(args), ...);

//         return result;
//     }
    
//   }//util
  
// }//rad
