#pragma once
#include <limits>

namespace rad{
  namespace consts{

    /**
     * Useful physics
     */
    constexpr double M_ele() { return 0.000510999;}
    constexpr double M_pro() { return 0.93827210;}
    constexpr double M_neu() { return 0.93956540;}
    constexpr double M_deut() { return 1.8756129;}
    constexpr double M_pi0() { return 0.13497680;}
    constexpr double M_pi() { return 0.13957040;}
    constexpr double M_K() { return 0.49367700;}
    constexpr double M_K0() { return 0.49761100;}
    constexpr double M_Jpsi() { return 3.0969000;}

    /**
     *    InvalidEntry
     */
    // Generic template declaration
    template<typename T>
    constexpr T InvalidEntry();

    template<typename T>
    constexpr T InvalidEntry() {
      return std::numeric_limits<T>::quiet_NaN();
    }
  
    // Specialization for int
    template<>
    constexpr int InvalidEntry<int>() {
      return std::numeric_limits<int>::max();
    }

    // Specialization for unsigned int
    template<>
    constexpr unsigned int InvalidEntry<unsigned int>() {
      return std::numeric_limits<unsigned int>::max();
    }

    template<typename T>
    inline bool IsInvalidEntry(const T& entry){
      return entry==InvalidEntry<T>();
    }
    // constexpr double InvalidEntry(){return std::numeric_limits<double>::quiet_NaN();};
    constexpr int InvalidIndex(){return -1;};


  }
}
