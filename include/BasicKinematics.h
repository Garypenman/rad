#pragma once

#include "Indicing.h"
#include "Constants.h"
#include "CommonDefines.h"

// Note: Assuming LorentzVector, PxPyPzMVector, XYZVector, ResultType_t, RVecResultType are defined in CommonDefines.h
// Note: Assuming Indices_t, RVecIndices are defined in Indicing.h or CommonDefines.h

void BasicKinematics(){}

namespace rad {

  ///\brief Helper functions and functors for RDF processing
  ///       combining momentum components into 4-vectors
  using ROOT::Math::VectorUtil::boost;

  //---------------------- 4-Vector Operations ----------------------

  /**
   * @brief Returns the 4-vector of a particle specified by a single index.
   * @tparam Tp Type of momentum components (e.g., RVec<float>).
   * @tparam Te Type of energy component (e.g., RVec<float>).
   * @param idx The index of the particle in the component vectors.
   * @param px The RVec of x-momentum components.
   * @param py The RVec of y-momentum components.
   * @param pz The RVec of z-momentum components.
   * @param e The RVec of energy components.
   * @return LorentzVector The four-momentum vector (Px, Py, Pz, E).
   */
  template<typename Tp, typename Te>
  inline LorentzVector FourVector(const uint idx, const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    return LorentzVector(px[idx], py[idx], pz[idx], e[idx]);
  }

  
  /**
   * @brief Adds the 4-vectors of multiple particles specified by indices to an existing 4-vector.
   * @tparam Tp Type of momentum components.
   * @tparam Te Type of energy component.
   * @param p4 Reference to the LorentzVector to which the sum is added.
   * @param ip The Indices_t (RVecI) containing the indices of particles to sum.
   * @param px The RVec of x-momentum components.
   * @param py The RVec of y-momentum components.
   * @param pz The RVec of z-momentum components.
   * @param e The RVec of energy components.
   */
  template<typename Tp, typename Te>
  inline void SumFourVector(LorentzVector& p4, const Indices_t &ip, const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    auto np = ip.size();
    for (size_t i = 0; i < np; ++i) {
      p4 += FourVector(ip[i], px, py, pz, e);
    }
  }

  /**
   * @brief Subtracts the 4-vectors of multiple particles specified by indices from an existing 4-vector.
   * @tparam Tp Type of momentum components.
   * @tparam Te Type of energy component.
   * @param p4 Reference to the LorentzVector from which the sum is subtracted.
   * @param ip The Indices_t (RVecI) containing the indices of particles to subtract.
   * @param px The RVec of x-momentum components.
   * @param py The RVec of y-momentum components.
   * @param pz The RVec of z-momentum components.
   * @param e The RVec of energy components.
   */
  template<typename Tp, typename Te>
  inline void SubtractFourVector(LorentzVector& p4, const Indices_t &ip, const Tp &px, const Tp &py, const Tp &pz, const Te &e) {// 
    auto np = ip.size();
    for (size_t i = 0; i <np ; ++i) {
      p4 -= FourVector(ip[i], px, py, pz, e);
    }
  }

  /**
   * @brief Returns the 4-vector resulting from the sum of particles specified by indices.
   * @tparam Tp Type of momentum components.
   * @tparam Te Type of energy component.
   * @param ipart The Indices_t (RVecI) containing the indices of particles to sum.
   * @param px The RVec of x-momentum components.
   * @param py The RVec of y-momentum components.
   * @param pz The RVec of z-momentum components.
   * @param e The RVec of energy components.
   * @return LorentzVector The summed four-momentum vector.
   */
  template<typename Tp, typename Te>
  inline LorentzVector FourVector(const Indices_t &ipart, const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    LorentzVector psum(0,0,0,0);
    SumFourVector(psum, ipart, px, py, pz, e);
    return psum;
  }

  /**
   * @brief Calculates the invariant mass of combined 4-vectors, allowing for both addition and subtraction of components.
   * * This function is the core single-combination kernel for calculating mass in combinatorial analysis.
   * * @tparam Tp Type of momentum components (e.g., RVec<double>).
   * @tparam Te Type of energy component.
   * @param indices The RVecIndices container holding all index sets. It is expected to contain at least two sets: 
   * indices[0] = Indices_t (RVecI) for particles to be ADDED (ipos).
   * indices[1] = Indices_t (RVecI) for particles to be SUBTRACTED (ineg).
   * @param px The RVec of x-momentum components.
   * @param py The RVec of y-momentum components.
   * @param pz The RVec of z-momentum components.
   * @param e The RVec of energy components.
   * @return ResultType_t The resulting mass. Returns M() of the resultant 4-vector.
   */
  template<typename Tp, typename Te>
  inline ResultType_t FourVectorMassCalc(const RVecIndices &indices, const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    const auto& ipos = indices[0];
    const auto& ineg = indices[1];
 
    
    LorentzVector psum(0,0,0,0);
    SumFourVector(psum, ipos, px, py, pz, e);
    SubtractFourVector(psum, ineg, px, py, pz, e);

    return psum.M();
  }

  /**
   * @brief Calculates the invariant mass squared of combined 4-vectors, allowing for both addition and subtraction of components.
   * * This function is the core single-combination kernel for calculating mass in combinatorial analysis.
   * * @tparam Tp Type of momentum components (e.g., RVec<double>).
   * @tparam Te Type of energy component.
   * @param indices The RVecIndices container holding all index sets. It is expected to contain at least two sets: 
   * indices[0] = Indices_t (RVecI) for particles to be ADDED (ipos).
   * indices[1] = Indices_t (RVecI) for particles to be SUBTRACTED (ineg).
   * @param px The RVec of x-momentum components.
   * @param py The RVec of y-momentum components.
   * @param pz The RVec of z-momentum components.
   * @param e The RVec of energy components.
   * @return ResultType_t The resulting mass. Returns M() of the resultant 4-vector.
   */
  template<typename Tp, typename Te>
  inline ResultType_t FourVectorMass2Calc(const RVecIndices &indices, const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    const auto& ipos = indices[0];
    const auto& ineg = indices[1];
    
    // Note: Assuming error checking for invalid indices is performed in the wrapper or upstream.
    
    LorentzVector psum(0,0,0,0);
    SumFourVector(psum, ipos, px, py, pz, e);
    SubtractFourVector(psum, ineg, px, py, pz, e);

    return psum.M2();
  }
  
  /**
   * @brief Calculates the transverse momentum of combined 4-vectors, allowing for both addition and subtraction of components.
   * * This function is the core single-combination kernel for calculating mass in combinatorial analysis.
   * * @tparam Tp Type of momentum components (e.g., RVec<double>).
   * @tparam Te Type of energy component.
   * @param indices The RVecIndices container holding all index sets. It is expected to contain at least two sets: 
   * indices[0] = Indices_t (RVecI) for particles to be ADDED (ipos).
   * indices[1] = Indices_t (RVecI) for particles to be SUBTRACTED (ineg).
   * @param px The RVec of x-momentum components.
   * @param py The RVec of y-momentum components.
   * @param pz The RVec of z-momentum components.
   * @param e The RVec of energy components.
   * @return ResultType_t The resulting mass. Returns M() of the resultant 4-vector.
   */
  template<typename Tp, typename Te>
  inline ResultType_t FourVectorPtCalc(const RVecIndices &indices, const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    const auto& ipos = indices[0];
    const auto& ineg = indices[1];
    
    // Note: Assuming error checking for invalid indices is performed in the wrapper or upstream.
    
    LorentzVector psum(0,0,0,0);
    SumFourVector(psum, ipos, px, py, pz, e);
    SubtractFourVector(psum, ineg, px, py, pz, e);

    return psum.Pt();
  }

  template<typename Tp, typename Te>
  inline ResultType_t FourVectorECalc(const RVecIndices &indices, const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    const auto& ipos = indices[0];
    const auto& ineg = indices[1];
    
    // Note: Assuming error checking for invalid indices is performed in the wrapper or upstream.
    
    LorentzVector psum(0,0,0,0);
    SumFourVector(psum, ipos, px, py, pz, e);
    SubtractFourVector(psum, ineg, px, py, pz, e);
    
    return psum.E();
    }

  //---------------------- 3-Vector Components (RVec Ops) ----------------------

  /**
   * @brief Calculates the magnitude (momentum) of 3-vectors.
   * @tparam T Type of momentum components (RVec<T>).
   * @param px The RVec of x-components.
   * @param py The RVec of y-components.
   * @param pz The RVec of z-components.
   * @return RVecResultType The RVec containing the magnitude of the 3-vector (p).
   */
  template<typename Tp>
  inline Tp ThreeVectorMag(const Tp &px, const Tp &py, const Tp &pz) {
    return sqrt(px * px + py * py + pz * pz);
  }

  /**
   * @brief prototype for a given indice, used by ApplyKinematics
   */
  template<typename Tp, typename Te>
  inline ResultType_t ThreeVectorMag(const RVecIndices &indices,const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    const Indice_t idx = indices[0][0];
    return ThreeVectorMag(px[idx],py[idx],pz[idx]);
  }
  /**
   * @brief Calculates the Theta angle of 3-vectors.
   * @tparam T Type of momentum components (RVec<T>).
   * @param px The RVec of x-components.
   * @param py The RVec of y-components.
   * @param pz The RVec of z-components.
   * @return RVecResultType The RVec containing the Theta angle.
   */
  template<typename Tp>
  inline Tp ThreeVectorTheta(const Tp &px, const Tp &py, const Tp &pz) {
    auto mag = ThreeVectorMag(px,py,pz);
    auto costh = pz/mag;
    return acos(costh);
  }
  
  /**
   * @brief prototype for a given indice, used by ApplyKinematics
   */
  template<typename Tp, typename Te>
  inline ResultType_t ThreeVectorTheta(const RVecIndices &indices,const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    const Indice_t idx = indices[0][0];
    return ThreeVectorTheta(px[idx],py[idx],pz[idx]);
  }
 

  /**
   * @brief Calculates the Phi angle of 3-vectors.
   * @tparam T Type of momentum components (RVec<T>).
   * @param px The RVec of x-components.
   * @param py The RVec of y-components.
   * @param pz The RVec of z-components.
   * @return RVecResultType The RVec containing the Phi angle.
   */
  template<typename Tp>
  inline Tp ThreeVectorPhi(const Tp &px, const Tp &py, const Tp &pz) {
    return atan2(py,px); 
  }
  /**
   * @brief prototype for a given indice, used by ApplyKinematics
   */
  template<typename Tp, typename Te>
  inline ResultType_t ThreeVectorPhi(const RVecIndices &indices,const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    const Indice_t idx = indices[0][0];
    return ThreeVectorPhi(px[idx],py[idx],pz[idx]);
  }

  /**
   * @brief Calculates the pseudorapidity (Eta) of 3-vectors.
   * @tparam T Type of momentum components (RVec<T>).
   * @param px The RVec of x-components.
   * @param py The RVec of y-components.
   * @param pz The RVec of z-components.
   * @return RVecResultType The RVec containing the pseudorapidity (Eta).
   */
  template<typename Tp>
  inline Tp ThreeVectorEta(const Tp &px, const Tp &py, const Tp &pz) {
    auto theta = ThreeVectorTheta(px,py,pz);
    return -log(tan(0.5 * theta));
  }
  /**
   * @brief prototype for a given indice, used by ApplyKinematics
   */
  template<typename Tp, typename Te>
  inline ResultType_t ThreeVectorEta(const RVecIndices &indices,const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    const Indice_t idx = indices[0][0];
    return ThreeVectorEta(px[idx],py[idx],pz[idx]);
  }

  /**
   * @brief Calculates the x-component of a vector given spherical coordinates (p, theta, phi).
   * @tparam T Type of magnitude/angle components (RVec<T>).
   * @return RVecResultType The RVec containing the x-component.
   */
  template<typename T>
  inline RVecResultType ThreeVectorX(const T &p, const T &theta, const T &phi) {
    return p*sin(theta)*cos(phi);
  }

  /**
   * @brief Calculates the y-component of a vector given spherical coordinates (p, theta, phi).
   * @tparam T Type of magnitude/angle components (RVec<T>).
   * @return RVecResultType The RVec containing the y-component.
   */
  template<typename T>
  inline RVecResultType ThreeVectorY(const T &p, const T &theta, const T &phi) {
    return p*sin(theta)*sin(phi);
  }

  /**
   * @brief Calculates the z-component of a vector given spherical coordinates (p, theta, phi).
   * @tparam T Type of magnitude/angle components (RVec<T>).
   * @return RVecResultType The RVec containing the z-component.
   */
  template<typename T>
  inline RVecResultType ThreeVectorZ(const T &p, const T &theta, const T &phi) {
    return p*cos(theta);
  }

  //---------------------- Delta Kinematics ----------------------

  /**
   * @brief Calculates the azimuthal angle difference (DeltaPhi) between two vectors.
   * * This is a single-combination kernel. It expects RVecIndices containing one set 
   * with exactly two particle indices (i0, i1).
   * * @tparam T Type of position/momentum components (RVec<T>).
   * @param indices The RVecIndices container. Expected structure: indices[0][0] = i0, indices[0][1] = i1.
   * @param px The RVec of x-components.
   * @param py The RVec of y-components.
   * @param pz The RVec of z-components.
   * @return double The DeltaPhi value. Returns InvalidEntry if indices are invalid or momentum is zero.
   */
  template<typename Tp, typename Te>
  inline ResultType_t DeltaPhi(const RVecIndices &indices, const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    // Needs single set of indices (indices[0]), containing exactly two particles
    const auto& i0 = indices[0][0];
    const auto& i1 = indices[0][1];
    
    if (i0 < 0 || i1 < 0) return rad::consts::InvalidEntry<double>();
    auto p0 = XYZVector(px[i0], py[i0], pz[i0]);
    auto p1 = XYZVector(px[i1], py[i1], pz[i1]);
    if (p0.Mag2() == 0 || p1.Mag2() == 0) return rad::consts::InvalidEntry<double>();
    
    return ROOT::Math::VectorUtil::DeltaPhi(p0, p1);
  }

  /**
   * @brief Calculates the angle difference (DeltaTheta) between two vectors.
   * * This is a single-combination kernel. It expects RVecIndices containing one set 
   * with exactly two particle indices (i0, i1).
   * * @tparam T Type of position/momentum components (RVec<T>).
   * @param indices The RVecIndices container. Expected structure: indices[0][0] = i0, indices[0][1] = i1.
   * @param px The RVec of x-components.
   * @param py The RVec of y-components.
   * @param pz The RVec of z-components.
   * @return double The DeltaTheta value. Returns InvalidEntry if indices are invalid or momentum is zero.
   */
  template<typename Tp, typename Te>
  inline ResultType_t DeltaTheta(const RVecIndices &indices, const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    const auto& i0 = indices[0][0];
    const auto& i1 = indices[0][1];
    
    if (i0 < 0 || i1 < 0) return rad::consts::InvalidEntry<double>();
    auto p0 = XYZVector(px[i0], py[i0], pz[i0]);
    auto p1 = XYZVector(px[i1], py[i1], pz[i1]);
    if (p0.Mag2() == 0 || p1.Mag2() == 0) return rad::consts::InvalidEntry<double>();
    
    return TMath::ACos(ROOT::Math::VectorUtil::CosTheta(p0, p1));
  }

  /**
   * @brief Calculates the magnitude of the momentum difference (DeltaP) between two vectors.
   * * This is a single-combination kernel. It expects RVecIndices containing one set 
   * with exactly two particle indices (i0, i1).
   * * @tparam T Type of position/momentum components (RVec<T>).
   * @param indices The RVecIndices container. Expected structure: indices[0][0] = i0, indices[0][1] = i1.
   * @param px The RVec of x-components.
   * @param py The RVec of y-components.
   * @param pz The RVec of z-components.
   * @return double The DeltaP value. Returns InvalidEntry if indices are invalid or momentum is zero.
   */
  template<typename Tp, typename Te>
  inline ResultType_t DeltaP(const RVecIndices &indices, const Tp &px, const Tp &py, const Tp &pz, const Te &e) {
    const auto& i0 = indices[0][0];
    const auto& i1 = indices[0][1];
    
    if (i0 < 0 || i1 < 0) return rad::consts::InvalidEntry<double>();
    auto p0 = XYZVector(px[i0], py[i0], pz[i0]);
    auto p1 = XYZVector(px[i1], py[i1], pz[i1]);
    if (p0.Mag2() == 0 || p1.Mag2() == 0) return rad::consts::InvalidEntry<double>();
    
    return TMath::Sqrt( (p0 - p1).Mag2() );
  }

  //---------------------- Debugging ----------------------

  /**
   * @brief Prints the 4-momentum and basic kinematic variables (p, theta) for all particles in an event.
   * @tparam Tpid Type of PID component (e.g., RVec<int>).
   * @tparam Tp Type of momentum components.
   * @tparam Te Type of energy component.
   * @param type String identifier for the data type (e.g., "tru" or "rec").
   * @param entry The event entry number (ULong64_t).
   * @param pid The RVec of PID values.
   * @param px The RVec of x-momentum components.
   * @param py The RVec of y-momentum components.
   * @param pz The RVec of z-momentum components.
   * @param e The RVec of energy components.
   * @return true Always returns true (useful for using this as an RDataFrame Filter/Define action).
   */
  template<typename Tpid, typename Tp, typename Tm>
  inline bool PrintParticles(const std::string& type, ULong64_t entry, const Tpid &pid, const Tp &px, const Tp &py, const Tp &pz, const Tm &m) {

    std::cout << setprecision(10);
    std::cout << type << " PrintParticles Event = " << entry << " number " <<px.size()<<std::endl;
    for (size_t idx = 0; idx < px.size(); ++idx) {
      auto p4 = PxPyPzMVector(px[idx], py[idx], pz[idx], m[idx]);
      std::cout << " " << pid[idx] << "\t" <<  p4
		<< " pmag " << ThreeVectorMag(px, py, pz)[idx] 
		<< " theta " << ThreeVectorTheta(px, py, pz)[idx]
		<< " mass " << p4.M()
		<< " Energy " << p4.E() << "\n";
      
    }
    return true;
  }

} // namespace rad
