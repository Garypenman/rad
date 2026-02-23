#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/FourVector.h"

namespace rad{
  namespace config{
    
    void BookHepMC3Snapshot(const string& filename){
      	try {
      	  ROOT::RDF::RSnapshotOptions opts;
      	  opts.fLazy = true;
      	  std::vector<std::string> cols = {
      	    "hepmc3_event",
      	    "GenRunInfo"
      	    //"hepmc3_event_modified"
      	    // add more as needed...
      	  };
	  
      	  SetupHepMC3SnapshotColumns(cols);
      	  auto snapshot_result = CurrFrame().Snapshot("hepmc3_tree",filename,cols,opts);
      	  _triggerSnapshots.emplace_back([snapshot = std::move(snapshot_result)]() mutable{});
      	} catch (const std::exception& ex) {
      	  std::cerr << "BookHepMC3Snapshot failed: " << ex.what() << std::endl;
      	  throw; // or handle gracefully
      	}
	
      }
      
      void HepMC3Snapshot(const string& filename){
      	try {
      	  std::vector<std::string> cols = {
      	    "hepmc3_event",
      	    "GenRunInfo"
      	  };
	  
      	  SetupHepMC3SnapshotColumns(cols);
      	  auto snapshot_result = CurrFrame().Snapshot("hepmc3_tree",filename,cols);
      	  _triggerSnapshots.emplace_back([snapshot = std::move(snapshot_result)]() mutable{});
      	} catch (const std::exception& ex) {
      	  std::cerr << "BookHepMC3Snapshot failed: " << ex.what() << std::endl;
      	  throw; // or handle gracefully
      	}
      }
      
      
      virtual void SetupHepMC3SnapshotColumns(std::vector<string>& cols){
	
      	setCurrFrame(CurrFrame().Redefine("hepmc3_event",[](const HepMC3::GenEventData &evt_data,const ROOT::RVec<double> &px,const ROOT::RVec<double> &py,const ROOT::RVec<double> &pz,const ROOT::RVec<double> &mass){
	      // Convert from data struct to GenEvent object
	      HepMC3::GenEvent new_evt(HepMC3::Units::GEV, HepMC3::Units::MM);
	      new_evt.read_data(evt_data);
	      
	      auto &particles = new_evt.particles();
	      size_t N = std::min({particles.size(), px.size(), py.size(), pz.size(), mass.size()});
	      
	      /* std::cout << "[DEBUG] new_evt.particles().size(): " << new_evt.particles().size() << std::endl; */
	      /* std::cout << "[DEBUG] px.size(): " << px.size() */
	      /* 	      << ", py.size(): " << py.size() */
	      /* 	      << ", pz.size(): " << pz.size() */
	      /* 	      << ", mass.size(): " << mass.size() << std::endl; */
	      
	      for (size_t i = 0; i < N; ++i) {
		HepMC3::FourVector new_mom(px[i], py[i], pz[i],std::sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i] + mass[i]*mass[i]));
		particles[i]->set_momentum(new_mom);
	      }
	      
	      // Convert back to data struct for ROOT
	      HepMC3::GenEventData modified_data;
	      new_evt.write_data(modified_data);
	      return modified_data;
	    },
	    {"hepmc3_event", MC()+"px", MC()+"py", MC()+"pz", MC()+"m"}
	    ));
      }
      
    
  }
}
