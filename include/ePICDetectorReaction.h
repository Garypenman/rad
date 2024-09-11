#pragma once

//!  Derived class to configure ePIC root files with associations made to detector and track info

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for ePIC files with fixed particle order. It also uses podio to link particles to their full detector reconstruction
*/
#include "ePICReaction.h"
#include "RVecHelpers.h"
#include "podio/CollectionIDTable.h"//podio
#include <TFile.h>
#include <TTree.h>


namespace rad{
  namespace config{
    
 
    //! Class definition

    class ePICDetectorReaction : public ePICReaction {


    public:

      ePICDetectorReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ePICReaction{treeName,fileNameGlob,columns} {
	//open file and copy podio table
	//auto file = std::unique_ptr<TFile>{TFile::Open(fileNameGlob.data())};
	TChain chain("podio_metadata");
	chain.Add(fileNameGlob.data());
	//	auto podio_metadata = file->Get<TTree>("podio_metadata");
	podio::CollectionIDTable* tab=nullptr;
	//	podio_metadata->SetBranchAddress( "events___idTable", &tab);
	//podio_metadata->GetEntry(0);
	chain.SetBranchAddress( "events___idTable", &tab);
	chain.GetEntry(0);

	//dont't want to do this with bare pointer, but no copy constructor
	_idTable =tab ;
	
      }
    ePICDetectorReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ePICReaction{treeName,filenames,columns} {
	//open file and copy podio table
	auto file = std::unique_ptr<TFile>{TFile::Open(filenames[0].data())};
	auto podio_metadata = file->Get<TTree>("podio_metadata");
	podio::CollectionIDTable* tab=nullptr;
	podio_metadata->SetBranchAddress( "events___idTable", &tab);
	podio_metadata->GetEntry(0);

	//dont't want to do this with bare pointer, but no copy constructor
	_idTable =tab ;

      }

      /**
       * Configure association of clusters data to particles
       */
      void AssociateClusters(const std::vector<std::string>& types,const std::vector<std::string>& members ){
	AssociateObjects("clusters",types,members);
      }
     /**
       * Configure association of tracks data to particles
       */
      void AssociateTracks(const std::vector<std::string>& types,const std::vector<std::string>& members ){
	AssociateObjects("tracks",types,members);
      }

      /**
       * Configure given associations to create objects which may 
       * be used to create columns matched to reconstructed particles
       */
      void AssociateObjects(const string& object,const std::vector<std::string>& types,const std::vector<std::string>& members){

	ROOT::RVecU collIndices; //create indices for the given collection types
      
	for(const auto& assoc_name:types){ //loop over the specified detector associations
	  // cout<<assoc_name<<" "<<_idTable->collectionID(assoc_name).has_value()<<" "<<endl;
	  if(_idTable->collectionID(assoc_name).has_value()==false){
	    std::cerr<<"Warning : ePICDetectorReaction::AssociateObjects, no detector object "
		     <<assoc_name<<" in podio_metadata. "<<std::endl;
	    continue;
	  }
	  
	  //save the collectionID value in indices
	  //this allows us to define a local index for the collection
	  collIndices.push_back(_idTable->collectionID(assoc_name).value());
	}
      
        //function to map hash collectionID vector to index in local defined vector.
	//copy collIndices
	auto LocalCollIdx= [collIndices](const ROOT::RVecU& collID){
	  ROOT::RVecI localID(collID.size());
	  uint i = 0;
	  for(const auto id:collID){
	    auto dist =rad::helpers::findIndex(collIndices,id);
	    if(dist==collIndices.size()) localID[i]=-1;
	    else localID[i]=dist;
	    
	    ++i;
	  }
	  return localID;
	};
	
	//define collection indices for this object association
	auto collIdxsName = object+"_idxs";
	Define(collIdxsName , LocalCollIdx, {"_ReconstructedParticles_"+object+".collectionID"});

	//now create a function to define the associated vector of member to ReconstructedParticle
	auto CreateAssocVector = [](const ROOT::RVec<ROOT::RVecF>& all,const ROOT::RVecI& localCollIdx,const ROOT::RVecI& idxs){
	  auto Nelements = localCollIdx.size();
    
	  ROOT::RVecF result(Nelements);
	  //loop over ReconstructedParticles and get the associated data
	  for(uint i=0;i<Nelements;++i){
	    if(localCollIdx[i]<0) result[i]=0; //return 0 if collection not requested
	    else result[i] = all[localCollIdx[i]][idxs[i]]; //get the value at specific collection and index
	  }
	  
	  return result;
	};
	
	//loop over given data members and create list for each association 
	std::string strNames;
	for(const auto& member:members){//loop over the data member we are interested in
	  for(const auto& assoc_name:types){ //loop over the specified detector associations
	    
	    strNames+=assoc_name; //association name  e.g. CentralCKFTracks
	    strNames+="."+member; //plus specific data member required e.g. momentum.x
	    strNames+=","; //seperate each member by a comma e.g. CentralCKFTracks.momentum.x,
	    
	  }
	  strNames.pop_back();//remove last ,
	  
	  //make a map of corresponding values
	  TString mapName = object+member+"_epicdet";
	  mapName.ReplaceAll(".","_");
	  Define(mapName,Form("ROOT::RVec<ROOT::RVecF>{%s}",strNames.data()));
	  //unravel map into synchronised (with ReconstructedParticles ordering)
	  TString assocName = object+member;
	  assocName.ReplaceAll(".","_");
	  Define(assocName,CreateAssocVector,{mapName.Data(), collIdxsName,"_ReconstructedParticles_"+object+".index"});
	  //reorder to match truth and rec if required
	  if(rad::config::ColumnExists("tru_match_id",CurrFrame()) == true ){
	    RedefineExpr(std::string(assocName), std::string(Form("rad::helpers::Reorder(%s,%s,%s,%s)",assocName.Data(),"rec_match_id","tru_match_id","tru_n")) );
	  }
	}
	
      }

      /**
       * do not write columns containing "_epicdet"
       */
      void RemoveSnapshotColumns(std::vector<string>& cols) override{

	ePICReaction::RemoveSnapshotColumns(cols);
	
	cols.erase( std::remove_if( cols.begin(), cols.end(),
			       []( const string& col ) -> bool
			       { return col.find("_epicdet") != std::string::npos; } ),
		    cols.end() );
	
      }
  
    private:
      
      //  podio::CollectionIDTable _idTable; //to be cloned from podio_metadata:events___idTable
      podio::CollectionIDTable *_idTable=nullptr; //to be cloned from podio_metadata:events___idTable, it has a deleted constructor so cannot use copy contructor
      //std::shared_ptr<podio::CollectionIDTable> _idTable; //to be cloned from podio_metadata:events___idTable
      /**
       * Local index for particle associated collection IDs
       * should be in range 0 - N declared trackers
       */
      /* ROOT::RVecU _clustersCollectionIdxs; */
      /* ROOT::RVecU _tracksCollectionIdxs; */
      /* ROOT::RVecU _particleIDsCollectionIdxs; */
      /* ROOT::RVecU _startVertexCollectionIdxs; */
      /* ROOT::RVecU _particleIDUsedCollectionIdxs; */
      
    };
  }
}
