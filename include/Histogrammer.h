#pragma once


/*!
 *  Class to configure histograms with data splits etc
 */

#include "ConfigReaction.h"
#include "DataSplitter.h"
#include "SplitHistoHelper.h"
#include <TCanvas.h>
#include <TFile.h>

namespace rad{
  
  namespace histo{


 
    using hist_ptr = std::shared_ptr<TH1>;
    using hists_splits_ptr =  ROOT::RDF::RResultPtr< std::vector< hist_ptr > >;
    using hists_results =  std::vector<hists_splits_ptr>;
      
       
    //! Class definition

    class Histogrammer {

    public:

      Histogrammer(config::ConfigReaction& rad):_rad{rad}{

      }
      
      Histogrammer(const std::string& name,config::ConfigReaction& rad) :_name{name},_rad{rad}{
      }

      /** 
       * Initialise the splitting scheme
       * requires type of variables to be histogrammed
       * e.g. rec_ tru_ mc_
       */
      void Init(const std::vector<string>& types={""}){
	//variable types
	_types = types;
	//create results vector for each type
	for(const auto& type:_types){
	  _typeResults[type]=hists_results();
	}
	
	//configure the bins
	_splits.Init(); 
	auto nbins = Splitter().NTotal();
	if(nbins==1){ //no splits, fix bin to 0
	  _rad.Define(_name, "static_cast<short>(0)" );
	  return;
	}//don't need to do more
	
	std::string cast_to_deque_double="std::deque<double>{";

	//dimension names must correspond to columns in the RDF
	auto dim_names = _splits.GetDimensionNames();
	for(const auto& dname:dim_names){
	  if(_verbose)std::cout<<"Histogrammer::Init() split variable  name  :"<<dname<<std::endl;
	  cast_to_deque_double+=Form("static_cast<double>(%s),",dname.data());
	}
	cast_to_deque_double+="}";
	if(_verbose)std::cout<<"Histogrammer::Init() variables cast"<<cast_to_deque_double<<std::endl;
	//define column of variables which required to split data
	_rad.Define(_name+"_vars", cast_to_deque_double.data());
	
	//define a column with index of bin for this event
	//copy splits for thread safety
	auto split = [splits=_splits](std::deque<double> vars){
	  auto res = splits.GetBin(vars);
	  return res;};
	
	if(_verbose)std::cout<<"Histogrammer::Init() Define splitter column : "<<_name<<std::endl;
	//Define a new column which will be
	//the split bin index for the current event
	_rad.Define(_name, split ,{_name+"_vars"});
      }
      
      /** 
       * Create a series of histograms split as defined in Splitter
       * Must define template histogram type, HISTOGRAM (TH1F, TH2D,...)
       *             and column types which will fill the histrogram, ColumTypes (double, float, int...)
       * Must give a histogram template of type HISTOGRAM
       *           a vector of column names to fill with, must work with ColumnTypes
       */
      template <typename HISTOGRAM,typename... ColumnTypes>
      void Create(const HISTOGRAM& thHist, const ROOT::RDF::ColumnNames_t&  columns){

	auto nbins = Splitter().NTotal();
	auto bin_names = Splitter().GetBinNames();
	//get hist name
	std::string hname = thHist.GetName();
	//create one for each split
	std::vector<HISTOGRAM> hists(nbins,thHist);
       
	uint ibin=0;
	for(auto& h:hists){
	  //give each hist name of bin/split
	  std::string name = hname +"_"+bin_names[ibin];
	  std::string title = std::string(h.GetTitle())+" "+bin_names[ibin];
	  h.SetNameTitle(name.data(),title.data());
	  ibin++;
	}
       
	//connect name to index in results vector
	_getIndexFromName[hname]=_nResults++;

	//use Helper Action class
	//https://root.cern/doc/master/classROOT_1_1RDF_1_1RInterface.html#a77b83f7955ca336487ce102e7d31e7d8
	//template types double and TH1D, the double and D should match. i.e. float, TH1F
	using Helper_t = SplitHistoHelper<double>;

	// loop over types of variables
	for(const auto& type:_types){
	  auto process  = Helper_t{hists};

	  // book my action. template types : short from splits.GetBin() aka _name, one double from 1D
	  // store ResultsPtr in vector datamember
	  ROOT::RDF::ColumnNames_t cols = {_name};
	  cols.insert(cols.end(), columns.begin(), columns.end());
	  auto badCol=false;
	  for(auto& col:cols){//prepend type
	    if(col==_name) continue;//dont prepend name
	    // col = type+col;
	    auto temp_col = type+col;
	    //sometimes not all variables are defined for all types 
	    if( CheckColumn(temp_col)==false){
	      //but may be defined without type
	      if( rad::config::ColumnExists( col,_rad.CurrFrame())==false){	      
		_typeResults[type].push_back(hists_splits_ptr()  );
		badCol=true;
		std::cout<<"Warning :: Histogrammer Column "<<col<<" "<<" does not exist "<<type<<std::endl ;
	      }
	      else{
		//column exists but with no type
		//just use no type column

	      }
	    }
	    else{ //good, type column exists
	      col = temp_col;
	    }
	  }
	  if(badCol==true){
	    continue; //column does not exist, ignore it
	  }
	  auto df = _rad.CurrFrame();

	  //Book the histogram action and store it as a result
	  // std::cout<<" book ";
	  //for(auto& col:cols){cout<<" "<<col;}cout<<" "<<std::endl;
	  auto result = df.Book<short, ColumnTypes...>(std::move(process), cols );
	  _typeResults[type].push_back( result );

	  _rad.setCurrFrame(df);
 	}
	
      
      }

      /**
       *  Check to see if column exists, or is part of an array
       *  If it is an array element define it as a new column for histogramming
       */
      bool CheckColumn(std::string& col){

	if(rad::config::ColumnExists(col,_rad.CurrFrame())==true) return true;
	
	// Check if we request an element of an array
	// in which case we need to define as a column
	if(col.find('[')!=string::npos){// [  indicates array element
	  //now check if array exists, i.e. string up to [

	  if(rad::config::ColumnExists(col.substr(0,col.find('[')),_rad.CurrFrame())==false)
	    return false;

	  //we have an array and an element so define it as a column
	  auto new_col = col;
	  new_col.replace(col.find('['), 1,1, '_');
	  new_col.replace(col.find(']'), 1,1, '_');

	  //Check if already defined new_col
	  if(rad::config::ColumnExists(new_col,_rad.CurrFrame())==false){
	    if(_verbose) std::cout<<"CheckColumn create new column from array element : "<<new_col<<endl;
	    _rad.Define(new_col,col);
	  }
	  //change col, so we use the new one
	  col = new_col;
	  return true;
	}
	//column or array do not exist
	return false;
      }
      /** 
       * Get DataSplitter to define splits etc
       */
      DataSplitter& Splitter(){return _splits;} //to configure splits

      /** 
       * Get histogram with name at split/bin index
       */
      hist_ptr GetResult(const std::string& type,const std::string& name, ushort index){
	if( index >= Splitter().NTotal() ){
	  std::cerr<< "Histogrammer::GetResult index out of range " <<index <<" >= "<<Splitter().N()<<std::endl;
	  return nullptr;
	}
	if(TypeResult(type).at(_getIndexFromName[name]).GetPtr()==nullptr) return hist_ptr();//nulltptr
	if(TypeResult(type).at(_getIndexFromName[name])->size()==0) return hist_ptr();//nulltptr
	return TypeResult(type).at(_getIndexFromName[name])->at(index);
      }
      
      //Draw 2D Histograms (default option colz)
      //need a clever way of auto drawing either
      //MC or Tru + Rec depending on what type of macro used
      //for now happy just auto drawing MC() OR Tru()
      void Draw2DMC(const std::string& name, const TString opt="colz", const TVirtualPad* pad = nullptr){
	using namespace rad::names::data_type;
	if(!pad)
	  new TCanvas();	
	for(const auto type:_types){
	  if (type==Rec()) continue;
	  if(GetResult(type,name,0).get()==nullptr){
	    continue;
	  }else{
	    auto h=GetResult(type,name,0);
	    h->DrawCopy(opt);
	  }
	}
      }
      
      
      /** 
       * Draw all histograms of type name  on a single canvas
       * optional pad argument if drawing on existing tcanvas/pad
       * useful for drawing on divided canvas, for example.
       */
      void DrawSame(const std::string& name, const TVirtualPad* pad = nullptr){
	using namespace rad::names::data_type;
	if(!pad)
	  new TCanvas();
	int iter=0;
	auto oldmax=-1;
	//loop over histograms to get the max y value
	for(const auto& type:_types){
	  //check if this histogram exists for this type
	  if(GetResult(type,name,0).get()==nullptr){
	    continue;
	  }
	  else{
	    auto newmax = GetResult(type,name,0)->GetMaximum();
	    if(newmax>oldmax)
	      oldmax=newmax;
	  }
	}
	//re-loop to draw them
	for(const auto& type:_types){
	  //check if this histogram exists for this type
	  if(GetResult(type,name,0).get()==nullptr){
	    continue;
	  }
	  else{
	    TString opt="hist";
	    if (iter>0) opt="hist same";
	    auto his = GetResult(type,name,0)->DrawCopy(opt);
	    if(oldmax!=-1)
	      his->SetMaximum(oldmax*1.2);
	    if(type==Rec())his->SetLineColor(kRed);
	    if(type==Truth())his->SetLineColor(kBlue);
	    if(type==MC())his->SetLineColor(kBlack);
	    iter++;
	    //std::cout << type << std::endl;
	  }
	}
      }
      
    
      void DrawAll(const std::string& name){
	//Loop over types
	for(const auto& type:_types){
	  //check if this histogram exists for this type
	  if(GetResult(type,name,0).get()==nullptr){
	    continue;
	  }
	  new TCanvas();
	  auto hmax =0.;
	  for(size_t i = 0; i < Splitter().NTotal(); ++i){
	    auto mymax  = GetResult(type,name,i)->GetMaximum();
	    if(mymax>hmax) hmax = mymax;
	  }
	  GetResult(type,name,0)->SetMaximum(hmax);
	  GetResult(type,name,0)->SetMinimum(0);
	  
	  for(size_t i = 0; i < Splitter().NTotal(); ++i){
	    TString opt="";
	    if(i>0) opt = "same";
	    auto his = GetResult(type,name,i)->DrawCopy(opt);
	    his->SetLineColor(i+1);
	  }
	}//type loop
      }
      /** 
       * Get histograms for a type
       */
      hists_results& TypeResult(const std::string& type){
	return _typeResults[type];
      }
      /** 
       * Write all historgams to file
       */
      void File(const string& filename){

	auto file = std::unique_ptr<TFile>{TFile::Open(filename.data(),"recreate")};

	for(auto result: _getIndexFromName){
	  //Loop over types
	  for(const auto& type:_types){
	    //check if this histogram exists for this type
	    if(GetResult(type,result.first,0).get()==nullptr){
	      continue;
	    }
	    
	    //make directory with given hist name and type
	    file->cd();
	    
	    auto dir = file->mkdir((type+result.first).data());
	    file->cd((type+result.first).data());
	    //create a summed histogram too
	    std::unique_ptr<TH1> htotal;
	    //loop over all splits and write histograms
	    for(size_t i=0;i<Splitter().NTotal(); ++i){
	      auto h = GetResult(type,result.first,i);
	      if(i==0) htotal.reset(static_cast<TH1*>(h->Clone((type+result.first).data())));
	      else htotal->Add(h.get());
	      h->SetDirectory(dir);
	      h->Write();
	    }//bins
	    //write summed histogram
	    htotal->Write();
	  }//type
	}//hist
      }

      void SetVerbose(int val){_verbose = val;}
    private:
     
      config::ConfigReaction& _rad;// = nullptr;
      DataSplitter _splits;
      //      hists_results _results;
      std::map<std::string, hists_results > _typeResults;
      
      std::map<std::string, ushort> _getIndexFromName;
      
      std::vector<string> _types;
      std::string _name;
      ushort _verbose=0;
      ushort _nResults=0;
    };

  }
}
