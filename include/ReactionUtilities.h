#pragma once

#include "DefineNames.h"

#include <TString.h>
#include <ROOT/RDFHelpers.hxx>

namespace rad{
  namespace reaction{
    namespace util{
      
    
      // void CountParticles(rad::config::ConfigReaction* rad, const std::string& type){
      template<typename T> //use template so can #include this in ConfigReaction
      void CountParticles(T* rdf, const std::string& type){
      
	rdf->Define(type+"Ngamma",Form("rad::helpers::Count(%spid,22)",type.data()) );
	rdf->Define(type+"Npip",Form("rad::helpers::Count(%spid,211)",type.data()) );
	rdf->Define(type+"Npim",Form("rad::helpers::Count(%spid,-211)",type.data()) );
	rdf->Define(type+"NKp",Form("rad::helpers::Count(%spid,321)",type.data()) );
	rdf->Define(type+"NKm",Form("rad::helpers::Count(%spid,-321)",type.data()) );
	rdf->Define(type+"Nele",Form("rad::helpers::Count(%spid,11)",type.data()) );
	rdf->Define(type+"Npos",Form("rad::helpers::Count(%spid,-11)",type.data()) );
	rdf->Define(type+"Npro",Form("rad::helpers::Count(%spid,2212)",type.data()) );
	rdf->Define(type+"Nneutron",Form("rad::helpers::Count(%spid,2112)",type.data()) );
      }
    
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
      std::string ColumnsToStringNoBraces(const ROOT::RDF::ColumnNames_t &cols) {
	if(cols.empty()==true) return "";
      
	string toString ="";
	for(const auto& p:cols){
	  toString=(toString + p + ",");
	}
	toString.pop_back(); //remove last ,
       
	return toString;
      }
    

      using rad::names::data_type::Rec;
      using rad::names::data_type::Truth;

      /**
       * calculate the difference in reconsutructed and truth variables
       * Case Reconstructed and truth synched via AliasColumnsAndMatchWithMC()
       */
      template<typename T> //use template so can #include this in ConfigReaction
      void Resolution(T* const  rdf,const string& var){
	// Define(string("res_")+var,[](const ROOT::RVec<T> &rec,const ROOT::RVec<T> &tru){
	//   return (rec - tru);
	// },{string(Rec())+var,string(Truth())+var});
	rdf->Define(string("res_")+var,Form("%s-%s",(Truth()+var).data(),(Rec()+var).data() ));
      }
      template<typename T> //use template so can #include this in ConfigReaction
      void ResolutionFraction(T* const rdf,const string& var){
	// rdf->Define(string("res_")+var,[](const ROOT::RVecD &rec,const ROOT::RVecD &tru){
	//   return ROOT::RVecD((rec - tru)/tru);
	// },{Rec()+var,Truth()+var});
	rdf->Define(string("res_")+var,Form("(%s-%s)/%s",(Truth()+var).data(),(Rec()+var).data(),(Truth()+var).data() ));
      }
      
     // template<typename T> //use template so can #include this in ConfigReaction
     // ColType DeduceColumnVectorType(T* const radf,const string& name){

     //   TString col_type = radf->CurrFrame().GetColumnType(name);
       
     // 	if(col_type.Contains("UInt_t")||col_type.Contains("uint")) return ColType::UInt;
     // 	if(col_type.Contains("Float_t")||col_type.Contains("float")) return ColType::Float;
     // 	if(col_type.Contains("Double_t")||col_type.Contains("double")) return ColType::Double;
     // 	if(col_type.Contains("Short_t")||col_type.Contains("short")) return ColType::Short;
     // 	if(col_type.Contains("Bool_t")||col_type.Contains("bool")) return ColType::Bool;
     // 	if(col_type.Contains("Long_t")||col_type.Contains("long")) return ColType::Long;
     // 	if(col_type.Contains("Int_t")||col_type.Contains("int")) return ColType::Int;
     // 	return ColType::Undef;
     //  }
      // std::string CreateDeducedFunctionString(config::ConfigReaction& cr,const std::string func,){
      // auto obj_types =cr.TypeColObjTypeString("mc_","components_p4");
      // auto function_expr = "rad::FourVectorMassCalc<" + obj_types + ">";

      // }
      //This should perhaps be more general and moved to REactionUtils or somewehre
      // template<typename T> //use template so can #include this in ConfigReaction
      // void RedefineFundamentalAliases(T* const radf){

      // 	auto alias_map = radf->AliasMap();
      // 	for(const auto& col :alias_map ){
      // 	  const auto& alias = col.first;
      // 	  switch(static_cast<int>(DeduceColumnVectorType(radf, col.second )) ) {
	    
      // 	  case static_cast<int>(ColType::Undef):
      // 	    break;
      // 	  case static_cast<int>(ColType::UInt):
      // 	    radf->template RedefineFundamental<UInt_t>(alias);
      // 	    break;
      // 	  case static_cast<int>(ColType::Int):
      // 	    radf->template RedefineFundamental<Int_t>(alias);
      // 	    break;
      // 	  case static_cast<int>(ColType::Float):
      // 	    radf->template RedefineFundamental<Float_t>(alias);
      // 	    break;
      // 	  case static_cast<int>(ColType::Double):
      // 	    radf->template RedefineFundamental<Double_t>(alias);
      // 	    break;
      // 	  case static_cast<int>(ColType::Short):
      // 	    radf->template RedefineFundamental<Short_t>(alias);
      // 	    break;
      // 	  case static_cast<int>(ColType::Bool):
      // 	    radf->template RedefineFundamental<Bool_t>(alias);
      // 	    break;
      // 	  case static_cast<int>(ColType::Long):
      // 	    radf->template RedefineFundamental<Long_t>(alias);
      // 	    break;
	    
      // 	  default:
      // 	    break;
      // 	  }
      // 	}
	
      // }
      
    }
  }//reaction
}//rad
