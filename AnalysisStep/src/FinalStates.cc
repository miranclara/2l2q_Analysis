#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

namespace {
  const unsigned nFS = 32;//changed
}


std::string finalState(int iFS) {
  if (iFS<0 || iFS>=(int)nFS) return "None";
  const std::string finalStates[nFS] = {"MMMM",      // 0
					"EEEE",      // 1
					"EEMM",      // 2
					"ZZ",        // 3
					"LLTT",      // 4
					"CRMMMMss",  // 5 
					"CRMMMMos",  // 6 
					"CREEEEss",  // 7 
					"CREEEEos",  //	8
					"CREEMMss",  // 9
					"CREEMMos",  // 10
					"CRMMEEss",  // 11
					"CRMMEEos",  // 12
					"ZL",        // 13
					"CRZLLHiSIP",   //14
					"CRZLLHiSIPMM", //15
					"CRZLLHiSIPKin",//16
					"CRZLL",     //17
					"ZLL",       //18
					"CRZ2mLL",   //19
					"CRZ2eLL",   //20
					"CRZLLss",   //21
					"CRZLLos_2P2F", //22
					"CRZLLos_3P1F", //23
					"CRZLLos_2P2F_ZZOnShell", //24
					"CRZLLos_3P1F_ZZOnShell", //25
					"ZZOnShell", //26
					"llTT",      //27
					"TTTT"       //28
					"EEQQ"       //29
					"MMQQ"       //30
					"TTQQ"       //30
  };
  return finalStates[iFS];			     	
}


std::string finalStateNiceName(int iFS) {
//  if (iFS<0||iFS>2) &&  ((iFS!=29) ||( iFS!=30)) return finalState(iFS);
  if ((iFS<0||iFS>2) && (iFS!=29 || iFS!=30)) return finalState(iFS);//changed
  int tempFS=0;
  const std::string finalStates[5] = {"4mu", "4e", "2e2mu","2e2q","2mu2q"};
 if(iFS>=0 && iFS<=2) tempFS=iFS;
 if (iFS==29 || iFS==30 ) tempFS=iFS-29+3;
//  return finalStates[iFS];
  return finalStates[tempFS];//changed
}



Channel finalState(std::string sFS) {
  for (unsigned i=0; i<nFS; ++i) {
    if (sFS==finalState(i)) return (Channel) i;
  }
  return NONE;
}
