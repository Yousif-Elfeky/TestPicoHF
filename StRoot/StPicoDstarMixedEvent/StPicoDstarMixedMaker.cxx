/* **************************************************
 *
 *  Authors: Yuanjing Ji
 Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StEvent/StDcaGeometry.h"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoETofPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StPicoEvent/StPicoMtdPidTraits.h"
#include "StPicoDstarMixedMaker.h"
#include "StAnaCuts.h"
#include "StMemStat.h"
#include "calmean.h"
#include "TLorentzVector.h"
#include "TVector3.h"
bool DEBUG = false;


ClassImp(StPicoDstarMixedMaker)
  StPicoDstarMixedMaker::StPicoDstarMixedMaker(char const * name, TString const inputFilesList, TString const outFileBaseName, StPicoDstMaker* picoDstMaker):
    StMaker(name), mPicoDstMaker(picoDstMaker),
    mInputFilesList(inputFilesList), mOutFileBaseName(outFileBaseName),mRunbyRunQA(true)
{}

Int_t StPicoDstarMixedMaker::Init()
{
  mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
  // -------------- USER VARIABLES -------------------------
  mFile = new TFile(mOutFileBaseName+".QA.root", "RECREATE");
  //mFile_RunID = new TFile(mOutFileBaseName+".RunID.root","RECREATE");
  //initialize trees
  initHists();

  return kStOK;
}
//-----------------------------------------------------------------------------
StPicoDstarMixedMaker::~StPicoDstarMixedMaker()
{}
//-----------------------------------------------------------------------------
void StPicoDstarMixedMaker::initHists(){
  //int totalNum = 170;//9.2GeV 
  //int totalNum = 661;//9.2GeV 2020b
  //int totalNum = 403;//9.2GeV 2020c
  //int totalNum = 335;//9.2GeV 2020c 20200807
  //int totalNum = 420;//9.2GeV 2020c 20200807
  //int totalNum = 508;//9.2GeV 2020c 20200828
  //int totalNum = 141;//7.7GeV 2020 20200911
  //int totalNum = 1127;//19.6GeV 2020c
  //int totalNum = 1950;//11.5GeV 
  //int totalNum = 29;//19.5GeV fixed target
  //int totalNum = 14;//54GeV 2020 reproduction test
  //int totalNum = 547;//54GeV 2020 reproduction 20201008
  //int totalNum = 53;//54GeV 2020 reproduction 20201008
  //int totalNum = 13;//200 GeV SL20d AuAu200
  //int totalNum = 567;//7.7GeV 2021
  //int totalNum = 562;//7.7GeV 2021
  //int totalNum = 62;//200GeV OO 2021 mb
  //int totalNum = 75;//200GeV OO 2021 ct
  //int totalNum = 170;//200GeV 17p3 2021
  //int totalNum = 37;//200GeV 17p3 20210628
  //int totalNum = 156;// 19p6 20210805
  //int totalNum = 30;// 19p6 20210805 SL21cB
  //int totalNum = 1145;// 19p6 20211221
  //int totalNum = 59;// AuAu200 2023
  //int totalNum = 1125; //19GeV 2019 
  //int totalNum = 32; //temp
  
  //char name_RunID[100];
ifstream readnum;
readnum.open(mRunNumList);

if (!readnum.is_open()) {
  cout << "Error: Could not open run number list file: " << mRunNumList << endl;
  return; 
}

int totalNum = 0; 

if (mRunbyRunQA) {
  cout << "Starting to initial run numbers..." << endl;

  int tmpRunNum; 
  int index = 0; 

  while (readnum >> tmpRunNum) {
    runnum.insert(pair<int, int>(tmpRunNum, index));
    if (DEBUG) cout << "Read run number: " << tmpRunNum << " -> assigned id: " << index << endl;
    index++;
  }
  readnum.close();

  // <<< UPDATE IT HERE. Now we are setting the value of the outer variable.
  totalNum = runnum.size();
  if(DEBUG) cout << "Successfully read " << totalNum << " run numbers in total." << endl;
}

if(QA){
  // event level QA
  hevt = new TH1D("hevt","hevt",totalNum,0,totalNum);
  hevtcut = new TH1D("hevtcut","hevtcut",totalNum,0,totalNum);
  hevtbadcut = new TH1D("hevtbadcut","Events after remove bad run;Run;Counts",totalNum,0,totalNum);
  hpassevtcut = new TH1D("hpassevtcut","pass event cut", 6, -0.5 , 5.5 );
  //run by run QA
  if (mRunbyRunQA){ 
    pVpdVz = new TProfile("VpdVz","VpdVz vs runId;runId;VpdVz(cm)",totalNum,0,totalNum);
    pVzVpdVz = new TProfile("VzVpdVz","VzVpdVz vs runId;runId;VpdVz-Vz(cm)",totalNum,0,totalNum);
    pRefmult = new TProfile("Refmult","Refmult vs runId;runId;Refmult",totalNum,0,totalNum);
    pVpdHitEast = new TProfile("pVpdHitEast","pVpdHitEast;runId;pVpdHitEast",totalNum,0,totalNum);
    pVpdHitWest = new TProfile("pVpdHitWest","pVpdHitWest;runId;pVpdHitWest",totalNum,0,totalNum);
    pVz  = new TProfile("pVz","pVz;runId;pVz(cm)",totalNum,0,totalNum);
    pVx = new TProfile("pVx","pVx;runId;pVx(cm)",totalNum,0,totalNum);
    pVy = new TProfile("pVy","pVy;runId;pVy(cm)",totalNum,0,totalNum);
    pVr = new TProfile("pVr","pVr;runId;pVr(cm)",totalNum,0,totalNum);
    
    //track level QA
    pTof = new TProfile("Tof","1/#beta vs runId;runId;1/#beta",totalNum,0,totalNum);
    pDedx = new TProfile("Dedx","dEdx vs runId;runId;dEdx",totalNum,0,totalNum);
    //pRMSDedx = new TProfile("RMSDedx","RMSdEdx vs runId;runId;RMSdEdx",totalNum,0,totalNum);
    pgDCA = new TProfile("gDCA","gDCA vs runId;runId;global DCA(cm)",totalNum,0,totalNum);
    ppDCA = new TProfile("pDCA","pDCA vs runId;runId;primary DCA(cm)",totalNum,0,totalNum);
    pgPt = new TProfile("gPt","global Pt vs runId;runId;global p_{T}(GeV/c)",totalNum,0,totalNum);
    pgPhi = new TProfile("gPhi","global Phi vs runId;runId;gPhi",totalNum,0,totalNum);
    pgEta = new TProfile("gEta","global Eta vs runId;runId;Eta",totalNum,0,totalNum);
    pNFits = new TProfile("NFits","NHitsFit vs runId;runId;nHitsFit",totalNum,0,totalNum);
    ppPt = new TProfile("pPt","primary Pt vs runId;runId;primary p_{T}(GeV/c)",totalNum,0,totalNum);
    ppEta = new TProfile("pEta","primary Eta vs runId;runId;pEta",totalNum,0,totalNum);
    ppPhi = new TProfile("pPhi","primary Phi vs runId;runId;pPhi",totalNum,0,totalNum);
    p_nSigmaE = new TProfile("p_nSigmaE","primary nSigmaE vs runId;runId;n#sigma_{e}",totalNum,0,totalNum);
    p_nSigmaPion = new TProfile("p_nSigmaPion","primary nSigmaPion vs runId;runId;n#sigma_{#pi}",totalNum,0,totalNum);
    p_nSigmaKaon = new TProfile("p_nSigmaKaon","primary nSigmaKaon vs runId;runId;n#sigma_{K}",totalNum,0,totalNum);
    p_nSigmaProton = new TProfile("p_nSigmaProton","primary nSigmaProton vs runId;runId;n#sigma_{p}",totalNum,0,totalNum);
    p_nTofMatch = new TProfile("p_nTofMatch","nTofMatch vs runId;runId;nTofMatch",totalNum,0,totalNum);
    
    p_nSigmaTofE = new TProfile("p_nSigmaTofE","primary nSigmaTofE vs runId;runId;nSigmaTof_{e}",totalNum,0,totalNum);
    p_nSigmaTofPion = new TProfile("p_nSigmaTofPion","primary nSigmaTofPion vs runId;runId;nSigmaTof_{#pi}",totalNum,0,totalNum);
    p_nSigmaTofKaon = new TProfile("p_nSigmaTofKaon","primary nSigmaTofKaon vs runId;runId;nSigmaTof_{K}",totalNum,0,totalNum);
    p_nSigmaTofProton = new TProfile("p_nSigmaTofProton","primary nSigmaTofProton vs runId;runId;nSigmaTof_{p}",totalNum,0,totalNum);
    p_btofYLocal = new TProfile("p_btofYLocal","btofYLocal vs runId;runId;btofYLocal",totalNum,0,totalNum);
    p_btofZLocal = new TProfile("p_btofZLocal","btofZLocal vs runId;runId;btofZLocal",totalNum,0,totalNum);

    p_nETofHits = new TProfile("p_nETofHits","ETof Hits vs runId;runId;nETofHits",totalNum,0,totalNum);
    p_nBemcPidTraits = new TProfile("p_nBemcPidTraits","nBemcPidTraits vs runId;runId;nBemcPidTraits",totalNum,0,totalNum);
    p_nMtdHits = new TProfile("p_nMtdHits","nMtdHits vs runId;runId;nMtdHits",totalNum,0,totalNum);
    p_nMtdPidTraits = new TProfile("p_nMtdPidTraits","nMtdPidTraits vs runId;runId;nMtdPidTraits",totalNum,0,totalNum);
    p_ETof_beta = new TProfile("p_ETof_beta","ETof beta vs runId;runId;1/beta(ETof)",totalNum,0,totalNum);
    p_ETof_deltaX = new TProfile("p_ETof_deltaX","ETof delta X vs runId;runId;#Delta X (ETof)",totalNum,0,totalNum);
    p_ETof_deltaY = new TProfile("p_ETof_deltaY","ETof delta Y vs runId;runId;#Delta Y (ETof)",totalNum,0,totalNum);
  }
   //event QA
    hVxVyVz = new TH3F("hVxVyVz","VxVyVz;Vx(cm);Vy(cm);Vz(cm)",100,-0.5,0.5,100,-0.5,0.5,240,-60,60);
    hVz = new TH1F("hVz","Vz;Vz(cm);Counts",800,-200,200);
    hVpdVz = new TH1F("hVpdVz","VpdVz;VpdVz(cm);Counts",800,-200,200);
    h_Vz_VpdVz = new TH2F("h_Vz_VpdVz","Vz vs VpdVz;Vz(cm);VpdVz(cm)",800,-200,200,800,-200,200);
    hVr = new TH1F("hVr","Vr;Vr(cm);Counts",100,0,1);
    hVzVpdVz = new TH1F("hVzVpdVz","Vz-VpdVz(cm)",2000,-100,100);
    h_Vx_Vy = new TH2F("h_Vx_Vy","Vertex XY",600,-3,3,600,-3,3);
    hnEvsEtavsVz = new TH2F("hnEvsEtavsVz","nElectron;#eta;Vz(cm)",40,-1.5,1.5,240,-60,60);
    hnEvsPhivsVz = new TH2F("hnEvsPhivsVz","nElectron;#phi;Vz(cm)",100,-1*TMath::Pi(),TMath::Pi(),240,-60,60);
    hnTofMult = new TH1F("hnTofMult","TOF Multiplicity",2000,0,2000); 
    hnTofMulvsRef = new TH2F("hnTofMulvsRef","RefMult vs nTofMult;RefMult;btofTrayMultiplicity",900,0,900,2000,0,2000); 
    hnTofMatch = new TH1F("hnTofMatch","BTOF matched tracks",900,0,900);
    hnTofMatvsRef= new TH2F("hnTofMatvsRef","RefMult VS nTofmatch;RefMult;nTofMatch",900,0,900,900,0,900);
    hnTofHitvsRef= new TH2F("hnTofHitvsRef","hnTofHit vs RefMult;nTofHits;refMult",900,0,900,900,0,900);
    hrefmult = new TH1F("hrefmult","RefMult for all tracks",700,0,700);
    hrefmult_Pos = new TH1F("hrefmult_Pos","RefMult for positive tracks",700,0,700);
    hrefmult_Neg = new TH1F("hrefmult_Neg","RefMult for negative tracks",700,0,700);

    h_Vz_btofZLocal = new TH2F("h_Vz_btofZLocal"," ;Vz;bTofZLocal",400,-200,200,1000,-5,5);
    h_Vz_btofYLocal = new TH2F("h_Vz_btofYLocal"," ;Vz;bTofYLocal",400,-200,200,1000,-5,5);

    //tracl level QA
    hnHitsFit = new TH1F("hnHitsFit","nHitsFit before cut;nHitsFit",180,-90,90);
    hnHitsFit_cut = new TH1F("hnHitsFit_cut","nHitsFit after cut;nHitsFit",180,-90,90);
    hnHitsPoss = new TH1F("hnHitsPoss","nHitsPoss before cut;nHitsPoss",180,-90,90);
    hnHitsPoss_cut = new TH1F("hnHitsPoss_cut","nHitsPoss after cut;nHitsPoss",180,-90,90);
    hnHitsDedx = new TH1F("hnHitsDedx","nHitsDedx before cut;nHitsDedx",180,-90,90);
    hnHitsDedx_cut = new TH1F("hnHitsDedx_cut","nHitsDedx after cut;nHitsDedx",180,-90,90);
    h_nHitsDedx_p = new TH2F("h_nHitsDedx_p","nHitsDedx vs p*charge after cut;p*charge (GeV/c);nHitsDedx",400,-10,10,90,0,90);
    h_p_E0 = new TH1F("h_p_E0","#frac{p}{E0};#frac{p}{E0}",40,0,4);
    ph_p_E0 = new TProfile("ph_p_E0","#frac{p}{E0};RunId;#frac{p}{E0}",totalNum,0,totalNum);
    h_bemcdz = new TH1F("h_bemcdz","bemcdz;bemcdz",40,-20,20);
    ph_bemcdz = new TProfile("ph_bemcdz","bemcdz;RunId;bemcdz",totalNum,0,totalNum);
    h_bemcDphi = new TH1F("h_bemcDphi","bemcDphi;bemcDphi",100,-0.5,0.5);
    ph_bemcDphi = new TProfile("ph_bemcDphi","bemcDphi;RunId;bemcDphi",totalNum,0,totalNum);
    h_nSMDphi = new TH1F("h_nSMDphi","nSMDphi;nSMDphi",20,-10,10);
    ph_nSMDphi = new TProfile("ph_nSMDphi","nSMDphi;RunId;nSMDphi",totalNum,0,totalNum);
    h_nSMDeta = new TH1F("h_nSMDeta","nSMDeta;nSMDeta",20,-10,10);
    ph_nSMDeta = new TProfile("ph_nSMDeta","nSMDeta;RunId;nSMDeta",totalNum,0,totalNum);
    h_btowPhiDz = new TH1F("h_btowPhiDz","btowPhiDz;btowPhiDz",1000,-50,50);
    ph_btowPhiDz = new TProfile("ph_btowPhiDz","btowPhiDz;RunId;btowPhiDz",totalNum,0,totalNum);
    h_btowEtaDz = new TH1F("h_btowEtaDz","btowEtaDz;btowEtaDz",1000,-50,50);
    ph_btowEtaDz = new TProfile("ph_btowEtaDz","btowEtaDz;RunId;btowEtaDz",totalNum,0,totalNum);
    
    h_mtdDeltaY = new TH1F("h_mtdDeltaY","#Delta Y(mtd);#Delta Y",1000,-50,50);
    ph_mtdDeltaY = new TProfile("ph_mtdDeltaY","#Delta Y(mtd);RunId;#Delta Y",totalNum,0,totalNum);
    h_mtdDeltaZ = new TH1F("h_mtdDeltaZ","#Delta Z(mtd);#Delta Z",1000,-50,50);
    ph_mtdDeltaZ = new TProfile("ph_mtdDeltaZ","#Delta Z(mtd);RunId;#Delta Z",totalNum,0,totalNum);
    h_mtdDeltaTOF = new TH1F("h_mtdDeltaTOF","#Delta TOF(mtd);#Delta TOF",20000,-2000,2000);
    ph_mtdDeltaTOF = new TProfile("ph_mtdDeltaTOF","#Delta TOF(mtd);RunId;#Delta TOF",totalNum,0,totalNum);

    hgDca = new TH1F("hgDca","gDca",50,0,5);
    hpDca = new TH1F("hpDca","pDca",50,0,5);
    hinvBetavsP = new TH2F("hinvBetavsP","#frac{1}{#beta} vs p;p(GeV/c);#frac{1}{#beta}",300,0,3,200,0.5,2.5);
    hdEdx = new TH2F("hdEdx","dEdx vs p*charge;p*charge(GeV/c);dEdx",200,-2,2,400,0,25);
    h_mTpc = new TH2F("h_mTpc","mass^{2} vs p*charge;p*charge(GeV/c);mass square (GeV/c^{2})^{2}",200,-2,2,120,0,1.2);
    hNsigEvsinvBeta = new TH3F("hNsigEvsinvBeta","nSigmaE vs 1/#beta;nSigmaE;1/#beta;p",200,-10,10,100,0.8,1.2,100,0.15,2.5);
    hpt_Pos = new TH1F("hpt_Pos","Positive track p_{T} before cut;p_{T}(GeV/c)",300,0,30);
    hpt_Pos_cut = new TH1F("hpt_Pos_cut","Positive track p_{T} after cut;p_{T}(GeV/c)",300,0,30);
    hpt_Neg = new TH1F("hpt_Neg","Negative track p_{T} before cut;p_{T}(GeV/c)",300,0,30);
    hpt_Neg_cut = new TH1F("hpt_Neg_cut","Negative track p_{T} after cut;p_{T}(GeV/c)",300,0,30);
    hGpt_Pos = new TH1F("hGpt_Pos","Positive track global p_{T} before cut;global p_{T}(GeV/c)",300,0,30);
    hGpt_Pos_cut = new TH1F("hGpt_Pos_cut","Positive track global p_{T} after cut;global p_{T}(GeV/c)",300,0,30);
    hGpt_Neg = new TH1F("hGpt_Neg","Negative track global p_{T} before cut;global p_{T}(GeV/c)",300,0,30);
    hGpt_Neg_cut = new TH1F("hGpt_Neg_cut","Negative track global p_{T} after cut;global p_{T}(GeV/c)",300,0,30);
    hEta = new TH1F("hEta","#eta before cut;#eta",120,-3.0,3.0);
    hEta_cut = new TH1F("hEta_cut","#eta after cut;#eta",80,-2.0,2.0);
    hPhi = new TH1F("hPhi","#phi before cut;#phi",80,-4,4);
    hPhi_cut = new TH1F("hPhi_cut","#phi after cut;#phi",80,-4,4);
    h_pT_Eta = new TH2F("h_pT_Eta","p_{T}*charge vs #eta;p_{T}*charge (GeV/c);#eta",400,-10,10,60,-1.5,1.5);
    h_pT_Phi = new TH2F("h_pT_Phi","p_{T}*charge vs #phi;p_{T}*charge (GeV/c);#phi",400,-10,10,160,-4,4);
    h_EtavsPhi_Pos = new TH2F("h_EtavsPhi_Pos","Positive Tracks #eta vs #phi;#eta;#phi",60,-1.5,1.5,160,-4,4);
    h_EtavsPhi_Neg = new TH2F("h_EtavsPhi_Neg","Negative Tracks #eta vs #phi;#eta;#phi",60,-1.5,1.5,160,-4,4);
    hBadTofId = new TH3F("hBadTofId","hBadTofId;tray;module;cell",125,-0.5,124.5,33,-0.5,32.5,7,-0.5,6.5);
    TofId = new TH1F("TofId","1.13<1/#beta<1.24 0.3<p<0.5;TofId",23100,0,23100);
    TofId_nSigmaPi = new TH1F("TofId_nSigmaPi","1.13<1/#beta<1.24 0.3<p<0.5;nSigmaPi",200,-10,10);
    h_ETof_beta = new TH1F("h_ETof_beta","ETof beta;#frac{1}{#beta}(ETof)",250,0,2.5);
    h_ETof_betavsP = new TH2F("h_ETof_betavsP","ETof beta;p(GeV/c);#frac{1}{#beta}(ETof)",300,0,3,250,0,2.5);
    h_ETof_deltaX = new TH1F("h_ETof_deltaX","ETof delta X;#Delta X (ETof);Counts",2000,-100,100);
    h_ETof_deltaY = new TH1F("h_ETof_deltaY","ETof delta Y;#Delta Y (ETof);Counts",2000,-100,100);

  }
    //invariant mass electron
    //hMeeCount = new TH1F("hMee","hMee;Count;Mee(GeV/c^{2})",400,0,4);
    hMeeCount = new TH1F("hMee","hMee;Mee(GeV/c^{2})",40000,0,4);
    hMeeCount_like1 = new TH1F("hMee_like1","hMee like sign electron;Mee(GeV/c^{2})",40000,0,4);
    hMeeCount_like2 = new TH1F("hMee_like2","hMee like sign positron;Mee(GeV/c^{2})",40000,0,4);
    hMeeCountPt = new TH2F("hMeePt","hMee vs p_{T};Mee(GeV/c^{2});p_{T}",40000,0,4,200,0,10);
    hMeeCountPt_like1 = new TH2F("hMeePt_like1","hMee vs p_{T} like sign electron;Mee(GeV/c^{2});p_{T}",40000,0,4,200,0,10);
    hMeeCountPt_like2 = new TH2F("hMeePt_like2","hMee vs p_{T} like sign positron;Mee(GeV/c^{2});p_{T}",40000,0,4,200,0,10);

    // Mkpi histograms for D^{0} analysis
    hMkpiCount = new TH1F("hMkpi","hMkpi;M_{K#pi}(GeV/c^{2})",40000,0,4);
    hMkpiCount_like1 = new TH1F("hMkpi_like1","hMkpi like sign negative;M_{K#pi}(GeV/c^{2})",40000,0,4);
    hMkpiCount_like2 = new TH1F("hMkpi_like2","hMkpi like sign positive;M_{K#pi}(GeV/c^{2})",40000,0,4);
    hMkpiCountPt = new TH2F("hMkpiPt","hMkpi vs p_{T};M_{K#pi}(GeV/c^{2});p_{T}",40000,0,4,200,0,10);
    hMkpiCountPt_like1 = new TH2F("hMkpiPt_like1","hMkpi vs p_{T} like sign negative;M_{K#pi}(GeV/c^{2});p_{T}",40000,0,4,200,0,10);
    hMkpiCountPt_like2 = new TH2F("hMkpiPt_like2","hMkpi vs p_{T} like sign positive;M_{K#pi}(GeV/c^{2});p_{T}",40000,0,4,200,0,10);
    h_nSigmaElectron_P = new TH2F("h_nSigmaElectron_P","nSigmaElectron_P;P(GeV/c);nSigmaElectron",300,0,3,150,-10,5);
    h_nSigmaElectron_P_tpc = new TH2F("h_nSigmaElectron_P_tpc","nSigmaElectron_P;P(GeV/c);nSigmaElectron",300,0,3,150,-10,5);

    // =================================================================
    // INITIALIZE THE NEW HISTOGRAMS HERE
    // =================================================================

    // --- 1D Inclusive nSigma Distributions ---
    h_nSigmaElectron_Inclusive = new TH1F("h_nSigmaElectron_Inclusive", "Inclusive TPC n#sigma_{Electron};n#sigma_{e};Counts", 300, -15, 15);
    h_nSigmaPion_Inclusive     = new TH1F("h_nSigmaPion_Inclusive", "Inclusive TPC n#sigma_{Pion};n#sigma_{#pi};Counts", 300, -15, 15);
    h_nSigmaKaon_Inclusive     = new TH1F("h_nSigmaKaon_Inclusive", "Inclusive TPC n#sigma_{Kaon};n#sigma_{K};Counts", 300, -15, 15);
    h_nSigmaProton_Inclusive   = new TH1F("h_nSigmaProton_Inclusive", "Inclusive TPC n#sigma_{Proton};n#sigma_{p};Counts", 300, -15, 15);
    
    // --- 2D Inclusive nSigma vs. p*charge ---
    h_nSigmaVsPcharge_Electron = new TH2F("h_nSigmaVsPcharge_Electron", "n#sigma_{Electron} vs p*charge;p*charge (GeV/c);n#sigma_{e}", 1000, -5, 5, 300, -15, 15);
    h_nSigmaVsPcharge_Pion     = new TH2F("h_nSigmaVsPcharge_Pion", "n#sigma_{Pion} vs p*charge;p*charge (GeV/c);n#sigma_{#pi}", 1000, -5, 5, 300, -15, 15);
    h_nSigmaVsPcharge_Kaon     = new TH2F("h_nSigmaVsPcharge_Kaon", "n#sigma_{Kaon} vs p*charge;p*charge (GeV/c);n#sigma_{K}", 1000, -5, 5, 300, -15, 15);
    h_nSigmaVsPcharge_Proton   = new TH2F("h_nSigmaVsPcharge_Proton", "n#sigma_{Proton} vs p*charge;p*charge (GeV/c);n#sigma_{p}", 1000, -5, 5, 300, -15, 15);
    
    // --- TPC nSigma distributions for tracks passing TOF cuts ---
    h_TpcNsigmaE_afterTofCut = new TH1F("h_TpcNsigmaE_afterTofCut", "TPC n#sigma_{e} (for TOF Electrons);n#sigma_{e};Counts", 200, -10, 10);
    h_TpcNsigmaPi_afterTofCut = new TH1F("h_TpcNsigmaPi_afterTofCut", "TPC n#sigma_{#pi} (for TOF Pions);n#sigma_{#pi};Counts", 200, -10, 10);
    h_TpcNsigmaK_afterTofCut = new TH1F("h_TpcNsigmaK_afterTofCut", "TPC n#sigma_{K} (for TOF Kaons);n#sigma_{K};Counts", 200, -10, 10);
    h_TpcNsigmaP_afterTofCut = new TH1F("h_TpcNsigmaP_afterTofCut", "TPC n#sigma_{p} (for TOF Protons);n#sigma_{p};Counts", 200, -10, 10);
    
    // --- TOF nSigma distributions for tracks passing TPC cuts ---
    h_TofNsigmaE_afterTpcCut = new TH1F("h_TofNsigmaE_afterTpcCut", "TOF n#sigma_{e} (for TPC Electrons);n#sigma_{e};Counts", 200, -10, 10);
    h_TofNsigmaPi_afterTpcCut = new TH1F("h_TofNsigmaPi_afterTpcCut", "TOF n#sigma_{#pi} (for TPC Pions);n#sigma_{#pi};Counts", 200, -10, 10);
    h_TofNsigmaK_afterTpcCut = new TH1F("h_TofNsigmaK_afterTpcCut", "TOF n#sigma_{K} (for TPC Kaons);n#sigma_{K};Counts", 200, -10, 10);
    h_TofNsigmaP_afterTpcCut = new TH1F("h_TofNsigmaP_afterTpcCut", "TOF n#sigma_{p} (for TPC Protons);n#sigma_{p};Counts", 200, -10, 10);


    const int nBinsP = 1000, nBinsNsigma = 300;
    const float pMin = -5.0, pMax = 5.0;
    const float nSigmaMin = -15.0, nSigmaMax = 15.0;

    // --- 2D TPC nSigma vs p*charge for tracks passing TOF cuts ---
    h2_TpcNsigmaE_vs_p_afterTofCut = new TH2F("h2_TpcNsigmaE_vs_p_afterTofCut", "TPC n#sigma_{e} vs p*charge (for TOF Electrons);p*charge (GeV/c);n#sigma_{e}", nBinsP, pMin, pMax, nBinsNsigma, nSigmaMin, nSigmaMax);
    h2_TpcNsigmaPi_vs_p_afterTofCut = new TH2F("h2_TpcNsigmaPi_vs_p_afterTofCut", "TPC n#sigma_{#pi} vs p*charge (for TOF Pions);p*charge (GeV/c);n#sigma_{#pi}", nBinsP, pMin, pMax, nBinsNsigma, nSigmaMin, nSigmaMax);
    h2_TpcNsigmaK_vs_p_afterTofCut = new TH2F("h2_TpcNsigmaK_vs_p_afterTofCut", "TPC n#sigma_{K} vs p*charge (for TOF Kaons);p*charge (GeV/c);n#sigma_{K}", nBinsP, pMin, pMax, nBinsNsigma, nSigmaMin, nSigmaMax);
    h2_TpcNsigmaP_vs_p_afterTofCut = new TH2F("h2_TpcNsigmaP_vs_p_afterTofCut", "TPC n#sigma_{p} vs p*charge (for TOF Protons);p*charge (GeV/c);n#sigma_{p}", nBinsP, pMin, pMax, nBinsNsigma, nSigmaMin, nSigmaMax);

    // --- 2D TOF nSigma vs p*charge for tracks passing TPC cuts ---
    h2_TofNsigmaE_vs_p_afterTpcCut = new TH2F("h2_TofNsigmaE_vs_p_afterTpcCut", "TOF n#sigma_{e} vs p*charge (for TPC Electrons);p*charge (GeV/c);n#sigma_{e}", nBinsP, pMin, pMax, nBinsNsigma, nSigmaMin, nSigmaMax);
    h2_TofNsigmaPi_vs_p_afterTpcCut = new TH2F("h2_TofNsigmaPi_vs_p_afterTpcCut", "TOF n#sigma_{#pi} vs p*charge (for TPC Pions);p*charge (GeV/c);n#sigma_{#pi}", nBinsP, pMin, pMax, nBinsNsigma, nSigmaMin, nSigmaMax);
    h2_TofNsigmaK_vs_p_afterTpcCut = new TH2F("h2_TofNsigmaK_vs_p_afterTpcCut", "TOF n#sigma_{K} vs p*charge (for TPC Kaons);p*charge (GeV/c);n#sigma_{K}", nBinsP, pMin, pMax, nBinsNsigma, nSigmaMin, nSigmaMax);
    h2_TofNsigmaP_vs_p_afterTpcCut = new TH2F("h2_TofNsigmaP_vs_p_afterTpcCut", "TOF n#sigma_{p} vs p*charge (for TPC Protons);p*charge (GeV/c);n#sigma_{p}", nBinsP, pMin, pMax, nBinsNsigma, nSigmaMin, nSigmaMax);



    // =================================================================


    //tof module id
    /*ModuleId_1 = new TH1F("ModuleId 1","0.8<1/#beta<0.9 0.4<P;ModuleId",40,0,40);
    TofId_1 = new TH1F("TofId 1","0.8<1/#beta<0.9 0.4<P;TofId",23100,0,23100);
    TrayId_1 = new TH1F("TrayId 1","0.8<1/#beta<0.9 0.4<P;TrayId",130,0,130);
    ModuleId_2 = new TH1F("ModuleId 2","0.82<1/#beta<0.9 0.4<P;ModuleId",40,0,40);
    TofId_2 = new TH1F("TofId 2","0.8<1/#beta<0.9 0.4<P;TofId",23100,0,23100);
    ModuleId_3 = new TH1F("ModuleId 3","0.82<1/#beta<0.88 0.4<P;ModuleId",40,0,40);
    TofId_3 = new TH1F("TofId 3","0.8<1/#beta<0.9 0.4<P;TofId",23100,0,23100);
    ModuleId_4 = new TH1F("ModuleId 4","0.84<1/#beta<0.88 0.4<P;ModuleId",40,0,40);
    TofId_4 = new TH1F("TofId 4","0.8<1/#beta<0.9 0.4<P;TofId",23100,0,23100);
    ModuleId_5 = new TH1F("ModuleId 5","0.84<1/#beta<0.9 & 0.4<P;ModuleId",40,0,40);
    TofId_5 = new TH1F("TofId 5","0.8<1/#beta<0.9 0.4<P;TofId",23100,0,23100);*/
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarMixedMaker::Finish()
{
  mFile->cd();
  if(QA){
    hVxVyVz->Write();
    hVz->Write();
    hVpdVz->Write();
    hVr->Write();
    hVzVpdVz->Write();
    hnEvsEtavsVz->Write();
    hnEvsPhivsVz->Write();
    hnTofMulvsRef->Write(); 
    hnTofMatvsRef->Write();
    hnTofHitvsRef->Write();
    hnTofMult->Write();
    hnTofMatch->Write();
    hevt->Write();
    h_Vx_Vy->Write();
    h_Vz_VpdVz->Write();
    hevtcut->Write();
    hevtbadcut->Write();
    hpassevtcut->Write(); 
    hrefmult->Write();
    hrefmult_Pos->Write();
    hrefmult_Neg->Write();
    hNsigEvsinvBeta->Write();
    //tracl level QA
    hnHitsFit->Write();
    hnHitsFit_cut->Write();
    hgDca->Write();
    hpDca->Write();
    hinvBetavsP->Write();
    TofId->Write();
    TofId_nSigmaPi->Write();
    // hinvBetavsY->Write();
    hdEdx->Write();
    h_mTpc->Write();
    hpt_Pos->Write();
    hpt_Pos_cut->Write();
    hpt_Neg->Write();
    hpt_Neg_cut->Write();
    hGpt_Pos->Write();
    hGpt_Pos_cut->Write();
    hGpt_Neg->Write();
    hGpt_Neg_cut->Write();
    hEta->Write();
    hEta_cut->Write();
    hPhi->Write();
    hPhi_cut->Write();
    hnHitsPoss->Write();
    hnHitsPoss_cut->Write();
    hnHitsDedx->Write();
    hnHitsDedx_cut->Write();
    h_nHitsDedx_p->Write();
    h_pT_Eta->Write();
    h_pT_Phi->Write();
    h_EtavsPhi_Pos->Write();
    h_EtavsPhi_Neg->Write();
    hBadTofId->Write();
    h_Vz_btofZLocal->Write();
    h_Vz_btofYLocal->Write();
    h_bemcdz->Write();
    h_bemcDphi->Write();
    h_nSMDphi->Write();
    h_nSMDeta->Write();
    h_btowPhiDz->Write();
    h_btowEtaDz->Write();
    h_mtdDeltaY->Write();
    h_mtdDeltaZ->Write();
    h_mtdDeltaTOF->Write();
    /*ModuleId_1->Write();
    TofId_1->Write();
    TrayId_1->Write();
    ModuleId_2->Write();
    TofId_2->Write();
    ModuleId_3->Write();
    TofId_3->Write();
    ModuleId_4->Write();
    TofId_4->Write();
    ModuleId_5->Write();
    TofId_5->Write();*/
    h_ETof_beta->Write();
    h_ETof_betavsP->Write();
    h_ETof_deltaX->Write();
    h_ETof_deltaY->Write();
  
    if (mRunbyRunQA) {
      pVpdVz->Write();
      pVzVpdVz->Write();
      pRefmult->Write();
      pVpdHitEast->Write();
      pVpdHitWest->Write();
      pVz->Write();
      pVx->Write();
      pVy->Write();
      pVr->Write();
      pTof->Write(); 
      pDedx->Write();
      pgDCA->Write();
      ppDCA->Write();
      pgPt->Write();
      pgPhi->Write();
      pgEta->Write();
      pNFits->Write();
      ppPt->Write();
      ppEta->Write();
      ppPhi->Write();
      p_nSigmaE->Write();
      p_nSigmaPion->Write();
      p_nSigmaKaon->Write();
      p_nSigmaProton->Write();
      p_nTofMatch->Write();
      p_nSigmaTofE->Write();
      p_nSigmaTofPion->Write();
      p_nSigmaTofKaon->Write();
      p_nSigmaTofProton->Write();
      p_btofYLocal->Write();
      p_btofZLocal->Write();
      p_nETofHits->Write();
      p_nBemcPidTraits->Write();
      p_nMtdHits->Write();
      p_nMtdPidTraits->Write();
      p_ETof_beta->Write();
      p_ETof_deltaX->Write();
      p_ETof_deltaY->Write();
      ph_p_E0->Write();
      ph_bemcdz->Write();
      ph_bemcDphi->Write();
      ph_nSMDphi->Write();
      ph_nSMDeta->Write();
      ph_btowPhiDz->Write();
      ph_btowEtaDz->Write();
      ph_mtdDeltaY->Write();
      ph_mtdDeltaZ->Write();
      ph_mtdDeltaTOF->Write();
    }
  }
  hMeeCount->Write();
  hMeeCount_like1->Write();
  hMeeCount_like2->Write();
  hMeeCountPt->Write();
  hMeeCountPt_like1->Write();
  hMeeCountPt_like2->Write();
  hMkpiCount->Write();
  hMkpiCount_like1->Write();
  hMkpiCount_like2->Write();
  hMkpiCountPt->Write();
  hMkpiCountPt_like1->Write();
  hMkpiCountPt_like2->Write();
  h_nSigmaElectron_P->Write();
  h_nSigmaElectron_P_tpc->Write();

  h_nSigmaElectron_Inclusive->Write();
  h_nSigmaPion_Inclusive->Write();
  h_nSigmaKaon_Inclusive->Write();
  h_nSigmaProton_Inclusive->Write();
  
  h_nSigmaVsPcharge_Electron->Write();
  h_nSigmaVsPcharge_Pion->Write();
  h_nSigmaVsPcharge_Kaon->Write();
  h_nSigmaVsPcharge_Proton->Write();
  
  // TPC after TOF cuts
  h_TpcNsigmaE_afterTofCut->Write();
  h_TpcNsigmaPi_afterTofCut->Write();
  h_TpcNsigmaK_afterTofCut->Write();
  h_TpcNsigmaP_afterTofCut->Write();

  // TOF after TPC cuts
  h_TofNsigmaE_afterTpcCut->Write();
  h_TofNsigmaPi_afterTpcCut->Write();
  h_TofNsigmaK_afterTpcCut->Write();
  h_TofNsigmaP_afterTpcCut->Write();

  // TPC after TOF cuts (2D)
  h2_TpcNsigmaE_vs_p_afterTofCut->Write();
  h2_TpcNsigmaPi_vs_p_afterTofCut->Write();
  h2_TpcNsigmaK_vs_p_afterTofCut->Write();
  h2_TpcNsigmaP_vs_p_afterTofCut->Write();

  // TOF after TPC cuts (2D)
  h2_TofNsigmaE_vs_p_afterTpcCut->Write();
  h2_TofNsigmaPi_vs_p_afterTpcCut->Write();
  h2_TofNsigmaK_vs_p_afterTpcCut->Write();
  h2_TofNsigmaP_vs_p_afterTpcCut->Write();  

  // =================================================================

  mFile->Close();

  /*mFile_RunID->cd();
  int x_RunID=0;
  for(x_RunID=0;x_RunID<78;x_RunID++)
  {
    hinvBetavsP_RunID[x_RunID]->Write();
  }
  mFile_RunID->Close();*/

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarMixedMaker::Make()
{
  if(DEBUG) cout<<"star make"<<endl;
  ParticleInfo particleinfo;
  vector<ParticleInfo> electroninfo;
  vector<ParticleInfo> positroninfo;
  vector<ParticleInfo> kaoninfo_pos;
  vector<ParticleInfo> kaoninfo_neg;
  vector<ParticleInfo> pioninfo_pos;
  vector<ParticleInfo> pioninfo_neg;
  // StMemStat mem;
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoDstarMixedMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }
  StPicoDst const* picoDst = mPicoDstMaker->picoDst();
  if (!picoDst)
  {
    LOG_WARN << "StPicoDstarMixedMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  // -------------- USER ANALYSIS -------------------------
  StPicoEvent const * picoEvent = picoDst->event();
  //trigger
  

  //if (!isGoodTrigger(picoEvent)) return 0;    
  mRunId = picoEvent->runId();
  
  //if(mRunId < 22057001) return 0;   

  if(QA) hevt->Fill(runnum[mRunId]);

  //if(mRunId < 22043001) return 0;// form 20210212

    TVector3 pVtx = picoEvent->primaryVertex();
  if (mRunbyRunQA && isGoodQaEvent(picoEvent) && QA ){ 
    //primary vertex
    // StThreeVectorF pVtx = picoEvent->primaryVertex();
    if(DEBUG) cout<<"star runbyrun QA"<<endl;
    mVx = pVtx.x();
    mVy = pVtx.y();
    mVz = pVtx.z();
    mVr = sqrt(mVx*mVx+mVy*mVy);
    
    mVpdVz = picoEvent->vzVpd();
    if(DEBUG) cout<<"runbyrun VpdVz: "<<mVpdVz<<endl;
    if(DEBUG) cout<<"runbyrun Vz: "<<mVz<<endl;
    mRefmult = picoEvent->refMult();
    mVpdHitEast = picoEvent->nVpdHitsEast();
    mVpdHitWest = picoEvent->nVpdHitsWest();
    if (DEBUG) cout<<"start filling "<<mRunId<<" "<<runnum[mRunId]<<endl; 
    if(fabs(mVpdVz) < 200){//fill event level profile
    pVpdVz->Fill(runnum[mRunId],mVpdVz);
    pVzVpdVz->Fill(runnum[mRunId],mVpdVz-mVz);
    }
    pRefmult->Fill(runnum[mRunId],mRefmult);
    pVpdHitEast->Fill(runnum[mRunId],mVpdHitEast);
    pVpdHitWest->Fill(runnum[mRunId],mVpdHitWest);
    pVx->Fill(runnum[mRunId],mVx);
    pVy->Fill(runnum[mRunId],mVy);
    pVz->Fill(runnum[mRunId],mVz);
    pVr->Fill(runnum[mRunId],mVr);

    p_nTofMatch->Fill(runnum[mRunId],picoEvent->nBTOFMatch());
    p_nETofHits->Fill(runnum[mRunId],picoDst->numberOfETofHits());
    p_nBemcPidTraits->Fill(runnum[mRunId],picoDst->numberOfBEmcPidTraits());
    p_nMtdPidTraits->Fill(runnum[mRunId],picoDst->numberOfMtdPidTraits());
    p_nMtdHits->Fill(runnum[mRunId],picoDst->numberOfMtdHits());
    //track level 

    int nTracks = picoDst->numberOfTracks();
    //if (DEBUG)  cout << nTracks <<endl;
    //global
    for (int itrack=0;itrack<nTracks;itrack++){
      StPicoTrack* trk = picoDst->track(itrack);
      bool goodQAtrack = isGoodQaTrack(trk); 
      if (!goodQAtrack) continue;
      if (!(fabs(trk->gDCA(pVtx.x(),pVtx.y(),pVtx.z()))<anaCuts::qaDca)) continue;
      bool isprimary = trk->isPrimary();
      double ptot = trk->gMom(pVtx, picoEvent->bField()).Mag();
      float beta = getTofBeta(trk);
      bool tofmatch =( beta!=std::numeric_limits<float>::quiet_NaN() )&& beta>0;
      bool tofQaPion = false;
      if (tofmatch) {
        float beta_pi = ptot / sqrt(ptot * ptot + M_PION_PLUS * M_PION_PLUS);
        tofQaPion = fabs(1. / beta - 1. / beta_pi) < anaCuts::qaTofPion;
      }
      bool tpcQaPion = fabs(trk->nSigmaPion()) < anaCuts::qaTpcPion;
      //global
      //runbyrunQA
      //pid performance 
      //if (tofQaPion && ptot<1) pTof->Fill(runnum[mRunId],(1./beta)); 
      //if (tpcQaPion && ptot<0.5) {
      //  pDedx->Fill(runnum[mRunId],trk->dEdx());
     // }
      //global tracking performance
      pgPt->Fill(runnum[mRunId],trk->gMom().Perp());
      pgDCA->Fill(runnum[mRunId],trk->gDCA(pVtx.x(),pVtx.y(),pVtx.z()));
      pgPhi->Fill(runnum[mRunId],trk->gMom().Phi());
      pgEta->Fill(runnum[mRunId],trk->gMom().Eta());

      //primary 
      if (isprimary){
        ppDCA->Fill(runnum[mRunId],trk->gDCA(pVtx.x(),pVtx.y(),pVtx.z()));
        if(tofmatch){pTof->Fill(runnum[mRunId],(1./beta));} 
        pDedx->Fill(runnum[mRunId],trk->dEdx());
        pNFits->Fill(runnum[mRunId],trk->nHitsFit()); 
      
        ppPt->Fill(runnum[mRunId],trk->pMom().Perp());
        ppEta->Fill(runnum[mRunId],trk->pMom().Eta());
        ppPhi->Fill(runnum[mRunId],trk->pMom().Phi());
        if(fabs(pow(trk->pMom().Mag()*sqrt(1-beta*beta)*1.0/beta,2) - 0.019)<0.003 && fabs(trk->nSigmaPion())<5){
           p_nSigmaPion->Fill(runnum[mRunId],trk->nSigmaPion());
        }
        if(fabs(pow(trk->pMom().Mag()*sqrt(1-beta*beta)*1.0/beta,2) - 0.243)<0.005 && fabs(trk->nSigmaKaon())<5){
           p_nSigmaKaon->Fill(runnum[mRunId],trk->nSigmaKaon());
        }
        if(fabs(pow(trk->pMom().Mag()*sqrt(1-beta*beta)*1.0/beta,2) - 0.879)<0.02 && fabs(trk->nSigmaProton())<5){
           p_nSigmaProton->Fill(runnum[mRunId],trk->nSigmaProton());
        }
        


      int pbemcId = trk->bemcPidTraitsIndex();
      if(pbemcId>=0){
         StPicoBEmcPidTraits * pbemctrait = picoDst->bemcPidTraits(pbemcId);
         if(pbemctrait){
            float pnSMDphi = pbemctrait->bemcSmdNPhi();
            ph_nSMDphi->Fill(runnum[mRunId],pnSMDphi);
            float pnSMDeta = pbemctrait->bemcSmdNEta();
            ph_nSMDeta->Fill(runnum[mRunId],pnSMDeta);
            float pbtowPhiDz = pbemctrait->btowPhiDist();
            ph_btowPhiDz->Fill(runnum[mRunId],pbtowPhiDz);
            float pbtowEtaDz = pbemctrait->btowEtaDist();
            ph_btowEtaDz->Fill(runnum[mRunId],pbtowEtaDz);
            float pbemcdz = pbemctrait->bemcZDist();
            ph_bemcdz->Fill(runnum[mRunId],pbemcdz);
            float pbemcDphi = pbemctrait->bemcPhiDist();
            ph_bemcDphi->Fill(runnum[mRunId],pbemcDphi);
            float pE0 = pbemctrait->bemcE0(); 
            ph_p_E0->Fill(runnum[mRunId],trk->gMom().Mag()*1.0/pE0); 
           }
      }
    int pmtdId = trk->mtdPidTraitsIndex();
    if(pmtdId>=0){
       StPicoMtdPidTraits * pmtdtrait = picoDst->mtdPidTraits(pmtdId);
       if(pmtdtrait){
          float pmtdDeltaY = pmtdtrait->deltaY();
          ph_mtdDeltaY->Fill(runnum[mRunId],pmtdDeltaY);
          float pmtdDeltaZ = pmtdtrait->deltaZ();
          ph_mtdDeltaZ->Fill(runnum[mRunId],pmtdDeltaZ);
          float pmtdDeltaTOF = pmtdtrait->deltaTimeOfFlight();
          ph_mtdDeltaTOF->Fill(runnum[mRunId],pmtdDeltaTOF);
       }
    }

    int pbtofId = trk->bTofPidTraitsIndex();
    if(pbtofId>=0 && tofmatch){
       StPicoBTofPidTraits * pbtoftrait = picoDst->btofPidTraits(pbtofId);
       if(pbtoftrait){
         // float pbtofnSigmaE = pbtoftrait->nSigmaElectron();
          float pbtofnSigmaPion = pbtoftrait->nSigmaPion();
          float pbtofnSigmaKaon = pbtoftrait->nSigmaKaon();
          float pbtofnSigmaProton = pbtoftrait->nSigmaProton();
          float pbtofYLocal = pbtoftrait->btofYLocal();
          float pbtofZLocal = pbtoftrait->btofZLocal();
          //float pbtofdeltaY = pbtoftrait->deltaY(); 
          //cout<<"btofdeltaY: "<<pbtofdeltaY<<endl;

          p_btofYLocal->Fill(runnum[mRunId],pbtofYLocal);  
          p_btofZLocal->Fill(runnum[mRunId],pbtofZLocal);  

          h_Vz_btofZLocal->Fill(mVz,pbtofZLocal);
          h_Vz_btofYLocal->Fill(mVz,pbtofYLocal);

        if(fabs(pow(trk->pMom().Mag()*sqrt(1-beta*beta)*1.0/beta,2) - 0.019)<0.003 && fabs(trk->nSigmaPion())<5){
           p_nSigmaTofPion->Fill(runnum[mRunId],pbtofnSigmaPion);
        }
        if(fabs(pow(trk->pMom().Mag()*sqrt(1-beta*beta)*1.0/beta,2) - 0.243)<0.005 && fabs(trk->nSigmaKaon())<5){
           p_nSigmaTofKaon->Fill(runnum[mRunId],pbtofnSigmaKaon);
        }
        if(fabs(pow(trk->pMom().Mag()*sqrt(1-beta*beta)*1.0/beta,2) - 0.879)<0.02 && fabs(trk->nSigmaProton())<5){
           p_nSigmaTofProton->Fill(runnum[mRunId],pbtofnSigmaProton);
        }
       }
    }
     
    int petofId = trk->eTofPidTraitsIndex();
    if(petofId>=0){
      StPicoETofPidTraits * petoftrait = picoDst->etofPidTraits(petofId);
      //int etofmatch = petoftrait->matchFlag();         
      if(petoftrait){
         float petofbeta = 0;
         float petofdeltaX = 0;
         float petofdeltaY = 0; 
         petofbeta = petoftrait->beta(); 
         petofdeltaX = petoftrait->deltaX(); 
         petofdeltaY = petoftrait->deltaY(); 
         int etofmatch = petoftrait->matchFlag();         
     
         if(etofmatch == 1 || etofmatch == 2){
            if(petofbeta != 0){p_ETof_beta->Fill(runnum[mRunId],1.0/petofbeta);}
            p_ETof_deltaX->Fill(runnum[mRunId],petofdeltaX);
            p_ETof_deltaY->Fill(runnum[mRunId],petofdeltaY);
            h_ETof_beta->Fill(1.0/petofbeta);
            h_ETof_betavsP->Fill((trk->pMom()).Mag(),1.0/petofbeta);
            h_ETof_deltaX->Fill(petofdeltaX);
            h_ETof_deltaY->Fill(petofdeltaY);
         }
      
      }
    }


      }
    }
  } //runbyrun QA 


  electroninfo.clear();
  positroninfo.clear();  
  kaoninfo_pos.clear();
  kaoninfo_neg.clear();
  pioninfo_pos.clear();
  pioninfo_neg.clear();

  //event and track level QA
  if(QA) hpassevtcut->Fill(0);
  if (!isBadrun(mRunId)){
  if(QA) hpassevtcut->Fill(1); //bad run list
  // bool vzcut = fabs(pVtx.z()) < 30;
  bool vzcut = fabs(pVtx.z()) < 60;
  bool verrcut = !(fabs(pVtx.x()) < anaCuts::Verror && fabs(pVtx.y()) < anaCuts::Verror && fabs(pVtx.z()) < anaCuts::Verror);
  bool vrcut =  sqrt(TMath::Power(pVtx.x(), 2) + TMath::Power(pVtx.y(), 2)) <=  anaCuts::Vr ;
  // bool vpdvzcut = fabs(pVtx.z() - picoEvent->vzVpd()) < 3;
  bool vpdvzcut = true;
  if (QA && vzcut) hpassevtcut->Fill(2);
  if (QA && vzcut &&  vrcut) hpassevtcut->Fill(3);
  // if (vzcut && vrcut  &&  vpdvzcut ) hpassevtcut->Fill(4);
  if (QA && vzcut && vrcut  &&  vpdvzcut && verrcut ) hpassevtcut->Fill(4);
  bool refusepileup = picoEvent->refMult()<picoEvent->btofTrayMultiplicity()*0.36+45;
  bool refusebadtof = picoEvent->refMult()>picoEvent->btofTrayMultiplicity()*0.28-115;
  bool passCentralityCut = refusepileup && refusebadtof  && verrcut && vrcut && fabs(pVtx.z()) < 10; 
  if(QA)if (passCentralityCut) hrefmult->Fill(picoEvent->refMult());
  if(QA)if (passCentralityCut) hrefmult_Pos->Fill(picoEvent->refMultPos());
  if(QA)if (passCentralityCut) hrefmult_Neg->Fill(picoEvent->refMultNeg());
  if (isGoodEvent(picoEvent)){
    // StThreeVectorF pVtx = picoEvent->primaryVertex();
    
    
   // if(mRunId < 22106001) return 0;   

    if(DEBUG) cout<<"star event QA"<<endl;
    TVector3 pVtx = picoEvent->primaryVertex();
    mVx = pVtx.x();
    mVy = pVtx.y();
    mVz = pVtx.z();
    mVpdVz = picoEvent->vzVpd();
    if(QA){
    h_Vx_Vy->Fill(mVx,mVy);
      hevtcut->Fill(runnum[mRunId]);
      hVz->Fill(mVz);
      hVpdVz->Fill(mVpdVz);
      hVxVyVz->Fill(mVx,mVy,mVz);
      hVr->Fill(sqrt(mVy*mVy+mVx*mVx));
      hVzVpdVz->Fill(mVpdVz-mVz);
      h_Vz_VpdVz->Fill(mVz,mVpdVz);  
      hnTofMult->Fill(picoEvent->btofTrayMultiplicity());  
      hnTofMulvsRef->Fill(picoEvent->refMult(),picoEvent->btofTrayMultiplicity());  
      hnTofMatch->Fill(picoEvent->nBTOFMatch());  
      hnTofMatvsRef->Fill(picoEvent->refMult(),picoEvent->nBTOFMatch());  
    }
    double ntofhits = 0;
    //    int ntrack_tof_hits =0; 
    int nTracks = picoDst->numberOfTracks();
    for (int itrack=0;itrack<nTracks;itrack++){
      StPicoTrack* trk = picoDst->track(itrack);

      TVector3 mom = trk->pMom();

      if(QA){
        hgDca->Fill(trk->gDCA(mVx,mVy,mVz));
        if(trk->charge()>0){
          hpt_Pos->Fill(mom.Perp());
          hGpt_Pos->Fill(trk->gMom().Perp());
        }else{
          hpt_Neg->Fill(mom.Perp());
          hGpt_Neg->Fill(trk->gMom().Perp());
        }
        hEta->Fill(mom.Eta());
        hPhi->Fill(mom.Phi());
        hnHitsFit->Fill(trk->nHitsFit()*trk->charge());
        hnHitsPoss->Fill(trk->nHitsPoss()*trk->charge());
        hnHitsDedx->Fill(trk->nHitsDedx()*trk->charge());
      }
      bool isprimary = trk->isPrimary();
      bool goodtrack = isGoodTrack(trk,trk->gDCA(mVx,mVy,mVz));
      if (!goodtrack) continue;
      if (!isprimary) continue;
      
      if(QA)hpDca->Fill(trk->gDCA(mVx,mVy,mVz));

      int bemcId = trk->bemcPidTraitsIndex();
      if(bemcId>=0 && QA){
         StPicoBEmcPidTraits * bemctrait = picoDst->bemcPidTraits(bemcId);
         if(bemctrait){
            float nSMDphi = bemctrait->bemcSmdNPhi();
            h_nSMDphi->Fill(nSMDphi);
            float nSMDeta = bemctrait->bemcSmdNEta();
            h_nSMDeta->Fill(nSMDeta);
            float btowPhiDz = bemctrait->btowPhiDist();
            h_btowPhiDz->Fill(btowPhiDz);
            float btowEtaDz = bemctrait->btowEtaDist();
            h_btowEtaDz->Fill(btowEtaDz);
            float bemcdz = bemctrait->bemcZDist();
            h_bemcdz->Fill(bemcdz);
            float bemcDphi = bemctrait->bemcPhiDist();
            h_bemcDphi->Fill(bemcDphi);
            float E0 = bemctrait->bemcE0(); 
            h_p_E0->Fill(trk->gMom().Mag()*1.0/E0); 
          }
      }
     
      int mtdId = trk->mtdPidTraitsIndex();
      if(mtdId>=0 && QA){
        StPicoMtdPidTraits * mtdtrait = picoDst->mtdPidTraits(mtdId);
        if(mtdtrait){
          float mtdDeltaY = mtdtrait->deltaY();
          h_mtdDeltaY->Fill(mtdDeltaY);
          float mtdDeltaZ = mtdtrait->deltaZ();
          h_mtdDeltaZ->Fill(mtdDeltaZ);
          float mtdDeltaTOF = mtdtrait->deltaTimeOfFlight();
          h_mtdDeltaTOF->Fill(mtdDeltaTOF);
          //cout<<"mtdDeltaTOF: "<<mtdDeltaTOF<<endl;
        }
      }

      if(QA){
        // StThreeVectorF mom = trk->pMom();  
        if(trk->charge()>0){
          hpt_Pos_cut->Fill(mom.Perp());
          hGpt_Pos_cut->Fill(trk->gMom().Perp());
          h_EtavsPhi_Pos->Fill(mom.Eta(),mom.Phi());
        }else{
          hpt_Neg_cut->Fill(mom.Perp());
          hGpt_Neg_cut->Fill(trk->gMom().Perp());
          h_EtavsPhi_Neg->Fill(mom.Eta(),mom.Phi());
        }
        // hpDca->Fill(trk->pDca(mVx,mVy,mVz));
        hPhi_cut->Fill(mom.Phi());
        hEta_cut->Fill(mom.Eta());
        h_pT_Eta->Fill(mom.Perp()*trk->charge(),mom.Eta());
        h_pT_Phi->Fill(mom.Perp()*trk->charge(),mom.Phi());
        hnHitsFit_cut->Fill(trk->nHitsFit()*trk->charge());
        hnHitsPoss_cut->Fill(trk->nHitsPoss()*trk->charge());
        hnHitsDedx_cut->Fill(trk->nHitsDedx()*trk->charge());
        h_nHitsDedx_p->Fill(mom.Mag()*trk->charge(),trk->nHitsDedx());
      }
      
      double beta = getTofBeta(trk);
      bool tofmatch = (beta!=std::numeric_limits<float>::quiet_NaN()) && beta>0;
      
      double p = mom.Mag();
      int charge = trk->charge();
      double nSigmaE = trk->nSigmaElectron();
      double nSigmaPi = trk->nSigmaPion();
      double nSigmaK = trk->nSigmaKaon();
      double nSigmaP = trk->nSigmaProton();
      // Fill 1D inclusive histograms
      h_nSigmaElectron_Inclusive->Fill(nSigmaE);
      h_nSigmaPion_Inclusive->Fill(nSigmaPi);
      h_nSigmaKaon_Inclusive->Fill(nSigmaK);
      h_nSigmaProton_Inclusive->Fill(nSigmaP);
      
      // Fill 2D inclusive histograms vs p*charge
      h_nSigmaVsPcharge_Electron->Fill(p * charge, nSigmaE);
      h_nSigmaVsPcharge_Pion->Fill(p * charge, nSigmaPi);
      h_nSigmaVsPcharge_Kaon->Fill(p * charge, nSigmaK);
      h_nSigmaVsPcharge_Proton->Fill(p * charge, nSigmaP);

	//choose inclusive electron
      // bool isTPCElectron =  trk->nSigmaElectron()<2 && trk->nSigmaElectron()>0.75;
      bool isTPCElectron=0;
      if (mom.Mag()>0.8) isTPCElectron =  trk->nSigmaElectron()<2 && trk->nSigmaElectron()>-0.75;
      else isTPCElectron = trk->nSigmaElectron()<2 && trk->nSigmaElectron()>(3*mom.Mag()-3.15);
      bool isTOFElectron = tofmatch?fabs(1./beta-1.)<0.025:false;
  
      h_nSigmaElectron_P_tpc->Fill(mom.Mag(),trk->nSigmaElectron());      

      // --- Pion ---
      bool isTpcPion = fabs(nSigmaPi) < 2.0;
      bool isTofPion = false;
      if (tofmatch && p > 0.1) {
          float oneOverBeta_pion_expected = sqrt(p * p + M_PION_PLUS * M_PION_PLUS) / p;
          if (fabs(1.0 / beta - oneOverBeta_pion_expected) < 0.03) {
              isTofPion = true;
          }
      }

      // --- Kaon ---
      bool isTpcKaon = fabs(nSigmaK) < 2.0;
      bool isTofKaon = false;
      if (tofmatch && p > 0.1) {
          float oneOverBeta_kaon_expected = sqrt(p * p + M_KAON_PLUS * M_KAON_PLUS) / p;
          if (fabs(1.0 / beta - oneOverBeta_kaon_expected) < 0.03) {
              isTofKaon = true;
          }
      }

      // --- Proton ---
      bool isTpcProton = fabs(nSigmaP) < 2.0;
      bool isTofProton = false;
      if (tofmatch && p > 0.1) {
          float oneOverBeta_proton_expected = sqrt(p * p + M_PROTON * M_PROTON) / p;
          if (fabs(1.0 / beta - oneOverBeta_proton_expected) < 0.03) {
              isTofProton = true;
          }
      }

      // --- Fill TPC nSigma if TOF cut is passed ---
      if (isTOFElectron) {h_TpcNsigmaE_afterTofCut->Fill(nSigmaE);h2_TpcNsigmaE_vs_p_afterTofCut->Fill(p * charge, nSigmaE);}
      if (isTofPion) {h_TpcNsigmaPi_afterTofCut->Fill(nSigmaPi);h2_TpcNsigmaPi_vs_p_afterTofCut->Fill(p * charge, nSigmaPi);}
      if (isTofKaon) {h_TpcNsigmaK_afterTofCut->Fill(nSigmaK);h2_TpcNsigmaK_vs_p_afterTofCut->Fill(p * charge, nSigmaK);}
      if (isTofProton) {h_TpcNsigmaP_afterTofCut->Fill(nSigmaP);h2_TpcNsigmaP_vs_p_afterTofCut->Fill(p * charge, nSigmaP);}


      if (isTOFElectron && isTPCElectron) {
        if(QA) hnEvsEtavsVz->Fill(mom.Eta(),mVz); 
        if(QA) hnEvsPhivsVz->Fill(mom.Phi(),mVz);
        
        if(QA) p_nSigmaE->Fill(runnum[mRunId],trk->nSigmaElectron());
        //p_nSigmaTofE->Fill(runnum[mRunId],pbtofnSigmaE);
        
        if(trk->charge()<0 && tofmatch)
         {
            particleinfo.charge = trk->charge();
            //cout<<"charge: "<<trk->charge()<<endl;
            particleinfo.pt = mom.Perp();
            particleinfo.eta = mom.Eta();
            particleinfo.phi = mom.Phi();
            particleinfo.p = mom.Mag();
            particleinfo.nSigmaPi = trk->nSigmaPion();
            particleinfo.beta = beta;
            particleinfo.energy = sqrt(pow(M_electron,2.0)+pow(mom.Mag(),2.0));
            particleinfo.p1 = mom.X();
            particleinfo.p2 = mom.Y();
            particleinfo.p3 = mom.Z();
            electroninfo.push_back(particleinfo);
            //            cout<<"debug01"<<endl;
            /* current_eMinus[current_nEMinus].SetPx(mom.x());
             current_eMinus[current_nEMinus].SetPy(mom.y());
             current_eMinus[current_nEMinus].SetPz(mom.z());
             current_eMinus[current_nEMinus].SetE(sqrt(pow(M_electron,2.0)+pow(mom.Mag(),2.0)));
             current_nEMinus++;*/
        }     

        if(trk->charge()>0 && tofmatch){

            particleinfo.charge = trk->charge();
            // cout<<"charge: "<<trk->charge()<<endl;
            particleinfo.pt = mom.Perp();
            particleinfo.eta = mom.Eta();
            particleinfo.phi = mom.Phi();
            particleinfo.p = mom.Mag();
            particleinfo.nSigmaPi = trk->nSigmaPion();
            particleinfo.beta = beta;
            particleinfo.energy = sqrt(pow(M_electron,2.0)+pow(mom.Mag(),2.0));
            particleinfo.p1 = mom.X();
            particleinfo.p2 = mom.Y();
            particleinfo.p3 = mom.Z();
            positroninfo.push_back(particleinfo);
            // cout<<"debug02"<<endl;
             /*current_ePlus[current_nEPlus].SetPx(mom.x());
             current_ePlus[current_nEPlus].SetPy(mom.y());
             current_ePlus[current_nEPlus].SetPz(mom.z());
             current_ePlus[current_nEPlus].SetE(sqrt(pow(M_electron,2.0)+pow(mom.Mag(),2.0)));
             current_nEPlus++;*/
        }
      }
         //current_nE++;

      // ---- Store Kaon candidates for D0 analysis ----
      if (isTofKaon && isTpcKaon && tofmatch) {
        if (charge < 0) {
          ParticleInfo kineg;
          kineg.charge = charge;
          kineg.pt = mom.Perp();
          kineg.eta = mom.Eta();
          kineg.phi = mom.Phi();
          kineg.p = mom.Mag();
          kineg.nSigmaPi = nSigmaPi;
          kineg.beta = beta;
          kineg.energy = sqrt(pow(M_kaon,2.0)+pow(mom.Mag(),2.0));
          kineg.p1 = mom.X();
          kineg.p2 = mom.Y();
          kineg.p3 = mom.Z();
          kaoninfo_neg.push_back(kineg);
        } else if (charge > 0) {
          ParticleInfo kipos;
          kipos.charge = charge;
          kipos.pt = mom.Perp();
          kipos.eta = mom.Eta();
          kipos.phi = mom.Phi();
          kipos.p = mom.Mag();
          kipos.nSigmaPi = nSigmaPi;
          kipos.beta = beta;
          kipos.energy = sqrt(pow(M_kaon,2.0)+pow(mom.Mag(),2.0));
          kipos.p1 = mom.X();
          kipos.p2 = mom.Y();
          kipos.p3 = mom.Z();
          kaoninfo_pos.push_back(kipos);
        }
      }

      // ---- Store Pion candidates for D0 analysis ----
      if (isTofPion && isTpcPion) {
        if (charge < 0) {
          ParticleInfo pineg;
          pineg.charge = charge;
          pineg.pt = mom.Perp();
          pineg.eta = mom.Eta();
          pineg.phi = mom.Phi();
          pineg.p = mom.Mag();
          pineg.nSigmaPi = nSigmaPi;
          pineg.beta = beta;
          pineg.energy = sqrt(pow(M_pion,2.0)+pow(mom.Mag(),2.0));
          pineg.p1 = mom.X();
          pineg.p2 = mom.Y();
          pineg.p3 = mom.Z();
          pioninfo_neg.push_back(pineg);
        } else if (charge > 0) {
          ParticleInfo pipos;
          pipos.charge = charge;
          pipos.pt = mom.Perp();
          pipos.eta = mom.Eta();
          pipos.phi = mom.Phi();
          pipos.p = mom.Mag();
          pipos.nSigmaPi = nSigmaPi;
          pipos.beta = beta;
          pipos.energy = sqrt(pow(M_pion,2.0)+pow(mom.Mag(),2.0));
          pipos.p1 = mom.X();
          pipos.p2 = mom.Y();
          pipos.p3 = mom.Z();
          pioninfo_pos.push_back(pipos);
        }
      } 

      if (tofmatch) {
        ntofhits++;
        if(QA) hinvBetavsP->Fill(mom.Mag(),1./beta);
        
        if(fabs(1.0/beta - 1) < 0.025)
        {
         h_nSigmaElectron_P->Fill(mom.Mag(),trk->nSigmaElectron()); 
        }       

       // hinvBetavsP_RunID[runnum[mRunId]]->Fill(mom.Mag(),1./beta);

    //    bool istoftrack = trk->isTofTrack();
    //    if(istoftrack && 0.8<1.0/beta && 1.0/beta<0.9 && mom.Mag()>0.4){
    //        ntrack_tof_hits = trk->numberOfBTofHits();
    //        for(int itrack_tof_hits=0;itrack_tof_hits<ntrack_tof_hits;itrack_tof_hits++)
    //           {
    //            StPicoBTofHit* tofhit = trk->btofHit(itrack_tof_hits);
    //            ModuleId_1->Fill(tofhit->module()); 
    //           }
    //   }

        int index2tof = trk->bTofPidTraitsIndex();
        StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
        if (1./beta<0.88 && mom.Mag()>0.8 && mom.Mag()<1.5) {
          int tofid = tofPid->btofCellId();
        //  TofId_1->Fill(tofid);
          if(QA) hBadTofId->Fill(tofid/192+1,(tofid%192)/6+1,tofid%6+1);
        }

        if(tofPid){
          if (isTPCElectron) {h_TofNsigmaE_afterTpcCut->Fill(tofPid->nSigmaElectron());h2_TofNsigmaE_vs_p_afterTpcCut->Fill(p * charge, tofPid->nSigmaElectron());}
          if (isTpcPion) {h_TofNsigmaPi_afterTpcCut->Fill(tofPid->nSigmaPion());h2_TofNsigmaPi_vs_p_afterTpcCut->Fill(p * charge, tofPid->nSigmaPion());}
          if (isTpcKaon) {h_TofNsigmaK_afterTpcCut->Fill(tofPid->nSigmaKaon());h2_TofNsigmaK_vs_p_afterTpcCut->Fill(p * charge, tofPid->nSigmaKaon());}
          if (isTpcProton) {h_TofNsigmaP_afterTpcCut->Fill(tofPid->nSigmaProton());h2_TofNsigmaP_vs_p_afterTpcCut->Fill(p * charge, tofPid->nSigmaProton());}
        }


        if(1.0/beta>1.13 && 1.0/beta<1.24 && mom.Mag()>0.3 && mom.Mag()<0.5)
         {
            int tofid = tofPid->btofCellId();
            if(QA) TofId->Fill(tofid);
            if(QA) TofId_nSigmaPi->Fill(trk->nSigmaPion());
            //TrayId_1->Fill(tofid/192+1);
            //ModuleId_1->Fill((tofid%192)/6+1);
         }

         
        if (isTOFElectron && isTPCElectron) {

        if(QA) p_nSigmaTofE->Fill(runnum[mRunId],tofPid->nSigmaElectron());
        }
       /* if(1.0/beta>0.8 && 1.0/beta<0.9 && mom.Mag()>0.4 && mom.Mag()<3)
         {
            int tofid = tofPid->btofCellId();
            TofId_1->Fill(tofid);
            TrayId_1->Fill(tofid/192+1);
            ModuleId_1->Fill((tofid%192)/6+1);
         }

        if(1.0/beta>0.82 && 1.0/beta<0.9 && mom.Mag()>0.4 && mom.Mag()<3)
         {
           int  tofid = tofPid->btofCellId();
          TofId_2->Fill(tofid);
            ModuleId_2->Fill((tofid%192)/6+1);
         }
        if(1.0/beta>0.82 && 1.0/beta<0.88 && mom.Mag()>0.4 && mom.Mag()<3)
         {
           int  tofid = tofPid->btofCellId();
          TofId_3->Fill(tofid);
            ModuleId_3->Fill((tofid%192)/6+1);
         }
        if(1.0/beta>0.84 && 1.0/beta<0.88 && mom.Mag()>0.4 && mom.Mag()<3)
         {
           int   tofid = tofPid->btofCellId();
          TofId_4->Fill(tofid);
            ModuleId_4->Fill((tofid%192)/6+1);
         }
        if(1.0/beta>0.84 && 1.0/beta<0.9 && mom.Mag()>0.4 && mom.Mag()<3)
         {
           int  tofid = tofPid->btofCellId();
          TofId_5->Fill(tofid);
            ModuleId_5->Fill((tofid%192)/6+1);
         }*/
        if(QA) hNsigEvsinvBeta->Fill(trk->nSigmaElectron(),1./beta,mom.Mag());
      }
      if(QA) hdEdx->Fill(mom.Mag()*trk->charge(),trk->dEdx());
      if(QA) h_mTpc->Fill(mom.Mag()*trk->charge(),pow(mom.Mag()*sqrt(1-beta*beta)*1.0/beta,2));
    }
    if(QA) hnTofHitvsRef->Fill(ntofhits,picoEvent->refMult());
    
      int x=0;
      int y=0;
      int num_electron = electroninfo.size();
      int num_positron = positroninfo.size();
      float inv_mass=0;
      TVector3 momentum_particle;
      TLorentzVector eepair(0,0,0,0);
      TLorentzVector particle1_4V(0,0,0,0);
      TLorentzVector particle2_4V(0,0,0,0);
      for(x=0;x<num_electron;x++)
         {
             particle1_4V.SetPx(electroninfo[x].p1);
             particle1_4V.SetPy(electroninfo[x].p2);
             particle1_4V.SetPz(electroninfo[x].p3);
             particle1_4V.SetE(electroninfo[x].energy);
               for(y=x+1;y<num_electron;y++)
                  {
                    particle2_4V.SetPx(electroninfo[y].p1);
                    particle2_4V.SetPy(electroninfo[y].p2);
                    particle2_4V.SetPz(electroninfo[y].p3);
                    particle2_4V.SetE(electroninfo[y].energy);
                    eepair = particle1_4V + particle2_4V;
                    
                    //if(eepair.Perp()<0.2){hMeeCount_like1->Fill(eepair.M());}
                    hMeeCount_like1->Fill(eepair.M());
                    hMeeCountPt_like1->Fill(eepair.M(),eepair.Perp());
//                    cout<<"debug03"<<endl;
                    //if(eepair.Perp()<=10){hMeelike1_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M());}
                  }
         }
     
      for(x=0;x<num_positron;x++)
         {
             particle1_4V.SetPx(positroninfo[x].p1);
             particle1_4V.SetPy(positroninfo[x].p2);
             particle1_4V.SetPz(positroninfo[x].p3);
             particle1_4V.SetE(positroninfo[x].energy);
               for(y=x+1;y<num_positron;y++)
                  {
                    particle2_4V.SetPx(positroninfo[y].p1);
                    particle2_4V.SetPy(positroninfo[y].p2);
                    particle2_4V.SetPz(positroninfo[y].p3);
                    particle2_4V.SetE(positroninfo[y].energy);
                    eepair = particle1_4V + particle2_4V;
                    //if(eepair.Perp()<0.2){hMeeCount_like2->Fill(eepair.M());}
                    hMeeCount_like2->Fill(eepair.M());
                    hMeeCountPt_like2->Fill(eepair.M(),eepair.Perp());
                    //if(eepair.Perp()<=10){hMeelike2_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M());}
                  }
         }
      for(x=0;x<num_positron;x++)
         {
             particle1_4V.SetPx(positroninfo[x].p1);
             particle1_4V.SetPy(positroninfo[x].p2);
             particle1_4V.SetPz(positroninfo[x].p3);
             particle1_4V.SetE(positroninfo[x].energy);
               for(y=0;y<num_electron;y++)
                  {
                    particle2_4V.SetPx(electroninfo[y].p1);
                    particle2_4V.SetPy(electroninfo[y].p2);
                    particle2_4V.SetPz(electroninfo[y].p3);
                    particle2_4V.SetE(electroninfo[y].energy);
                    eepair = particle1_4V + particle2_4V;
                    //if(eepair.Perp()<0.2){hMeeCount->Fill(eepair.M());}
                    hMeeCount->Fill(eepair.M());
                    hMeeCountPt->Fill(eepair.M(),eepair.Perp());
                    //if(eepair.Perp()<=10){hMee_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M());}
                  }
         }

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++ ADD THE NEW D0 RECONSTRUCTION LOGIC HERE +++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TVector3 pVtx = picoEvent->primaryVertex();
double bField = picoEvent->bField();
int nTracks = picoDst->numberOfTracks();

  // --- Loop over all track pairs to form D0 candidates ---
  for (int i = 0; i < nTracks; ++i) {
      StPicoTrack* trk1 = picoDst->track(i);

      // Apply daughter track cuts
      if (!isGoodTrack(trk1, trk1->gDCA(pVtx))) continue;
      if (trk1->pMom().Perp() < 0.4) continue; // Daughter pT > 400 MeV

      for (int j = i + 1; j < nTracks; ++j) { // Start from i+1 to avoid double counting
          StPicoTrack* trk2 = picoDst->track(j);

          // Apply daughter track cuts
          if (!isGoodTrack(trk2, trk2->gDCA(pVtx))) continue;
          if (trk2->pMom().Perp() < 0.4) continue; // Daughter pT > 400 MeV

          // --- Check charge sign ---
          bool isOppositeSign = (trk1->charge() * trk2->charge() < 0);
          bool isLikeSign = (trk1->charge() * trk2->charge() > 0);

          if (!isOppositeSign && !isLikeSign) continue; // Skip neutral tracks

          // --- Particle Identification ---
          StPicoTrack *trk_K = nullptr, *trk_pi = nullptr;

          if (isKaon(trk1) && isPion(trk2)) {
              trk_K = trk1;
              trk_pi = trk2;
          } else if (isPion(trk1) && isKaon(trk2)) {
              trk_K = trk2;
              trk_pi = trk1;
          } else {
              continue; // Pair is not a K-pi candidate
          }

          // --- Topological Cuts (Applied to both Signal and Background) ---
          StPicoPhysicalHelix kaonHelix = trk_K->helix(bField);
          StPicoPhysicalHelix pionHelix = trk_pi->helix(bField);

          pair<double, double> s = kaonHelix.pathLengths(pionHelix);
          TVector3 dcaVtx = (kaonHelix.at(s.first) + pionHelix.at(s.second)) * 0.5;
          double dcaDaughters = (kaonHelix.at(s.first) - pionHelix.at(s.second)).Mag();
          
          // These cuts are crucial for signal
          if (dcaDaughters > 0.01) continue; // DCA between daughters < 100 um

          TVector3 d0Mom = trk_K->pMom() + trk_pi->pMom();
          double decayLength = (dcaVtx - pVtx).Mag();
          if (decayLength < 0.02) continue; // Decay length > 200 um

          double cosPointingAngle = d0Mom.Dot(dcaVtx - pVtx) / (d0Mom.Mag() * (dcaVtx - pVtx).Mag());
          if (cosPointingAngle < 0.98) continue;

          if (trk_K->gDCA(pVtx) < 0.01) continue; // Kaon DCA to PV > 100 um
          if (trk_pi->gDCA(pVtx) < 0.01) continue; // Pion DCA to PV > 100 um

          // --- Calculate Invariant Mass and Fill Histograms ---
          TLorentzVector kaon4V, pion4V;
          kaon4V.SetVectM(trk_K->pMom(), M_KAON_PLUS);
          pion4V.SetVectM(trk_pi->pMom(), M_PION_PLUS);
          TLorentzVector d0pair = kaon4V + pion4V;

          if (d0pair.Perp() < 1.5) continue; // D0 candidate pT > 1.5 GeV/c

          // Fill Signal or Background histograms based on charge
          if (isOppositeSign) {
              hMkpiCount->Fill(d0pair.M());
              hMkpiCountPt->Fill(d0pair.M(), d0pair.Perp());
          } else if (isLikeSign) {
              // Distinguish between -- and ++ pairs if needed
              if (trk_K->charge() < 0) { // K-pi-
                  hMkpiCount_like1->Fill(d0pair.M());
                  hMkpiCountPt_like1->Fill(d0pair.M(), d0pair.Perp());
              } else { // K+pi+
                  hMkpiCount_like2->Fill(d0pair.M());
                  hMkpiCountPt_like2->Fill(d0pair.M(), d0pair.Perp());
              }
          }
      } // end inner track loop
    } // end outer track loop
  } //Good Event
}
  if(DEBUG) cout<<"end make"<<endl;
  return kStOK;
}
bool StPicoDstarMixedMaker::isGoodTrigger(StPicoEvent const* const picoEvent) const
{
  for (auto trg : anaCuts::triggers)
  {
    if (picoEvent->isTrigger(trg)) return true;
  }

  return false;
}
bool StPicoDstarMixedMaker::isGoodTrack(StPicoTrack const* trk, float dca) const
{
  // StThreeVectorF const vtx = mPicoDstMaker->picoDst()->event()->primaryVertex();
  return trk->gPt() > anaCuts::GPt && fabs(trk->nHitsFit()) >= anaCuts::NHitsFit && 
    fabs(trk->gMom().Eta())<anaCuts::Eta &&
    fabs(trk->nHitsFit()*1.0/trk->nHitsMax())>=anaCuts::NHitsFit2Poss &&
    fabs(trk->nHitsDedx())>=anaCuts::NHitsDedx && fabs(dca)<=anaCuts::Dca;
    // fabs(trk->nHitsDedx())>=anaCuts::NHitsDedx &&
     //fabs( trk->gDCA(vtx.x() , vtx.y(), vtx.z() )) <= anaCuts::Dca;
}
bool StPicoDstarMixedMaker::isGoodQaTrack(StPicoTrack const* const trk) const
{
  // StThreeVectorF vtx = mPicoDstMaker->picoDst()->event()->primaryVertex();
  return trk->gPt() > anaCuts::qaGPt && fabs(trk->nHitsFit()) >= anaCuts::qaNHitsFit && 
    fabs(trk->gMom().Eta())<anaCuts::qaEta &&
    fabs(trk->nHitsFit()*1.0/trk->nHitsMax())>=anaCuts::NHitsFit2Poss &&
    // fabs(trk->nHitsDedx())>=anaCuts::qaNHitsDedx && fabs(trk->gDCA(vtx.x(),vtx.y(),vtx.z()))<=anaCuts::qaDca;
    fabs(trk->nHitsDedx())>=anaCuts::qaNHitsDedx ;
}
bool StPicoDstarMixedMaker::isGoodQaEvent(StPicoEvent const* const picoEvent) const
{
  TVector3 pVtx = picoEvent->primaryVertex();
  // StThreeVectorF pVtx = picoEvent->primaryVertex();
  return fabs(pVtx.z()) < anaCuts::qavz &&
     fabs(pVtx.z() - picoEvent->vzVpd()) < anaCuts::qavzVpdVz &&
    !(fabs(pVtx.x()) < anaCuts::qaVerror && fabs(pVtx.y()) < anaCuts::qaVerror &&
        fabs(pVtx.z()) < anaCuts::qaVerror) &&
    sqrt(TMath::Power(pVtx.x(), 2) + TMath::Power(pVtx.y(), 2)) <=  anaCuts::qaVr;
}
bool StPicoDstarMixedMaker::isGoodEvent(StPicoEvent const* const picoEvent) const
{
  TVector3 pVtx = picoEvent->primaryVertex();
  // StThreeVectorF pVtx = picoEvent->primaryVertex();
  return fabs(pVtx.z()) < anaCuts::vz &&
     fabs(pVtx.z() - picoEvent->vzVpd()) < anaCuts::vzVpdVz &&
    !(fabs(pVtx.x()) < anaCuts::Verror && fabs(pVtx.y()) < anaCuts::Verror &&
        fabs(pVtx.z()) < anaCuts::Verror) &&
    sqrt(TMath::Power(pVtx.x(), 2) + TMath::Power(pVtx.y(), 2)) <=  anaCuts::Vr;
}
float StPicoDstarMixedMaker::getTofBeta(StPicoTrack const* const trk) const
{
  int index2tof = trk->bTofPidTraitsIndex();
  float beta = std::numeric_limits<float>::quiet_NaN();
  if (index2tof >= 0)
  {
    StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
    if (tofPid)
    {
      beta = tofPid->btofBeta();
      if (beta < 1e-4)
      {
        // StThreeVectorF const vtx = mPicoDstMaker->picoDst()->event()->primaryVertex();
        TVector3 const vtx3 = mPicoDstMaker->picoDst()->event()->primaryVertex();
        StThreeVectorF vtx(vtx3.x(),vtx3.y(),vtx3.z());
        TVector3 const btofHitPos3 = tofPid->btofHitPos();
        // StThreeVectorF const btofHitPos = tofPid->btofHitPos();
        StThreeVectorF btofHitPos(btofHitPos3.x(),btofHitPos3.y(),btofHitPos3.z());
        StPicoPhysicalHelix helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());
        // StPhysicalHelixD helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());
        float L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
      }
    }
  } 
  return beta;
}

bool StPicoDstarMixedMaker::isPion(StPicoTrack const* const trk) const
{
    double p = trk->pMom().Mag();
    double beta = getTofBeta(trk);
    bool tofmatch = (beta != std::numeric_limits<float>::quiet_NaN()) && beta > 0;

    bool isTpcPion = fabs(trk->nSigmaPion()) < 2.0;

    bool isTofPion = false;
    if (tofmatch) {
        float beta_expected = sqrt(p*p + M_PION_PLUS*M_PION_PLUS) / p;
        if (fabs(1.0/beta - beta_expected) < 0.03) {
            isTofPion = true;
        }
    }
    return isTpcPion && isTofPion;
}

bool StPicoDstarMixedMaker::isKaon(StPicoTrack const* const trk) const
{
    double p = trk->pMom().Mag();
    double beta = getTofBeta(trk);
    bool tofmatch = (beta != std::numeric_limits<float>::quiet_NaN()) && beta > 0;

    bool isTpcKaon = fabs(trk->nSigmaKaon()) < 2.0;
    
    bool isTofKaon = false;
    if (tofmatch) {
        float beta_expected = sqrt(p*p + M_KAON_PLUS*M_KAON_PLUS) / p;
        if (fabs(1.0/beta - beta_expected) < 0.03) {
            isTofKaon = true;
        }
    }
    return isTpcKaon && isTofKaon;
}