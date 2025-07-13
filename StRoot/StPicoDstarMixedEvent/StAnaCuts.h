#ifndef StAnaCuts_H
#define StAnaCuts_H

/* **************************************************
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include "Rtypes.h"
#include <string>
#include <array>

namespace anaCuts
{
   /*std::array<unsigned int, 3> const triggers = {
   710000,710010,710020
   };*/ // 11p5 2020 MB
   /*std::array<unsigned int, 1> const triggers = {
   780020
   };*///9p2 2020 MB
  /* std::array<unsigned int, 1> const triggers = {
   800010
   };*///9p2 2020 MB
   /*std::array<unsigned int, 13> const triggers = {
   600003,600001,600011,600021,600031,600009,600019,600029,600002,600012,600022,600032,600042
   };//9p2 2020 MB*/
   /*std::array<unsigned int, 10> const triggers = {
   450008,450018,450009,450050,450060,450005,450015,450025,450014,450024
   };//AuAu200 2014 MB*/
   /*std::array<unsigned int, 7> const triggers = {
   810010,810004,819002,810003,810002,810009,810008
   };//7p7 2021 MB */
   /*std::array<unsigned int, 4> const triggers = {
   810010,810020,810030,810040
   };*///7p7 2021 MB
   //std::array<unsigned int, 4> const triggers = {900001, 900011, 900021, 900031}; //200GeV AuAu 2023 ZDC-MB
   std::array<unsigned int, 3> const triggers = {640041,640029,640032,}; //19GeV AuAu 2019

   //cut before QA
   float const qavz = 100.0;// < cm.
   float const qaVerror = 1.0e-5; //
   float const qaVr = 2.0; //cm
   //float const qavzVpdVz = 3; //cm
   float const qavzVpdVz = 6; //cm

   // QA tracks cuts
   float const qaGPt = 0.2;
   int const qaNHitsFit = 15;
   int const qaNHitsDedx = 5;
   float const qaDca = 3;// < cm
   float const qaEta = 1.5;
   float const qaTofPion=4;
   float const qaTpcPion=4;

   //cut
   float const vz = 100.0;// < cm.
   float const Verror = 1.0e-5; //
   float const Vr = 2.0; //cm
   //float const vzVpdVz = 3; //cm
   float const vzVpdVz = 6; //cm

   // QA tracks cuts
   float const GPt = 0.2;
   int const NHitsFit = 15;
   int const NHitsDedx = 5;
   float const NHitsFit2Poss = 0.51;
   float const Dca = 3;// < cm
   float const Eta = 1.5;
   float const TofPion=4;
   float const TpcPion=4;

   namespace pid {
      float const nSigmaPion = 2.0;
      float const nSigmaKaon = 2.0;
      float const tofBetaPion = 0.03;
      float const tofBetaKaon = 0.03;
   }

   namespace d0 {
      float const daughterPtMin = 0.4; // GeV/c
      float const dcaDaughtersMax = 0.01; // cm
      float const decayLengthMin = 0.02; // cm
      float const cosPointingAngleMin = 0.98;
      float const daughterDCAPVMin = 0.01; // cm
      float const d0PtMin = 1.5; // GeV/c
   }
}
#endif