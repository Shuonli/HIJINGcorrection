// these include guards are not really needed, but if we ever include this
// file somewhere they would be missed and we will have to refurbish all macros
#ifndef MACRO_FUN4ALLISOLATEDTRACKS_C
#define MACRO_FUN4ALLISOLATEDTRACKS_C

#include <G4Setup_sPHENIX.C>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <g4centrality/PHG4CentralityReco.h>

#include </sphenix/user/shuhangli/HIJINGcorrection/HIJINGcorrection/source/EnergyCorrection.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libEnergyCorrection.so)
// example code for running the upweighting afterburner, modified from Emma's isotrack analysis macro(thanks ;) )
void Fun4All_IsolatedTrackAnalysis(
    const string &trackFile = "dst_tracks.list",
    const string &clusterFile = "dst_calo_g4hit.list",
    const string &truthFile = "g4hits.list",
    const string &bbcFile = "dst_bbc_g4hit.list",
    const string &trackClusters = "dst_trkr_cluster.list",
    const string &trackG4hits = "dst_trkr_g4hit.list",
    const string &outputFile = "isotrack_output.root") {
  // this convenience library knows all our i/o objects so you don't
  // have to figure out what is in each dst type
  gSystem->Load("libg4dst.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1); // set it to 1 if you want event printouts
  /*
  Fun4AllInputManager *inTrack = new Fun4AllDstInputManager("DSTTracks");
  std::cout << "Adding file list " << trackFile << std::endl;
  inTrack->AddListFile(trackFile,1);
  se->registerInputManager(inTrack);
  */
  /*
  Fun4AllInputManager *inCluster = new Fun4AllDstInputManager("DSTClusters");
  std::cout << "Adding file list " << clusterFile << std::endl;
  inCluster->AddListFile(clusterFile,1);
  se->registerInputManager(inCluster);
  */
  //-----------------------------
  //you only need a list that have all G4hits and primary truth to run the correction
  //is is safe to just have the g4hit list
  //-----------------------------
  Fun4AllInputManager *truthTracks = new Fun4AllDstInputManager("DSTTruth");
  std::cout << "Adding file list " << truthFile << std::endl;
  truthTracks->AddListFile(truthFile, 1);
  se->registerInputManager(truthTracks);
  /*
  Fun4AllInputManager *inCentrality = new Fun4AllDstInputManager("DSTCentral");
  std::cout << "Adding file list " << bbcFile << std::endl;
  inCentrality->AddListFile(bbcFile,1);
  se->registerInputManager(inCentrality);
  */
  /*
  Fun4AllInputManager *inTrackCluster = new
  Fun4AllDstInputManager("DSTTrackClusters"); std::cout << "Adding file list "
  << trackClusters << std::endl; inTrackCluster->AddListFile(trackClusters,1);
  se->registerInputManager(inTrackCluster);

  Fun4AllInputManager *inTrackG4hits = new
  Fun4AllDstInputManager("DSTTrackG4hits"); std::cout << "Adding file list " <<
  trackG4hits << std::endl; inTrackG4hits->AddListFile(trackG4hits,1);
  se->registerInputManager(inTrackG4hits);
  */
  EnergyCorrection *energycorrect = new EnergyCorrection();
  energycorrect->SetHitNodeName("G4HIT_CEMC");
  //only want to upweight the truth once
  energycorrect->SetUpWeightTruth(true);
  se->registerSubsystem(energycorrect);
  energycorrect = new EnergyCorrection();
  energycorrect->SetHitNodeName("G4HIT_HCALIN");
  se->registerSubsystem(energycorrect);
  energycorrect->SetUpWeightTruth(false);
  energycorrect = new EnergyCorrection();
  energycorrect->SetHitNodeName("G4HIT_HCALOUT");
  energycorrect->SetUpWeightTruth(false);
  se->registerSubsystem(energycorrect);
  /*
    PHG4CylinderCellReco *cemc_cells =
        new PHG4CylinderCellReco("CEMCCYLCELLRECO");
    cemc_cells->Detector("CEMC");
    for (int i = G4CEMC::Min_cemc_layer; i <= G4CEMC::Max_cemc_layer; i++) {
      const double radius = 95;
      cemc_cells->cellsize(i, 2 * M_PI / 256. * radius, 2 * M_PI / 256. *
    radius);
    }
    cemc_cells->Verbosity(1);
    se->registerSubsystem(cemc_cells);
  */

  PHG4FullProjSpacalCellReco *cemc_cells =
      new PHG4FullProjSpacalCellReco("CEMCCYLCELLRECO");
  cemc_cells->Detector("CEMC");
  cemc_cells->Verbosity(0);
  cemc_cells->get_light_collection_model().load_data_file(
      string(getenv("CALIBRATIONROOT")) +
          string("/CEMC/LightCollection/Prototype3Module.xml"),
      "data_grid_light_guide_efficiency", "data_grid_fiber_trans");
  se->registerSubsystem(cemc_cells);

  PHG4HcalCellReco *ihc = new PHG4HcalCellReco("HCALIN_CELLRECO");
  ihc->Detector("HCALIN");
  ihc->Verbosity(0);
  se->registerSubsystem(ihc);

  PHG4HcalCellReco *ohc = new PHG4HcalCellReco("HCALOUT_CELLRECO");
  ohc->Detector("HCALOUT");
  ohc->Verbosity(1);
  se->registerSubsystem(ohc);

  RawTowerBuilder *ETowerBuilder = new RawTowerBuilder("EmcRawTowerBuilder");
  ETowerBuilder->Detector("CEMC");
  ETowerBuilder->set_sim_tower_node_prefix("SIM");
  ETowerBuilder->Verbosity(0);
  se->registerSubsystem(ETowerBuilder);

  HcalRawTowerBuilder *ITowerBuilder =
      new HcalRawTowerBuilder("HcalInRawTowerBuilder");
  ITowerBuilder->Detector("HCALIN");
  ITowerBuilder->set_sim_tower_node_prefix("SIM");
  if (!isfinite(G4HCALIN::phistart)) {
    G4HCALIN::phistart =
        0.0445549893; // offet in phi (from zero) extracted from geantinos
  }
  ITowerBuilder->set_double_param("phistart", G4HCALIN::phistart);
  if (isfinite(G4HCALIN::tower_emin)) {
    ITowerBuilder->set_double_param("emin", G4HCALIN::tower_emin);
  }
  if (G4HCALIN::tower_energy_source >= 0) {
    ITowerBuilder->set_int_param("tower_energy_source",
                                 G4HCALIN::tower_energy_source);
  }
  ITowerBuilder->Verbosity(0);
  se->registerSubsystem(ITowerBuilder);

  HcalRawTowerBuilder *OTowerBuilder =
      new HcalRawTowerBuilder("HcalOutRawTowerBuilder");
  OTowerBuilder->Detector("HCALOUT");
  OTowerBuilder->set_sim_tower_node_prefix("SIM");
  if (!isfinite(G4HCALOUT::phistart)) {
    G4HCALOUT::phistart =
        0.0240615415; // offet in phi (from zero) extracted from geantinos
  }
  OTowerBuilder->set_double_param("phistart", G4HCALOUT::phistart);
  if (isfinite(G4HCALOUT::tower_emin)) {
    OTowerBuilder->set_double_param("emin", G4HCALOUT::tower_emin);
  }
  if (G4HCALOUT::tower_energy_source >= 0) {
    OTowerBuilder->set_int_param("tower_energy_source",
                                 G4HCALOUT::tower_energy_source);
  }
  OTowerBuilder->Verbosity(0);
  se->registerSubsystem(OTowerBuilder);

  /*
  PHG4CentralityReco *cent = new PHG4CentralityReco();
  cent->Verbosity(0);
  cent->GetCalibrationParameters().ReadFromFile("centrality", "xml", 0, 0,
  string(getenv("CALIBRATIONROOT")) + string("/Centrality/"));
  se->registerSubsystem(cent);
 */
  //------------------------------analysis
  // here------------------------------------------ 

  se->run(300);
  se->PrintTimer();
  se->End();

  delete se;
  cout << "Analysis Completed" << endl;

  gSystem->Exit(0);
}

#endif // MACRO_FUN4ALLG4SLOPECAL_C
