//____________________________________________________________________________..
//
//____________________________________________________________________________..

#include "EnergyCorrection.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h> // for GenVertex, GenVertex::part...
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h> // for GenParticle
#include <HepMC/GenRanges.h>
#include <HepMC/HeavyIon.h> // for HeavyIon

#include <TFile.h>
#include <TH1.h>

//____________________________________________________________________________..
EnergyCorrection::EnergyCorrection(const std::string &name) : SubsysReco(name) {
  std::cout << "EnergyCorrection::EnergyCorrection(const std::string &name) "
               "Calling ctor"
            << std::endl;
}

//____________________________________________________________________________..
EnergyCorrection::~EnergyCorrection() {
  std::cout << "EnergyCorrection::~EnergyCorrection() Calling dtor"
            << std::endl;
}

//____________________________________________________________________________..
int EnergyCorrection::Init(PHCompositeNode *topNode) {
  std::cout << "EnergyCorrection::Init(PHCompositeNode *topNode) Initializing"
            << std::endl;
  // read correction histogram from file
  std::string filename = "/sphenix/user/shuhangli/dETdeta/macro/fitratio" + m_generatortype + ".root";
  TFile *f_upweight =
      new TFile(filename.c_str());
  std::string postfix[5] = {"0010_ratio", "1020_ratio", "2040_ratio",
                            "4060_ratio", "6092_ratio"};
  for (int i = 0; i < 5; i++) {
    h_pimi[i] = (TH1F *)f_upweight->Get(Form("h_pimi%s", postfix[i].c_str()));
    h_pip[i] = (TH1F *)f_upweight->Get(Form("h_pip%s", postfix[i].c_str()));
    h_p[i] = (TH1F *)f_upweight->Get(Form("h_p%s", postfix[i].c_str()));
    h_pbar[i] = (TH1F *)f_upweight->Get(Form("h_pbar%s", postfix[i].c_str()));
    h_kp[i] = (TH1F *)f_upweight->Get(Form("h_kp%s", postfix[i].c_str()));
    h_kmi[i] = (TH1F *)f_upweight->Get(Form("h_kmi%s", postfix[i].c_str()));
  }

  //set centralities average
  if(m_generatortype == "HIJING") {
    avgcent[0] = 329.815;
    avgcent[1] = 238.602;
    avgcent[2] = 144.272;
    avgcent[3] = 64.9728;
    avgcent[4] = 17.8337;
  }
  else if(m_generatortype == "AMPT") {
    avgcent[0] = 340.095;
    avgcent[1] = 254.515;
    avgcent[2] = 160.105;
    avgcent[3] = 76.1558;
    avgcent[4] = 22.783;
  }
  else if (m_generatortype == "EPOS") {
    avgcent[0] = 325.8;
    avgcent[1] = 236.1;
    avgcent[2] = 141.5;
    avgcent[3] = 61.6;
    avgcent[4] = 14.7;
  }
  else {
    std::cout << "EnergyCorrection::Init(PHCompositeNode *topNode) "
                 "generator type not supported"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EnergyCorrection::process_event(PHCompositeNode *topNode) {
  if (Verbosity() > 0)
    std::cout << "EnergyCorrection::process_event(PHCompositeNode *topNode) "
                 "Processing Event"
              << std::endl;

  PHHepMCGenEventMap *genevtmap =
      findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!genevtmap) {
    std::cout << "no genevtmap" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  for (PHHepMCGenEventMap::Iter iter = genevtmap->begin();
       iter != genevtmap->end(); ++iter) {
    PHHepMCGenEvent *genevt = iter->second;
    // check if embedded
    if (genevt->get_embedding_id() != 0)
      continue;
    HepMC::GenEvent *event = genevt->getEvent();
    if (!event) {
      std::cout << PHWHERE << " no evt pointer under HEPMC Node found"
                << std::endl;
    } else {
      HepMC::HeavyIon *hi = event->heavy_ion();
      if (!hi) {
        std::cout << PHWHERE
                  << ": Fermi Motion Afterburner needs the Heavy Ion Event "
                     "Info, GenEvent::heavy_ion() returns NULL"
                  << std::endl;
        exit(1);
      }
      m_npart = hi->Npart_proj() + hi->Npart_targ();
    }
  }
  if (m_npart < 0) {
    std::cout << "cant find npart" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get truthinfo
  PHG4TruthInfoContainer *truthinfo =
      findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  // get primary particles
  if (!truthinfo) {
    std::cout << "EnergyCorrection::process_event(PHCompositeNode *topNode) "
                 "Could not locate G4TruthInfo node"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (m_upweighttruth) {
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();

    for (PHG4TruthInfoContainer::Iterator iter = range.first;
         iter != range.second; ++iter) {
      PHG4Particle *particle = iter->second;
      int pid = particle->get_pid();
      float pt = sqrt(particle->get_px() * particle->get_px() +
                      particle->get_py() * particle->get_py());
      float scale = findcorrection(m_npart, pid, pt);
      particle->set_e(particle->get_e() * scale);
      particle->set_px(particle->get_px() * scale);
      particle->set_py(particle->get_py() * scale);
      particle->set_pz(particle->get_pz() * scale);
    }
  }

  // get hits
  PHG4HitContainer *hits =
      findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  if (!hits) {
    std::cout << "EnergyCorrection::process_event(PHCompositeNode *topNode) "
                 "Could not locate g4 hit node "
              << m_HitNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHG4HitContainer::ConstRange hit_range = hits->getHits();
  for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first;
       hit_iter != hit_range.second; hit_iter++) {
    PHG4Hit *hit = hit_iter->second;
    int showerid = hit_iter->second->get_shower_id();
    PHG4Shower *shower = truthinfo->GetPrimaryShower(showerid);
    if (!shower) {
      if (Verbosity() > 0)
        std::cout << "EnergyCorrection::process_event(PHCompositeNode "
                     "*topNode) No shower found showerid: "
                  << showerid << std::endl;
      continue;
    }
    int trkid = shower->get_parent_particle_id();
    PHG4Particle *part = truthinfo->GetParticle(trkid);
    if (!part) {
      std::cout << "EnergyCorrection::process_event(PHCompositeNode *topNode) "
                   "No parent particle found, track id: "
                << trkid << std::endl;
    }
    // get particle pid and pt
    int pid = part->get_pid();
    float pt =
        sqrt(part->get_px() * part->get_px() + part->get_py() * part->get_py());
    // find correction factor for G4Hits
    float scale = findcorrection(m_npart, pid, pt);

    // apply correction
    hit->set_edep(hit->get_edep() * scale);
    hit->set_light_yield(hit->get_light_yield() * scale);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EnergyCorrection::End(PHCompositeNode *topNode) {
  std::cout
      << "EnergyCorrection::End(PHCompositeNode *topNode) This is the End..."
      << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}