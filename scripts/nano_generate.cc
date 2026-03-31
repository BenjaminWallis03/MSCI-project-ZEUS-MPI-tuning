// Custom Pythia8 + Rivet event generation wrapper
// Used for ZEUS photoproduction MPI tuning
// Reads run.cmnd files and produces Rivet YODA output
//
// Based on standard nano.generate structure with HepMC3 support

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/FourVector.h"

#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Math/Vector4.hh"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// ---------- helpers ----------
static inline std::string trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  const auto e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}
static inline bool starts_with(const std::string& s, const std::string& pfx) {
  return s.rfind(pfx, 0) == 0;
}
static inline std::string to_lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  return s;
}
static inline bool parse_bool(const std::string& rhs, bool def=false) {
  const std::string v = to_lower(trim(rhs));
  if (v == "on" || v == "true" || v == "1" || v == "yes") return true;
  if (v == "off" || v == "false" || v == "0" || v == "no") return false;
  return def;
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " run.cmnd\n";
    return 1;
  }
  const std::string cmndPath = argv[1];

  // ---- Defaults (overridden by Main:* in run.cmnd) ----
  long long numberOfEvents   = 1000000;
  long long timesAllowErrors = 100000;
  bool runRivet              = true;
  std::string analysisName   = "ZEUS_2007_I753991";
  bool outputLog             = true;

  std::ofstream log;
  auto logprint = [&](const std::string& msg) {
    if (outputLog) {
      if (!log.is_open()) log.open("nano.log");
      log << msg << "\n";
    }
    std::cout << msg << "\n";
  };

  // ---- Parse run.cmnd into Pythia settings + Main:* wrapper settings ----
  std::vector<std::string> pythiaLines;
  {
    std::ifstream in(cmndPath);
    if (!in) {
      std::cerr << "Could not open " << cmndPath << "\n";
      return 1;
    }
    std::string line;
    while (std::getline(in, line)) {
      line = trim(line);
      if (line.empty()) continue;
      if (starts_with(line, "!")) continue;

      // strip inline comments after '!'
      const auto excl = line.find('!');
      if (excl != std::string::npos) line = trim(line.substr(0, excl));
      if (line.empty()) continue;

      const auto eq = line.find('=');
      if (eq == std::string::npos) continue;

      const std::string key = trim(line.substr(0, eq));
      const std::string rhs = trim(line.substr(eq + 1));

      if (starts_with(key, "Main:")) {
        if (key == "Main:numberOfEvents") numberOfEvents = std::stoll(rhs);
        else if (key == "Main:timesAllowErrors") timesAllowErrors = std::stoll(rhs);
        else if (key == "Main:runRivet") runRivet = parse_bool(rhs, true);
        else if (key == "Main:analyses") analysisName = rhs;
        else if (key == "Main:outputLog") outputLog = parse_bool(rhs, true);
        continue;
      }

      pythiaLines.push_back(key + " = " + rhs);
    }
  }

  // ---- Configure + init Pythia ----
  Pythia8::Pythia pythia;
  for (const auto& s : pythiaLines) {
    if (!pythia.readString(s)) {
      logprint("Pythia rejected setting: " + s);
    }
  }

  if (!pythia.init()) {
    logprint("ERROR: Pythia failed to initialise.");
    return 1;
  }

  // ---- Beam IDs & energies from Pythia settings (i.e. from run.cmnd) ----
  const int idA = pythia.settings.mode("Beams:idA");
  const int idB = pythia.settings.mode("Beams:idB");
  const double eA = pythia.settings.parm("Beams:eA");
  const double eB = pythia.settings.parm("Beams:eB");

  // Build approximate on-axis 4-momenta in the lab/collider frame.
  // For ultra-relativistic beams, pz ~ +/-E is fine for Rivet beam bookkeeping.
  const double pxA = 0.0, pyA = 0.0, pzA = +eA, EA = eA;
  const double pxB = 0.0, pyB = 0.0, pzB = -eB, EB = eB;

  // ---- Configure Rivet ----
  Rivet::AnalysisHandler rivet;
  bool rivetInitialised = false;

  if (runRivet) {
    // Split comma-separated analyses and add each one
  {
    std::stringstream ss(analysisName);
    std::string token;
    while (std::getline(ss, token, ',')) {
      while (!token.empty() && token.front() == ' ') token.erase(token.begin());
      while (!token.empty() && token.back() == ' ') token.pop_back();
      if (!token.empty()) rivet.addAnalysis(token);
    }
  }

  // Set the Rivet "run beams"
  const Rivet::Particle beamA(idA, Rivet::FourMomentum(pxA, pyA, pzA, EA));
  const Rivet::Particle beamB(idB, Rivet::FourMomentum(pxB, pyB, pzB, EB));
  rivet.setRunBeams(std::make_pair(beamA, beamB));

    logprint("Rivet enabled with analysis: " + analysisName);
  } else {
    logprint("Rivet disabled (Main:runRivet=off).");
  }

  // ---- HepMC3 converter ----
  HepMC3::Pythia8ToHepMC3 toHepMC;
  toHepMC.set_print_inconsistency(false);
  toHepMC.set_free_parton_warnings(false);

  long long nAccepted = 0;
  long long nAborted  = 0;

  // ---- Event loop ----
  for (long long i = 0; i < numberOfEvents; ++i) {
    if (!pythia.next()) {
      ++nAborted;
      if (nAborted > timesAllowErrors) {
        logprint("ERROR: too many aborted events (" + std::to_string(nAborted) +
                 " > " + std::to_string(timesAllowErrors) + "). Stopping.");
        break;
      }
      continue;
    }
    ++nAccepted;

    if (!runRivet) continue;

    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
    toHepMC.fill_next_event(pythia, &evt);

  // ---- IMPORTANT FIX: create a vertex and attach beams as incoming particles ----
  auto beamA_evt = std::make_shared<HepMC3::GenParticle>(
    HepMC3::FourVector(pxA, pyA, pzA, EA), idA, 4
);
  auto beamB_evt = std::make_shared<HepMC3::GenParticle>(
    HepMC3::FourVector(pxB, pyB, pzB, EB), idB, 4
);

  auto vtx = std::make_shared<HepMC3::GenVertex>();
  vtx->add_particle_in(beamA_evt);
  vtx->add_particle_in(beamB_evt);

  // Add the vertex (and hence beams) to the event graph
  evt.add_vertex(vtx);

  // Also explicitly mark them as beams
  evt.set_beam_particles(beamA_evt, beamB_evt);

    // Init Rivet on first event (your API requires this)
    if (!rivetInitialised) {
      rivet.init(evt);
      rivetInitialised = true;
    }

    rivet.analyze(evt);
  }

  if (runRivet) {
    if (!rivetInitialised) {
      logprint("ERROR: Rivet never initialised (no accepted events). Check your process settings.");
      return 1;
    }
    rivet.finalize();
    rivet.writeData("rivet.yoda");
    logprint("Wrote: rivet.yoda");
  }

  logprint("Requested events: " + std::to_string(numberOfEvents));
  logprint("Accepted events : " + std::to_string(nAccepted));
  logprint("Aborted events  : " + std::to_string(nAborted));

  pythia.stat();
  return 0;
}

