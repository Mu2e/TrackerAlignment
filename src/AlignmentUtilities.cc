// Ryunosuke O'Neil, 2020
// roneil@fnal.gov
// ryunoneil@gmail.com

#include "TrackerAlignment/inc/AlignmentUtilities.hh"
#include "TrackerAlignment/inc/MilleDataWriter.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "Minuit2/MnUserCovariance.h"

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <vector>
using namespace mu2e;

namespace AlignmentUtilities {
  using xyzVec = CLHEP::Hep3Vector; // switch to XYZVec TODO

int hitAmbiguity(CosmicTimeTrack const& track, Hep3Vector const& straw_mp,
                 Hep3Vector const& straw_dir) {
  Hep3Vector sep = track.intercept() - straw_mp;
  Hep3Vector perp = (track.direction().cross(straw_dir)).unit();
  double dperp = perp.dot(sep);

  return (dperp > 0 ? -1 : 1);
}

// straw alignment function
// carbon copy of AlignedTrackerMaker::fromDb
std::pair<Hep3Vector, Hep3Vector> alignStraw(Tracker const& tracker,
                                             StrawId const& strawId,
                                             TrkAlignParams const& align_tracker,
                                             TrkAlignParams const& align_plane,
                                             TrkAlignParams const& align_panel,
                                             TrkStrawEndAlign const& align_straw) {

  auto const& plane = tracker.getPlane(strawId);
  auto const& panel = tracker.getPanel(strawId);
  auto const& straw = tracker.getStraw(strawId);

  // nominal place of plane in the tracker
  auto plane_to_ds = plane.planeToDS();
  // nominal panel in the plane
  auto panel_to_plane = plane.dsToPlane()*panel.panelToDS();
  // chained alignment including panel
  auto aligned_panel_to_ds = align_tracker.transform() * (plane_to_ds * align_plane.transform()) * (panel_to_plane * align_panel.transform());
  // Move wire back to the nominal panel frame (UVW)
  auto ds_to_panel = panel.dsToPanel();
  std::array<xyzVec,2> wireends;
  for(int iend=0;iend < StrawEnd::nends; iend++){
    auto end = static_cast<StrawEnd::End>(iend);
    // include the straw end alignment
    auto wireend_UVW = ds_to_panel*straw.wireEnd(end) + align_straw.wireDeltaUVW(end);
    wireends[iend] = aligned_panel_to_ds*wireend_UVW;
  }
// compute the net position and direction from these
  auto aligned_strawmid = 0.5*(wireends[0] + wireends[1]);
  auto aligned_strawdir = (wireends[StrawEnd::cal]- wireends[StrawEnd::hv]).unit(); // convention  is direction points from HV to Cal
  return {aligned_strawmid, aligned_strawdir};
}



/* Numerical partial DOCA/TOCA derivatives
 *
 */


namespace {

// do the alignment parameters really need to be double?
  double tocaGlobalDep(CosmicTimeTrack const& track, StrawId const& strawId,
                                std::vector<double> const& globals, Tracker const& nominalTracker, StrawResponse const& strawRes) {


    TrkAlignParams align_tracker{strawId, StrawIdMask::tracker, 0, 0, 0, 0, 0, 0};
    TrkAlignParams align_plane{strawId, StrawIdMask::plane, globals[0], globals[1], globals[2], globals[3], globals[4], globals[5]};
    TrkAlignParams align_panel{strawId, StrawIdMask::uniquepanel, globals[6], globals[7],  globals[8], globals[9], globals[10], globals[11]};
    TrkStrawEndAlign align_straw{strawId.uniqueStraw(), strawId, 0,0,0,0,0,0,0,0}; // not sure how to really initialize this FIXME!

    // returns pair of vectors { straw_pos, straw_dir }
    Hep3Vector straw_pos, straw_dir;
    std::tie(straw_pos, straw_dir) = alignStraw(nominalTracker,
                                                strawId, align_tracker, align_plane, align_panel, align_straw);

    TwoLinePCA pca(track.intercept(), track.direction(), straw_pos, straw_dir);

    double traj_time = (pca.point1() - track.intercept()).dot(track.direction()) / 299.9;
    double d2t_doca = strawRes.driftDistanceToTime(strawId, pca.dca(), 0);
    double t_offset = strawRes.driftTimeOffset(strawId, pca.dca(), 0);

    double predictedTime = traj_time + t_offset + track.params[CosmicTimeTrack::t0] + d2t_doca;

    return predictedTime;
  }

  double docaGlobalDep(CosmicTimeTrack const& track, StrawId const& strawId,
                                std::vector<double> const& globals, Tracker const& nominalTracker, StrawResponse const& strawRes) {
// this function should be consolidated with the above FIXME!
    TrkAlignParams align_tracker{strawId, StrawIdMask::tracker, 0, 0, 0, 0, 0, 0};
    TrkAlignParams align_plane{strawId, StrawIdMask::plane, globals[0], globals[1], globals[2], globals[3], globals[4], globals[5]};
    TrkAlignParams align_panel{strawId, StrawIdMask::uniquepanel, globals[6], globals[7],  globals[8], globals[9], globals[10], globals[11]};
    TrkStrawEndAlign align_straw{strawId.uniqueStraw(), strawId, 0,0,0,0,0,0,0,0}; // not sure how to really initialize this FIXME!
    // returns pair of vectors { straw_pos, straw_dir }
    Hep3Vector straw_pos, straw_dir;
    std::tie(straw_pos, straw_dir) = alignStraw(nominalTracker,
                                                strawId, align_tracker, align_plane, align_panel,align_straw);

    TwoLinePCA pca(track.intercept(), track.direction(), straw_pos, straw_dir);

    int ambig = hitAmbiguity(track, straw_pos, straw_dir);
    double doca = ambig * pca.dca();

    return doca;
  }


  // not meant to be called from outside of this namespace
  double _numericalDerivative(StrawId const& straw, CosmicTimeTrack& track,
                            std::vector<double>& globals, Tracker const& nominalTracker,
                            StrawResponse const& strawRes, bool isGlobalParam,
                            size_t const& paramIdx, double step_size, bool useTimeDomain) {
    // calculate numerical partial derivative wrt param at paramIdx in either
    // local, or global param array

    double x;

    if (isGlobalParam) {
      x = globals[paramIdx];
      globals[paramIdx] = x + step_size;
    } else {
      x = track.params[paramIdx];
      track.params[paramIdx] = x + step_size;
    }

    double pdiff;

    if (useTimeDomain) {
      pdiff = tocaGlobalDep(track, straw, globals, nominalTracker, strawRes);
    }
    else {
      pdiff = docaGlobalDep(track, straw, globals, nominalTracker, strawRes);
    }

    if (isGlobalParam) {
      globals[paramIdx] = x - step_size;
    } else {
      track.params[paramIdx] = x - step_size;
    }

    if (useTimeDomain) {
      pdiff -= tocaGlobalDep(track, straw, globals, nominalTracker, strawRes);
    }
    else {
      pdiff -= docaGlobalDep(track, straw, globals, nominalTracker, strawRes);
    }

    pdiff /= (2.0 * step_size);

    if (isGlobalParam) {
      globals[paramIdx] = x;
    } else {
      track.params[paramIdx] = x;
    }

    Straw const& nominalStraw = nominalTracker.getStraw(straw);
    Hep3Vector const& nominalStraw_mp = nominalStraw.wirePosition();
    Hep3Vector const& nominalStraw_dir = nominalStraw.wireDirection();

    return pdiff;
  }
}

std::pair<std::vector<double>, std::vector<double>>
numericalDerivatives(CosmicTimeTrack const& _track, StrawId const& straw,
                         TrkAlignParams const& alignPlane,
                         TrkAlignParams const& alignPanel,
                         Tracker const& nominalTracker,
                         StrawResponse const& strawRes,
                         bool useTimeDomain) {

  std::vector<double> result_locals;
  std::vector<double> result_globals;

  // copy track into this scope - we will modify the params in place a lot
  CosmicTimeTrack track(_track);

  // same here also
  std::vector<double> globals{alignPlane.dx(), alignPlane.dy(), alignPlane.dz(),
                              alignPlane.rx(), alignPlane.ry(), alignPlane.rz(),

                              alignPanel.dx(), alignPanel.dy(), alignPanel.dz(),
                              alignPanel.rx(), alignPanel.ry(), alignPanel.rz()};


  // locals first
  for (size_t paramIdx = 0; paramIdx < track.npars(); ++paramIdx) {
    result_locals.emplace_back(
        _numericalDerivative(straw, track, globals, nominalTracker, strawRes,
                              false, paramIdx, 1e-5, useTimeDomain));
  }

  for (size_t paramIdx = 0; paramIdx < globals.size(); ++paramIdx) {
    result_globals.emplace_back(
        _numericalDerivative(straw, track, globals, nominalTracker, strawRes,
                              true, paramIdx, 1e-5, useTimeDomain));
  }


  return {result_locals, result_globals};
}

std::pair<std::vector<double>, std::vector<double>>
analyticDerivatives(StrawId const& straw_id,
    TwoLinePCA const& pca,
    Tracker const& nominalTracker,
    StrawResponse const& strawRes,
    bool useTimeDomain) {

  std::vector<double> derivativesLocal;
  std::vector<double> derivativesGlobal;

  auto dca_dir = (pca.point2()-pca.point1()).unit();

  //FIXME
  double drift_time_p = strawRes.driftDistanceToTime(straw_id, pca.dca()+1e-5, 0) +
    strawRes.driftTimeOffset(straw_id, pca.dca()+1e-5, 0);
  double drift_time_n = strawRes.driftDistanceToTime(straw_id, pca.dca()-1e-5, 0) +
    strawRes.driftTimeOffset(straw_id, pca.dca()-1e-5, 0);
  auto drdx = (drift_time_p-drift_time_n)/(2*1e-5) * dca_dir;

  double drdt = -1;
  double drda0 = -drdx.x();
  double drdb0 = -drdx.z();
  double drda1 = pca.point2().y()*drdx.x();
  double drdb1 = pca.point2().y()*drdx.z();

  derivativesLocal.push_back(drda0);
  derivativesLocal.push_back(drdb0);
  derivativesLocal.push_back(drda1);
  derivativesLocal.push_back(drdb1);
  derivativesLocal.push_back(drdt);

  const Panel& panel = nominalTracker.getPanel(straw_id);
  const Plane& plane = nominalTracker.getPlane(straw_id);
  std::vector<double> panelDerivs(6,0);
  std::vector<double> planeDerivs(6,0);
  // calculate derivative wrt alignment directions
  panelDerivs[0] = panel.uDirection().dot(drdx);
  panelDerivs[1] = panel.vDirection().dot(drdx);
  panelDerivs[2] = panel.wDirection().dot(drdx);

  auto panel_delta = pca.point2() - panel.origin();
  panelDerivs[3] = panel.uDirection().cross(panel_delta).dot(drdx);
  panelDerivs[4] = panel.vDirection().cross(panel_delta).dot(drdx);
  panelDerivs[5] = panel.wDirection().cross(panel_delta).dot(drdx);

  CLHEP::Hep3Vector plane_uDirection(1,0,0);
  CLHEP::Hep3Vector plane_vDirection(0,1,0);
  CLHEP::Hep3Vector plane_wDirection(0,0,1);
  if (plane.planeToDS().rotation().getTheta() != 0){
    plane_uDirection = CLHEP::Hep3Vector(-1,0,0);
    plane_wDirection = CLHEP::Hep3Vector(0,0,-1);
  }

  planeDerivs[0] = plane_uDirection.dot(drdx); //plane.uDirection().dot(drdx);
  planeDerivs[1] = plane_vDirection.dot(drdx); //plane.vDirection().dot(drdx);
  planeDerivs[2] = plane_wDirection.dot(drdx); //plane.wDirection().dot(drdx);

  auto plane_delta = pca.point2() - plane.origin();
  planeDerivs[3] = plane_uDirection.cross(plane_delta).dot(drdx); //plane.uDirection().cross(plane_delta).dot(drdx);
  planeDerivs[4] = plane_vDirection.cross(plane_delta).dot(drdx); //plane.vDirection().cross(plane_delta).dot(drdx);
  planeDerivs[5] = plane_wDirection.cross(plane_delta).dot(drdx); //plane.wDirection().cross(plane_delta).dot(drdx);

  for (size_t i=0;i<6;i++)
    derivativesGlobal.push_back(planeDerivs[i]);
  for (size_t i=0;i<6;i++)
    derivativesGlobal.push_back(panelDerivs[i]);

  return {derivativesLocal, derivativesGlobal};
}
} // namespace AlignmentUtilities
