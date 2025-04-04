// Ryunosuke O'Neil, 2019
// roneil@fnal.gov
// ryunoneil@gmail.com

// A module to collect Cosmic NoField tracks and write out 'Mille' data files used as input to a
// Millepede-II alignment fit.

// Consult README.md for more information
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdint.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Offline/CosmicReco/inc/MinuitDriftFitter.hh"
#include "Offline/CosmicReco/inc/PDFFit.hh"
#include "Offline/DbTables/inc/TrkAlignElement.hh"
#include "Offline/DbTables/inc/TrkAlignStraw.hh"
#include "Offline/GeneralUtilities/inc/BitMap.hh"
#include "Offline/GeneralUtilities/inc/HepTransform.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Minuit2/MnUserCovariance.h"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"

#include "boost/math/distributions/chi_squared.hpp"
#include "boost/math/distributions/normal.hpp"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/ProducerTable.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"
#include "Offline/TrackerGeom/inc/Plane.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "Offline/DbService/inc/DbHandle.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"

#include "Offline/DataProducts/inc/StrawId.hh"

#include "TAxis.h"
#include "TH1F.h"
#include "TMatrixDSym.h"
#include "TTree.h"

#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "cetlib_except/exception.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/detail/validationException.h"

#include "TrackerAlignment/inc/AlignmentDerivatives.hh"
#include "TrackerAlignment/inc/MilleDataWriter.hh"
#include "TrackerAlignment/inc/AlignmentUtilities.hh"



namespace mu2e {

class AlignTrackCollector : public art::EDAnalyzer {
public:
  size_t _dof_per_plane = 6; // dx, dy, dz, a, b, g (translation, rotation)
  size_t _dof_per_panel = 6; // dx, dy, dz, a, b, g (translation, rotation)

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<int> diaglvl{Name("diagLevel"), Comment("diagnostic level")};

    fhicl::Atom<art::InputTag> costag{    Name("CosmicTrackSeedCollection"), Comment("tag for cosmic track seed collection")};
    fhicl::Atom<bool> millefilegzip{      Name("GzipCompression"), Comment("Enable gzip compression for millepede output file")};
    fhicl::Atom<std::string> millefile{   Name("MilleFile"), Comment("Output filename for Millepede track data file")};
    fhicl::Atom<std::string> steerfile{   Name("SteerFile"), Comment("Output filename for Millepede steering file")};
    fhicl::Atom<std::string> paramfile{   Name("ParamFile"), Comment("Output filename for Millepede parameters file")};
    fhicl::Atom<std::string> extrafile{   Name("ExtraFile"), Comment("Output filename for Millepede parameters file")};
    fhicl::Atom<std::string> constrfile{  Name("ConstrFile"), Comment("Output filename for Millepede constraints file")};
    fhicl::Sequence<std::string> mpsteers{Name("SteeringOpts"), Comment("Additional configuration options to add to generated steering file.")};

    fhicl::Atom<bool> usenumerical{       Name("UseNumericalDiffn"), Comment("Whether or not to use numerical derivatives. Default is false.")};

    fhicl::Atom<int> minplanetraverse{Name("MinTraversedPlanes"), Comment("How many planes must be traversed for a track to be accepted.")};
    fhicl::Atom<double> maxtimeres{Name("MaxTimeRes"), Comment("Require that the maximum ABSOLUTE time residual over all track hits < MaxTimeRes.")};
    fhicl::Atom<int> mintrackhits{Name("MinTrackSH"), Comment("Require that the minimum straw hits in a track > MinTrackSH.")};

    fhicl::Sequence<int> fixplane{        Name("FixPlane"), Comment("Planes to fix. The parameters are fixed to the proditions value.")};
    fhicl::Sequence<int> fixpanel{        Name("FixPanel"), Comment("Panels to fix. The parameters are fixed to the proditions value.")};
    fhicl::Atom<bool> enableplanetranslation{ Name("EnablePlaneTranslationDOF"), Comment("Fit for plane translations")};
    fhicl::Atom<bool> enableplanerotation{ Name("EnablePlaneRotationDOF"), Comment("Fit for plane rotations")};
    fhicl::Atom<bool> enablepaneltranslation{ Name("EnablePanelTranslationDOF"), Comment("Fit for panel translations")};
    fhicl::Atom<bool> enablepanelrotation{ Name("EnablePanelRotationDOF"), Comment("Fit for panel rotations")};

    fhicl::Atom<std::string> weakconstraints{ Name("WeakConstraints"), Comment("Which weak constraint strategy to use. Either 'None', 'Fix', or 'Measurement'")};
    fhicl::Atom<bool> panelconstraints{       Name("PanelConstraints"), Comment("Whether to constrain each panel DOF")};
 
    fhicl::Sequence<float> weakvalues{Name("WeakValues"), Comment("Constrain weak modes at these values"), std::vector<float>{}};
    fhicl::Sequence<float> weaksigmas{Name("WeakSigmas"), Comment("Constrain weak modes with these sigmas"), std::vector<float>{}};
    fhicl::Sequence<float> panelvaluesx{Name("PanelValuesX"), Comment("Constrain panel x at these values"), std::vector<float>{}};
    fhicl::Sequence<float> panelvaluesy{Name("PanelValuesY"), Comment("Constrain panel y at these values"), std::vector<float>{}};
    fhicl::Sequence<float> panelvaluesz{Name("PanelValuesZ"), Comment("Constrain panel z at these values"), std::vector<float>{}};
    fhicl::Sequence<float> panelvaluesrx{Name("PanelValuesRX"), Comment("Constrain panel rx at these values"), std::vector<float>{}};
    fhicl::Sequence<float> panelvaluesry{Name("PanelValuesRY"), Comment("Constrain panel ry at these values"), std::vector<float>{}};
    fhicl::Sequence<float> panelvaluesrz{Name("PanelValuesRZ"), Comment("Constrain panel rz at these values"), std::vector<float>{}};
    fhicl::Atom<float> panelsigmax{Name("PanelSigmaX"), Comment("Constrain panel x with this sigma"),0};
    fhicl::Atom<float> panelsigmay{Name("PanelSigmaY"), Comment("Constrain panel y with this sigma"),0};
    fhicl::Atom<float> panelsigmaz{Name("PanelSigmaZ"), Comment("Constrain panel z with this sigma"),0};
    fhicl::Atom<float> panelsigmarx{Name("PanelSigmaRX"), Comment("Constrain panel rx with this sigma"),0};
    fhicl::Atom<float> panelsigmary{Name("PanelSigmaRY"), Comment("Constrain panel ry with this sigma"),0};
    fhicl::Atom<float> panelsigmarz{Name("PanelSigmaRZ"), Comment("Constrain panel rz with this sigma"),0};
  };

  typedef art::EDAnalyzer::Table<Config> Parameters;

  explicit AlignTrackCollector(const Parameters& settings);
  virtual ~AlignTrackCollector() {}

  void beginJob();
  void endJob();
  void beginRun(art::Run const&);
  void analyze(art::Event const&);

private:

  int _diag;
  art::ProductToken<CosmicTrackSeedCollection> _costag;
  bool _gzipCompress;
  std::string _milleFilename;
  std::string _steerFilename;
  std::string _paramFilename;
  std::string _extraFilename;
  std::string _constrFilename;
  std::vector<std::string> _steerLines;
  bool _useNumericalDerivs;
  int _minPlaneTraverse;
  double _maxTimeRes;
  int _minTrackHits;

  std::vector<bool> _fixedPlanes;
  std::vector<bool> _fixedPanels;
  bool _enablePlaneTranslationDOF, _enablePlaneRotationDOF, _enablePanelTranslationDOF, _enablePanelRotationDOF;

  std::string _weakConstraints;
  bool _panelConstraints;
  std::vector<float> _weakValues;
  std::vector<float> _weakSigmas;
  std::vector<std::vector<float>> _panelValues;
  std::vector<float> _panelSigmas;


  bool _wroteMillepedeParams;
  MilleDataWriter<double> _milleFile;
  std::map<int, int> _DOFCounts;
  int _tracksWritten;
  std::vector<float> _startingAlignPlanes; 
  std::vector<float> _startingAlignPanels; 

  const CosmicTrackSeedCollection* _cscol;
  const Tracker* _nominalTracker;

  ProditionsHandle<Tracker> _alignedTracker_h;
  ProditionsHandle<StrawResponse> srep_h;

  std::unique_ptr<DbHandle<TrkAlignPlane>> _trkAlignPlane_h;
  std::unique_ptr<DbHandle<TrkAlignPanel>> _trkAlignPanel_h;

  // Tree and tree fill members
  static constexpr int MAX_NHITS = 100;
  TTree* _diagtree;
  Int_t _nHits;
  Float_t _doca_residual[MAX_NHITS];
  Float_t _time_residual[MAX_NHITS];
  Float_t _residual_err[MAX_NHITS];
  Float_t _ddx[MAX_NHITS];
  Float_t _ddy[MAX_NHITS];
  Float_t _ddz[MAX_NHITS];
  Float_t _ddrx[MAX_NHITS];
  Float_t _ddry[MAX_NHITS];
  Float_t _ddrz[MAX_NHITS];
  Float_t _ddxpa[MAX_NHITS];
  Float_t _ddypa[MAX_NHITS];
  Float_t _ddzpa[MAX_NHITS];
  Float_t _ddrxpa[MAX_NHITS];
  Float_t _ddrypa[MAX_NHITS];
  Float_t _ddrzpa[MAX_NHITS];
  Float_t _pull_hittime[MAX_NHITS];
  Float_t _doca[MAX_NHITS];
  Float_t _deltax[MAX_NHITS];
  Float_t _deltay[MAX_NHITS];
  Float_t _deltaz[MAX_NHITS];
  Float_t _time[MAX_NHITS];
  Int_t _plane_uid[MAX_NHITS];
  Int_t _panel_uid[MAX_NHITS];
  Double_t _A0;
  Double_t _A1;
  Double_t _B0;
  Double_t _B1;
  Double_t _T0;
  Double_t _chisq;
  Int_t _ndof;
  Double_t _pvalue;
  Int_t _planes_trav;

  int getLabel(int const&, int const&, int const&);
  std::vector<int> generateDOFLabels(uint16_t plane, uint16_t panel);
  std::vector<int> generateDOFLabels(StrawId const& strw);
  std::vector<double> fixDerivativesGlobal(uint16_t plane, uint16_t panel, std::vector<double> &derivativesGlobal);
  bool isDOFenabled(int object_class, int object_id, int dof_n);
  bool isDOFfixed(int object_class, int object_id, int dof_n);
  void cacheStartingParams(TrkAlignPlane const& alignConstPlanes, TrkAlignPanel const& alignConstPanels);
  void writeMillepedeConstraints(Tracker const& tracker);
  void writeMillepedeSteering();
  void writeMillepedeParams();
  void writeMillepedeExtras();

  bool goodTrack(CosmicTrackSeed const& sts);
};

AlignTrackCollector::AlignTrackCollector(const Parameters& conf) :
  art::EDAnalyzer(conf), 
  _diag(conf().diaglvl()),
  _costag(consumes<CosmicTrackSeedCollection>(conf().costag())),
  _gzipCompress(conf().millefilegzip()), 
  _milleFilename(conf().millefile()),
  _steerFilename(conf().steerfile()), 
  _paramFilename(conf().paramfile()),
  _extraFilename(conf().extrafile()),
  _constrFilename(conf().constrfile()),
  _steerLines(conf().mpsteers()),
  _useNumericalDerivs(conf().usenumerical()),
  _minPlaneTraverse(conf().minplanetraverse()),
  _maxTimeRes(conf().maxtimeres()),
  _minTrackHits(conf().mintrackhits()),
      
  _enablePlaneTranslationDOF(conf().enableplanetranslation()),
  _enablePlaneRotationDOF(conf().enableplanerotation()),
  _enablePanelTranslationDOF(conf().enablepaneltranslation()),
  _enablePanelRotationDOF(conf().enablepanelrotation()),
  _weakConstraints(conf().weakconstraints()),
  _panelConstraints(conf().panelconstraints()),
  _weakValues(conf().weakvalues()),
  _weakSigmas(conf().weaksigmas()),

  _wroteMillepedeParams(false), 
  _milleFile(_milleFilename, _gzipCompress),
  _tracksWritten(0)
{
  _fixedPlanes = std::vector<bool>(StrawId::_nplanes,false);
  for (size_t i=0;i<conf().fixplane().size();i++)
    _fixedPlanes[conf().fixplane()[i]] = true;
  _fixedPanels = std::vector<bool>(StrawId::_nupanels,false);
  for (size_t i=0;i<conf().fixpanel().size();i++)
    _fixedPanels[conf().fixpanel()[i]] = true;

  if (_panelConstraints){
    _panelValues.push_back(conf().panelvaluesx());
    _panelValues.push_back(conf().panelvaluesy());
    _panelValues.push_back(conf().panelvaluesz());
    _panelValues.push_back(conf().panelvaluesrx());
    _panelValues.push_back(conf().panelvaluesry());
    _panelValues.push_back(conf().panelvaluesrz());
    _panelSigmas.push_back(conf().panelsigmax());
    _panelSigmas.push_back(conf().panelsigmay());
    _panelSigmas.push_back(conf().panelsigmaz());
    _panelSigmas.push_back(conf().panelsigmarx());
    _panelSigmas.push_back(conf().panelsigmary());
    _panelSigmas.push_back(conf().panelsigmarz());
    for (size_t i=0;i<6;i++){
      if (_panelValues.size() != StrawId::_nupanels){
        if (_panelValues[i].size() == 0){
          std::cout << "PanelValues " << i << " not set! using 0 as default" << std::endl;
        } else{
          std::cout << "Warning! Incorrect number of panelValues " << i << std::endl;
        }
        _panelValues[i] = std::vector<float>(StrawId::_nupanels,0);
      }
    }
  }

  for (size_t p=0;p<StrawId::_nplanes;p++){
    for (size_t j=0;j<_dof_per_plane;j++)
      _DOFCounts[getLabel(1,p,j)] = 0;
  }
  for (size_t p=0;p<StrawId::_nupanels;p++){
    for (size_t j=0;j<_dof_per_panel;j++)
      _DOFCounts[getLabel(2,p,j)] = 0;
  }
}



void AlignTrackCollector::beginJob() {
  _trkAlignPlane_h = std::make_unique<DbHandle<TrkAlignPlane>>();
  _trkAlignPanel_h = std::make_unique<DbHandle<TrkAlignPanel>>();

  if (_diag > 0) {
    art::ServiceHandle<art::TFileService> tfs;

    _diagtree = tfs->make<TTree>("tracks", "Tracks collected for an alignment iteration");
    _diagtree->Branch("nHits", &_nHits, "nHits/I");
    _diagtree->Branch("doca_resid", &_doca_residual, "doca_resid[nHits]/F");
    _diagtree->Branch("time_resid", &_time_residual, "time_resid[nHits]/F");

    _diagtree->Branch("ddx",&_ddx,"ddx[nHits]/F");
    _diagtree->Branch("ddy",&_ddy,"ddy[nHits]/F");
    _diagtree->Branch("ddz",&_ddz,"ddz[nHits]/F");
    _diagtree->Branch("ddrx",&_ddrx,"ddrx[nHits]/F");
    _diagtree->Branch("ddry",&_ddry,"ddry[nHits]/F");
    _diagtree->Branch("ddrz",&_ddrz,"ddrz[nHits]/F");
    _diagtree->Branch("ddxpa",&_ddxpa,"ddxpa[nHits]/F");
    _diagtree->Branch("ddypa",&_ddypa,"ddypa[nHits]/F");
    _diagtree->Branch("ddzpa",&_ddzpa,"ddzpa[nHits]/F");
    _diagtree->Branch("ddrxpa",&_ddrxpa,"ddrxpa[nHits]/F");
    _diagtree->Branch("ddrypa",&_ddrypa,"ddrypa[nHits]/F");
    _diagtree->Branch("ddrzpa",&_ddrzpa,"ddrzpa[nHits]/F");

    _diagtree->Branch("resid_err", &_residual_err, "resid_err[nHits]/F");

    _diagtree->Branch("pull_hittime", &_pull_hittime, "pull_hittime[nHits]/F");

    _diagtree->Branch("doca", &_doca, "doca[nHits]/F");
    _diagtree->Branch("deltax", &_deltax, "deltax[nHits]/F");
    _diagtree->Branch("deltay", &_deltay, "deltay[nHits]/F");
    _diagtree->Branch("deltaz", &_deltaz, "deltaz[nHits]/F");
    _diagtree->Branch("time", &_time, "time[nHits]/F");
    _diagtree->Branch("plane", &_plane_uid, "plane[nHits]/I");
    _diagtree->Branch("panel", &_panel_uid, "panel[nHits]/I");

    _diagtree->Branch("A0", &_A0, "A0/D");
    _diagtree->Branch("A1", &_A1, "A1/D");
    _diagtree->Branch("B0", &_B0, "B0/D");
    _diagtree->Branch("B1", &_B1, "B1/D");
    _diagtree->Branch("T0", &_T0, "T0/D");

    _diagtree->Branch("chisq", &_chisq, "chisq/D");

    _diagtree->Branch("ndof", &_ndof, "ndof/I");
    _diagtree->Branch("pvalue", &_pvalue, "pvalue/D");

    _diagtree->Branch("planes_trav", &_planes_trav, "planes_trav/I");
  }

}

void AlignTrackCollector::endJob() {
  if (_diag > 0) {
    std::cout << "AlignTrackCollector: wrote " 
              << _tracksWritten << " tracks to " 
              << _milleFilename
              << std::endl;
  }

  writeMillepedeConstraints(*_nominalTracker);
  writeMillepedeParams();
  writeMillepedeSteering();
  writeMillepedeExtras();
}

void AlignTrackCollector::beginRun(art::Run const& run) {
  GeomHandle<Tracker> track;
  _nominalTracker = track.get();
}

void AlignTrackCollector::analyze(art::Event const& event) {
  StrawResponse const& _srep = srep_h.get(event.id());
  Tracker const& alignedTracker = _alignedTracker_h.get(event.id());

  auto alignConsts_planes = _trkAlignPlane_h->get(event.id());
  auto alignConsts_panels = _trkAlignPanel_h->get(event.id());

  if (!_wroteMillepedeParams) {
    cacheStartingParams(alignConsts_planes, alignConsts_panels);
    _wroteMillepedeParams = true;
  }

  auto stH = event.getValidHandle<CosmicTrackSeedCollection>(_costag);

  if (stH.product() == 0) {
    return;
  }

  CosmicTrackSeedCollection const& coscol = *stH.product();

  // dedicated to CosmicTrackSeedCollection
  for (CosmicTrackSeed const& sts : coscol) {
    if (!goodTrack(sts))
      continue;

    CosmicTrack const& st = sts._track;

    std::set<uint16_t> planes_traversed;

    _A0 = st.MinuitParams.A0;
    _A1 = st.MinuitParams.A1;
    _B0 = st.MinuitParams.B0;
    _B1 = st.MinuitParams.B1;
    _T0 = st.MinuitParams.T0;

    GaussianDriftFit fit_object(sts._straw_chits, _srep, &alignedTracker);


    // for the max timeresidual track quality cut
    double max_time_res_track = -1;

    int ngood_hits = 0;
    _nHits = 0;
    _pvalue = 0;
    _chisq = 0;
    bool bad_track = false;

    std::vector<double> residuals;
    std::vector<double> residual_errs;
    std::vector<std::vector<double>> global_derivs_temp;
    std::vector<std::vector<double>> local_derivs_temp;
    std::vector<std::vector<int>> labels_temp;


    // get residuals and their derivatives with respect
    // to all local and global parameters

    for (ComboHit const& straw_hit : sts._straw_chits) {
      // straw and plane info
      StrawId const& straw_id = straw_hit.strawId();
      Straw const& alignedStraw = alignedTracker.getStraw(straw_id);

      auto plane_id = straw_id.getPlane();
      auto panel_id = straw_id.uniquePanel();

      auto fitpos = GenVector::Hep3Vec(st.MinuitEquation.Pos);
      auto fitdir = GenVector::Hep3Vec(st.MinuitEquation.Dir.unit());

      // hit DOCA, DOCA/TOCA residuals, errors
      TwoLinePCA pca(fitpos, fitdir, 
          alignedStraw.wirePosition(), alignedStraw.wireDirection());

      double time_resid = fit_object.TimeResidual(straw_hit, sts);
      double drift_res = _srep.driftTimeError(straw_hit.strawId(), pca.dca(), 0);

      if (isnan(time_resid) || isnan(drift_res)) {
        bad_track = true;
        break;
      }

      // alignment constants for derivatives
      auto const& rowpl = alignConsts_planes.rowAt(plane_id);
      auto const& rowpa = alignConsts_panels.rowAt(panel_id);

      // now calculate the derivatives
      std::vector<double> derivativesLocal, derivativesGlobal;

      AlignmentUtilities::CosmicTimeTrack track { _A0, _B0, _A1, _B1, _T0 };

      if (_useNumericalDerivs){
        std::tie(derivativesLocal, derivativesGlobal) = AlignmentUtilities::numericalDerivatives(
            track, straw_id, rowpl, rowpa, *_nominalTracker, _srep);
      }else{
        auto dca_dir = (pca.point2()-pca.point1()).unit();

        //FIXME
        double drift_time_p = _srep.driftDistanceToTime(straw_id, pca.dca()+1e-5, 0) +
          _srep.driftTimeOffset(straw_id, pca.dca()+1e-5, 0);
        double drift_time_n = _srep.driftDistanceToTime(straw_id, pca.dca()-1e-5, 0) +
          _srep.driftTimeOffset(straw_id, pca.dca()-1e-5, 0);
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

        const Panel& panel = _nominalTracker->getPanel(straw_id);
        const Plane& plane = _nominalTracker->getPlane(straw_id);
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
      }

      for (size_t i=0;i<derivativesGlobal.size();i++)
        derivativesGlobal[i] *= -1;
      std::vector<double> derivativesGlobalFixed = fixDerivativesGlobal(straw_id.getPlane(), straw_id.uniquePanel(), derivativesGlobal);

      planes_traversed.insert(plane_id);

      // Time residual cut
      // avoid outlier hits when applying this cut
      // FIXME
      if (!straw_hit._flag.hasAnyProperty(StrawHitFlag::outlier)) {
        if (abs(time_resid) > max_time_res_track) {
          max_time_res_track = abs(time_resid);
        }
      }

      std::vector<int> labels = generateDOFLabels(straw_id);

      if (labels.size() != derivativesGlobalFixed.size()) {
        throw cet::exception("RECO") << "N global derivatives != N labels"
          << " ... Something is wrong!";
      }

      // write the hit to the track buffer
      residuals.emplace_back(time_resid);
      residual_errs.emplace_back(drift_res);

      global_derivs_temp.emplace_back(derivativesGlobalFixed);
      local_derivs_temp.emplace_back(derivativesLocal);
      labels_temp.push_back(labels);

      ngood_hits++;

      if (_nHits < MAX_NHITS){
        _doca_residual[_nHits] = fit_object.DOCAresidual(straw_hit, track.as_vector());
        _time_residual[_nHits] = time_resid;
        _residual_err[_nHits] = drift_res;
        _pull_hittime[_nHits] = time_resid/drift_res;
        _doca[_nHits] = pca.dca();
        _deltax[_nHits] = (pca.point1()-pca.point2()).x();
        _deltay[_nHits] = (pca.point1()-pca.point2()).y();
        _deltaz[_nHits] = (pca.point1()-pca.point2()).z();
        _time[_nHits] = straw_hit.time();
        _plane_uid[_nHits] = straw_id.plane();
        _panel_uid[_nHits] = straw_id.uniquePanel();
        _ddx[_nHits] = derivativesGlobal[0];
        _ddy[_nHits] = derivativesGlobal[1];
        _ddz[_nHits] = derivativesGlobal[2];
        _ddrx[_nHits] = derivativesGlobal[3];
        _ddry[_nHits] = derivativesGlobal[4];
        _ddrz[_nHits] = derivativesGlobal[5];
        _ddxpa[_nHits] = derivativesGlobal[6];
        _ddypa[_nHits] = derivativesGlobal[7];
        _ddzpa[_nHits] = derivativesGlobal[8];
        _ddrxpa[_nHits] = derivativesGlobal[9];
        _ddrypa[_nHits] = derivativesGlobal[10];
        _ddrzpa[_nHits] = derivativesGlobal[11];
        _nHits++;
      }
    }

    if (bad_track){
      if (_diag > 1)
        std::cout << "AlignTrackCollector: track failed quality cuts" << std::endl;
      if (_diag > 2) 
        std::cout << "  reason: bad track" << std::endl;
      continue;
    }


    _ndof = ngood_hits - 5;

    // track acceptance cuts
    if (_ndof <= 0){
      if (_diag > 1)
        std::cout << "AlignTrackCollector: track failed quality cuts" << std::endl;
      if (_diag > 2) 
        std::cout << "  reason: insufficient dof " << _ndof << std::endl;
      continue;
    }

    _pvalue = 1.0-boost::math::cdf(boost::math::chi_squared(_ndof), _chisq);
    _chisq /= _ndof;
    _planes_trav = planes_traversed.size();

    if ((int) planes_traversed.size() < _minPlaneTraverse){
      if (_diag > 1)
        std::cout << "AlignTrackCollector: track failed quality cuts" << std::endl;
      if (_diag > 2) 
        std::cout << "  reason: planes traversed " << planes_traversed.size() << std::endl;
      continue;
    }
    if (max_time_res_track > _maxTimeRes){
      if (_diag > 1)
        std::cout << "AlignTrackCollector: track failed quality cuts" << std::endl;
      if (_diag > 2) 
        std::cout << "  reason: max time residual " << max_time_res_track << std::endl;
      continue;
    }
    if (ngood_hits < _minTrackHits){
      if (_diag > 1)
        std::cout << "AlignTrackCollector: track failed quality cuts" << std::endl;
      if (_diag > 2) 
        std::cout << "  reason: insufficient hits " << ngood_hits << std::endl;
      continue;
    }

    // write hits to buffer
    for (size_t i=0;i<residuals.size();i++){
      _milleFile.pushHit(local_derivs_temp[i], global_derivs_temp[i], labels_temp[i],
          residuals[i], residual_errs[i]);

      for (size_t j=0;j<labels_temp[i].size();j++)
        _DOFCounts[labels_temp[i][j]]++;
    }
    if (_diag > 0){
      _diagtree->Fill();
    }
    
    // Write the track buffer to file
    _milleFile.flushTrack();

    _tracksWritten++;

    if (_diag > 1) {
      std::cout << "wrote track " << _tracksWritten << std::endl;
    }
  }
}

bool AlignTrackCollector::goodTrack(CosmicTrackSeed const& sts){
    CosmicTrack const& st = sts._track;

    if (!sts._status.hasAllProperties(TrkFitFlag::helixOK)) {
      return false;
    }

    if (!st.converged || !st.minuit_converged) {
      return false;
    }

    if (isnan(st.MinuitParams.A0)) {
      return false;
    }
    return true;
}


int AlignTrackCollector::getLabel(int const& object_cls, int const& obj_uid, int const& dof_id) {
  // object class: 0 - 9 - i.e. 1 for planes, 2 for panels
  // object unique id: 0 - 999 supports up to 999 unique objects which is fine for this level of
  // alignment
  // object dof id: 0 - 9

  // (future) a straw label might be:
  // id  plane panel straw dof
  // 3   00    0     00    0
  // min
  // 3 00 0 00 0
  // max
  // 3,356,959

  // 1 000 0

  return object_cls * 10000 + obj_uid * 10 + dof_id;
}

std::vector<int> AlignTrackCollector::generateDOFLabels(StrawId const& strw) {
  return generateDOFLabels(strw.getPlane(), strw.uniquePanel());
}

std::vector<int> AlignTrackCollector::generateDOFLabels(uint16_t plane, uint16_t panel) {
  std::vector<int> labels;
  labels.reserve(_dof_per_plane + _dof_per_panel);

  for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
    if (!isDOFenabled(1, plane, dof_n)) {
      continue;
    }
    labels.push_back(getLabel(1, plane, dof_n));
  }
  for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
    if (!isDOFenabled(2, panel, dof_n)) {
      continue;
    }
    labels.push_back(getLabel(2, panel, dof_n));
  }
  
  return labels;
}

std::vector<double> AlignTrackCollector::fixDerivativesGlobal(uint16_t plane, uint16_t panel, std::vector<double> &derivativesGlobal) {
  std::vector<double> fixedDerivativesGlobal;
  size_t index = 0;
  for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
    if (isDOFenabled(1, plane, dof_n)) {
      fixedDerivativesGlobal.push_back(derivativesGlobal[index]);
    }
    index++;
  }
  for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
    if (isDOFenabled(2, panel, dof_n)) {
      fixedDerivativesGlobal.push_back(derivativesGlobal[index]);
    }
    index++;
  }
  
  return fixedDerivativesGlobal;
}

bool AlignTrackCollector::isDOFenabled(int object_class, int object_id, int dof_n) {
  if (object_class == 2 && !_enablePanelTranslationDOF && dof_n <= 2) {
    return false;
  }
  if (object_class == 2 && !_enablePanelRotationDOF && dof_n > 2) {
    return false;
  }
  if (object_class == 2 && dof_n == 0)
    return false;
  if (object_class == 1 && !_enablePlaneTranslationDOF && dof_n <= 2) {
    return false;
  }
  if (object_class == 1 && !_enablePlaneRotationDOF && dof_n > 2) {
    return false;
  }
  return true;
}

bool AlignTrackCollector::isDOFfixed(int object_class, int object_id, int dof_n) {
  if (object_class == 1) 
    return _fixedPlanes[object_id];
  if (object_class == 2) 
    return _fixedPanels[object_id];
  return false;
}

void AlignTrackCollector::writeMillepedeConstraints(Tracker const& nominalTracker){
  std::ofstream output_file(_constrFilename);
  output_file << "! Generated by AlignTrackCollector" << std::endl;
  output_file << "! Weak constraint strategy: " << _weakConstraints << std::endl;

  // plane overall translations and rotations are fixed
  output_file << "! planes" << std::endl;
  for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {

    // x and y parallelogramming not constrained here
    if (dof_n == 3 || dof_n == 4)
      continue;
    // check if any of these DOF are enabled
    bool has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, dof_n)) {
        has_enabled = true;
        break;
      }
    }
    if (!has_enabled)
      continue;
    output_file << "Constraint   0" << std::endl;
    for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
      if (!isDOFenabled(1, p, dof_n))
        continue;
      if (dof_n == 0 || dof_n == 2 || dof_n == 3 || dof_n == 5){
        // if this plane is rotated, need to invert constraint
        if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0)
          output_file << getLabel(1, p, dof_n) << "    1" << std::endl;
        else
          output_file << getLabel(1, p, dof_n) << "    -1" << std::endl;
      }else{
        output_file << getLabel(1, p, dof_n) << "    1" << std::endl;
      }
    }
  }

  // panel overall z translation is fixed (relative to plane)
  output_file << "! panels follow" << std::endl;
  // one constraint per plane
  for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
    // check if any of these DOF are enabled
    bool has_enabled = false;
    for (uint16_t pa=0;pa<StrawId::_npanels;++pa){
      if (isDOFenabled(2, p*StrawId::_npanels+pa, 2)) {
        has_enabled = true;
        break;
      }
    }
    if (!has_enabled)
      continue;

    output_file << "Constraint   0" << std::endl;
    for (uint16_t pa=0; pa < StrawId::_npanels; ++ pa){
      if (!isDOFenabled(2, p*StrawId::_npanels+pa, 2))
        continue;
      StrawId tempid(p,pa,0);
      if (nominalTracker.getPanel(tempid).panelToDS().rotation().getTheta() == 0)
        output_file << getLabel(2, p*StrawId::_npanels+pa, 2) << "    1" << std::endl;
      else
        output_file << getLabel(2, p*StrawId::_npanels+pa, 2) << "    -1" << std::endl;
    }
  }

  // panel v translations dot xhat and yhat are fixed (equal to overall plane translation)
  CLHEP::Hep3Vector xhat(1,0,0);
  CLHEP::Hep3Vector yhat(0,1,0);
  for (size_t pl=0;pl<StrawId::_nplanes;pl++){
    bool has_enabled = false;
    for (uint16_t pa=0;pa<StrawId::_npanels;++pa){
      if (isDOFenabled(2, pl*StrawId::_npanels+pa, 1)) {
        has_enabled = true;
        break;
      }
    }
    if (!has_enabled)
      continue;
    output_file << "Constraint   0" << std::endl;
    for (size_t pa=0;pa<StrawId::_npanels;pa++){
      if (!isDOFenabled(2, pl*StrawId::_npanels+pa, 1))
        continue;
      StrawId tempid(pl,pa,0);
      output_file << getLabel(2, pl*StrawId::_npanels+pa, 1) << "    " << nominalTracker.getPanel(tempid).vDirection().dot(xhat) << std::endl;
    }
    output_file << "Constraint   0" << std::endl;
    for (size_t pa=0;pa<StrawId::_npanels;pa++){
      if (!isDOFenabled(2, pl*StrawId::_npanels+pa, 1))
        continue;
      StrawId tempid(pl,pa,0);
      output_file << getLabel(2, pl*StrawId::_npanels+pa, 1) << "    " << nominalTracker.getPanel(tempid).vDirection().dot(yhat) << std::endl;
    }
  }

  if (_panelConstraints){
    for (uint16_t p=0; p< StrawId::_nupanels; ++p){
      for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
        if (!isDOFenabled(2, p, dof_n)) {
          continue;
        }
        output_file << "Measurement " << _panelValues[dof_n][p]-_startingAlignPanels[p*6+dof_n] << " " << _panelSigmas[dof_n] << std::endl;
        output_file << getLabel(2, p, dof_n) << "  1" << std::endl;
      }
    }
  }

  output_file << "! weak modes follow" << std::endl;
  if (_weakConstraints == "None"){
    return;
    // if fixing specific panels then no further overall constraints
  }else{

    // calculate current values
    double xskew = 0;
    double yskew = 0;
    double xparallel = 0;
    double yparallel = 0;
    double zsqueeze = 0;
    double ztwist = 0;
    for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
      // FIXME check if enabled
      if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0){
        xskew += _startingAlignPlanes[p*6+0] * (p-17.5);
        yskew += _startingAlignPlanes[p*6+1] * (p-17.5);
        zsqueeze += _startingAlignPlanes[p*6+2] * (p-17.5);
        ztwist += _startingAlignPlanes[p*6+5] * (p-17.5);
        xparallel += _startingAlignPlanes[p*6+3];
        yparallel += _startingAlignPlanes[p*6+4];
      }else{
        xskew -= _startingAlignPlanes[p*6+0] * (p-17.5);
        yskew += _startingAlignPlanes[p*6+1] * (p-17.5);
        zsqueeze -= _startingAlignPlanes[p*6+2] * (p-17.5);
        ztwist -= _startingAlignPlanes[p*6+5] * (p-17.5);
        xparallel -= _startingAlignPlanes[p*6+3];
        yparallel += _startingAlignPlanes[p*6+4];
      }
    }
    std::vector<double> vpairoffsets(6,0);
    std::vector<double> vpairslopes(6,0);
    for (uint16_t pl = 0; pl < StrawId::_nplanes; ++pl) {
      for (uint16_t pa = 0; pa < StrawId::_npanels; ++pa) {
        uint16_t p = pl*6+pa;
        // FIXME check if enabled
        if (nominalTracker.getPlane(pl).planeToDS().rotation().getTheta() == 0){
          vpairoffsets[pa] += _startingAlignPanels[p*6+1];
          vpairslopes[pa]  += _startingAlignPanels[p*6+1] * (pl-17.5);
        }else{
          vpairoffsets[5-pa] -= _startingAlignPanels[p*6+1];
          vpairslopes[5-pa] -= _startingAlignPanels[p*6+1] * (pl-17.5);
        }
      }
    }

    // Fix all weak modes to 0
    // or Constrain weak modes by measurements

    // X skew
    bool has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, 0)) {
        has_enabled = true;
        break;
      }
    }
    if (has_enabled){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << -1*xskew << std::endl;
      }else{
        output_file << "Measurement " << _weakValues[0]-xskew << " " << _weakSigmas[0] << std::endl;
      }
      for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        if (!isDOFenabled(1, p, 0))
          continue;
        // if this plane is rotated, need to invert constraint
        if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0)
          output_file << getLabel(1, p, 0) << "    " << (p-17.5) << std::endl;
        else
          output_file << getLabel(1, p, 0) << "    " << -1*(p-17.5) << std::endl;
      }
    }
    
    // Y skew
    has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, 1)) {
        has_enabled = true;
        break;
      }
    }
    if (has_enabled){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << -1*yskew << std::endl;
      }else{
        output_file << "Measurement " << _weakValues[1]-yskew << " " << _weakSigmas[1] << std::endl;
      }
      for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        if (!isDOFenabled(1, p, 1))
          continue;
        output_file << getLabel(1, p, 1) << "    " << (p-17.5) << std::endl;
      }
    }

    // Z squeeze
    has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, 2)) {
        has_enabled = true;
        break;
      }
    }
    if (has_enabled){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << -1*zsqueeze << std::endl;
      }else{
        output_file << "Measurement " << _weakValues[2]-zsqueeze << " " << _weakSigmas[2] << std::endl;
      }
      for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        if (!isDOFenabled(1, p, 2))
          continue;
        // if this plane is rotated, need to invert constraint
        if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0)
          output_file << getLabel(1, p, 2) << "    " << (p-17.5) << std::endl;
        else
          output_file << getLabel(1, p, 2) << "    " << -1*(p-17.5) << std::endl;
      }
    }
    // X parallel
    has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, 3)) {
        has_enabled = true;
        break;
      }
    }
    if (has_enabled){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << -1*xparallel << std::endl;
      }else{
        output_file << "Measurement " << _weakValues[3]-xparallel << " " << _weakSigmas[3] << std::endl;
      }
      for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        if (!isDOFenabled(1, p, 3))
          continue;
        // if this plane is rotated, need to invert constraint
        if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0)
          output_file << getLabel(1, p, 3) << "    " << 1 << std::endl;
        else
          output_file << getLabel(1, p, 3) << "    " << -1 << std::endl;
      }
    }
    // Y parallel
    has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, 4)) {
        has_enabled = true;
        break;
      }
    }
    if (has_enabled){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << -1*yparallel << std::endl;
      }else{
        output_file << "Measurement " << _weakValues[4]-yparallel << " " << _weakSigmas[4] << std::endl;
      }
      for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        if (!isDOFenabled(1, p, 4))
          continue;
        output_file << getLabel(1, p, 4) << "    " << 1 << std::endl;
      }
    }
    // z twist
    has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, 5)) {
        has_enabled = true;
        break;
      }
    }
    if (has_enabled){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << -1*ztwist << std::endl;
      }else{
        output_file << "Measurement " << _weakValues[5]-ztwist << " " << _weakSigmas[5] << std::endl;
      }
      for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        if (!isDOFenabled(1, p, 5))
          continue;
        // if this plane is rotated, need to invert constraint
        if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0)
          output_file << getLabel(1, p, 5) << "    " << (p-17.5) << std::endl;
        else
          output_file << getLabel(1, p, 5) << "    " << -1*(p-17.5) << std::endl;
      }
    }
    //FIXME starting value
    //FIXME has enabled
    // panel radial out pair
    for (size_t i=0;i<6;i++){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << -1*vpairoffsets[i] << std::endl;
      }else{
        output_file << "Measurement " << _weakValues[6+i]-vpairoffsets[i] << " " << _weakSigmas[6+i] << std::endl;
      }
      for (uint16_t pl = 0; pl < StrawId::_nplanes; ++pl) {
        for (uint16_t pa = 0; pa < StrawId::_npanels; ++pa) {
          auto p = pl*6 + pa;
          if (!isDOFenabled(2, p, 1))
            continue;
          if (nominalTracker.getPlane(pl).planeToDS().rotation().getTheta() == 0){
            if (pa == i)
              output_file << getLabel(2, p, 1) << "    " << 1 << std::endl;
          }else{
            if (pa == 5-i)
              output_file << getLabel(2, p, 1) << "    " << -1 << std::endl;
          }
        }
      }
    }
    for (size_t i=0;i<6;i++){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << -1*vpairslopes[i] << std::endl;
      }else{
        output_file << "Measurement " << _weakValues[7+i]-vpairslopes[i] << " " << _weakSigmas[7+i] << std::endl;
      }
      for (uint16_t pl = 0; pl < StrawId::_nplanes; ++pl) {
        for (uint16_t pa = 0; pa < StrawId::_npanels; ++pa) {
          auto p = pl*6 + pa;
          if (!isDOFenabled(2, p, 1))
            continue;
          if (nominalTracker.getPlane(pl).planeToDS().rotation().getTheta() == 0){
            if (pa == i)
              output_file << getLabel(2, p, 1) << "    " << (pl-17.5) << std::endl;
          }else{
            if (pa == 5-i)
              output_file << getLabel(2, p, 1) << "    " << -1*(pl-17.5) << std::endl;
          }
        }
      }
    }
  }
}

void AlignTrackCollector::writeMillepedeSteering() {
  std::ofstream output_file(_steerFilename);

  output_file << "! Steering file generated by AlignTrackCollector" << std::endl
              << "Cfiles" << std::endl
              << _paramFilename << std::endl
              << _constrFilename << std::endl
              << _milleFilename << std::endl
              << std::endl;

  for (std::string const& line : _steerLines) {
    output_file << line << std::endl;
  }
  
  //FIXME
  output_file << std::endl
              << "method inversion 10 0.001" << std::endl
              << "end" << std::endl;
}

void AlignTrackCollector::writeMillepedeParams() {
  // write a params.txt telling millepede what the alignment constants
  // were for this track collection iteration
  std::ofstream output_file(_paramFilename);
  output_file << "! Generated by AlignTrackCollector" << std::endl;
  output_file << "! Weak constraint strategy: " << _weakConstraints << std::endl;
  output_file << "! columns: label, alignment parameter start value, "
                 "presigma (-ve: fixed, 0: variable)" << std::endl;

  output_file << "Parameter" << std::endl;
  output_file << "! Plane DOFs:" << std::endl;

  for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
    for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
      if (!isDOFenabled(1, p, dof_n)) {
        continue;
      }

      // label(int)   initial value (float)   presigma (-ve: fixed, 0: variable)
      output_file << getLabel(1, p, dof_n) << "  " << 0 << "  " 
                  << (isDOFfixed(1, p, dof_n) ? "-1" : "0") << std::endl;
    }
  }
  output_file << std::endl << "! end plane DOFs" << std::endl;

  output_file << std::endl << "! Panel DOFs:" << std::endl;

  for (uint16_t p = 0; p < StrawId::_nupanels; ++p) {
    for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
      if (!isDOFenabled(2, p, dof_n)) {
        continue;
      }

      // label(int)   initial value (float)   presigma (-ve: fixed, 0: variable)
      output_file << getLabel(2, p, dof_n) << "  " << 0 << "  " << (isDOFfixed(2, p, dof_n) ? "-1" : "0") << std::endl;
    }
  }
  output_file << std::endl << "! end panel DOFs" << std::endl;

  output_file.close();
}

void AlignTrackCollector::writeMillepedeExtras() {
  std::ofstream output_file2(_extraFilename);
  for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
    for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
      if (!isDOFenabled(1, p, dof_n)) {
        continue;
      }

      // label(int)   initial value (float)   presigma (-ve: fixed, 0: variable)
      output_file2 << getLabel(1, p, dof_n) << "  " << _startingAlignPlanes[p*6+dof_n] << "  " << _DOFCounts[getLabel(1,p,dof_n)] << std::endl;
    }
  }
  for (uint16_t p = 0; p < StrawId::_nupanels; ++p) {
    for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
      if (!isDOFenabled(2, p, dof_n)) {
        continue;
      }

      // label(int)   initial value (float)   presigma (-ve: fixed, 0: variable)
      output_file2 << getLabel(2, p, dof_n) << "  " << _startingAlignPanels[p*6+dof_n] << "  " << _DOFCounts[getLabel(2,p,dof_n)] << std::endl;
    }
  }

  output_file2.close();
}

void AlignTrackCollector::cacheStartingParams(TrkAlignPlane const& alignConstPlanes,
                                               TrkAlignPanel const& alignConstPanels) {
  _startingAlignPlanes = std::vector<float>(StrawId::_nplanes*6,0);
  _startingAlignPanels = std::vector<float>(StrawId::_nupanels*6,0);
  for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
    auto const& rowpl = alignConstPlanes.rowAt(p);

    std::vector<float> pl_consts{rowpl.dx(), rowpl.dy(), rowpl.dz(),
                                 rowpl.rx(), rowpl.ry(), rowpl.rz()};

    for (size_t i=0;i<6;i++) 
      _startingAlignPlanes[p*6+i] = pl_consts[i];
  }
  for (uint16_t p = 0; p < StrawId::_nupanels; ++p) {
    auto const& rowpa = alignConstPanels.rowAt(p);

    std::vector<float> pa_consts{rowpa.dx(), rowpa.dy(), rowpa.dz(),
                                 rowpa.rx(), rowpa.ry(), rowpa.rz()};

    for (size_t i=0;i<6;i++) 
      _startingAlignPanels[p*6+i] = pa_consts[i];
  }
}

} // namespace mu2e

using mu2e::AlignTrackCollector;
DEFINE_ART_MODULE(AlignTrackCollector)
