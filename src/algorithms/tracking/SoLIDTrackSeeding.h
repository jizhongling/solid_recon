#pragma once

#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4hep/Vector3f.h>
#include <spdlog/logger.h>
#include <memory>
#include <vector>
#include <tuple>

#include "ActsGeometryProvider.h"
#include "DD4hepBField.h"
#include "SoLIDTrackSeedingConfig.h"
#include "algorithms/interfaces/WithPodConfig.h"

namespace eicrecon {

class SoLIDTrackSeeding : public WithPodConfig<SoLIDTrackSeedingConfig> {
public:
    void init(std::shared_ptr<const ActsGeometryProvider> geo_svc, 
              std::shared_ptr<spdlog::logger> log);
    
    std::unique_ptr<edm4eic::TrackParametersCollection> 
    produce(const edm4eic::TrackerHitCollection& trk_hits);

private:
    struct Seed {
        std::vector<edm4eic::TrackerHit> hits;
        float momentum;
        float theta;
        float phi;
        edm4hep::Vector3f position;
        int planeCombo; // bit mask: bit0=plane3, bit1=plane4, bit2=plane5
    };

    void configure();
    
    // Find seeds using pairs of hits
    std::vector<Seed> findDoubleSeeds(
        const std::vector<edm4eic::TrackerHit>& plane1_hits,
        const std::vector<edm4eic::TrackerHit>& plane2_hits,
        int plane1_idx, int plane2_idx);
    
    // Calculate initial parameters for a hit pair
    bool calcInitParForPair(
        const edm4eic::TrackerHit& hit1,
        const edm4eic::TrackerHit& hit2,
        float& momentum,
        float& theta,
        float& phi);

    // Get layer index from hit
    int getLayer(const edm4eic::TrackerHit& hit);
    
    // Propagate track to target and check validity
    bool propagateToTarget(Seed& seed);
    
    // Propagate track to ECal and check matching
    bool matchECal(const Seed& seed);
    
    // Merge seeds that share common hits
    std::vector<Seed> mergeSeeds(const std::vector<Seed>& seeds);
    
    // Convert seeds to track parameters
    std::unique_ptr<edm4eic::TrackParametersCollection> 
    makeTrackParams(const std::vector<Seed>& seeds);
    
    // Helper: group hits by plane
    std::map<int, std::vector<edm4eic::TrackerHit>> 
    groupHitsByPlane(const edm4eic::TrackerHitCollection& hits);

    std::shared_ptr<spdlog::logger> m_log;
    std::shared_ptr<const ActsGeometryProvider> m_geoSvc;
    std::shared_ptr<const eicrecon::BField::DD4hepBField> m_BField;
    Acts::MagneticFieldContext m_fieldctx;
};

} // namespace eicrecon