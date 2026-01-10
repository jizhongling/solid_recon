#include "SoLIDTrackSeeding.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>
#include <Acts/EventData/GenericBoundTrackParameters.hpp>
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <DD4hep/BitFieldCoder.h>
#include <edm4eic/Cov6f.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4hep/Vector2f.h>
#include <edm4hep/Vector3f.h>
#include <spdlog/logger.h>
#include <memory>
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>

namespace eicrecon {

void SoLIDTrackSeeding::init(
    std::shared_ptr<const ActsGeometryProvider> geo_svc,
    std::shared_ptr<spdlog::logger> log) {
    
    m_log = log;
    m_geoSvc = geo_svc;
    m_BField = std::dynamic_pointer_cast<const eicrecon::BField::DD4hepBField>(
        m_geoSvc->getFieldProvider());
    m_fieldctx = eicrecon::BField::BFieldVariant(m_BField);
    
    configure();
}

void SoLIDTrackSeeding::configure() {
    // Any additional configuration can go here
    m_log->info("SoLIDTrackSeeding configured with planes {}, {}, {}",
                m_cfg.seedPlane1, m_cfg.seedPlane2, m_cfg.seedPlane3);
}

std::unique_ptr<edm4eic::TrackParametersCollection> 
SoLIDTrackSeeding::produce(const edm4eic::TrackerHitCollection& trk_hits) {
    
    // Group hits by plane (layer)
    auto hits_by_plane = groupHitsByPlane(trk_hits);
    
    std::vector<Seed> all_seeds;
    
    // Find seeds using different plane combinations
    // Combo 1: Planes 4 & 5
    if (hits_by_plane.count(m_cfg.seedPlane1) && 
        hits_by_plane.count(m_cfg.seedPlane2)) {
        auto seeds = findDoubleSeeds(
            hits_by_plane[m_cfg.seedPlane1],
            hits_by_plane[m_cfg.seedPlane2],
            m_cfg.seedPlane1, m_cfg.seedPlane2);
        all_seeds.insert(all_seeds.end(), seeds.begin(), seeds.end());
    }
    
    // Combo 2: Planes 4 & 6
    if (hits_by_plane.count(m_cfg.seedPlane1) && 
        hits_by_plane.count(m_cfg.seedPlane3)) {
        auto seeds = findDoubleSeeds(
            hits_by_plane[m_cfg.seedPlane1],
            hits_by_plane[m_cfg.seedPlane3],
            m_cfg.seedPlane1, m_cfg.seedPlane3);
        all_seeds.insert(all_seeds.end(), seeds.begin(), seeds.end());
    }
    
    // Combo 3: Planes 5 & 6
    if (hits_by_plane.count(m_cfg.seedPlane2) && 
        hits_by_plane.count(m_cfg.seedPlane3)) {
        auto seeds = findDoubleSeeds(
            hits_by_plane[m_cfg.seedPlane2],
            hits_by_plane[m_cfg.seedPlane3],
            m_cfg.seedPlane2, m_cfg.seedPlane3);
        all_seeds.insert(all_seeds.end(), seeds.begin(), seeds.end());
    }
    
    m_log->debug("Found {} raw seeds before merging", all_seeds.size());
    
    // Merge seeds that share hits
    auto merged_seeds = mergeSeeds(all_seeds);
    
    m_log->debug("After merging: {} seeds", merged_seeds.size());
    
    // Convert to track parameters
    return makeTrackParams(merged_seeds);
}

std::map<int, std::vector<edm4eic::TrackerHit>> 
SoLIDTrackSeeding::groupHitsByPlane(const edm4eic::TrackerHitCollection& hits) {
    std::map<int, std::vector<edm4eic::TrackerHit>> grouped;
    
    for (const auto& hit : hits) {
        int layer = getLayer(hit);
        grouped[layer].push_back(hit);
    }
    
    return grouped;
}

std::vector<SoLIDTrackSeeding::Seed> SoLIDTrackSeeding::findDoubleSeeds(
    const std::vector<edm4eic::TrackerHit>& plane1_hits,
    const std::vector<edm4eic::TrackerHit>& plane2_hits,
    int plane1_idx, int plane2_idx) {
    
    std::vector<Seed> seeds;
    
    for (const auto& hit1 : plane1_hits) {
        const auto pos1 = hit1.getPosition();
        float r1 = std::hypot(pos1.x, pos1.y);
        float phi1 = std::atan2(pos1.y, pos1.x);
        
        for (const auto& hit2 : plane2_hits) {
            const auto pos2 = hit2.getPosition();
            float r2 = std::hypot(pos2.x, pos2.y);
            float phi2 = std::atan2(pos2.y, pos2.x);
            
            // Apply cuts on deltaR and deltaPhi
            float deltaR = std::abs(r2 - r1);
            float deltaPhi = std::remainder(phi2 - phi1, 2.0 * M_PI);
            
            m_log->debug("Hit pair deltaR: {}, deltaPhi: {}", deltaR, deltaPhi);
            if (deltaR > m_cfg.maxDeltaR) continue;
            if (std::abs(deltaPhi) > m_cfg.maxDeltaPhi) continue;
            
            // Calculate initial parameters
            float momentum, theta, phi;
            if (!calcInitParForPair(hit1, hit2, momentum, theta, phi)) continue;
            
            // Create seed
            Seed seed;
            seed.hits.push_back(hit1);
            seed.hits.push_back(hit2);
            seed.momentum = momentum;
            seed.theta = theta;
            seed.phi = phi;
            seed.position = pos1; // Use first hit position
            
            // Set plane combination bit mask
            seed.planeCombo = 0;
            if (plane1_idx == m_cfg.seedPlane1) seed.planeCombo |= (1 << 0);
            if (plane1_idx == m_cfg.seedPlane2) seed.planeCombo |= (1 << 1);
            if (plane1_idx == m_cfg.seedPlane3) seed.planeCombo |= (1 << 2);
            if (plane2_idx == m_cfg.seedPlane1) seed.planeCombo |= (1 << 0);
            if (plane2_idx == m_cfg.seedPlane2) seed.planeCombo |= (1 << 1);
            if (plane2_idx == m_cfg.seedPlane3) seed.planeCombo |= (1 << 2);
            
            // Check if track propagates to target region
            if (!propagateToTarget(seed)) continue;
            
            // Check if track matches ECal (if you have ECal hits available)
            // if (!matchECal(seed)) continue;
            
            seeds.push_back(seed);
        }
    }
    
    return seeds;
}

bool SoLIDTrackSeeding::calcInitParForPair(
    const edm4eic::TrackerHit& hit1,
    const edm4eic::TrackerHit& hit2,
    float& momentum,
    float& theta,
    float& phi) {
    
    const auto pos1 = hit1.getPosition();
    const auto pos2 = hit2.getPosition();

    float deltaX = pos2.x - pos1.x;
    float deltaY = pos2.y - pos1.y;
    float deltaZ = pos2.z - pos1.z;
    
    float deltaR = std::hypot(pos2.x, pos2.y) - std::hypot(pos1.x, pos1.y);

    float phi1 = std::atan2(pos1.y, pos1.x);
    float phi2 = std::atan2(pos2.y, pos2.x);
    float deltaPhi = std::remainder(phi2 - phi1, 2.0 * M_PI);
    deltaPhi = std::abs(deltaPhi);
    
    // Estimate momentum from deltaPhi
    momentum = 0;
    
    // Estimate theta from deltaPhi
    theta = atan(deltaR/deltaZ);
    
    // Estimate phi from deltaPhi
    phi = atan2(deltaY, deltaX);

    // Calculate charge (assume electron for DIS)
    int charge = -1;

    if (m_cfg.detConf == 3){ //SIDIS NH3
        if (m_cfg.type == SoLIDTrackSeedingConfig::kFAEC && getLayer(hit1) == 5 && getLayer(hit2) == 6){
            if (deltaPhi < -3.54398e-04 + 1e-6) return false;
            momentum = 1.82273e-01/(deltaPhi - -3.54398e-04);
            theta += -3.15393e-05 + -0.00194545*deltaPhi + 1.92086*deltaPhi*deltaPhi;
            phi += 0.000503764 + 0.843494*deltaPhi + 0.613864*deltaPhi*deltaPhi;
        }
        else if (m_cfg.type == SoLIDTrackSeedingConfig::kFAEC && getLayer(hit1) == 4 && getLayer(hit2) == 5){
            if (deltaPhi < -1.44365e-04 + 1e-6) return false;
            momentum = 1.52902e-01/(deltaPhi - -1.44365e-04);
            theta += -6.27864e-05 + 0.00150745*deltaPhi + 1.87546*deltaPhi*deltaPhi;
            phi +=  0.00283427+ 0.925668*deltaPhi + 1.29002*deltaPhi*deltaPhi;
        }
        else if (m_cfg.type == SoLIDTrackSeedingConfig::kFAEC && getLayer(hit1) == 4 && getLayer(hit2) == 6){
            if (deltaPhi < -2.67014e-04 + 1e-6) return false;
            momentum = 3.40491e-01/(deltaPhi - -2.67014e-04);
            theta +=  -0.000104712+ 0.00152549*deltaPhi + 0.451974*deltaPhi*deltaPhi;
            phi += 0.00410935 + 0.832924*deltaPhi + 0.607294*deltaPhi*deltaPhi;
        }
        else if (m_cfg.type == SoLIDTrackSeedingConfig::kLAEC && getLayer(hit1) == 3 && getLayer(hit2) == 4){
            if (deltaPhi < -1.18851e-03 + 1e-6) return false;
            momentum = 1.12414e-01/(deltaPhi - -1.18851e-03);
            theta +=  0.000239174+ -0.0296545*deltaPhi + 5.32953*deltaPhi*deltaPhi;
            phi += 0.00472525+ 0.88319*deltaPhi + 5.14816*deltaPhi*deltaPhi;
        }
        else if (m_cfg.type == SoLIDTrackSeedingConfig::kLAEC && getLayer(hit1) == 3 && getLayer(hit2) == 4){
            if (deltaPhi < -9.32190e-04 + 1e-6) return false;
            momentum = 6.42127e-02/(deltaPhi - -9.32190e-04);
            theta += 4.88206e-05 + -0.0184821*deltaPhi + 9.93942*deltaPhi*deltaPhi;
            phi += 0.00560415+ 0.608656*deltaPhi + 20.6661*deltaPhi*deltaPhi;
        }
        else if (m_cfg.type == SoLIDTrackSeedingConfig::kLAEC && getLayer(hit1) == 2 && getLayer(hit2) == 4){
            if (deltaPhi < -2.10386e-03 + 1e-6) return false;
            momentum = 1.76726e-01/(deltaPhi - -2.10386e-03);
            theta += 2.66966e-05 + -0.00712552*deltaPhi + 1.88514*deltaPhi*deltaPhi;
            phi +=0.0103316 + 0.798952*deltaPhi + 4.67825*deltaPhi*deltaPhi;
        }
    }
    else{ //SIDIS He3 and JPsi
        if (m_cfg.type == SoLIDTrackSeedingConfig::kFAEC && getLayer(hit1) == 5 && getLayer(hit2) == 6){
            if (deltaPhi < -0.000260251 + 1e-6) return false;
            momentum = 0.184566/(deltaPhi - -0.000260251);
            theta += 0.00117211 + -0.0519514*deltaPhi + 2.25303*deltaPhi*deltaPhi;
            phi +=(-1.*charge)*( 0.000617591 + 0.827572*deltaPhi + 0.53233*deltaPhi*deltaPhi );
        }
        else if (m_cfg.type == SoLIDTrackSeedingConfig::kFAEC && getLayer(hit1) == 4 && getLayer(hit2) == 5){
            if (deltaPhi < -0.000456619 + 1e-6) return false;
            momentum = 0.156196/(deltaPhi - -0.000456619);
            theta += 0.000614522 + -0.0317596*deltaPhi + 2.14209*deltaPhi*deltaPhi;
            phi += (-1*charge)*(-0.000117883 + 1.0794*deltaPhi + -0.178718*deltaPhi*deltaPhi);
        }
        else if (m_cfg.type == SoLIDTrackSeedingConfig::kFAEC && getLayer(hit1) == 4 && getLayer(hit2) == 6){
            if (deltaPhi < -0.00068094 + 1e-6) return false;
            momentum = 0.340493/(deltaPhi - -0.00068094);
            theta += 0.000992593 + -0.0239743*deltaPhi + 0.549664*deltaPhi*deltaPhi;
            phi += (-1*charge)*(0.000672179 + 0.908331*deltaPhi + 0.166936*deltaPhi*deltaPhi);
        }
        else if (m_cfg.type == SoLIDTrackSeedingConfig::kLAEC && getLayer(hit1) == 3 && getLayer(hit2) == 4){
            if (deltaPhi < -0.000396581 + 1e-6) return false;
            momentum = 0.109309/(deltaPhi - -0.000396581);
            theta += 3.61953e-05 + -1.06941e-05*deltaPhi + 4.23056*deltaPhi*deltaPhi;
            phi += (-1*charge)*(-0.000283361 + 1.2515*deltaPhi + -0.637376*deltaPhi*deltaPhi);
        }
        else if (m_cfg.type == SoLIDTrackSeedingConfig::kLAEC && getLayer(hit1) == 2 && getLayer(hit2) == 3){
            if (deltaPhi < -0.000371282 + 1e-6) return false;
            momentum = 0.0610608/(deltaPhi - -0.000371282);
            theta += 5.30086e-05 + -0.00477822*deltaPhi + 8.29537*deltaPhi*deltaPhi;
            phi += (-1*charge)*(-0.00018356 + 1.39156*deltaPhi + -1.29865*deltaPhi*deltaPhi);
        }
        else if (m_cfg.type == SoLIDTrackSeedingConfig::kLAEC && getLayer(hit1) == 2 && getLayer(hit2) == 4){
            if (deltaPhi < -0.000743132 + 1e-6) return false;
            momentum = 0.17019/(deltaPhi - -0.000743132);
            theta += -5.87987e-06 + 0.00195822*deltaPhi + 1.58111*deltaPhi*deltaPhi;
            phi += (-1*charge)*(-0.000508897 + 1.30589*deltaPhi + -0.468172*deltaPhi*deltaPhi);
        }
    }
    
    phi = std::remainder(phi, 2.0 * M_PI);
    m_log->debug("Initial parameters: p = {} GeV, theta = {}, phi = {}", 
                 momentum / Acts::UnitConstants::GeV, theta, phi);

    return true;
}

int SoLIDTrackSeeding::getLayer(const edm4eic::TrackerHit& hit) {
    dd4hep::BitFieldCoder coder("system:8,layer:5,module:16,sensor:3,x:32:-16,z:-16");
    int layer = coder.get(hit.getCellID(), "layer");
    return layer;
}

bool SoLIDTrackSeeding::propagateToTarget(Seed& seed) {
    // Create target surface at target position
    Acts::Transform3 transform = Acts::Transform3::Identity();
    transform.translation() = Acts::Vector3(m_cfg.targetX, m_cfg.targetY, m_cfg.targetZ);

    auto targetSurface =
      Acts::Surface::makeShared<Acts::PlaneSurface>(
        transform,
        std::make_shared<Acts::RectangleBounds>(1.0 * Acts::UnitConstants::m, 1.0 * Acts::UnitConstants::m)
      );

    // Create initial track parameters from seed
    Acts::Vector3 position(seed.position.x, seed.position.y, seed.position.z);
    Acts::Vector3 direction(
        std::sin(seed.theta) * std::cos(seed.phi),
        std::sin(seed.theta) * std::sin(seed.phi),
        std::cos(seed.theta)
    );
    
    // Create a reference surface at the seed position (plane perpendicular to z-axis)
    auto seedSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
        Acts::Vector3(seed.position.x, seed.position.y, seed.position.z),
        Acts::Vector3(0, 0, 1)); // normal along z
    
    // Convert global position to local coordinates on seed surface
    auto localResult = seedSurface->globalToLocal(
        m_geoSvc->getActsGeometryContext(),
        position, direction);
    
    if (!localResult.ok()) {
        m_log->trace("Failed to convert seed position to local coordinates");
        return false;
    }
    
    Acts::Vector2 localPos = localResult.value();
    
    // Assume electron charge for DIS
    int charge = -1;
    
    // Create bound track parameters at seed position
    Acts::BoundVector boundParams;
    boundParams << localPos(0), localPos(1), 
                   seed.phi, seed.theta, 
                   charge / seed.momentum, 
                   0.0; // time
    
    // Create covariance matrix (use large errors for seed)
    Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Identity();
    cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = 100.0 * Acts::UnitConstants::mm2;
    cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = 100.0 * Acts::UnitConstants::mm2;
    cov(Acts::eBoundPhi, Acts::eBoundPhi) = 0.1; // ~18 degrees
    cov(Acts::eBoundTheta, Acts::eBoundTheta) = 0.1;
    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = 0.01 / (Acts::UnitConstants::GeV * Acts::UnitConstants::GeV);
    cov(Acts::eBoundTime, Acts::eBoundTime) = 1.0;
    
    Acts::BoundTrackParameters initialParams(
        seedSurface,
        boundParams,
        cov,
        Acts::ParticleHypothesis::electron());
    
    // Set up the propagator
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField = m_geoSvc->getFieldProvider();
    using Stepper = Acts::EigenStepper<>;
    using Propagator = Acts::Propagator<Stepper>;
    
    Stepper stepper(magneticField);
    Propagator propagator(stepper);
    
    // Set up propagation options
    Acts::PropagatorOptions<> options(
        m_geoSvc->getActsGeometryContext(),
        m_fieldctx);
    options.direction = Acts::Direction::Backward;
    options.maxStepSize = m_cfg.propagationStepSize;
    options.maxSteps = 10000;
    
    // Propagate to target surface
    auto result = propagator.propagate(initialParams, *targetSurface, options);
    
    // Check if propagation succeeded
    if (!result.ok()) {
        m_log->trace("Propagation to target failed: {}", result.error().message());
        return false;
    }
    
    // Get the propagated parameters
    const auto& endParams = *result->endParameters;
    
    // Check if track reaches target within reasonable distance
    // Get position at target
    Acts::Vector3 targetPos = endParams.position(m_geoSvc->getActsGeometryContext());
    
    // Check distance from nominal target position
    Acts::Vector3 targetNominal(m_cfg.targetX, m_cfg.targetY, m_cfg.targetZ);
    float distance = (targetPos - targetNominal).norm();
    
    // Accept if within reasonable range (e.g., 200 mm from target)
    const float maxTargetDist = 200.0 * Acts::UnitConstants::mm;
    
    if (distance > maxTargetDist) {
        m_log->trace("Track misses target by {} mm", distance / Acts::UnitConstants::mm);
        return false;
    }
    
    m_log->trace("Track propagated to target successfully, distance = {} mm", 
                 distance / Acts::UnitConstants::mm);
    
    // Update seed with propagated parameters
    seed.position.x = targetPos.x();
    seed.position.y = targetPos.y();
    seed.position.z = targetPos.z();
    
    // Get momentum and direction
    auto propagatedMomentum = endParams.absoluteMomentum();
    auto propagatedDirection = endParams.direction();
    
    // Extract theta and phi from direction vector
    seed.theta = std::acos(propagatedDirection.z());
    seed.phi = std::atan2(propagatedDirection.y(), propagatedDirection.x());
    seed.momentum = propagatedMomentum;
    
    return true;
}

bool SoLIDTrackSeeding::matchECal(const Seed& seed) {
    // Placeholder for ECal matching
    // You would propagate to ECal surface and check distance to ECal hits
    return true;
}

std::vector<SoLIDTrackSeeding::Seed> 
SoLIDTrackSeeding::mergeSeeds(const std::vector<Seed>& seeds) {
    
    std::vector<Seed> merged;
    std::vector<bool> used(seeds.size(), false);
    
    for (size_t i = 0; i < seeds.size(); ++i) {
        if (used[i]) continue;
        
        Seed merged_seed = seeds[i];
        used[i] = true;
        
        // Look for seeds to merge
        for (size_t j = i + 1; j < seeds.size(); ++j) {
            if (used[j]) continue;
            
            // Check if seeds share a common hit on middle plane
            bool share_hit = false;
            for (const auto& hit1 : seeds[i].hits) {
                for (const auto& hit2 : seeds[j].hits) {
                    if (hit1.getCellID() == hit2.getCellID()) {
                        share_hit = true;
                        break;
                    }
                }
                if (share_hit) break;
            }
            
            if (share_hit) {
                // Merge the seeds
                for (const auto& hit : seeds[j].hits) {
                    bool already_has = false;
                    for (const auto& existing_hit : merged_seed.hits) {
                        if (hit.getCellID() == existing_hit.getCellID()) {
                            already_has = true;
                            break;
                        }
                    }
                    if (!already_has) {
                        merged_seed.hits.push_back(hit);
                    }
                }
                merged_seed.planeCombo |= seeds[j].planeCombo;
                used[j] = true;
            }
        }
        
        merged.push_back(merged_seed);
    }
    
    return merged;
}

std::unique_ptr<edm4eic::TrackParametersCollection> 
SoLIDTrackSeeding::makeTrackParams(const std::vector<Seed>& seeds) {
    
    auto trackparams = std::make_unique<edm4eic::TrackParametersCollection>();
    
    for (const auto& seed : seeds) {
        // Create target surface at target position
        Acts::Transform3 transform = Acts::Transform3::Identity();
        transform.translation() = Acts::Vector3(m_cfg.targetX, m_cfg.targetY, m_cfg.targetZ);

        auto targetSurface =
          Acts::Surface::makeShared<Acts::PlaneSurface>(
            transform,
            std::make_shared<Acts::RectangleBounds>(1.0 * Acts::UnitConstants::m, 1.0 * Acts::UnitConstants::m)
          );
        
        // Calculate charge (assume electron for DIS)
        int charge = -1;
        
        // Convert to Acts parameters
        Acts::Vector3 global(seed.position.x, seed.position.y, seed.position.z);
        Acts::Vector3 direction(
            std::sin(seed.theta) * std::cos(seed.phi),
            std::sin(seed.theta) * std::sin(seed.phi),
            std::cos(seed.theta)
        );
        
        auto local = targetSurface->globalToLocal(
            m_geoSvc->getActsGeometryContext(),
            global, direction);
        
        if (!local.ok()) continue;
        
        Acts::Vector2 localpos = local.value();
        localpos = 0.0 * localpos; // Set to zero
        
        auto trackparam = trackparams->create();
        trackparam.setType(-1); // seed type
        trackparam.setLoc({static_cast<float>(localpos(0)), 
                          static_cast<float>(localpos(1))});
        trackparam.setPhi(seed.phi);
        trackparam.setTheta(seed.theta);
        trackparam.setQOverP(charge / seed.momentum);
        trackparam.setTime(0);
        
        // Set covariance
        edm4eic::Cov6f cov;
        cov(0,0) = m_cfg.locaError / Acts::UnitConstants::mm; // loc0
        cov(1,1) = m_cfg.locbError / Acts::UnitConstants::mm; // loc1
        cov(2,2) = m_cfg.phiError / Acts::UnitConstants::rad; // phi
        cov(3,3) = m_cfg.thetaError / Acts::UnitConstants::rad; // theta
        cov(4,4) = m_cfg.qOverPError * Acts::UnitConstants::GeV; // qOverP
        cov(5,5) = m_cfg.timeError / Acts::UnitConstants::ns; // time
        trackparam.setCovariance(cov);

        m_log->trace("Created track parameter: p = {} GeV, theta = {}, phi = {}", 
                     seed.momentum / Acts::UnitConstants::GeV, seed.theta, seed.phi);
    }
    
    return trackparams;
}

} // namespace eicrecon