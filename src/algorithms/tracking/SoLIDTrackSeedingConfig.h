#pragma once

#include <Acts/Definitions/Units.hpp>

namespace eicrecon {

struct SoLIDTrackSeedingConfig {
    // Target position
    float targetX = 0.0 * Acts::UnitConstants::mm;
    float targetY = 0.0 * Acts::UnitConstants::mm;
    float targetZ = -3500.0 * Acts::UnitConstants::mm;

    int detConf = 3; // Detector configuration
    enum ECType {kLAEC = 0, kFAEC};
    ECType type = kFAEC; // ECal type
    
    // GEM plane indices for seeding (most downstream planes)
    int seedPlane1 = 4;  // GEM 4
    int seedPlane2 = 5;  // GEM 5
    int seedPlane3 = 6;  // GEM 6
    
    // Cut parameters for seed finding (from your simulation)
    float maxDeltaR = 100.0 * Acts::UnitConstants::mm;
    float maxDeltaPhi = 0.1;  // radians
    
    // ECal matching parameters
    float maxECalMatchDist = 50.0 * Acts::UnitConstants::mm;
    
    // Propagation step size
    float propagationStepSize = 10.0 * Acts::UnitConstants::mm;

    // Seed Covariance Error Matrix
    float locaError   = 1.5 * Acts::UnitConstants::mm;     //Error on Loc a
    float locbError   = 1.5 * Acts::UnitConstants::mm;     //Error on Loc b
    float phiError    = 0.02 * Acts::UnitConstants::rad;     //Error on phi
    float thetaError  = 0.002 * Acts::UnitConstants::rad;  //Error on theta
    float qOverPError = 0.025 / Acts::UnitConstants::GeV; //Error on q over p
    float timeError   = 0.1 * Acts::UnitConstants::mm;      //Error on time
};

} // namespace eicrecon