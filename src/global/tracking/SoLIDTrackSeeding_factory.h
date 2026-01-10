#pragma once

#include <JANA/JEvent.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "algorithms/tracking/SoLIDTrackSeeding.h"
#include "algorithms/tracking/SoLIDTrackSeedingConfig.h"
#include "extensions/jana/JOmniFactory.h"
#include "services/algorithms_init/AlgorithmsInit_service.h"
#include "services/geometry/acts/ACTSGeo_service.h"

namespace eicrecon {

class SoLIDTrackSeeding_factory :
    public JOmniFactory<SoLIDTrackSeeding_factory, SoLIDTrackSeedingConfig> {

private:
    using AlgoT = eicrecon::SoLIDTrackSeeding;
    std::unique_ptr<AlgoT> m_algo;

    PodioInput<edm4eic::TrackerHit> m_hits_input {this};
    PodioOutput<edm4eic::TrackParameters> m_parameters_output {this};

    ParameterRef<float> m_targetZ {this, "targetZ", config().targetZ};
    ParameterRef<int> m_detConf {this, "detConf", config().detConf};
    ParameterRef<float> m_maxDeltaR {this, "maxDeltaR", config().maxDeltaR};
    ParameterRef<float> m_maxDeltaPhi {this, "maxDeltaPhi", config().maxDeltaPhi};
    
    Service<ACTSGeo_service> m_ACTSGeoSvc {this};

public:
    void Configure() {
        m_algo = std::make_unique<AlgoT>();
        m_algo->applyConfig(config());
        m_algo->init(m_ACTSGeoSvc().actsGeoProvider(), logger());
    }

    void ChangeRun(int64_t run_number) {
    }

    void Process(int64_t run_number, uint64_t event_number) {
        m_parameters_output() = m_algo->produce(*m_hits_input());
    }
};

} // namespace eicrecon