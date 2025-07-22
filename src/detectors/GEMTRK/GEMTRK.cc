// Copyright 2025, Zhongling Ji
// Subject to the terms in the LICENSE file found in the top-level directory.
//
//

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplication.h>
#include <string>

#include "algorithms/interfaces/WithPodConfig.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/digi/GEMTrackerDigi_factory.h"
#include "factories/tracking/TrackerHitReconstruction_factory.h"

extern "C" {
void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);

    using namespace eicrecon;

    // Digitization
    app->Add(new JOmniFactoryGeneratorT<GEMTrackerDigi_factory>(
        "GEMTrackerRawHits",
        {
          "GEMTrackerHits"
        },
        {
          "GEMTrackerRawHits",
          "GEMTrackerRawHitAssociations"
        },
        {
            .threshold = 100 * dd4hep::eV,
            .timeResolution = 10,
        },
        app
    ));


    // Convert raw digitized hits into hits with geometry info (ready for tracking)
    app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
        "GEMTrackerRecHits",
        {"GEMTrackerRawHits"},
        {"GEMTrackerRecHits"},
        {
            .timeResolution = 10,
        },
        app
    ));

}
} // extern "C"
