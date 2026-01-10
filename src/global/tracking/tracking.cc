// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2024, Dmitry Romanov, Tyler Kutz, Wouter Deconinck, Dmitry Kalinkin

#include <DD4hep/Detector.h>
#include <JANA/JApplication.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/TrackCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <fmt/core.h>
#include <algorithm>
#include <gsl/pointers>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "ActsToTracks.h"
#include "ActsToTracks_factory.h"
#include "AmbiguitySolver_factory.h"
#include "CKFTracking_factory.h"
#include "IterativeVertexFinder_factory.h"
#include "TrackParamTruthInit_factory.h"
#include "TrackProjector_factory.h"
#include "TrackPropagationConfig.h"
#include "TrackPropagation_factory.h"
#include "TrackSeeding_factory.h"
#include "TrackerMeasurementFromHits_factory.h"
#include "TracksToParticles_factory.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/meta/CollectionCollector_factory.h"
#include "services/geometry/dd4hep/DD4hep_service.h"
#include "SoLIDTrackSeeding_factory.h"

//
extern "C" {
void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);

    using namespace eicrecon;

    app->Add(new JOmniFactoryGeneratorT<TrackParamTruthInit_factory>(
            "CentralTrackTruthSeeds",
            {"MCParticles"},
            {"CentralTrackTruthSeeds"},
            {},
            app
            ));

    // Possible collections from arches, brycecanyon and craterlake configurations
    std::vector<std::tuple<std::string, std::string, std::string, std::string>> possible_collections = {
        {"GEMTrackerHits", "GEMTrackerRawHits", "GEMTrackerRawHitAssociations", "GEMTrackerRecHits"}
    };

    // Filter out collections that are not present in the current configuration
    std::vector<std::string> input_rec_collections;
    std::vector<std::string> input_raw_assoc_collections;
    auto readouts = app->GetService<DD4hep_service>()->detector()->readouts();
    for (const auto& [hit_collection, raw_collection, raw_assoc_collection, rec_collection] : possible_collections) {
      if (readouts.find(hit_collection) != readouts.end()) {
        // Add the collection to the list of input collections
        input_rec_collections.push_back(rec_collection);
        input_raw_assoc_collections.push_back(raw_assoc_collection);
      }
    }

    // Tracker hits collector
    app->Add(new JOmniFactoryGeneratorT<CollectionCollector_factory<edm4eic::TrackerHit>>(
        "CentralTrackingRecHits",
        input_rec_collections,
        {"CentralTrackingRecHits"}, // Output collection name
        app));

    // Tracker hit associations collector
    app->Add(new JOmniFactoryGeneratorT<CollectionCollector_factory<edm4eic::MCRecoTrackerHitAssociation>>(
        "CentralTrackingRawHitAssociations",
        input_raw_assoc_collections,
        {"CentralTrackingRawHitAssociations"}, // Output collection name
        app));

    app->Add(new JOmniFactoryGeneratorT<TrackerMeasurementFromHits_factory>(
            "CentralTrackerMeasurements",
            {"CentralTrackingRecHits"},
            {"CentralTrackerMeasurements"},
            app
            ));

    app->Add(new JOmniFactoryGeneratorT<CKFTracking_factory>(
        "CentralCKFTruthSeededTrajectories",
        {
            "CentralTrackTruthSeeds",
            "CentralTrackerMeasurements"
        },
        {
            "CentralCKFTruthSeededActsTrajectoriesUnfiltered",
            "CentralCKFTruthSeededActsTracksUnfiltered",
        },
        app
    ));

    app->Add(new JOmniFactoryGeneratorT<ActsToTracks_factory>(
        "CentralCKFTruthSeededTracksUnfiltered",
        {
            "CentralTrackerMeasurements",
            "CentralCKFTruthSeededActsTrajectoriesUnfiltered",
            "CentralTrackingRawHitAssociations",
        },
        {
            "CentralCKFTruthSeededTrajectoriesUnfiltered",
            "CentralCKFTruthSeededTrackParametersUnfiltered",
            "CentralCKFTruthSeededTracksUnfiltered",
            "CentralCKFTruthSeededTrackUnfilteredAssociations",
        },
        app
    ));

    app->Add(new JOmniFactoryGeneratorT<AmbiguitySolver_factory>(
        "TruthSeededAmbiguityResolutionSolver",
        {
             "CentralCKFTruthSeededActsTracksUnfiltered",
             "CentralTrackerMeasurements"
        },
        {
             "CentralCKFTruthSeededActsTracks",
             "CentralCKFTruthSeededActsTrajectories",
        },
        app
    ));

    app->Add(new JOmniFactoryGeneratorT<ActsToTracks_factory>(
        "CentralCKFTruthSeededTracks",
        {
            "CentralTrackerMeasurements",
            "CentralCKFTruthSeededActsTrajectories",
            "CentralTrackingRawHitAssociations",
        },
        {
            "CentralCKFTruthSeededTrajectories",
            "CentralCKFTruthSeededTrackParameters",
            "CentralCKFTruthSeededTracks",
            "CentralCKFTruthSeededTrackAssociations",
        },
        app
    ));

    // Replace the standard TrackSeeding with SoLIDTrackSeeding for GEM hits
    app->Add(new JOmniFactoryGeneratorT<SoLIDTrackSeeding_factory>(
        "CentralTrackSeedingResults",
        {"CentralTrackingRecHits"},
        {"CentralTrackSeedingResults"},
        {},
        app
    ));

    // Continue with CKF tracking using the SoLID seeds
    app->Add(new JOmniFactoryGeneratorT<CKFTracking_factory>(
        "CentralCKFTrajectories",
        {
            "CentralTrackSeedingResults",
            "CentralTrackerMeasurements"
        },
        {
            "CentralCKFActsTrajectoriesUnfiltered",
            "CentralCKFActsTracksUnfiltered",
        },
        app
    ));

    app->Add(new JOmniFactoryGeneratorT<ActsToTracks_factory>(
        "CentralCKFTracksUnfiltered",
        {
            "CentralTrackerMeasurements",
            "CentralCKFActsTrajectoriesUnfiltered",
            "CentralTrackingRawHitAssociations",
        },
        {
            "CentralCKFTrajectoriesUnfiltered",
            "CentralCKFTrackParametersUnfiltered",
            "CentralCKFTracksUnfiltered",
            "CentralCKFTrackUnfilteredAssociations",
        },
        app
    ));

    app->Add(new JOmniFactoryGeneratorT<AmbiguitySolver_factory>(
        "AmbiguityResolutionSolver",
        {
             "CentralCKFActsTracksUnfiltered",
             "CentralTrackerMeasurements"
        },
        {
             "CentralCKFActsTracks",
             "CentralCKFActsTrajectories",
        },
        app
    ));

    app->Add(new JOmniFactoryGeneratorT<ActsToTracks_factory>(
        "CentralCKFTracks",
        {
            "CentralTrackerMeasurements",
            "CentralCKFActsTrajectories",
            "CentralTrackingRawHitAssociations",
        },
        {
            "CentralCKFTrajectories",
            "CentralCKFTrackParameters",
            "CentralCKFTracks",
            "CentralCKFTrackAssociations",
        },
        app
    ));

    app->Add(new JOmniFactoryGeneratorT<TrackProjector_factory>(
        "CentralTrackSegments",
        {
            "CentralCKFActsTrajectories",
            "CentralCKFTracks",
        },
        {
            "CentralTrackSegments",
        },
        app
    ));

    app->Add(new JOmniFactoryGeneratorT<IterativeVertexFinder_factory>(
            "CentralTrackVertices",
            {"CentralCKFActsTrajectories","ReconstructedChargedParticles"},
            {"CentralTrackVertices"},
            {},
            app
            ));

    app->Add(new JOmniFactoryGeneratorT<TrackPropagation_factory>(
            "CalorimeterTrackPropagator",
            {"CentralCKFTracks", "CentralCKFActsTrajectories", "CentralCKFActsTracks"},
            {"CalorimeterTrackProjections"},
            {
                .target_surfaces{
                    // DiscSurfaceConfig{id, zmin, rmin, rmax}
                    // EICrecon uses id name, but solid_recon uses id number.
                    // Need to modify algorithms/tracking/TrackPropagation.cc
                    // LAECPreShower
                    eicrecon::DiscSurfaceConfig{"3", "-65*cm",         "0.9*83*cm", "1.1*140*cm"},
                    eicrecon::DiscSurfaceConfig{"3", "-65*cm + 50*mm", "0.9*83*cm", "1.1*140*cm"},
                    // LAECShower
                    eicrecon::DiscSurfaceConfig{"4", "-57*cm",         "0.9*83*cm", "1.1*140*cm"},
                    eicrecon::DiscSurfaceConfig{"4", "-57*cm + 50*mm", "0.9*83*cm", "1.1*140*cm"},
                    // FAECPreShower
                    eicrecon::DiscSurfaceConfig{"5", "425*cm",         "0.9*98*cm", "1.1*230*cm"},
                    eicrecon::DiscSurfaceConfig{"5", "425*cm + 50*mm", "0.9*98*cm", "1.1*230*cm"},
                    // FAECShower
                    eicrecon::DiscSurfaceConfig{"6", "433*cm",         "0.9*98*cm", "1.1*230*cm"},
                    eicrecon::DiscSurfaceConfig{"6", "433*cm + 50*mm", "0.9*98*cm", "1.1*230*cm"},
                }
            },
            app
            ));



    std::vector<std::string> input_track_collections;
    //Check size of input_rec_collections to determine if CentralCKFTracks should be added to the input_track_collections
    if (input_rec_collections.size() > 0) {
        input_track_collections.push_back("CentralCKFTracks");
    }
    //Check if the TaggerTracker readout is present in the current configuration
    if (readouts.find("TaggerTrackerHits") != readouts.end()) {
        input_track_collections.push_back("TaggerTrackerTracks");
    }

    // Add central and other tracks
    app->Add(new JOmniFactoryGeneratorT<CollectionCollector_factory<edm4eic::Track>>(
            "CombinedTracks",
            input_track_collections,
            {"CombinedTracks"},
            app
            ));

    app->Add(new JOmniFactoryGeneratorT<TracksToParticles_factory>(
            "ChargedTruthSeededParticlesWithAssociations",
            {
              "CentralCKFTruthSeededTracks",
              "CentralCKFTruthSeededTrackAssociations",
            },
            {"ReconstructedTruthSeededChargedWithoutPIDParticles",
             "ReconstructedTruthSeededChargedWithoutPIDParticleAssociations"
            },
            {},
            app
            ));

    app->Add(new JOmniFactoryGeneratorT<TracksToParticles_factory>(
            "ChargedParticlesWithAssociations",
            {
              "CombinedTracks",
              "CentralCKFTrackAssociations",
            },
            {
              "ReconstructedChargedWithoutPIDParticles",
              "ReconstructedChargedWithoutPIDParticleAssociations"
            },
            {},
            app
            ));
}
} // extern "C"
