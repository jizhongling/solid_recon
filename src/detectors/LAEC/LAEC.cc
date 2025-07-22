// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2021 - 2024, Chao Peng, Sylvester Joosten, Whitney Armstrong, David Lawrence, Friederike Bock, Wouter Deconinck, Kolja Kauder, Sebouh Paul

#include <edm4eic/EDM4eicVersion.h>
#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplication.h>
#include <math.h>
#include <string>

#include "algorithms/calorimetry/CalorimeterHitDigiConfig.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/calorimetry/CalorimeterClusterRecoCoG_factory.h"
#include "factories/calorimetry/CalorimeterHitDigi_factory.h"
#include "factories/calorimetry/CalorimeterHitReco_factory.h"
#include "factories/calorimetry/CalorimeterIslandCluster_factory.h"
#include "factories/calorimetry/CalorimeterTruthClustering_factory.h"
#include "factories/calorimetry/CalorimeterClusterShape_factory.h"
#include "factories/calorimetry/TrackClusterMergeSplitter_factory.h"

extern "C" {
void InitPlugin(JApplication* app) {

  using namespace eicrecon;

  InitJANAPlugin(app);
  // Make sure digi and reco use the same value
  decltype(CalorimeterHitDigiConfig::capADC) LAEC_Sh_capADC =
      1048576; //1048576, assuming 20 bits. For approximate HGCROC resolution use 65536
  decltype(CalorimeterHitDigiConfig::dyRangeADC) LAEC_Sh_dyRangeADC   = 100 * dd4hep::GeV;
  decltype(CalorimeterHitDigiConfig::pedMeanADC) LAEC_Sh_pedMeanADC   = 200;
  decltype(CalorimeterHitDigiConfig::pedSigmaADC) LAEC_Sh_pedSigmaADC = 2.4576;
  decltype(CalorimeterHitDigiConfig::resolutionTDC) LAEC_Sh_resolutionTDC =
      10 * dd4hep::picosecond;
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "LAEC_ShRawHits", {"LAEC_ShHits"},
#if EDM4EIC_VERSION_MAJOR >= 7
      {"LAEC_ShRawHits", "LAEC_ShRawHitAssociations"},
#else
      {"LAEC_ShRawHits"},
#endif
      {
          .eRes      = {0.1 * sqrt(dd4hep::GeV), 0.01,
                        0.0 * dd4hep::GeV}, // (10% / sqrt(E)) \oplus 1%
          .tRes      = 0.0,
          .threshold = 0.0,
          // .threshold = 15 * dd4hep::MeV for a single tower, applied on ADC level
          .capADC        = LAEC_Sh_capADC,
          .capTime       = 100, // given in ns, 4 samples in HGCROC
          .dyRangeADC    = LAEC_Sh_dyRangeADC,
          .pedMeanADC    = LAEC_Sh_pedMeanADC,
          .pedSigmaADC   = LAEC_Sh_pedSigmaADC,
          .resolutionTDC = LAEC_Sh_resolutionTDC,
          .corrMeanScale = "1.0",
          .readout       = "LAEC_ShHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "LAEC_ShRecHits", {"LAEC_ShRawHits"}, {"LAEC_ShRecHits"},
      {
          .capADC          = LAEC_Sh_capADC,
          .dyRangeADC      = LAEC_Sh_dyRangeADC,
          .pedMeanADC      = LAEC_Sh_pedMeanADC,
          .pedSigmaADC     = LAEC_Sh_pedSigmaADC,
          .resolutionTDC   = LAEC_Sh_resolutionTDC,
          .thresholdFactor = 0.0,
          .thresholdValue =
              157, // The ADC of a 15 MeV particle is adc = 200 + 15 * ( 1.0 + 0) / 100000 * 1048576 = 200 + 157.2864
          .sampFrac = "1.0",
          .readout  = "LAEC_ShHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "LAEC_ShTruthProtoClusters", {"LAEC_ShRecHits", "LAEC_ShHits"},
      {"LAEC_ShTruthProtoClusters"},
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "LAEC_ShIslandProtoClusters", {"LAEC_ShRecHits"}, {"LAEC_ShIslandProtoClusters"},
      {
          .sectorDist                    = 6.25 * dd4hep::cm,
          .dimScaledLocalDistXY          = {1.5, 1.5},
          .splitCluster                  = false,
          .minClusterHitEdep             = 0.0 * dd4hep::MeV,
          .minClusterCenterEdep          = 60.0 * dd4hep::MeV,
          .transverseEnergyProfileMetric = "dimScaledLocalDistXY",
          .transverseEnergyProfileScale  = 1.,
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "LAEC_ShTruthClustersWithoutShapes",
      {
        "LAEC_ShTruthProtoClusters", // edm4eic::ProtoClusterCollection
#if EDM4EIC_VERSION_MAJOR >= 7
            "LAEC_ShRawHitAssociations"
      }, // edm4eic::MCRecoCalorimeterHitAssociationCollection
#else
            "LAEC_ShHits"
      }, // edm4hep::SimCalorimeterHitCollection
#endif
      {"LAEC_ShTruthClustersWithoutShapes",             // edm4eic::Cluster
       "LAEC_ShTruthClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = true},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "LAEC_ShTruthClusters",
      {"LAEC_ShTruthClustersWithoutShapes", "LAEC_ShTruthClusterAssociationsWithoutShapes"},
      {"LAEC_ShTruthClusters", "LAEC_ShTruthClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "LAEC_ShClustersWithoutShapes",
      {
        "LAEC_ShIslandProtoClusters", // edm4eic::ProtoClusterCollection
#if EDM4EIC_VERSION_MAJOR >= 7
            "LAEC_ShRawHitAssociations"
      }, // edm4eic::MCRecoCalorimeterHitAssociationCollection
#else
            "LAEC_ShHits"
      }, // edm4hep::SimCalorimeterHitCollection
#endif
      {"LAEC_ShClustersWithoutShapes",             // edm4eic::Cluster
       "LAEC_ShClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {
          .energyWeight    = "log",
          .sampFrac        = 1.0,
          .logWeightBase   = 3.6,
          .enableEtaBounds = false,
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "LAEC_ShClusters",
      {"LAEC_ShClustersWithoutShapes", "LAEC_ShClusterAssociationsWithoutShapes"},
      {"LAEC_ShClusters", "LAEC_ShClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));

  app->Add(new JOmniFactoryGeneratorT<TrackClusterMergeSplitter_factory>(
      "LAEC_ShSplitMergeProtoClusters",
      {"LAEC_ShIslandProtoClusters", "CalorimeterTrackProjections"},
      {"LAEC_ShSplitMergeProtoClusters"},
      {.idCalo                       = "LAEC_Sh_ID",
       .minSigCut                    = -2.0,
       .avgEP                        = 1.0,
       .sigEP                        = 0.10,
       .drAdd                        = 0.30,
       .sampFrac                     = 1.0,
       .transverseEnergyProfileScale = 1.0},
      app // TODO: remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "LAEC_ShSplitMergeClustersWithoutShapes",
      {
        "LAEC_ShSplitMergeProtoClusters", // edm4eic::ProtoClusterCollection
#if EDM4EIC_VERSION_MAJOR >= 7
            "LAEC_ShRawHitAssociations"
      }, // edm4hep::MCRecoCalorimeterHitAssociationCollection
#else
            "LAEC_ShHits"
      }, // edm4hep::SimCalorimeterHitCollection
#endif
      {"LAEC_ShSplitMergeClustersWithoutShapes",             // edm4eic::Cluster
       "LAEC_ShSplitMergeClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 3.6, .enableEtaBounds = false},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "LAEC_ShSplitMergeClusters",
      {"LAEC_ShSplitMergeClustersWithoutShapes",
       "LAEC_ShSplitMergeClusterAssociationsWithoutShapes"},
      {"LAEC_ShSplitMergeClusters", "LAEC_ShSplitMergeClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));

  // Preshower is identical to regular Ecal
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "LAEC_PrShRawHits", {"LAEC_PrShHits"},
#if EDM4EIC_VERSION_MAJOR >= 7
      {"LAEC_PrShRawHits", "LAEC_PrShRawHitAssociations"},
#else
      {"LAEC_PrShRawHits"},
#endif
      {
          .eRes      = {0.1 * sqrt(dd4hep::GeV), 0.01,
                        0.0 * dd4hep::GeV}, // (10% / sqrt(E)) \oplus 1%
          .tRes      = 0.0,
          .threshold = 0.0,
          // .threshold = 15 * dd4hep::MeV for a single tower, applied on ADC level
          .capADC        = LAEC_Sh_capADC,
          .capTime       = 100, // given in ns, 4 samples in HGCROC
          .dyRangeADC    = LAEC_Sh_dyRangeADC,
          .pedMeanADC    = LAEC_Sh_pedMeanADC,
          .pedSigmaADC   = LAEC_Sh_pedSigmaADC,
          .resolutionTDC = LAEC_Sh_resolutionTDC,
          .corrMeanScale = "1.0",
          .readout       = "LAEC_PrShHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "LAEC_PrShRecHits", {"LAEC_PrShRawHits"}, {"LAEC_PrShRecHits"},
      {
          .capADC          = LAEC_Sh_capADC,
          .dyRangeADC      = LAEC_Sh_dyRangeADC,
          .pedMeanADC      = LAEC_Sh_pedMeanADC,
          .pedSigmaADC     = LAEC_Sh_pedSigmaADC,
          .resolutionTDC   = LAEC_Sh_resolutionTDC,
          .thresholdFactor = 0.0,
          .thresholdValue =
              157, // The ADC of a 15 MeV particle is adc = 200 + 15 * ( 1.0 + 0) / 100000 * 1048576 = 200 + 157.2864
          .sampFrac = "1.0",
          .readout  = "LAEC_PrShHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "LAEC_PrShTruthProtoClusters", {"LAEC_PrShRecHits", "LAEC_PrShHits"},
      {"LAEC_PrShTruthProtoClusters"},
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "LAEC_PrShIslandProtoClusters", {"LAEC_PrShRecHits"},
      {"LAEC_PrShIslandProtoClusters"},
      {
          .sectorDist                    = 6.25 * dd4hep::cm,
          .dimScaledLocalDistXY          = {1.5, 1.5},
          .splitCluster                  = false,
          .minClusterHitEdep             = 0.0 * dd4hep::MeV,
          .minClusterCenterEdep          = 60.0 * dd4hep::MeV,
          .transverseEnergyProfileMetric = "dimScaledLocalDistXY",
          .transverseEnergyProfileScale  = 1.,
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "LAEC_PrShTruthClustersWithoutShapes",
      {
        "LAEC_PrShTruthProtoClusters", // edm4eic::ProtoClusterCollection
#if EDM4EIC_VERSION_MAJOR >= 7
            "LAEC_PrShRawHitAssociations"
      }, // edm4eic::MCRecoCalorimeterHitCollection
#else
            "LAEC_PrShHits"
      }, // edm4hep::SimCalorimeterHitCollection
#endif
      {"LAEC_PrShTruthClustersWithoutShapes",             // edm4eic::Cluster
       "LAEC_PrShTruthClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = true},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "LAEC_PrShTruthClusters",
      {"LAEC_PrShTruthClustersWithoutShapes",
       "LAEC_PrShTruthClusterAssociationsWithoutShapes"},
      {"LAEC_PrShTruthClusters", "LAEC_PrShTruthClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "LAEC_PrShClustersWithoutShapes",
      {
        "LAEC_PrShIslandProtoClusters", // edm4eic::ProtoClusterCollection
#if EDM4EIC_VERSION_MAJOR >= 7
            "LAEC_PrShRawHitAssociations"
      }, // edm4eic::MCRecoCalorimeterHitCollection
#else
            "LAEC_PrShHits"
      }, // edm4hep::SimCalorimeterHitCollection
#endif
      {"LAEC_PrShClustersWithoutShapes",             // edm4eic::Cluster
       "LAEC_PrShClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {
          .energyWeight    = "log",
          .sampFrac        = 1.0,
          .logWeightBase   = 3.6,
          .enableEtaBounds = false,
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "LAEC_PrShClusters",
      {"LAEC_PrShClustersWithoutShapes",
       "LAEC_PrShClusterAssociationsWithoutShapes"},
      {"LAEC_PrShClusters", "LAEC_PrShClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));
}
}
