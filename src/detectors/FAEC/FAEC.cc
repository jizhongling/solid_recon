// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2021 - 2024, Chao Peng, Sylvester Joosten, Whitney Armstrong, David Lawrence, Friederike Bock, Wouter Deconinck, Kolja Kauder, Sebouh Paul

#include <DD4hep/Detector.h>
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
#include "services/geometry/dd4hep/DD4hep_service.h"

extern "C" {
void InitPlugin(JApplication* app) {

  using namespace eicrecon;

  InitJANAPlugin(app);
  // Make sure digi and reco use the same value
  decltype(CalorimeterHitDigiConfig::capADC) FAEC_Sh_capADC =
      1048576; //1048576, assuming 20 bits. For approximate HGCROC resolution use 65536
  decltype(CalorimeterHitDigiConfig::dyRangeADC) FAEC_Sh_dyRangeADC   = 100 * dd4hep::GeV;
  decltype(CalorimeterHitDigiConfig::pedMeanADC) FAEC_Sh_pedMeanADC   = 200;
  decltype(CalorimeterHitDigiConfig::pedSigmaADC) FAEC_Sh_pedSigmaADC = 2.4576;
  decltype(CalorimeterHitDigiConfig::resolutionTDC) FAEC_Sh_resolutionTDC =
      10 * dd4hep::picosecond;
  auto detector = app->GetService<DD4hep_service>()->detector();
  double rmod = detector->constant<double>("FAEC_rmod");
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "FAEC_ShRawHits", {"FAEC_ShHits"},
#if EDM4EIC_VERSION_MAJOR >= 7
      {"FAEC_ShRawHits", "FAEC_ShRawHitAssociations"},
#else
      {"FAEC_ShRawHits"},
#endif
      {
          .eRes      = {0.1 * sqrt(dd4hep::GeV), 0.01,
                        0.0 * dd4hep::GeV}, // (10% / sqrt(E)) \oplus 1%
          .tRes      = 0.0,
          .threshold = 0.0,
          // .threshold = 15 * dd4hep::MeV for a single tower, applied on ADC level
          .capADC        = FAEC_Sh_capADC,
          .capTime       = 100, // given in ns, 4 samples in HGCROC
          .dyRangeADC    = FAEC_Sh_dyRangeADC,
          .pedMeanADC    = FAEC_Sh_pedMeanADC,
          .pedSigmaADC   = FAEC_Sh_pedSigmaADC,
          .resolutionTDC = FAEC_Sh_resolutionTDC,
          .corrMeanScale = "1.0",
          .readout       = "FAEC_ShHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "FAEC_ShRecHits", {"FAEC_ShRawHits"}, {"FAEC_ShRecHits"},
      {
          .capADC          = FAEC_Sh_capADC,
          .dyRangeADC      = FAEC_Sh_dyRangeADC,
          .pedMeanADC      = FAEC_Sh_pedMeanADC,
          .pedSigmaADC     = FAEC_Sh_pedSigmaADC,
          .resolutionTDC   = FAEC_Sh_resolutionTDC,
          .thresholdFactor = 0.0,
          .thresholdValue =
              157, // The ADC of a 15 MeV particle is adc = 200 + 15 * ( 1.0 + 0) / 100000 * 1048576 = 200 + 157.2864
          .sampFrac = "1.0",
          .readout  = "FAEC_ShHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "FAEC_ShTruthProtoClusters", {"FAEC_ShRecHits", "FAEC_ShHits"},
      {"FAEC_ShTruthProtoClusters"},
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "FAEC_ShIslandProtoClusters", {"FAEC_ShRecHits"}, {"FAEC_ShIslandProtoClusters"},
      {
          .sectorDist                    = rmod * 2 * 1.5,
          .dimScaledLocalDistXY          = {1.5, 1.5},
          .splitCluster                  = false,
          .minClusterHitEdep             = 0.0 * dd4hep::MeV,
          .minClusterCenterEdep          = 20.0 * dd4hep::MeV,
          .transverseEnergyProfileMetric = "dimScaledLocalDistXY",
          .transverseEnergyProfileScale  = 1.,
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "FAEC_ShTruthClustersWithoutShapes",
      {
        "FAEC_ShTruthProtoClusters", // edm4eic::ProtoClusterCollection
#if EDM4EIC_VERSION_MAJOR >= 7
            "FAEC_ShRawHitAssociations"
      }, // edm4eic::MCRecoCalorimeterHitAssociationCollection
#else
            "FAEC_ShHits"
      }, // edm4hep::SimCalorimeterHitCollection
#endif
      {"FAEC_ShTruthClustersWithoutShapes",             // edm4eic::Cluster
       "FAEC_ShTruthClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = true},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "FAEC_ShTruthClusters",
      {"FAEC_ShTruthClustersWithoutShapes", "FAEC_ShTruthClusterAssociationsWithoutShapes"},
      {"FAEC_ShTruthClusters", "FAEC_ShTruthClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "FAEC_ShClustersWithoutShapes",
      {
        "FAEC_ShIslandProtoClusters", // edm4eic::ProtoClusterCollection
#if EDM4EIC_VERSION_MAJOR >= 7
            "FAEC_ShRawHitAssociations"
      }, // edm4eic::MCRecoCalorimeterHitAssociationCollection
#else
            "FAEC_ShHits"
      }, // edm4hep::SimCalorimeterHitCollection
#endif
      {"FAEC_ShClustersWithoutShapes",             // edm4eic::Cluster
       "FAEC_ShClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {
          .energyWeight    = "log",
          .sampFrac        = 1.0,
          .logWeightBase   = 3.6,
          .enableEtaBounds = false,
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "FAEC_ShClusters",
      {"FAEC_ShClustersWithoutShapes", "FAEC_ShClusterAssociationsWithoutShapes"},
      {"FAEC_ShClusters", "FAEC_ShClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));

  app->Add(new JOmniFactoryGeneratorT<TrackClusterMergeSplitter_factory>(
      "FAEC_ShSplitMergeProtoClusters",
      {"FAEC_ShIslandProtoClusters", "CalorimeterTrackProjections"},
      {"FAEC_ShSplitMergeProtoClusters"},
      {.idCalo                       = "FAEC_Sh_ID",
       .minSigCut                    = -2.0,
       .avgEP                        = 1.0,
       .sigEP                        = 0.10,
       .drAdd                        = 0.30,
       .sampFrac                     = 1.0,
       .transverseEnergyProfileScale = 1.0},
      app // TODO: remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "FAEC_ShSplitMergeClustersWithoutShapes",
      {
        "FAEC_ShSplitMergeProtoClusters", // edm4eic::ProtoClusterCollection
#if EDM4EIC_VERSION_MAJOR >= 7
            "FAEC_ShRawHitAssociations"
      }, // edm4hep::MCRecoCalorimeterHitAssociationCollection
#else
            "FAEC_ShHits"
      }, // edm4hep::SimCalorimeterHitCollection
#endif
      {"FAEC_ShSplitMergeClustersWithoutShapes",             // edm4eic::Cluster
       "FAEC_ShSplitMergeClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 3.6, .enableEtaBounds = false},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "FAEC_ShSplitMergeClusters",
      {"FAEC_ShSplitMergeClustersWithoutShapes",
       "FAEC_ShSplitMergeClusterAssociationsWithoutShapes"},
      {"FAEC_ShSplitMergeClusters", "FAEC_ShSplitMergeClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));

  // Preshower is identical to regular Ecal
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "FAEC_PrShRawHits", {"FAEC_PrShHits"},
#if EDM4EIC_VERSION_MAJOR >= 7
      {"FAEC_PrShRawHits", "FAEC_PrShRawHitAssociations"},
#else
      {"FAEC_PrShRawHits"},
#endif
      {
          .eRes      = {0.1 * sqrt(dd4hep::GeV), 0.01,
                        0.0 * dd4hep::GeV}, // (10% / sqrt(E)) \oplus 1%
          .tRes      = 0.0,
          .threshold = 0.0,
          // .threshold = 15 * dd4hep::MeV for a single tower, applied on ADC level
          .capADC        = FAEC_Sh_capADC,
          .capTime       = 100, // given in ns, 4 samples in HGCROC
          .dyRangeADC    = FAEC_Sh_dyRangeADC,
          .pedMeanADC    = FAEC_Sh_pedMeanADC,
          .pedSigmaADC   = FAEC_Sh_pedSigmaADC,
          .resolutionTDC = FAEC_Sh_resolutionTDC,
          .corrMeanScale = "1.0",
          .readout       = "FAEC_PrShHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "FAEC_PrShRecHits", {"FAEC_PrShRawHits"}, {"FAEC_PrShRecHits"},
      {
          .capADC          = FAEC_Sh_capADC,
          .dyRangeADC      = FAEC_Sh_dyRangeADC,
          .pedMeanADC      = FAEC_Sh_pedMeanADC,
          .pedSigmaADC     = FAEC_Sh_pedSigmaADC,
          .resolutionTDC   = FAEC_Sh_resolutionTDC,
          .thresholdFactor = 0.0,
          .thresholdValue =
              157, // The ADC of a 15 MeV particle is adc = 200 + 15 * ( 1.0 + 0) / 100000 * 1048576 = 200 + 157.2864
          .sampFrac = "1.0",
          .readout  = "FAEC_PrShHits",
      },
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "FAEC_PrShTruthProtoClusters", {"FAEC_PrShRecHits", "FAEC_PrShHits"},
      {"FAEC_PrShTruthProtoClusters"},
      app // TODO: Remove me once fixed
      ));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "FAEC_PrShIslandProtoClusters", {"FAEC_PrShRecHits"},
      {"FAEC_PrShIslandProtoClusters"},
      {
          .sectorDist                    = rmod * 2 * 1.5,
          .dimScaledLocalDistXY          = {1.5, 1.5},
          .splitCluster                  = false,
          .minClusterHitEdep             = 0.0 * dd4hep::MeV,
          .minClusterCenterEdep          = 20.0 * dd4hep::MeV,
          .transverseEnergyProfileMetric = "dimScaledLocalDistXY",
          .transverseEnergyProfileScale  = 1.,
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "FAEC_PrShTruthClustersWithoutShapes",
      {
        "FAEC_PrShTruthProtoClusters", // edm4eic::ProtoClusterCollection
#if EDM4EIC_VERSION_MAJOR >= 7
            "FAEC_PrShRawHitAssociations"
      }, // edm4eic::MCRecoCalorimeterHitCollection
#else
            "FAEC_PrShHits"
      }, // edm4hep::SimCalorimeterHitCollection
#endif
      {"FAEC_PrShTruthClustersWithoutShapes",             // edm4eic::Cluster
       "FAEC_PrShTruthClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = true},
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "FAEC_PrShTruthClusters",
      {"FAEC_PrShTruthClustersWithoutShapes",
       "FAEC_PrShTruthClusterAssociationsWithoutShapes"},
      {"FAEC_PrShTruthClusters", "FAEC_PrShTruthClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "FAEC_PrShClustersWithoutShapes",
      {
        "FAEC_PrShIslandProtoClusters", // edm4eic::ProtoClusterCollection
#if EDM4EIC_VERSION_MAJOR >= 7
            "FAEC_PrShRawHitAssociations"
      }, // edm4eic::MCRecoCalorimeterHitCollection
#else
            "FAEC_PrShHits"
      }, // edm4hep::SimCalorimeterHitCollection
#endif
      {"FAEC_PrShClustersWithoutShapes",             // edm4eic::Cluster
       "FAEC_PrShClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {
          .energyWeight    = "log",
          .sampFrac        = 1.0,
          .logWeightBase   = 3.6,
          .enableEtaBounds = false,
      },
      app // TODO: Remove me once fixed
      ));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "FAEC_PrShClusters",
      {"FAEC_PrShClustersWithoutShapes",
       "FAEC_PrShClusterAssociationsWithoutShapes"},
      {"FAEC_PrShClusters", "FAEC_PrShClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));
}
}
