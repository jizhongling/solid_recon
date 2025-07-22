// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2025 Zhongling Ji

#pragma once

#include <TRandomGen.h>
#include <algorithms/algorithm.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <functional>
#include <string>
#include <string_view>

#include "GEMTrackerDigiConfig.h"
#include "algorithms/interfaces/WithPodConfig.h"

namespace eicrecon {

using GEMTrackerDigiAlgorithm =
    algorithms::Algorithm<algorithms::Input<edm4hep::SimTrackerHitCollection>,
                          algorithms::Output<edm4eic::RawTrackerHitCollection,
                                             edm4eic::MCRecoTrackerHitAssociationCollection>>;

class GEMTrackerDigi : public GEMTrackerDigiAlgorithm,
                           public WithPodConfig<GEMTrackerDigiConfig> {

public:
  GEMTrackerDigi(std::string_view name)
      : GEMTrackerDigiAlgorithm{name,
                                    {"inputHitCollection"},
                                    {"outputRawHitCollection", "outputHitAssociations"},
                                    "Apply threshold, digitize within ADC range, "
                                    "convert time with smearing resolution."} {}

  void init() final;
  void process(const Input&, const Output&) const final;

private:
  /** Random number generation*/
  TRandomMixMax m_random;
  std::function<double()> m_gauss;

  // FIXME replace with standard random engine
  // std::default_random_engine generator; // TODO: need something more appropriate here
  // std::normal_distribution<double> m_normDist; // defaults to mean=0, sigma=1

  // algorithms::Generator m_rng = algorithms::RandomSvc::instance().generator();
};

} // namespace eicrecon
