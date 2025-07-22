#include <iostream>
#include <string>
#include <vector>

#include <JANA/JApplication.h>
#include <JANA/JVersion.h>

#include "solidrecon_cli.h"

/// The default plugins
/// Add new default plugin names here and the main() will do JApplication::AddPlugin() for you.
std::vector<std::string> SOLIDRECON_DEFAULT_PLUGINS = {
    // clang-format off
    "log",
    "dd4hep",
    "evaluator",
    "acts",
    "algorithms_init",
    "reco",
    "tracking",
    "LAEC",
    "FAEC",
    "GEMTRK",
    "podio",
    // clang-format on
};

int main(int narg, char **argv) {
  std::cout << "========================" << std::endl;
  std::cout << "  SoLIDrecon v0.0.1   " << std::endl;
  std::cout << "    JANA2 v" << JVersion::GetVersion() << std::endl;
  std::cout << "    EICrecon v1.20.0-2 " << std::endl;
  std::cout << "========================" << std::endl;

  std::vector<std::string> default_plugins = SOLIDRECON_DEFAULT_PLUGINS;

  auto options = jana::GetCliOptions(narg, argv, false);

  if (jana::HasPrintOnlyCliOptions(options, default_plugins))
    return -1;

  AddAvailablePluginsToOptionParams(options, default_plugins);

  japp = jana::CreateJApplication(options);

  auto exit_code = jana::Execute(japp, options);

  delete japp;
  return exit_code;
}
