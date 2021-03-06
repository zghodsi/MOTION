// MIT License
//
// Copyright (c) 2019 Oleksandr Tkachenko
// Cryptography and Privacy Engineering Group (ENCRYPTO)
// TU Darmstadt, Germany
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <regex>

#include <fmt/format.h>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include "base/party.h"
#include "common/mytest.h"
#include "communication/communication_layer.h"
#include "communication/tcp_transport.h"
#include "statistics/analysis.h"
#include "utility/typedefs.h"

#include <sys/resource.h>

namespace program_options = boost::program_options;

bool CheckPartyArgumentSyntax(const std::string& party_argument);

std::pair<program_options::variables_map, bool> ParseProgramOptions(int ac, char* av[]);

encrypto::motion::PartyPointer CreateParty(const program_options::variables_map& user_options);

constexpr std::size_t kIllegalProtocol{100}, kIllegalOperationType{100};

class memMeter{	
	struct rusage usage;
	double memory_usage;
public: 
	/*https://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-runtime-using-c
	https://elinux.org/Runtime_Memory_Measurement*/
	void process_mem_usage(double& resident_set, double& resident_set_max){		
		getrusage(RUSAGE_SELF, &usage);
		resident_set_max = (double)usage.ru_maxrss/1024;
		
		using std::ios_base;
		using std::ifstream;
		using std::string;		
		
		long program_size, resident_set_size, shared_pages, text, stack , library , dirty_pages; 
		ifstream stat_stream("/proc/self/statm",ios_base::in);		
		stat_stream >> program_size >> resident_set_size >> shared_pages >> text >> stack  >> library  >> dirty_pages;		
		stat_stream.close();
		
		long page_size = sysconf(_SC_PAGE_SIZE); // in Bytes
		resident_set = (double)(resident_set_size*page_size)/1024/1024;		
	}
	
	void print(std::string message){
		double resident_set, resident_set_max;
		process_mem_usage(resident_set, resident_set_max);
		
		std::cout << message << ": ";
		if(resident_set < 1024) std::cout << resident_set << "MB, " << resident_set_max << "MB" << std::endl;
		else std::cout << resident_set/1024 << "GB, " << resident_set_max/1024 << "GB" << std::endl;
	}
};

struct Combination {
  Combination(std::size_t bit_size, encrypto::motion::MpcProtocol protocol,
              encrypto::motion::IntegerOperationType operation_type, std::size_t number_of_simd)
      : bit_size_(bit_size),
        protocol_(protocol),
        operation_type_(operation_type),
        number_of_simd_(number_of_simd) {}

  std::size_t bit_size_{0};
  encrypto::motion::MpcProtocol protocol_{kIllegalProtocol};
  encrypto::motion::IntegerOperationType operation_type_{kIllegalOperationType};
  std::size_t number_of_simd_{0};
};


std::vector<Combination> GenerateAllCombinations() {
  using T = encrypto::motion::IntegerOperationType;

  const std::array kArithmeticBitSizes = {8, 16, 32, 64};
  const std::array kNumbersOfSimd = {1};
  const std::array kOperationTypes = {T::kAdd}; //{T::kAdd, T::kMul, T::kDiv, T::kEq, T::kGt, T::kSub};

  std::vector<Combination> combinations;

  for (const auto bit_size : kArithmeticBitSizes) {
    for (const auto operation_type : kOperationTypes) {
      for (const auto number_of_simd : kNumbersOfSimd) {
        combinations.emplace_back(bit_size, encrypto::motion::MpcProtocol::kArithmeticGmw,
                                  operation_type, number_of_simd);
        combinations.emplace_back(bit_size, encrypto::motion::MpcProtocol::kBooleanGmw,
                                  operation_type, number_of_simd);
        combinations.emplace_back(bit_size, encrypto::motion::MpcProtocol::kBmr, operation_type,
                                  number_of_simd);
      }
    }
  }
  return combinations;
}



int main(int ac, char* av[]) {
  auto [user_options, help_flag] = ParseProgramOptions(ac, av);
  // if help flag is set - print allowed command line arguments and exit
  if (help_flag) return EXIT_SUCCESS;

  const auto number_of_repetitions{user_options["repetitions"].as<std::size_t>()};

  std::vector<Combination> combinations;

  // TODO: add custom combination instead of generating all of them if needed

  combinations = GenerateAllCombinations();

  for (const auto combination : combinations) {
    encrypto::motion::AccumulatedRunTimeStatistics accumulated_statistics;
    encrypto::motion::AccumulatedCommunicationStatistics accumulated_communication_statistics;
    memMeter M;
    for (std::size_t i = 0; i < number_of_repetitions; ++i) {
      encrypto::motion::PartyPointer party{CreateParty(user_options)};
      // establish communication channels with other parties
      auto statistics = EvaluateProtocol(party, combination.number_of_simd_, combination.bit_size_,
                                         combination.protocol_, combination.operation_type_);
      accumulated_statistics.Add(statistics);
      auto communication_statistics = party->GetCommunicationLayer().GetTransportStatistics();
      accumulated_communication_statistics.Add(communication_statistics);
    }
    M.print("usage");
    std::cout << fmt::format(encrypto::motion::to_string(combination.protocol_),
                             encrypto::motion::to_string(combination.operation_type_),
                             combination.bit_size_, combination.number_of_simd_);
    std::cout << encrypto::motion::PrintStatistics(
        fmt::format("Protocol {} operation {} bit size {} SIMD {}",
                    encrypto::motion::to_string(combination.protocol_),
                    encrypto::motion::to_string(combination.operation_type_), combination.bit_size_,
                    combination.number_of_simd_),
        accumulated_statistics, accumulated_communication_statistics);
    std::cout << "my test ended" << std::endl;
  }
  return EXIT_SUCCESS;
}

const std::regex kPartyArgumentRegex(
    "(\\d+),(\\d{1,3}\\.\\d{1,3}\\.\\d{1,3}\\.\\d{1,3}),(\\d{1,5})");

bool CheckPartyArgumentSyntax(const std::string& party_argument) {
  // other party's id, IP address, and port
  return std::regex_match(party_argument, kPartyArgumentRegex);
}

std::tuple<std::size_t, std::string, std::uint16_t> ParsePartyArgument(
    const std::string& party_argument) {
  std::smatch match;
  std::regex_match(party_argument, match, kPartyArgumentRegex);
  auto id = boost::lexical_cast<std::size_t>(match[1]);
  auto host = match[2];
  auto port = boost::lexical_cast<std::uint16_t>(match[3]);
  return {id, host, port};
}

// <variables map, help flag>
std::pair<program_options::variables_map, bool> ParseProgramOptions(int ac, char* av[]) {
  using namespace std::string_view_literals;
  constexpr std::string_view kConfigFileMessage =
      "configuration file, other arguments will overwrite the parameters read from the configuration file"sv;
  bool print, help;
  boost::program_options::options_description description("Allowed options");
  // clang-format off
  description.add_options()
      ("help,h", program_options::bool_switch(&help)->default_value(false),"produce help message")
      ("disable-logging,l","disable logging to file")
      ("print-configuration,p", program_options::bool_switch(&print)->default_value(false), "print configuration")
      ("configuration-file,f", program_options::value<std::string>(), kConfigFileMessage.data())
      ("my-id", program_options::value<std::size_t>(), "my party id")
      ("other-parties", program_options::value<std::vector<std::string>>()->multitoken(), "(other party id, IP, port, my role), e.g., --other-parties 1,127.0.0.1,7777")
      ("online-after-setup", program_options::value<bool>()->default_value(true), "compute the online phase of the gate evaluations after the setup phase for all of them is completed (true/1 or false/0)")
      ("repetitions", program_options::value<std::size_t>()->default_value(1), "number of repetitions");
  // clang-format on

  program_options::variables_map user_options;

  program_options::store(program_options::parse_command_line(ac, av, description), user_options);
  program_options::notify(user_options);

  // argument help or no arguments (at least a configuration file is expected)
  if (help) {
    std::cout << description << "\n";
    return std::make_pair<program_options::variables_map, bool>({}, true);
  }

  // read configuration file
  if (user_options.count("configuration-file")) {
    std::ifstream ifs(user_options["configuration-file"].as<std::string>().c_str());
    program_options::variables_map user_option_config_file;
    program_options::store(program_options::parse_config_file(ifs, description), user_options);
    program_options::notify(user_options);
  }

  // print parsed parameters
  if (user_options.count("my-id")) {
    if (print) std::cout << "My id " << user_options["my-id"].as<std::size_t>() << std::endl;
  } else
    throw std::runtime_error("My id is not set but required");

  if (user_options.count("other-parties")) {
    const std::vector<std::string> other_parties{
        user_options["other-parties"].as<std::vector<std::string>>()};
    std::string parties("Other parties: ");
    for (auto& party : other_parties) {
      if (CheckPartyArgumentSyntax(party)) {
        if (print) parties.append(" " + party);
      } else {
        throw std::runtime_error("Incorrect party argument syntax " + party);
      }
    }
    if (print) std::cout << parties << std::endl;
  } else
    throw std::runtime_error("Other parties' information is not set but required");

  if (print) {
    std::cout << "Number of SIMD AES evaluations: " << user_options["num-simd"].as<std::size_t>()
              << std::endl;

    std::cout << "MPC Protocol: " << user_options["protocol"].as<std::string>() << std::endl;
  }
  return std::make_pair(user_options, help);
}

encrypto::motion::PartyPointer CreateParty(const program_options::variables_map& user_options) {
  const auto parties_string{user_options["other-parties"].as<const std::vector<std::string>>()};
  const auto number_of_parties{parties_string.size()};
  const auto my_id{user_options["my-id"].as<std::size_t>()};
  if (my_id >= number_of_parties) {
    throw std::runtime_error(fmt::format(
        "My id needs to be in the range [0, #parties - 1], current my id is {} and #parties is {}",
        my_id, number_of_parties));
  }

  encrypto::motion::communication::TcpPartiesConfiguration parties_configuration(number_of_parties);

  for (const auto& party_string : parties_string) {
    const auto [party_id, host, port] = ParsePartyArgument(party_string);
    if (party_id >= number_of_parties) {
      throw std::runtime_error(
          fmt::format("Party's id needs to be in the range [0, #parties - 1], current id "
                      "is {} and #parties is {}",
                      party_id, number_of_parties));
    }
    parties_configuration.at(party_id) = std::make_pair(host, port);
  }
  encrypto::motion::communication::TcpSetupHelper helper(my_id, parties_configuration);
  auto communication_layer = std::make_unique<encrypto::motion::communication::CommunicationLayer>(
      my_id, helper.SetupConnections());
  auto party = std::make_unique<encrypto::motion::Party>(std::move(communication_layer));
  auto configuration = party->GetConfiguration();
  // disable logging if the corresponding flag was set
  const auto logging{!user_options.count("disable-logging")};
  configuration->SetLoggingEnabled(logging);
  configuration->SetOnlineAfterSetup(user_options["online-after-setup"].as<bool>());
  return party;
}
