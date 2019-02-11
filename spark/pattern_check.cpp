#include <iostream>
#include <regex>

int main(int argc, char **argv) {
  std::string cur_line;
  std::string pattern;
  while (1) {
    std::cout << "Print expression or quit for exit:" << std::endl;
    std::getline(std::cin, cur_line);
    if (!cur_line.compare("quit"))
      break;
    std::cout << "Write the pattern:" << std::endl;
    std::getline(std::cin, pattern);
    if (pattern.empty()) {
      pattern =
          "^\\s*(\\w+)\\s*:\\s*(\\d+(\\.\\d*([eE]?[+-]?\\d+)?)?)\\s*(//.*)?$";
      std::cout << "Setting default pattern: " << pattern << std::endl;
    }
    std::regex rgx(pattern.c_str());
    std::smatch matches;
    if (std::regex_search(cur_line, matches, rgx)) {
      for (auto m : matches) {
        std::cout << m.str() << std::endl;
      }
    } else {
      std::cout << "No matches" << std::endl;
    }
  }
  return 0;
}