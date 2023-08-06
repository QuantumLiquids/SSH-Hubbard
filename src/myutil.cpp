#include "gqmps2/gqmps2.h"
using gqmps2::kMpsPath;
using gqmps2::kMpsTenBaseName;
using gqmps2::kGQTenFileSuffix;

//number of mps file in default mps path("./mps")
size_t GetNumofMps() {
  size_t NumberOfMpsFile = 0;
  for (NumberOfMpsFile = 0; NumberOfMpsFile < 1e5; NumberOfMpsFile++) {
    std::string file;
    file = kMpsPath + "/" + kMpsTenBaseName + std::to_string(NumberOfMpsFile) + "." + kGQTenFileSuffix;
    std::ifstream ifs(file, std::ifstream::binary);
    if (ifs.good()) {
      ifs.close();
    } else {
      break;
    }
  }
  return NumberOfMpsFile;
}

//If site `i` is electron site, for SSH-Hubbard model
bool IsElectron(size_t i, size_t Ly, size_t Np) {
  size_t residue = i % ((2 * Np + 1) * Ly);
  if (residue < (Np + 1) * Ly && residue % (Np + 1) == 0) {
    return true;
  } else return false;
}

void Show(std::vector<size_t> v) {
  for (auto iter = v.begin(); iter < v.end(); iter++) {
    std::cout << *iter << ",";
  }
  std::cout << '\b' << std::endl;
}

// When used to measure, note if should not set start too small to exceed canonical center.
bool Parser(const int argc, char *argv[],
            size_t &start,
            size_t &end) {
  int nOptionIndex = 1;

  std::string arguement1 = "--start=";
  std::string arguement2 = "--end=";
  bool start_argument_has(false), end_argument_has(false);
  while (nOptionIndex < argc) {
    if (strncmp(argv[nOptionIndex], arguement1.c_str(), arguement1.size()) == 0) {
      std::string para_string = &argv[nOptionIndex][arguement1.size()];
      start = atoi(para_string.c_str());
      start_argument_has = true;
    } else if (strncmp(argv[nOptionIndex], arguement2.c_str(), arguement2.size()) == 0) {
      std::string para_string = &argv[nOptionIndex][arguement2.size()];
      end = atoi(para_string.c_str());
      end_argument_has = true;
    }
    nOptionIndex++;
  }

  if (start_argument_has != end_argument_has) {
    std::cout << "Only setting one start/end argument, exit(1)." << std::endl;
    exit(1);
  }

  if (!start_argument_has) {
    std::cout << "Note: no start/end argument, set it by default (L/4, 3*L/4)." << std::endl;
  }

  return start_argument_has;
}

bool ParserBondDimension(int argc, char *argv[],
                         std::vector<size_t> &D_set) {
  int nOptionIndex = 1;
  std::string D_string;
  std::string arguement1 = "--D=";
  bool has_D_parameter(false);
  while (nOptionIndex < argc) {
    if (strncmp(argv[nOptionIndex], arguement1.c_str(), arguement1.size()) == 0) {
      D_string = &argv[nOptionIndex][arguement1.size()];
      has_D_parameter = true;
    }
    nOptionIndex++;
  }

  //split thread num list
  const char *split = ",";
  char *p;
  const size_t MAX_CHAR_LENTH = 1000;
  char D_char[MAX_CHAR_LENTH];
  for (size_t i = 0; i < MAX_CHAR_LENTH; i++) {
    D_char[i] = 0;
  }

  strcpy(D_char, D_string.c_str());

  p = strtok(D_char, split);
  while (p != nullptr) {
    D_set.push_back(atoi(p));
    p = strtok(nullptr, split);
  }

  return has_D_parameter;
}
