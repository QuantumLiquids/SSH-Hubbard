#pragma once
#include <stdlib.h>

bool IsElectron(size_t i, size_t Ly, size_t Np);
size_t GetNumofMps();
void Show(std::vector<size_t> v);
bool Parser(const int argc, char *argv[],
            size_t& start,
            size_t& end);
