#pragma once

#include <iostream>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>
#include <set>
#include <regex>

class HelperMethods {
public:
    static std::vector<std::vector<int>> parseMonomersFile(std::string filename);
    static std::vector<std::vector<int>> add_unit_monomers(std::vector<std::vector<int>> monomers);
    static std::vector<int> vectorAdd(std::vector<int> v1, std::vector<int> v2);
    static std::vector<int> vectorSub(std::vector<int> v1, std::vector<int> v2);
    static std::vector<int> vectorNegative (std::vector<int> v);
    static bool is_lex_leq(std::vector<int> a, std::vector<int> b);
    static void printVector(std::vector<int>& v);
    static std::vector<int> coeffToVector(std::vector<std::vector<int>> monomers, std::vector<int> coeff);
    static std::vector<std::vector<int>> remove_unit_monomers(std::vector<std::vector<int>> basis, int n_dim);
};