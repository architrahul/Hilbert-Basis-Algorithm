#pragma once

#include <vector>
#include <numeric>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>
#include <set>
#include <regex>
#include "HelperMethods.hxx"

#define MAX_NORM 100
#define DEBUG 0
#define INPUT_FILE "vectors.txt"
class dualAlgorithm {
public:

    static int polymers;
    static int unsplittablePolymers;
    
    static std::vector<std::vector<int>> formF(std::vector<std::vector<int>> monomers) {
        std::vector<std::vector<int>> F = monomers;
        for (size_t i = 0; i < monomers.size(); i++) {
            F.push_back(HelperMethods::vectorNegative(monomers[i]));
        }
        return F;
    }

    static std::vector<std::vector<int>> formSD(std::vector<std::vector<int>> monomers) {
        std::vector<std::vector<int>> SD;
        for (size_t i = 0; i < monomers.size(); i++) {
            for (size_t j = i + 1; j < monomers.size(); j++) {
                std::vector<int> v = HelperMethods::vectorAdd(monomers[i], monomers[j]);
                if (v != std::vector<int>(monomers.size(), 0)) {
                    SD.push_back(v);
                }
            }
        }
        return SD;
    }

    static std::vector<int> remainder(std::vector<int> v, std::vector<std::vector<int>> divisors) {
        std::vector<int> R(v.size(), 0);
        for (size_t i = 0; i < divisors.size(); i++) {
            bool sameOrthant = true;
            std::vector<int> d = divisors[i];
            for (size_t j = 0; j < v.size(); j++) {
                if (v[j] * d[j] < 0) {
                    sameOrthant = false;
                }
            }

            if (!sameOrthant) {
                continue;
            }

            std::vector<int> q(v.size());
            for (size_t j = 0; j < v.size(); j++) {
                if (d[j] != 0) {
                    q[j] = int(std::abs(v[j]) / std::abs(d[j]));
                } else {
                    q[j] = std::numeric_limits<int>::max();
                }
            }

            int Q = *std::min_element(q.begin(), q.end());
            for (size_t j = 0; j < v.size(); j++) {
                R[j] = v[j] - (d[j] * Q);
            }
            v = R;
        }
        return R;
    }

};

int main(int argc, char* argv[]) {

    std::vector<std::vector<int>> monomers;

    if (argc < 2) {
        std::cerr << "no input file found, defaulting to vector.txt" << std::endl;
        try {
            monomers = HelperMethods::parseMonomersFile(INPUT_FILE);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return 1;
        }
    } else {
        try {
            monomers = HelperMethods::parseMonomersFile(argv[1]);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return 1;
        }
    }

    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    // Implement procedure completion algorithm
    std::vector<std::vector<int>> F = dualAlgorithm::formF(monomers);
    std::vector<std::vector<int>> SD = dualAlgorithm::formSD(F);

    while (!SD.empty()) {
        std::vector<int> v = SD.back(); // Choose the last element from SD
        SD.pop_back(); // Remove v from SD
        SD.erase(std::remove(SD.begin(), SD.end(), HelperMethods::vectorNegative(v)), SD.end()); // Remove -v from SD
        v = dualAlgorithm::remainder(v, F);
        if (v != std::vector<int>(monomers.size(), 0)) {
            for (size_t i = 0; i < F.size(); i++) {
                std::vector<int> v1 = HelperMethods::vectorAdd(v, F[i]);
                if (std::find(SD.begin(), SD.end(), v1) == SD.end()) {
                    SD.push_back(v1);
                }
                std::vector<int> v2 = HelperMethods::vectorSub(v, F[i]);
                if (std::find(SD.begin(), SD.end(), v2) == SD.end()) {
                    SD.push_back(v2);
                }
            }
            if (std::find(F.begin(), F.end(), v) == F.end()) {
                F.push_back(v);
            }
            if (std::find(F.begin(), F.end(), HelperMethods::vectorNegative(v)) == F.end()) {
                F.push_back(HelperMethods::vectorNegative(v));
            }
        }
    }
    // Implement Procedure reduction algorithm
    // while there exist v, v' in F such that v != v' and v' divides v
    // F = F - {v, -v} union {dualAlgorithm::remainder(v, {v'}), -dualAlgorithm::remainder(v, {v'})}
    
    std::vector<std::vector<int>> solution = F;
    // End timing
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);



    // Print execution time
    std::cout << "\nExecution time: " << duration.count() << " microseconds";
    std::cout << " (" << duration.count() / 1000.0 << " milliseconds)" << std::endl;

    return 0;
}