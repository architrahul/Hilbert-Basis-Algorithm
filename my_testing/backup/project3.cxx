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
class naiveAlgorithm {
public:

    static int polymers;
    static int unsplittablePolymers;
    

    static bool isComplementary(std::vector<std::vector<int>> monomers, std::vector<int> coeff1, std::vector<int> coeff2) {
        if (coeff1.size() != coeff2.size()) {
            throw std::invalid_argument("Vectors must be of the same size for complement check.");
        }
        // Check if the coefficients are complementary
        std::vector<int> v1 = HelperMethods::coeffToVector(monomers, coeff1);
        std::vector<int> v2 = HelperMethods::coeffToVector(monomers, coeff2);
        for (size_t i = 0; i < v1.size(); ++i) {
            if (v1[i]*v2[i] < 0) {
                return true;
            }
        }
        return false;
    }

    static bool isUnsplittable(std::vector<int> v, std::vector<std::vector<int>> monomers) {
        int N = v.size();
        std::vector<int> b(N, 0);

        if (v.size() != monomers.size()) {
            throw std::invalid_argument("Input vector size does not match monomers size.");
        }

        // print the input vector
        if (DEBUG) {
            std::cout << "Checking unsplittability for vector: ";
            HelperMethods::printVector(v);
        }
    
        while (true) {
            // Skip the iteration if b is a null vector
            if (!(std::all_of(b.begin(), b.end(), [](int x) { return x == 0; }))) {
                // Compute c = a - b
                std::vector<int> c(N);
                c = HelperMethods::vectorSub(v, b);
    
                // Only do if b is lex ≤ c
                if (HelperMethods::is_lex_leq(b, c)) {
                    std::vector<int> iVector = HelperMethods::coeffToVector(monomers, b);
                    std::vector<int> iComplementVector = HelperMethods::coeffToVector(monomers, c);
                    if (!isComplementary(monomers, b, c)) {
                        if (DEBUG) {
                            std::cout << "Found uncomplementary pair. Polymer is splittable." << std::endl;
                            std::cout << "b: ";
                            HelperMethods::printVector(b);
                            std::cout << "c: ";
                            HelperMethods::printVector(c);
                        }
                        
                        return false; // Found a uncomplementary pair
                    }
                }
            }
    
            // Update b
            int i = N - 1;
            while (i >= 0) {
                b[i]++;
                if (b[i] <= v[i]) break;
                b[i] = 0;
                i--;
            }
            if (i < 0) break; // Done
        }
        unsplittablePolymers++;
        if (unsplittablePolymers % 1000 == 0) {
            std::cout << "Found " << unsplittablePolymers << " unsplittable polymers so far." << std::endl;
        }
        if (DEBUG) {
        std::cout << "Polymer is unsplittable." << std::endl;
        }
        return true;
    }
};

int naiveAlgorithm::polymers = 0;
int naiveAlgorithm::unsplittablePolymers = 0;

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

    std::set<std::vector<int>> S;
    for (int i = 0; i < monomers.size(); i++) {
        std::vector<int> unitVector(monomers.size(), 0);
        unitVector[i] = 1;

        S.insert(unitVector);
    }

    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    for (auto it1 = S.begin(); it1 != S.end(); ++it1) {
        for (auto it2 = S.begin(); it2 != it1; ++it2) {
            std::vector<int> p = HelperMethods::vectorAdd(*it1, *it2);
            naiveAlgorithm::polymers++;
            if (naiveAlgorithm::polymers % 10000 == 0) {
                std::cout << "Processed " << naiveAlgorithm::polymers << " polymers so far." << std::endl;
            }
    
            // Skip the unsplittable check if polymer is already in S
            if (S.find(p) == S.end()) { // Only check if not already in S
                if ((std::accumulate(p.begin(), p.end(), 0) <= MAX_NORM) && naiveAlgorithm::isUnsplittable(p, monomers)) {
                    std::cout << "Adding polymer to S: ";
                    // Print polymer
                    std::cout << "(";
                    for (size_t i = 0; i < p.size(); i++) {
                        std::cout << p[i];
                        if (i < p.size() - 1) std::cout << ", ";
                    }
                    std::cout << ")";
                    std::cout << std::endl;
                    S.insert(p); // Add polymer to S
                }
            }
        }
    }

    // End timing
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Print S
    std::cout << "Final set of unsplittable polymers (S):" << std::endl;
    for (const auto& solution : S) {
        std::cout << "(";
        for (size_t i = 0; i < solution.size(); i++) {
            std::cout << solution[i];
            if (i < solution.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;
    }

    // Print execution time
    std::cout << "\nExecution time: " << duration.count() << " microseconds";
    std::cout << " (" << duration.count() / 1000.0 << " milliseconds)" << std::endl;

    return 0;
}