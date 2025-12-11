#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <string>
#include <numeric>
#include <sstream>
#include <fstream>
#include "HelperMethods.hxx"

// Set to 1 to enable debug output, 0 to disable
#define DEBUG 0
#define level_limit 200
#define mode 1 // 0 for Hilbert Basis, 1 for naive algorithm

class HilbertBasis {
public:
    const std::vector<std::vector<int>>& monomers;
    const int numVars;
    const int nummonomers;

    // Check if a vector is a solution vector (all zeros)
    bool isSolutionVector(const std::vector<int>& vec) const {
        return std::all_of(vec.begin(), vec.end(), [](int x) { return x == 0; });
    }

    // Check if two vectors have a negative dot product
    bool hasNegativeDotProduct(const std::vector<int>& v1, const std::vector<int>& v2) const {
        return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0) < 0;
    }

    // Check if a vector is dominated by any basis vector
    bool isGreaterThanAnyBasis(const std::vector<int>& vec, 
                              const std::vector<std::vector<int>>& basis) const {
        return std::any_of(basis.begin(), basis.end(),
            [&vec](const std::vector<int>& basisVec) {
                return std::equal(vec.begin(), vec.end(), basisVec.begin(),
                    [](int a, int b) { return a >= b; });
            });
    }

public:
    HilbertBasis(const std::vector<std::vector<int>>& eqs) 
        : monomers(eqs), nummonomers(eqs.size()), numVars(eqs[0].size()) {}

    std::vector<std::vector<int>> compute() {
        std::vector<std::vector<int>> basis;
        std::vector<std::pair<std::vector<int>, std::vector<bool>>> currentLevelPairs;
        
        basis.reserve(nummonomers);
        currentLevelPairs.reserve(nummonomers);
        // Start at level 1 with unit vectors and their initial frozen states
        for (int i = 0; i < nummonomers; i++) {
            std::vector<int> unitVector(nummonomers, 0);
            unitVector[i] = 1;

            std::vector<bool> initialFrozenStatus(nummonomers, false);
            for (int j = i + 1; j < nummonomers; j++) {
                initialFrozenStatus[j] = true;
            }
            currentLevelPairs.push_back({unitVector, initialFrozenStatus});
        }
        
        int levelCount = 1;
        
        while (!currentLevelPairs.empty() && levelCount <= level_limit) {
            if (DEBUG) {
                std::cout << "\nProcessing level " << levelCount << " with " 
                      << currentLevelPairs.size() << " pairs." << std::endl;
            }
            
            std::vector<std::pair<std::vector<int>, std::vector<bool>>> nextLevelPairs;
            nextLevelPairs.reserve(currentLevelPairs.size() * nummonomers);
    
            for (const auto& currentPair : currentLevelPairs) {
                if (DEBUG) {
                    std::cout << "Current combination: ";
                    for (const auto& val : currentPair.first) {
                        std::cout << val << " ";
                    }
                    std::cout << "\n";
                }
                
                const std::vector<int>& currentCombination = currentPair.first;
                std::vector<bool> currentFrozenStatus = currentPair.second;
                
                auto actualVector = HelperMethods::coeffToVector(monomers, currentCombination);
                
                if (isSolutionVector(actualVector)) {
                    basis.push_back(currentCombination);
                    std::cout << "Added to basis (solution vector): ";
                    for (const auto& val : currentCombination) {
                        std::cout << val << " ";
                    }
                    std::cout << std::endl;
                    continue;
                }
                
                std::vector<bool> possiblePaths = currentFrozenStatus;
                int prevPathIdx = -1;

                for (int path_taken_idx = nummonomers - 1; path_taken_idx >= 0; path_taken_idx--) {
                    if (currentFrozenStatus[path_taken_idx]) {
                        continue;
                    }
                    
                    if (hasNegativeDotProduct(monomers[path_taken_idx], actualVector)) {
                        auto newCombination = currentCombination;
                        newCombination[path_taken_idx]++;
                        
                        // Freeze positions after the path taken
                        if (prevPathIdx != -1) {
                            currentFrozenStatus[prevPathIdx] = true;
                        }
                        std::vector<bool> newFrozenStatus = currentFrozenStatus;
                        prevPathIdx = path_taken_idx;
                        
                        if (!isGreaterThanAnyBasis(newCombination, basis)) {
                            nextLevelPairs.push_back({newCombination, newFrozenStatus});
                        }
                    }
                }
            }
            levelCount++;
            currentLevelPairs = std::move(nextLevelPairs);
        }
        
        return basis;
    }
};

// Example usage
int main(int argc, char* argv[]) {
    std::vector<std::vector<int>> monomers;
    int og_monomers_size = 0;
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    try {
        monomers = HelperMethods::parseMonomersFile(argv[1]);
        if (mode == 1) {
            og_monomers_size = monomers.size();
            monomers = HelperMethods::add_unit_monomers(monomers);
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    
    auto start = std::chrono::high_resolution_clock::now();

    HilbertBasis hb(monomers);
    std::vector<std::vector<int>> basis = hb.compute();
    if (mode == 1) {
        std::cout << "Naive algorithm mode enabled. Removing unit monomers." << std::endl;
        basis = HelperMethods::remove_unit_monomers(basis, og_monomers_size);
    }

    // End timing
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // Print the final basis 
    std::cout << "\nHilbert Basis:" << std::endl;
    for (const auto& solution : basis) {
        std::cout << "(";
        for (size_t i = 0; i < solution.size(); i++) {
            std::cout << solution[i];
            if (i < solution.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;
    }


    std::cout << "\nExecution time: " << duration.count() << " microseconds";
    std::cout << " (" << duration.count() / 1000.0 << " milliseconds)" << std::endl;
    
    return 0;
}

