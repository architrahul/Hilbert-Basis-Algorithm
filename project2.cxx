#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <string>
#include <numeric>
#include <sstream>

// Set to 1 to enable debug output, 0 to disable
#define DEBUG 1

#define level_limit 10000

class HilbertBasis {
private:
    const std::vector<std::vector<int>>& equations;
    const int numVars;
    const int numEquations;

    // Calculate the actual vector by multiplying the combination with equations
    std::vector<int> calculateActualVector(const std::vector<int>& combination) {
        std::vector<int> result(numVars, 0);
        for (int i = 0; i < numEquations; i++) {
            if (combination[i] != 0) {
                for (int j = 0; j < numVars; j++) {
                    result[j] += combination[i] * equations[i][j];
                }
            }
        }
        return result;
    }

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
        : equations(eqs), numEquations(eqs.size()), numVars(eqs[0].size()) {}

    std::vector<std::vector<int>> compute() {
        std::vector<std::vector<int>> basis;
        std::vector<std::pair<std::vector<int>, std::vector<bool>>> currentLevelPairs;
        
        basis.reserve(numEquations); // As in your original code
        currentLevelPairs.reserve(numEquations); // As in your original code

        // Start at level 1 with unit vectors and their initial frozen states
        for (int i = 0; i < numEquations; i++) {
            std::vector<int> unitVector(numEquations, 0);
            unitVector[i] = 1;

            std::vector<bool> initialFrozenStatus(numEquations, false);
            for (int j = i + 1; j < numEquations; j++) {
                initialFrozenStatus[j] = true;
            }
            currentLevelPairs.push_back({unitVector, initialFrozenStatus});
        }
        
        int levelCount = 1;
#if DEBUG
        std::cout << "\n--- Level " << levelCount << " ---" << std::endl;
#endif
        
        while (!currentLevelPairs.empty() && levelCount < level_limit) {
            std::vector<std::pair<std::vector<int>, std::vector<bool>>> nextLevelPairs;
            nextLevelPairs.reserve(currentLevelPairs.size() * numEquations);
    
            for (const auto& currentPair : currentLevelPairs) {
                const std::vector<int>& currentCombination = currentPair.first;
                // currentFrozenStatus is a modifiable copy of the frozen status from the pair
                std::vector<bool> currentFrozenStatus_localCopy = currentPair.second; 
                
                auto actualVector = calculateActualVector(currentCombination);
#if DEBUG
                std::cout << "Current Combination: (";
                for (size_t i = 0; i < currentCombination.size(); i++) {
                    std::cout << currentCombination[i] << (i < currentCombination.size() - 1 ? ", " : "");
                }
                std::cout << ")" << std::endl;

                std::cout << "  Frozen States (from parent pair): (";
                for (size_t i = 0; i < currentPair.second.size(); i++) { 
                    std::cout << (currentPair.second[i] ? "T" : "F") << (i < currentPair.second.size() - 1 ? ", " : "");
                }
                std::cout << ")" << std::endl;
#endif
                
                if (isSolutionVector(actualVector)) {
                    basis.push_back(currentCombination);
#if DEBUG
                    std::cout << "  --> Added to Hilbert Basis." << std::endl;
#endif
                    continue;
                }
                
                // The `possiblePaths` variable was a copy and unused, so removed.
                // Debug for "possible paths" based on the frozen status *before* inner loop modifications
#if DEBUG
                std::cout << "  Possible Paths (evaluating from current combination):" << std::endl;
                for (int i = 0; i < numEquations; ++i) { 
                    std::cout << "    Path " << i << ": ";
                    if (currentPair.second[i]) { // Use the original frozen status for this display
                        std::cout << "FROZEN (from parent)";
                    } else if (hasNegativeDotProduct(equations[i], actualVector)) {
                        std::cout << "AVAILABLE";
                    } else {
                        std::cout << "NO NEG DOT PRODUCT";
                    }
                    std::cout << std::endl;
                }
#endif
                // This variable `possiblePaths` was unused in your original code after this line.
                // std::vector<bool> possiblePaths = currentFrozenStatus_localCopy; 
                int prevPathIdx = -1;

                // Inner loop iterates from high to low index.
                for (int path_taken_idx = numEquations - 1; path_taken_idx >= 0; path_taken_idx--) {
                    // Use the local copy `currentFrozenStatus_localCopy` which might be modified by prevPathIdx logic
                    if (currentFrozenStatus_localCopy[path_taken_idx]) { 
                        continue;
                    }
                    
                    if (hasNegativeDotProduct(equations[path_taken_idx], actualVector)) {
                        auto newCombination = currentCombination;
                        newCombination[path_taken_idx]++;
                        
                        // Your existing logic for modifying currentFrozenStatus_localCopy and creating newFrozenStatus
                        if (prevPathIdx != -1) {
                            currentFrozenStatus_localCopy[prevPathIdx] = true; 
                        }
                        std::vector<bool> newFrozenStatus = currentFrozenStatus_localCopy; 
                        prevPathIdx = path_taken_idx; 
                        
                        if (!isGreaterThanAnyBasis(newCombination, basis)) {
#if DEBUG
                            std::cout << "    Taking path " << path_taken_idx << ". New Combination: (";
                            for (size_t i = 0; i < newCombination.size(); i++) {
                                std::cout << newCombination[i] << (i < newCombination.size() - 1 ? ", " : "");
                            }
                            std::cout << ")" << std::endl;
                            std::cout << "      New Frozen States for next level: (";
                            for (size_t i = 0; i < newFrozenStatus.size(); i++) {
                                std::cout << (newFrozenStatus[i] ? "T" : "F") << (i < newFrozenStatus.size() - 1 ? ", " : "");
                            }
                            std::cout << ")" << std::endl;
#endif
                            nextLevelPairs.push_back({newCombination, newFrozenStatus});
                        }
                    }
                }
            }
            
            levelCount++;
            currentLevelPairs = std::move(nextLevelPairs);
#if DEBUG
            if (!currentLevelPairs.empty()) {
                std::cout << "\n--- Level " << levelCount << " ---" << std::endl;
            } else if (levelCount < level_limit) { 
                std::cout << "\n--- No more vectors to process. Algorithm finished. ---" << std::endl;
            }
#endif
        }
#if DEBUG
        if (levelCount >= level_limit && !currentLevelPairs.empty()) {
             std::cout << "\n--- Reached level_limit (" << level_limit << "). Stopping. ---" << std::endl;
        }
#endif
        return basis;
    }
};

// Example usage
int main() {
    std::vector<std::vector<int>> equations = {{-1, -1}, {1, 3}, {2, -2}, {-3, -1}};
    
    auto start = std::chrono::high_resolution_clock::now();

    HilbertBasis hb(equations);
    std::vector<std::vector<int>> basis = hb.compute();

    // End timing
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // Print the final basis (always show this regardless of DEBUG setting)
    std::cout << "\nHilbert Basis:" << std::endl;
    for (const auto& solution : basis) {
        std::cout << "(";
        for (size_t i = 0; i < solution.size(); i++) {
            std::cout << solution[i];
            if (i < solution.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;
    }

    // Print execution time (always show this regardless of DEBUG setting)
    std::cout << "\nExecution time: " << duration.count() << " microseconds";
    std::cout << " (" << duration.count() / 1000.0 << " milliseconds)" << std::endl;
    
    return 0;
}