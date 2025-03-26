#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <string>
#include <numeric>
#include <sstream>

// Set to 1 to enable debug output, 0 to disable
#define DEBUG 0

class HilbertBasis {
private:
    const std::vector<std::vector<int>>& equations;
    const int numVars;
    const int numEquations;

    // Calculate the actual vector by multiplying the combination with equations
    std::vector<int> calculateActualVector(const std::vector<int>& combination) {
        std::vector<int> result(numVars, 0);
        for (int i = 0; i < numEquations; i++) {
            if (combination[i] != 0) {  // Skip if coefficient is 0
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

    // Check if a position is frozen in the current vector
    bool isPositionFrozen(int position, const std::vector<int>& combination) const {
        // Find highest taken position (highest priority path taken)
        int highestTakenIndex = -1;
        for (int i = 0; i < numEquations; i++) {
            if (combination[i] > 0) highestTakenIndex = i;
        }
        
        // Positions after the highest taken index are frozen
        return position > highestTakenIndex;
    }

public:
    HilbertBasis(const std::vector<std::vector<int>>& eqs) 
        : equations(eqs), numEquations(eqs.size()), numVars(eqs[0].size()) {}

    std::vector<std::vector<int>> compute() {
        std::vector<std::vector<int>> basis;
        std::vector<std::vector<int>> currentLevel;
        basis.reserve(100);
        currentLevel.reserve(100);
        
        // Start at level 1 with unit vectors
        for (int i = 0; i < numEquations; i++) {
            std::vector<int> unitVector(numEquations, 0);
            unitVector[i] = 1;
            currentLevel.push_back(unitVector);
        }
        
        int levelCount = 1;  // Start at level 1
        
#if DEBUG
        std::cout << "\nLevel " << levelCount << ":" << std::endl;
        for (int i = 0; i < numEquations; i++) {
            std::cout << "path " << i << ": (";
            for (int j = 0; j < numEquations; j++) {
                if (j > 0) std::cout << ",";
                if (j == i) std::cout << "1";
                else std::cout << "0";
            }
            std::cout << ")" << std::endl;
        }
#endif
        
        while (!currentLevel.empty() && levelCount < 10) {
            std::vector<std::vector<int>> nextLevel;
            nextLevel.reserve(currentLevel.size() * numEquations);
    
            // Process each vector at the current level
            for (const auto& current : currentLevel) {
                auto actualVector = calculateActualVector(current);
                
#if DEBUG
                // Print the current vector
                std::cout << "(";
                for (size_t i = 0; i < current.size(); i++) {
                    std::cout << current[i];
                    if (i < current.size() - 1) std::cout << ",";
                }
                std::cout << ")";
                
                // Print which positions are frozen
                std::cout << " (";
                int highestTakenIndex = -1;
                for (int i = 0; i < numEquations; i++) {
                    if (current[i] > 0) highestTakenIndex = i;
                }
                
                if (highestTakenIndex != -1) {
                    bool firstFrozen = true;
                    for (int i = highestTakenIndex + 1; i < numEquations; i++) {
                        if (!firstFrozen) std::cout << ", ";
                        std::cout << "position " << i << " frozen";
                        firstFrozen = false;
                    }
                } else {
                    std::cout << "no positions frozen";
                }
                std::cout << ")" << std::endl;
#endif
                
                if (isSolutionVector(actualVector)) {
                    basis.push_back(current);
#if DEBUG
                    std::cout << "  → Added to basis (solution vector)" << std::endl;
#endif
                    continue;
                }
                
#if DEBUG
                // Print valid paths
                std::cout << "  Valid paths: ";
                bool anyValid = false;
                for (int i = 0; i < numEquations; i++) {
                    // Skip frozen positions
                    if (isPositionFrozen(i, current)) {
                        continue;
                    }
                    
                    if (hasNegativeDotProduct(equations[i], actualVector)) {
                        if (anyValid) std::cout << ", ";
                        std::cout << "path " << i << " (+";
                        for (int j = 0; j < numEquations; j++) {
                            if (j > 0) std::cout << ",";
                            if (j == i) std::cout << "1";
                            else std::cout << "0";
                        }
                        std::cout << ")";
                        anyValid = true;
                    }
                }
                if (!anyValid) std::cout << "none";
                std::cout << std::endl;
#endif

                // Process each valid path
                for (int i = 0; i < numEquations; i++) {
                    // Skip frozen positions
                    if (isPositionFrozen(i, current)) {
                        continue;
                    }
                    
                    if (hasNegativeDotProduct(equations[i], actualVector)) {
                        auto newCombination = current;
                        newCombination[i]++;
                        
                        // Only check for dominance, no need to check for duplicates
                        // since the freezing mechanism prevents duplicates
                        if (!isGreaterThanAnyBasis(newCombination, basis)) {
                            nextLevel.push_back(std::move(newCombination));
                        }
                    }
                }
            }
            
            levelCount++;
            if (!nextLevel.empty()) {
#if DEBUG
                std::cout << "\nLevel " << levelCount << ":" << std::endl;
#endif
            }
            currentLevel = std::move(nextLevel);
        }
        
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
    std::cout << "\nFinal Hilbert Basis:" << std::endl;
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