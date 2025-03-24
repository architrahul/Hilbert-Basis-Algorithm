#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <string>
#include <numeric>
#include <iomanip>
#include <sstream>

// Add a hash function for vectors to use with unordered_set
struct VectorHash {
    size_t operator()(const std::vector<int>& v) const {
        size_t hash = v.size();
        for (auto& i : v) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

class HilbertBasis {
private:
    const std::vector<std::vector<int>>& equations;
    const int numVars;
    const int numEquations;
    std::unordered_set<std::vector<int>, VectorHash> seenVectors;

    // Optimized actual vector calculation using vector operations
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

    // Optimized solution check
    bool isSolutionVector(const std::vector<int>& vec) const {
        return std::all_of(vec.begin(), vec.end(), [](int x) { return x == 0; });
    }

    // Optimized dot product using std::inner_product
    bool hasNegativeDotProduct(const std::vector<int>& v1, const std::vector<int>& v2) const {
        return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0) < 0;
    }

    bool isGreaterThanAnyBasis(const std::vector<int>& vec, 
                              const std::vector<std::vector<int>>& basis) const {
        return std::any_of(basis.begin(), basis.end(),
            [&vec](const std::vector<int>& basisVec) {
                return std::equal(vec.begin(), vec.end(), basisVec.begin(),
                    [](int a, int b) { return a >= b; });
            });
    }

    // New method to check if a position is frozen
    bool isPositionFrozen(int position, const std::vector<int>& combination, 
        const std::vector<int>& actualVector) const {
        // Check if any higher priority path (lower index) is already taken
        for (int i = 0; i < position; i++) {
            if (combination[i] > 0) {
                return true;  // Position is frozen because lower index has a value
            }
        }

        // Count available higher priority paths
        int availableHigherPriorityPaths = 0;
        for (int i = 0; i < position; i++) {
            if (hasNegativeDotProduct(equations[i], actualVector)) {
                availableHigherPriorityPaths++;
            }
        }

        // If there are available higher priority paths, then freeze this position
        return availableHigherPriorityPaths > 0;
    }

    // Helper method to get a string representation of a vector with underlined frozen positions
    std::string getVectorStringWithFrozenIndicators(const std::vector<int>& combination, 
                                                   const std::vector<int>& actualVector) const {
        std::stringstream ss;
        ss << "(";
        for (int i = 0; i < numEquations; i++) {
            // Check if this position is frozen
            bool frozen = isPositionFrozen(i, combination, actualVector);
            
            if (frozen) {
                // Underline for frozen positions
                ss << "\033[4m" << combination[i] << "\033[0m";
            } else {
                ss << combination[i];
            }
            
            if (i < numEquations - 1) {
                ss << ", ";
            }
        }
        ss << ")";
        return ss.str();
    }

public:
    HilbertBasis(const std::vector<std::vector<int>>& eqs) 
        : equations(eqs), numEquations(eqs.size()), numVars(eqs[0].size()) {
        seenVectors.reserve(1000);  // Prereserve space
    }

    std::vector<std::vector<int>> compute() {
        std::vector<std::vector<int>> basis;
        std::vector<std::vector<int>> currentLevel;
        basis.reserve(100);  // Prereserve space
        currentLevel.reserve(100);
        
        // Initialize with unit vectors
        for (int i = 0; i < numEquations; i++) {
            std::vector<int> unitVector(numEquations, 0);
            unitVector[i] = 1;
            currentLevel.push_back(std::move(unitVector));
        }

        int levelCount = 0;
        while (!currentLevel.empty() && levelCount++ < 10) {
            std::vector<std::vector<int>> nextLevel;
            nextLevel.reserve(currentLevel.size() * numEquations);

            std::cout << "\nLevel " << levelCount << ":" << std::endl;
            std::cout << "Current Level Vectors (underline = frozen position):" << std::endl;

            for (const auto& current : currentLevel) {
                auto actualVector = calculateActualVector(current);
                
                // Print the current vector with frozen positions underlined
                std::cout << getVectorStringWithFrozenIndicators(current, actualVector);
                // Add additional info about the actual vector
                std::cout << " -> Actual: (";
                for (size_t j = 0; j < actualVector.size(); j++) {
                    std::cout << actualVector[j];
                    if (j < actualVector.size() - 1) std::cout << ", ";
                }
                std::cout << ")" << std::endl;
                
                if (isSolutionVector(actualVector)) {
                    basis.push_back(current);
                    std::cout << "  → Added to basis (solution vector)" << std::endl;
                    continue;
                }
                
                bool anyPathAdded = false;
                for (int i = 0; i < numEquations; i++) {
                    // Check if this position is frozen
                    bool frozen = isPositionFrozen(i, current, actualVector);
                    bool negDotProduct = hasNegativeDotProduct(equations[i], actualVector);
                    
                    // Show detailed path analysis
                    std::cout << "  Path " << i << ": ";
                    if (frozen) std::cout << "FROZEN";
                    else if (!negDotProduct) std::cout << "NO NEG DOT PRODUCT";
                    else std::cout << "AVAILABLE";
                    std::cout << std::endl;
                    
                    // Skip frozen positions
                    if (frozen) {
                        continue;
                    }
            
                    if (negDotProduct) {
                        auto newCombination = current;
                        newCombination[i]++;
                        
                        // Use seenVectors cache to avoid duplicates
                        if (!isGreaterThanAnyBasis(newCombination, basis) && 
                            seenVectors.insert(newCombination).second) {
                            nextLevel.push_back(std::move(newCombination));
                            anyPathAdded = true;
                            std::cout << "    → Added to next level" << std::endl;
                        } else {
                            std::cout << "    → Skipped (already seen or dominated)" << std::endl;
                        }
                    }
                }
                
                if (!anyPathAdded) {
                    std::cout << "  No paths added from this vector" << std::endl;
                }
            }
            
            std::cout << "\nSummary:" << std::endl;
            std::cout << "Current Level Size: " << currentLevel.size() << std::endl;
            std::cout << "Next Level Size: " << nextLevel.size() << std::endl;
            std::cout << "Current Basis Size: " << basis.size() << std::endl;
            std::cout << "------------------------" << std::endl;

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
    
    // Print the final basis
    std::cout << "\nFinal Hilbert Basis:" << std::endl;
    for (const auto& solution : basis) {
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