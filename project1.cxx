#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <string>
#include <numeric>

// Set to 1 to enable debug output, 0 to disable
#define DEBUG 1

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
    const std::vector<std::vector<int>>& equations;  // Changed to reference
    const int numVars;
    const int numEquations;
    std::unordered_set<std::vector<int>, VectorHash> seenVectors;  // Cache seen vectors

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
            [&vec](const std::vector<int>& basisVec) {  // Explicitly specify type
                return std::equal(vec.begin(), vec.end(), basisVec.begin(),
                    [](int a, int b) { return a >= b; });
            });
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

#if DEBUG
            // Print current vectors being processed
            std::cout << "\nProcessing vectors at level " << levelCount << ":" << std::endl;
            for (const auto& vec : currentLevel) {
                std::cout << "  (";
                for (size_t i = 0; i < vec.size(); i++) {
                    std::cout << vec[i];
                    if (i < vec.size() - 1) std::cout << ",";
                }
                std::cout << ")" << std::endl;
            }
#endif

            for (const auto& current : currentLevel) {
                auto actualVector = calculateActualVector(current);
                
#if DEBUG
                // Print current vector and its actual vector
                std::cout << "(";
                for (size_t i = 0; i < current.size(); i++) {
                    std::cout << current[i];
                    if (i < current.size() - 1) std::cout << ",";
                }
                std::cout << ") -> Actual: (";
                for (size_t i = 0; i < actualVector.size(); i++) {
                    std::cout << actualVector[i];
                    if (i < actualVector.size() - 1) std::cout << ",";
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
#endif

                for (int i = 0; i < numEquations; i++) {
                    if (hasNegativeDotProduct(equations[i], actualVector)) {
#if DEBUG
                        if (anyValid) std::cout << ", ";
                        std::cout << "path " << i;
                        anyValid = true;
#endif
                        
                        auto newCombination = current;
                        newCombination[i]++;
                        
                        // Use seenVectors cache to avoid duplicates
                        if (!isGreaterThanAnyBasis(newCombination, basis) && 
                            seenVectors.insert(newCombination).second) {
                            nextLevel.push_back(std::move(newCombination));
#if DEBUG
                            std::cout << " (added)";
#endif
                        }
#if DEBUG
                        else {
                            std::cout << " (skipped - already seen or dominated)";
                        }
#endif
                    }
                }
#if DEBUG
                if (!anyValid) std::cout << "none";
                std::cout << std::endl;
#endif
            }

#if DEBUG
            std::cout << "\nLevel " << levelCount << " Summary:" << std::endl;
            std::cout << "Current Level Size: " << currentLevel.size() << std::endl;
            std::cout << "Next Level Size: " << nextLevel.size() << std::endl;
            std::cout << "Current Basis Size: " << basis.size() << std::endl;
            std::cout << "------------------------" << std::endl;
#endif

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
    
    // Print the basis (always show this regardless of DEBUG setting)
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