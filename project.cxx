#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <string>
#include <numeric>

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

            for (const auto& current : currentLevel) {
                auto actualVector = calculateActualVector(current);
                
                if (isSolutionVector(actualVector)) {
                    basis.push_back(current);
                    continue;
                }

                for (int i = 0; i < numEquations; i++) {
                    if (hasNegativeDotProduct(equations[i], actualVector)) {
                        auto newCombination = current;
                        newCombination[i]++;
                        
                        // Use seenVectors cache to avoid duplicates
                        if (!isGreaterThanAnyBasis(newCombination, basis) && 
                            seenVectors.insert(newCombination).second) {
                            nextLevel.push_back(std::move(newCombination));
                        }
                    }
                }
            }
            std::cout << "Level " << levelCount << ":" << std::endl;
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
    
    // Print the basis
    for (const auto& solution : basis) {
        for (int x : solution) {
            std::cout << x << " ";
        }
        std::cout << std::endl;
    }

    // Print execution time
    std::cout << "\nExecution time: " << duration.count() << " microseconds";
    std::cout << " (" << duration.count() / 1000.0 << " milliseconds)" << std::endl;
    
    return 0;
}