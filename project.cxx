#include <vector>
#include <stack>
#include <algorithm>
#include <ostream>
#include <iostream>

class HilbertBasis {
private:
    std::vector<std::vector<int>> equations; // Matrix representing a(x)
    int numVars;                            // Number of variables
    int numEquations;                       // Number of equations

    // Helper function to check if t1 <= t2 componentwise
    bool isGreaterThanAnyBasis(const std::vector<int>& vec, const std::vector<std::vector<int>>& basis) {
        for (const auto& basisVec : basis) {
            bool isGreater = true;
            for (size_t i = 0; i < vec.size(); i++) {
                if (vec[i] < basisVec[i]) {
                    isGreater = false;
                    break;
                }
            }
            if (isGreater) return true;
        }
        return false;
    }

    // Helper function to check dot product of two vectors
    bool dotProduct(const std::vector<int>& v1, const std::vector<int>& v2) {
        int sum = 0;
        for (int i = 0; i < numVars; i++) {
            sum += v1[i] * v2[i];
        }
        if (sum >= 0) return true;
        return false;
    }

public:
    HilbertBasis(const std::vector<std::vector<int>>& eqs) 
        : equations(eqs), numEquations(eqs.size()), numVars(eqs[0].size()) {}

        std::vector<std::vector<int>> compute() {
            std::vector<std::vector<int>> basis;
            std::vector<std::vector<int>> currentLevel;
            int levelCount = 0;
            
            // Initialize starting level with unit vectors
            for (int i = 0; i < numEquations; i++) {
                std::vector<int> unitVector(numEquations, 0);
                unitVector[i] = 1;
                currentLevel.push_back(unitVector);
            }
        
            while (!currentLevel.empty() && levelCount < 10) {
                levelCount++;
                std::vector<std::vector<int>> nextLevel;
                
                for (const auto& current : currentLevel) {
                    // Calculate the actual vector from linear combination
                    std::vector<int> actualVector(numVars, 0);
                    for (int i = 0; i < numEquations; i++) {
                        for (int j = 0; j < numVars; j++) {
                            actualVector[j] += current[i] * equations[i][j];
                        }
                    }
                    
                    // Check if this is a solution (sum = 0)
                    bool isSolution = true;
                    for (int j = 0; j < numVars; j++) {
                        if (actualVector[j] != 0) {
                            isSolution = false;
                            break;
                        }
                    }
                    if (isSolution) {
                        basis.push_back(current);  // Store the combination vector instead
                        continue;
                    }
                    
                    // Generate next level combinations
                    for (int i = 0; i < numEquations; i++) {
                        if (!dotProduct(equations[i], actualVector)) {
                            std::vector<int> newCombination = current;
                            newCombination[i]++;
                            
                            // Only add if not greater than any existing basis vector
                            if (!isGreaterThanAnyBasis(newCombination, basis)) {
                                nextLevel.push_back(newCombination);
                            }
                        }
                    }
                }
                
                // Debug output
                std::cout << "Level " << levelCount << ":" << std::endl;
                std::cout << "Current Level Size: " << currentLevel.size() << std::endl;
                std::cout << "Next Level Size: " << nextLevel.size() << std::endl;
                std::cout << "Current Basis Size: " << basis.size() << std::endl;
                std::cout << "------------------------" << std::endl;
                
                currentLevel = nextLevel;
            }
            
            if (levelCount >= 10) {
                std::cout << "Warning: Maximum level (10) reached." << std::endl;
            }
            
            return basis;
        }
};

// Example usage
int main() {
    std::vector<std::vector<int>> equations = {{-1, -1}, {1, 3}, {2, -2}, {-3, -1}};
    
    HilbertBasis hb(equations);
    std::vector<std::vector<int>> basis = hb.compute();
    
    // Print the basis
    for (const auto& solution : basis) {
        for (int x : solution) {
            std::cout << x << " ";
        }
        std::cout << std::endl;
    }
    
    return 0;
}