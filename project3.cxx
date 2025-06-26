#include <vector>
#include <numeric>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#define MAX_NORM 3
class naiveAlgorithm {
public:
    static std::vector<int> calculateWeights(std::vector<int> v) {
        int n = v.size();
        if (n == 0) {
            return {};
        }

        std::vector<int> weights(n);
        weights[n - 1] = 1;

        for (int i = n - 2; i >= 0; i--) {
            weights[i] = weights[i + 1] * ((v[i + 1]) + 1);
        }
        return weights;
    }

    static int vectorToInteger(std::vector<int> v, std::vector<int> weights) {
        if (v.empty()) return 0;

        int intValue = 0;
        for (size_t i = 0; i < v.size(); i++) {
            intValue += v[i] * weights[i];
        }
        return intValue;
    }

    static std::vector<int> integerToVector(int intValue, std::vector<int> weights) {
        std::vector<int> v(weights.size());
        int remainingValue = intValue;

        for (int i = 0; i < weights.size(); ++i) {
            if (weights[i] <= 0) {
                throw std::invalid_argument("Invalid weight (must be positive) at index " + std::to_string(i) + " during integer to vector conversion.");
            }
            v[i] = static_cast<int>(remainingValue / weights[i]);
            remainingValue %= weights[i];
        }

        if (remainingValue != 0) {
            throw std::overflow_error("Remainder is non-zero.");
        }

        return v;
    }

    static std::vector<int> coeffToVector(std::vector<std::vector<int>> monomers, std::vector<int> coeff) {
        int size = coeff.size();
        std::vector<int> actualVector(monomers[0].size(), 0);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < actualVector.size(); j++) {
                actualVector[j] += coeff[i] * monomers[i][j];
            }
        }
        return actualVector;
    }

    static bool isComplementary(std::vector<std::vector<int>> monomers, std::vector<int> coeff1, std::vector<int> coeff2) {
        if (coeff1.size() != coeff2.size()) {
            throw std::invalid_argument("Vectors must be of the same size for complement check.");
        }
        std::vector<int> v1 = coeffToVector(monomers, coeff1);
        std::vector<int> v2 = coeffToVector(monomers, coeff2);
        for (size_t i = 0; i < v1.size(); ++i) {
            if (v1[i]*v2[i] < 0) return true;
        }
        return false;
    }

    static bool isUnsplittable (std::vector<int> v, std::vector<std::vector<int>> monomers) {
        std::vector<int> weights = calculateWeights(v);
        int intValue = vectorToInteger(v, weights);
        int halfValue = intValue / 2;
        for (int  i = 1; i < halfValue; i++) {
            int iComplement = intValue - i;
            std::vector<int> iVectorCoeff = integerToVector(i, weights);
            std::vector<int> iComplementVectorCoeff = integerToVector(iComplement, weights);
            std::vector<int> iVector = coeffToVector(monomers, iVectorCoeff);
            std::vector<int> iComplementVector = coeffToVector(monomers, iComplementVectorCoeff);
            if (!isComplementary(monomers, iVector, iComplementVector)) {
                return false; // Found a uncomplementary pair
            }
        }
        return true;
    }

    static std::vector<int> vectorAdd (std::vector<int> v1, std::vector<int> v2) {
        std::vector<int> result(v1.size());
        for (size_t i = 0; i < v1.size(); ++i) {
            result[i] = v1[i] + v2[i];
        }
        return result;
    }
};

int main() {
    std::vector<std::vector<int>> monomers = {{2, 2, 2, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 2, 2, 2, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, -1, 2, 2, 2, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 2, 2, 2, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2},
    {0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-1, -1, -1, 0, 0, 0, -1, -1, -1, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    std::vector<std::vector<int>> basis;
    basis.reserve(monomers.size());

    std::vector<std::vector<int>> S;
    for (int i = 0; i < monomers.size(); i++) {
        std::vector<int> unitVector(monomers.size(), 0);
        unitVector[i] = 1;

        S.push_back(unitVector);
    }
    
    int i = 1;
    while (i < S.size()) {
        for (int j = 0; j < i; j++) {
            if(naiveAlgorithm::isComplementary(monomers, S[i], S[j])) {
                std::vector<int> p = naiveAlgorithm::vectorAdd(S[i], S[j]);
                if ((std::accumulate(p.begin(), p.end(), 0) <= MAX_NORM) && naiveAlgorithm::isUnsplittable(p, monomers)) {
                    S.push_back(p);
                }
            }
        }
        i++;
    }

    // Print S
    for (const auto& solution : S) {
        std::cout << "(";
        for (size_t i = 0; i < solution.size(); i++) {
            std::cout << solution[i];
            if (i < solution.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;
    }
}

