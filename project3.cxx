#include <vector>
#include <numeric>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <sstream>

#define MAX_NORM 3
class naiveAlgorithm {
public:
    static std::vector<std::vector<int>> parseMonomersFile(const std::string& filename) {
        std::vector<std::vector<int>> monomers;
        std::ifstream infile(filename);

        if (!infile) {
            throw std::runtime_error("Error: Unable to open file " + filename);
        }

        std::string line;
        while (std::getline(infile, line)) {
            std::vector<int> monomer;
            std::istringstream iss(line);
            int value;

            // Parse integers from the line
            while (iss >> value) {
                monomer.push_back(value);
            }

            if (!monomer.empty()) {
                monomers.push_back(monomer);
            }
        }

        infile.close();
        std::cout << "Parsed monomers from file: " << filename << std::endl;
        for (size_t i = 0; i < monomers.size(); ++i) {
            std::cout << "Monomer " << i << ": ";
            for (size_t j = 0; j < monomers[i].size(); ++j) {
                std::cout << monomers[i][j] << " ";
            }
            std::cout << std::endl;
        }
        return monomers;
    }

    static std::vector<int> coeffToVector(std::vector<std::vector<int>> monomers, std::vector<int> coeff) {
        int size = coeff.size();
        if (size != monomers.size()) {
            std::cerr << "Error: Coefficient size (" << size << ") does not match monomer size (" << monomers.size() << ")." << std::endl;
            std::cerr << "Coefficient vector: ";
            for (size_t i = 0; i < coeff.size(); ++i) {
                std::cerr << coeff[i] << " ";
            }
            std::cerr << std::endl;
    
            std::cerr << "Monomers vectors:" << std::endl;
            for (size_t i = 0; i < monomers.size(); ++i) {
                std::cerr << "Monomer " << i << ": ";
                for (size_t j = 0; j < monomers[i].size(); ++j) {
                    std::cerr << monomers[i][j] << " ";
                }
                std::cerr << std::endl;
            }
    
            throw std::invalid_argument("Coefficient size does not match monomer size.");
        }

        std::vector<int> actualVector(monomers.size(), 0);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < actualVector.size(); j++) {
                actualVector[j] += coeff[i] * monomers[i][j];
            }
        }
        std::cout << "Converted coefficients to vector: ";
        for (size_t i = 0; i < actualVector.size(); ++i) {
            std::cout << actualVector[i] << " ";
        }
        std::cout << std::endl;
        return actualVector;
    }

    static bool isComplementary(std::vector<std::vector<int>> monomers, std::vector<int> coeff1, std::vector<int> coeff2) {
        if (coeff1.size() != coeff2.size()) {
            throw std::invalid_argument("Vectors must be of the same size for complement check.");
        }
        std::vector<int> v1 = coeffToVector(monomers, coeff1);
        std::vector<int> v2 = coeffToVector(monomers, coeff2);
        for (size_t i = 0; i < v1.size(); ++i) {
            if (v1[i]*v2[i] < 0) {
                std::cout << "Vectors are complementary." << std::endl;
                return true;
            }
        }
        std::cout << "Vectors are not complementary." << std::endl;
        return false;
    }

    static std::vector<int> vectorAdd (std::vector<int> v1, std::vector<int> v2) {
        std::vector<int> result(v1.size());
        for (size_t i = 0; i < v1.size(); ++i) {
            result[i] = v1[i] + v2[i];
        }
        std::cout << "Added vectors: ";
        for (size_t i = 0; i < result.size(); ++i) {
            std::cout << result[i] << " ";
        }
        std::cout << std::endl;
        return result;
    }

    static std::vector<int> vectorSub (std::vector<int> v1, std::vector<int> v2) {
        std::vector<int> result(v1.size());
        for (size_t i = 0; i < v1.size(); ++i) {
            result[i] = v1[i] - v2[i];
        }
        std::cout << "Subtracted vectors: ";
        for (size_t i = 0; i < result.size(); ++i) {
            std::cout << result[i] << " ";
        }
        std::cout << std::endl;
        return result;
    }

    static bool is_lex_leq(std::vector<int> a, std::vector<int> b) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] < b[i]) return true;
            if (a[i] > b[i]) return false;
        }
        return true; // a == b
    }

    static bool isUnsplittable (std::vector<int> v, std::vector<std::vector<int>> monomers) {
        int N = v.size();
        std::vector<int> b(N, 0);

        while (true) {
            // Compute c = a - b
            std::vector<int> c(N);
            c = vectorSub(v, b);

            // Only do if b is lex ≤ c
            if (is_lex_leq(b, c)) {
                std::vector<int> iVector = coeffToVector(monomers, b);
                std::vector<int> iComplementVector = coeffToVector(monomers, c);
                if (!isComplementary(monomers, iVector, iComplementVector)) {
                    std::cout << "Found uncomplementary pair. Polymer is splittable." << std::endl;
                    return false; // Found a uncomplementary pair
                }
            }

            int i = N - 1;
            while (i >= 0) {
                b[i]++;
                if (b[i] <= v[i]) break;
                b[i] = 0;
                i--;
            }
            if (i < 0) break; // Done
        }
        std::cout << "Polymer is unsplittable." << std::endl;
        return true;
    }

};

int main(int argc, char* argv[]) {
    std::vector<std::vector<int>> monomers;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    try {
        monomers = naiveAlgorithm::parseMonomersFile(argv[1]);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

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
            if (naiveAlgorithm::isComplementary(monomers, S[i], S[j])) {
                std::vector<int> p = naiveAlgorithm::vectorAdd(S[i], S[j]);
                if ((std::accumulate(p.begin(), p.end(), 0) <= MAX_NORM) && naiveAlgorithm::isUnsplittable(p, monomers)) {
                    std::cout << "Adding polymer to S: ";
                    for (size_t k = 0; k < p.size(); ++k) {
                        std::cout << p[k] << " ";
                    }
                    std::cout << std::endl;
                    S.push_back(p);
                }
            }
        }
        i++;
    }

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

    return 0;
}