#include <vector>
#include <numeric>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>
#include <set>

#define MAX_NORM 100
#define DEBUG 0
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

    
    static std::vector<int> vectorAdd (std::vector<int> v1, std::vector<int> v2) {
        std::vector<int> result(v1.size());
        for (size_t i = 0; i < v1.size(); ++i) {
            result[i] = v1[i] + v2[i];
        }
        return result;
    }
    
    static std::vector<int> vectorSub (std::vector<int> v1, std::vector<int> v2) {
        std::vector<int> result(v1.size());
        for (size_t i = 0; i < v1.size(); ++i) {
            result[i] = v1[i] - v2[i];
        }
        return result;
    }
    
    static bool is_lex_leq(std::vector<int> a, std::vector<int> b) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] < b[i]) return true;
            if (a[i] > b[i]) return false;
        }
        return true; // a == b
    }
    
    static void printVector(const std::vector<int>& v) {
        if (DEBUG) {
            std::cout << "(";
            for (size_t i = 0; i < v.size(); ++i) {
                std::cout << v[i];
                if (i < v.size() - 1) std::cout << ", ";
            }
            std::cout << ")" << std::endl;
        }
    }
    
    static std::vector<int> coeffToVector(std::vector<std::vector<int>> monomers, std::vector<int> coeff) {
        if (monomers.size() != coeff.size()) {
            std::cout << "Monomers size: " << monomers.size() << ", Coefficients size: " << coeff.size() << std::endl;
            std::cout << "Monomers: ";
            for (const auto& monomer : monomers) {
                printVector(monomer);
            }
            std::cout << "Coefficients: ";
            printVector(coeff);
            throw std::invalid_argument("Monomers and coefficients must be of the same size.");
        }

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
        // Check if the coefficients are complementary
        std::vector<int> v1 = coeffToVector(monomers, coeff1);
        std::vector<int> v2 = coeffToVector(monomers, coeff2);
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
            printVector(v);
        }
    
        while (true) {
            // Skip the iteration if b is a null vector
            if (!(std::all_of(b.begin(), b.end(), [](int x) { return x == 0; }))) {
                // Compute c = a - b
                std::vector<int> c(N);
                c = vectorSub(v, b);
    
                // Only do if b is lex ≤ c
                if (is_lex_leq(b, c)) {
                    std::vector<int> iVector = coeffToVector(monomers, b);
                    std::vector<int> iComplementVector = coeffToVector(monomers, c);
                    if (!isComplementary(monomers, b, c)) {
                        if (DEBUG) {
                            std::cout << "Found uncomplementary pair. Polymer is splittable." << std::endl;
                            std::cout << "b: ";
                            printVector(b);
                            std::cout << "c: ";
                            printVector(c);
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
        if (DEBUG) {
        std::cout << "Polymer is unsplittable." << std::endl;
        }
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
            std::vector<int> p = naiveAlgorithm::vectorAdd(*it1, *it2);
            if ((std::accumulate(p.begin(), p.end(), 0) <= MAX_NORM) && naiveAlgorithm::isUnsplittable(p, monomers)) {
                if (DEBUG) {
                    std::cout << "Adding polymer to S: ";
                    naiveAlgorithm::printVector(p);
                    std::cout << std::endl;
                }
                S.insert(p);
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