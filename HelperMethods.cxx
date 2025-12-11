#include "HelperMethods.hxx"


std::vector<std::vector<int>> HelperMethods::parseMonomersFile(std::string filename) {
    std::vector<std::vector<int>> monomers;
    std::ifstream infile(filename);

    if (!infile) {
        throw std::runtime_error("Error: Unable to open file " + filename);
    }

    std::string line;
    std::regex numberRegex("-?\\d+");
    size_t expectedDimension = 0; // To store the dimension of the first vector

    while (std::getline(infile, line)) {
        std::vector<int> monomer;
        std::sregex_iterator it(line.begin(), line.end(), numberRegex);
        std::sregex_iterator end;

        // Extract integers from the line using regex
        while (it != end) {
            monomer.push_back(std::stoi(it->str()));
            ++it;
        }

        if (!monomer.empty()) {
            // Check dimensions
            if (monomers.empty()) {
                expectedDimension = monomer.size();
            } else if (monomer.size() != expectedDimension) {
                throw std::runtime_error("Error: Inconsistent vector dimensions in file " + filename +
                                         ". Expected dimension: " + std::to_string(expectedDimension) +
                                         ", but found: " + std::to_string(monomer.size()));
            }

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

std::vector<std::vector<int>> HelperMethods::add_unit_monomers(std::vector<std::vector<int>> monomers) {
    std::vector<std::vector<int>> result = monomers;

    int dimension = monomers.empty() ? 0 : monomers[0].size();

    for (int i = 0; i < dimension; ++i) {
        std::vector<int> unitVector(dimension, 0);
        unitVector[i] = 1;

        if (std::find(result.begin(), result.end(), unitVector) == result.end()) {
            result.push_back(unitVector);
        }

        unitVector[i] = -1;

        if (std::find(result.begin(), result.end(), unitVector) == result.end()) {
            result.push_back(unitVector);
        }
    }

    return result;
}

std::vector<std::vector<int>> HelperMethods::remove_unit_monomers(std::vector<std::vector<int>> basis, int n_dim) {
    std::vector<std::vector<int>> trimmedBasis;
    printf("dimension of basis: %zu\n", basis.size());
    printf("n_dim: %d\n", n_dim);

    // Trim each vector to n_dim
    for (auto& vec : basis) {
        if (vec.size() > n_dim) {
            vec.resize(n_dim);
        }
        trimmedBasis.push_back(vec);
    }

    // Remove duplicates
    std::set<std::vector<int>> uniqueBasis(trimmedBasis.begin(), trimmedBasis.end());
    trimmedBasis.assign(uniqueBasis.begin(), uniqueBasis.end());

    return trimmedBasis;
}

std::vector<int> HelperMethods::vectorAdd (std::vector<int> v1, std::vector<int> v2) {
    std::vector<int> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

std::vector<int> HelperMethods::vectorSub (std::vector<int> v1, std::vector<int> v2) {
    std::vector<int> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

std::vector<int> HelperMethods::vectorNegative (std::vector<int> v) {
    std::vector<int> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = -v[i];
    }
    return result;
}

bool HelperMethods::is_lex_leq(std::vector<int> a, std::vector<int> b) {
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] < b[i]) return true;
        if (a[i] > b[i]) return false;
    }
    return true; // a == b
}

void HelperMethods::printVector(std::vector<int>& v) {
        std::cout << "(";
        for (size_t i = 0; i < v.size(); ++i) {
            std::cout << v[i];
            if (i < v.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;
}

std::vector<int> HelperMethods::coeffToVector(std::vector<std::vector<int>> monomers, std::vector<int> coeff) {
    if (monomers.size() != coeff.size()) {
        std::cout << "Monomers size: " << monomers.size() << ", Coefficients size: " << coeff.size() << std::endl;
        std::cout << "Monomers: ";
        for (auto& monomer : monomers) {
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
