#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <numeric>
#include <set>

// Type definitions for clarity
using Monomer = std::vector<int>;        // Vector of binding sites
using Polymer = std::vector<int>;        // Vector of monomer counts

// Hash function for polymers
struct PolymerHash {
    size_t operator()(const Polymer& p) const {
        size_t hash = p.size();
        for (auto& i : p) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

class UnsplittablePolymers {
private:
    std::vector<Monomer> monomers;   // List of monomer types (binding site vectors)
    int numMonomers;                // Number of monomer types
    int numBindingSites;            // Number of binding site types

    // Calculate B(p) - the net binding site vector for polymer p
    std::vector<int> calculateNetBinding(const Polymer& p) {
        std::vector<int> result(numBindingSites, 0);
        
        // Sum the binding site vectors of all monomers
        for (int i = 0; i < numMonomers; i++) {
            if (p[i] > 0) {
                for (int j = 0; j < numBindingSites; j++) {
                    result[j] += p[i] * monomers[i][j];
                }
            }
        }
        return result;
    }
    
    // Check if two polymers are complementary
    bool areComplementary(const Polymer& p1, const Polymer& p2) {
        auto bp1 = calculateNetBinding(p1);
        auto bp2 = calculateNetBinding(p2);
        
        for (int i = 0; i < numBindingSites; i++) {
            if (bp1[i] * bp2[i] < 0) {
                return true;  // Found a complementary binding site
            }
        }
        return false;
    }
    
    // Add two polymers
    Polymer addPolymers(const Polymer& p1, const Polymer& p2) {
        Polymer result(numMonomers);
        for (int i = 0; i < numMonomers; i++) {
            result[i] = p1[i] + p2[i];
        }
        return result;
    }
    
    // Check if a polymer is splittable
    bool isSplittable(const Polymer& p) {
        // Generate all possible sub-multisets of p
        return generateAndCheckSubsets(p, Polymer(numMonomers, 0), 0);
    }
    
    // Recursive helper to generate all sub-multisets and check if any makes p splittable
    bool generateAndCheckSubsets(const Polymer& p, Polymer current, int index) {
        if (index == numMonomers) {
            // Check if current is a proper subset (not empty and not p itself)
            bool isProper = false;
            bool isNotEmpty = false;
            
            for (int i = 0; i < numMonomers; i++) {
                if (current[i] < p[i]) isProper = true;
                if (current[i] > 0) isNotEmpty = true;
            }
            
            if (isProper && isNotEmpty) {
                // Calculate p - current
                Polymer complement(numMonomers);
                for (int i = 0; i < numMonomers; i++) {
                    complement[i] = p[i] - current[i];
                }
                
                // If current and complement are not complementary, p is splittable
                if (!areComplementary(current, complement)) {
                    return true;  // Found a way to split p
                }
            }
            return false;
        }
        
        // Try all possible counts for this monomer type (0 to p[index])
        for (int count = 0; count <= p[index]; count++) {
            current[index] = count;
            if (generateAndCheckSubsets(p, current, index + 1)) {
                return true;
            }
        }
        
        return false;
    }
    
    // Print a polymer
    void printPolymer(const Polymer& p) {
        std::cout << "(";
        for (int i = 0; i < numMonomers; i++) {
            std::cout << p[i];
            if (i < numMonomers - 1) std::cout << ", ";
        }
        std::cout << ")";
    }

public:
    UnsplittablePolymers(const std::vector<Monomer>& m) 
        : monomers(m), numMonomers(m.size()), numBindingSites(m[0].size()) {}
    
    std::vector<Polymer> enumerate() {
        std::vector<Polymer> unsplittable;
        
        // Start with single monomers (trivially unsplittable)
        for (int i = 0; i < numMonomers; i++) {
            Polymer p(numMonomers, 0);
            p[i] = 1;
            unsplittable.push_back(p);
            std::cout << "Added single monomer: ";
            printPolymer(p);
            std::cout << std::endl;
        }
        
        // Iterate through the list, trying to merge polymers
        int i = 0;
        while (i < unsplittable.size()) {
            std::cout << "Processing polymer " << i << ": ";
            printPolymer(unsplittable[i]);
            std::cout << std::endl;
            
            for (int j = 0; j <= i; j++) {
                if (areComplementary(unsplittable[i], unsplittable[j])) {
                    // Create merged polymer
                    Polymer merged = addPolymers(unsplittable[i], unsplittable[j]);
                    
                    // Check if it's unsplittable
                    if (!isSplittable(merged)) {
                        // Check if we've already found this polymer
                        bool found = false;
                        for (const auto& p : unsplittable) {
                            if (p == merged) {
                                found = true;
                                break;
                            }
                        }
                        
                        if (!found) {
                            unsplittable.push_back(merged);
                            std::cout << "Found new unsplittable: ";
                            printPolymer(merged);
                            std::cout << " (from ";
                            printPolymer(unsplittable[i]);
                            std::cout << " + ";
                            printPolymer(unsplittable[j]);
                            std::cout << ")" << std::endl;
                        }
                    }
                }
            }
            
            i++;
            if (i % 10 == 0) {
                std::cout << "Progress: " << i << " / " << unsplittable.size() 
                          << " polymers processed" << std::endl;
            }
        }
        
        return unsplittable;
    }
    
    // Improved version with optimizations
    std::vector<Polymer> enumerateEfficient() {
        std::vector<Polymer> unsplittable;
        std::unordered_set<Polymer, PolymerHash> unsplittableSet;
        
        // Start with single monomers
        for (int i = 0; i < numMonomers; i++) {
            Polymer p(numMonomers, 0);
            p[i] = 1;
            unsplittable.push_back(p);
            unsplittableSet.insert(p);
        }
        
        int i = 0;
        int maxSize = 5; // Limit the size of polymers for efficiency
        
        while (i < unsplittable.size()) {
            int totalMonomers = std::accumulate(unsplittable[i].begin(), unsplittable[i].end(), 0);
            if (totalMonomers >= maxSize) {
                i++;
                continue; // Skip large polymers
            }
            
            for (int j = 0; j <= i; j++) {
                int j_size = std::accumulate(unsplittable[j].begin(), unsplittable[j].end(), 0);
                
                // Skip if combined size would exceed limit
                if (totalMonomers + j_size > maxSize) continue;
                
                if (areComplementary(unsplittable[i], unsplittable[j])) {
                    Polymer merged = addPolymers(unsplittable[i], unsplittable[j]);
                    
                    if (!isSplittable(merged) && unsplittableSet.find(merged) == unsplittableSet.end()) {
                        unsplittable.push_back(merged);
                        unsplittableSet.insert(merged);
                    }
                }
            }
            
            i++;
            if (i % 100 == 0) {
                std::cout << "Progress: " << i << " / " << unsplittable.size() << std::endl;
            }
        }
        
        return unsplittable;
    }
};

int main() {
    // Example: Define some monomers with binding sites
    // Monomer 1: {1 a, 0 b, 0 c}
    // Monomer 2: {0 a, 1 b, 0 c}
    // Monomer 3: {-1 a, 0 b, 1 c}
    std::vector<Monomer> monomers = {
        {-1, -1}, {1, 3}, {2, -2}, {-3, -1}
    };
    
    auto start = std::chrono::high_resolution_clock::now();
    
    UnsplittablePolymers solver(monomers);
    std::vector<Polymer> result = solver.enumerateEfficient();
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "\nUnsplittable Polymers:" << std::endl;
    for (const auto& p : result) {
        std::cout << "(";
        for (int i = 0; i < p.size(); i++) {
            std::cout << p[i];
            if (i < p.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;
    }
    
    std::cout << "\nTotal unsplittable polymers found: " << result.size() << std::endl;
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    
    return 0;
}