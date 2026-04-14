//
// Created by Jakob Goldmann on 14.04.25.
//

#ifndef GLOBAL_STORAGE_H
#define GLOBAL_STORAGE_H
#include <map>
#include <string>


// This serves as some "global" / static cache for various parts of the whole program. This implementation is a bit
// "dirty" but in this case since the mombf code is so convoluted and complex, making the caches available globally
// is the easiest option.

// This is a cache for storing the "final" th/thopt values calculated (for example in imomModeK) and re-use them for
// faster convergence.
thread_local inline static std::map<std::string, std::vector<double> > model_thopt_mapping;
thread_local inline std::string current_model;

enum PolyRootAlgo {
    JENKINS_TRAUB,
    NEWTON_MADSEN,
    ABERTH_EHRLICH
};

constexpr PolyRootAlgo root_finding_algo = PolyRootAlgo::NEWTON_MADSEN;

#endif //GLOBAL_STORAGE_H
