#ifndef SEQAN_ENABLE_PARALLELISM
#define SEQAN_ENABLE_PARALLELISM 1
#endif

#ifndef SEQAN_BGZF_NUM_THREADS
#define SEQAN_BGZF_NUM_THREADS 16
#endif

#include <eigen3/Eigen/Core>
#include "parser.hpp"
#include "options.hpp"
#include "genotyper.hpp"
#include "seqFileHandler.hpp"
#include <chrono>
#include <omp.h>
#include <unistd.h>
#include "variantProfile.hpp"

int main(int argc, char const ** argv) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    Eigen::initParallel();

    parser argParser(argc, argv);
    if(!argParser.wasSuccessful())
        return 0;

    ProgramOptions options = argParser.getOptions(); 
    omp_set_num_threads(options.getNumberOfThreads());

    Genotyper genotyper(options);
    
    genotyper.genotypeAllSamples();
    
    genotyper.writeResults();
    genotyper.writeStats();
    genotyper.writeToVCF();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    if (options.isOptionProfile())
        std::cout << "Finished after: " << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << "s" << std::endl;

    return 0;
}
