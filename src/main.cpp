#include <cstring>
#ifndef SEQAN_ENABLE_PARALLELISM
#define SEQAN_ENABLE_PARALLELISM 1
#endif

#ifndef SEQAN_BGZF_NUM_THREADS
#define SEQAN_BGZF_NUM_THREADS 16
#endif

#ifndef VERSION
#define VERSION "0.0"
#endif

#ifndef DATE
#define DATE "1.1.1970"
#endif

#include <eigen3/Eigen/Core>
#include "options.hpp"
#include <chrono>
#include <omp.h>
#include <unistd.h>
#include "profileSamples.hpp"
#include "profileVariants.hpp"
#include "genotype.hpp"

void printHelp()
{
    std::cerr << "GGTyper - Genotyping of complex structural variants" << std::endl;
    std::cerr << "===================================================" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mSYNOPSIS\033[0m" << std::endl;
    std::cerr << "    \033[1m" << "ggtyper" << " COMMAND\033[0m [\033[4mARGUMENTS\033[0m] [\033[4mOPTIONS\033[0m]" << std::endl;
    std::cerr << std::endl;

    std::cerr << "\033[1mCOMMANDS\033[0m" << std::endl;
    std::cerr << "    \033[1mprofile-samples\033[0m       Create profiles of given bam files and write them to disk." << std::endl;
    std::cerr << "    \033[1mprofile-variants\033[0m      Create profiles of given variants and write them to disk." << std::endl;
    std::cerr << "    \033[1mgenotype\033[0m              Genotype given variants (specified by profiles) in all samples (specified by profiles)." << std::endl;
    std::cerr << std::endl;

    std::cerr << "\033[1mVERSION\033[0m" << std::endl;
    std::cerr << "    " << "GGTyper" << " version: " << VERSION << std::endl;
    std::cerr << "    Last update " << DATE << std::endl;
    std::cerr << "    Contact: Tim Mirus (Tim.Mirus[at]ukr.de)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Try 'ggtyper COMMAND --help' for more information on each command." << std::endl;
    std::cerr << std::endl;
}

int main(int argc, char const ** argv) {
    
    Eigen::initParallel();

    if (argc < 2)
    {
        std::cerr << "Incorrect call. Try 'ggtyper -h' for help." << std::endl;
        return 1;
    }

    const char * command = argv[1];
    std::string call = std::string(argv[0]) + " " + std::string(argv[1]);
    argv[1] = call.c_str();
    ++argv;
    --argc;

    if (strcmp(command, "genotype") == 0)
    {
        return genotype (argc, argv);
    } 
    else if (strcmp(command, "profile-samples") == 0)
    {
        return profileSamples (argc, argv);
    }
    else if (strcmp(command, "profile-variants") == 0)
    {
       return profileVariants (argc, argv);
    }
    else if (strcmp(command, "--help") == 0 || strcmp(command, "-h") == 0)
    {
        printHelp();
        return 0;
    } else {
        std::cerr << "Incorrect call. Try 'ggtyper -h' for help." << std::endl;
        return 1; 
    }

    return 0;
}
