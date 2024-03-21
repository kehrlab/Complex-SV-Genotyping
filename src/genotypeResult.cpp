#include "genotypeResult.hpp"

GenotypeResult::GenotypeResult()
{
    this->filename = "";
    this->sampleName = "";
    this->observedReads = 0;
    this->callCertainty = 0;
    this->lowerBoundQuality = 0;
    this->upperBoundQuality = 0;
    this->minQuality = 0;
    this->maxQuality = 0;
    this->meanQuality = 0.0;
    this->outlierCount = 0;
    this->useQualities = false;
}

GenotypeResult::GenotypeResult(std::string filename, bool useQualities)
{
    this->filename = filename;
    this->sampleName = "unknown";
    this->observedReads = 0;
    this->callCertainty = 0;
    this->lowerBoundQuality = 0;
    this->upperBoundQuality = 0;
    this->minQuality = 60;
    this->maxQuality = 0;
    this->meanQuality = 0.0;
    this->outlierCount = 0;
    this->useQualities = useQualities;
}

GenotypeResult::GenotypeResult(std::string filename, std::string sampleName, bool useQualities)
{
    this->filename = filename;
    this->sampleName = sampleName;
    this->observedReads = 0;
    this->callCertainty = 0;
    this->lowerBoundQuality = 0;
    this->upperBoundQuality = 0;
    this->minQuality = 60;
    this->maxQuality = 0;
    this->meanQuality = 0.0;
    this->outlierCount = 0;
    this->useQualities = useQualities;
}

GenotypeResult::GenotypeResult(std::string filename, std::string sampleName, bool useQualities, std::unordered_map<std::string, float> priors)
{
    this->filename = filename;
    this->sampleName = sampleName;
    this->observedReads = 0;
    this->callCertainty = 0;
    this->lowerBoundQuality = 0;
    this->upperBoundQuality = 0;
    this->minQuality = 60;
    this->maxQuality = 0;
    this->meanQuality = 0.0;
    this->outlierCount = 0;
    this->useQualities = useQualities;
    this->genotypePriors = priors;
}


void GenotypeResult::initLikelihoods(std::vector<std::string> genotypeNames)
{
    this->genotypeNames = genotypeNames;
    for (auto gtName : genotypeNames) {
        std::vector<float> emptyVec;
        this->templateProbabilities.push_back(emptyVec);
    }
    
    // prepare genotype priors
    bool fittingPriors {true};
    if (this->genotypePriors.size() != this->genotypeNames.size())
    {
        fittingPriors = false; 
    } else {
        for (auto & gt : this->genotypeNames)
        {
            if (this->genotypePriors.find(gt) == this->genotypePriors.end())
            {
                fittingPriors = false;
                break;
            }
        }
    }
    // create equal priors if the given priors do not fit
    if(!fittingPriors)
    {
        this->genotypePriors = std::unordered_map<std::string, float>();
        for (auto & gt : this->genotypeNames)
            this->genotypePriors[gt] = 1.;   
    }

    // rescale to 1
    float sum = 0.;
    for (auto & p : this->genotypePriors)
        sum += p.second;
    for (auto & p : this->genotypePriors)
        p.second /= sum;

    return;
}

void GenotypeResult::callGenotype()
{
    // call genotype and quality based on R
    calculateLikelihoods();

    if (this->genotypeLikelihoods.size() < 3)
        throw std::runtime_error("ERROR: There must be at least 3 different genotypes.");

    std::vector<float> likelihoods = this->genotypeLikelihoods;
    std::sort(likelihoods.begin(), likelihoods.end());
    int minIdx = std::find(this->genotypeLikelihoods.begin(), this->genotypeLikelihoods.end(), likelihoods[0]) - this->genotypeLikelihoods.begin();
    int secondIdx = std::find(this->genotypeLikelihoods.begin(), this->genotypeLikelihoods.end(), likelihoods[1]) - this->genotypeLikelihoods.begin();
    float lMin = this->genotypeLikelihoods[minIdx];
    float lSecond = this->genotypeLikelihoods[secondIdx];
    this->calledGenotype = this->genotypeNames[minIdx];
    this->quality = lSecond - lMin;

    // create bootstrap samples and determine Q_k,i and Q(S) for each sample
    bootstrapLikelihoods();
    bootstrapQuality(minIdx, secondIdx);

    createOutputString();
}

std::string GenotypeResult::getCalledGenotype()
{
    return this->calledGenotype;
}

std::string GenotypeResult::getFilename()
{
    return this->filename;
}

void GenotypeResult::setFilename(std::string filename)
{
    this->filename = filename;
}

void GenotypeResult::addTemplateProbabilities(std::vector<std::string> genotypeNames, std::vector<float> probabilities, float weight)
{
    this->templateWeights.push_back(weight);
    for (uint32_t i = 0; i < genotypeNames.size(); ++i)
        addProbability(genotypeNames[i], probabilities[i]);
    return;
}

void GenotypeResult::addProbability(std::string genotypeName, float probability)
{
    int gtIdx = -1;
    for (uint32_t i = 0; i < this->genotypeNames.size(); ++i)
    {
        if (this->genotypeNames[i] == genotypeName)
        {
            gtIdx = (int) i;
            break;
        }
    }
    if (gtIdx != -1)
        this->templateProbabilities[gtIdx].push_back((-10) * std::log10(probability));
    else
        throw std::runtime_error("Did not find Genotype to adjust!");
    return;
}

void GenotypeResult::calculateLikelihoods()
{
    for (uint32_t i = 0; i < this->genotypeNames.size(); ++i)
    {
        float L = 0;
        for (uint32_t j = 0; j < this->templateProbabilities[i].size(); ++j)
            L += this->templateProbabilities[i][j];
        L += -10 * std::log10(this->genotypePriors[this->genotypeNames[i]]);
        this->genotypeLikelihoods.push_back(L);
    }
}

void GenotypeResult::bootstrapLikelihoods()
{
    std::random_device rd;
    std::mt19937 gen(rd());


    uint32_t n = this->templateProbabilities[0].size();
    uint32_t sampleSize = n;
    uint32_t nSamples = 10000;

    if (!this->useQualities)
        for (float & w : this->templateWeights)
            w = 1.0;
        
    std::discrete_distribution<> distrib(this->templateWeights.begin(), this->templateWeights.end());

    // init likelihood vector
    for (uint32_t i = 0; i < this->genotypeNames.size(); ++i)
    {
        std::vector<float> temp;
        this->bootstrappedLikelihoods.push_back(temp);   
    }
    
    float L = 0.0;
    for (uint32_t j = 0; j < nSamples; ++j)
    {
        // create bootstrap sample
        std::vector<int> indices(sampleSize, 0);
        for (uint32_t k = 0; k < sampleSize; ++k)
            indices[k] = distrib(gen);

        for (uint32_t i = 0; i < this->genotypeNames.size(); ++i)
        {
            L = 0.0;
            for (uint32_t k = 0; k < sampleSize; ++k) 
                L += this->templateProbabilities[i][indices[k]];
            this->bootstrappedLikelihoods[i].push_back(L);
        }
    }
}

void GenotypeResult::bootstrapQuality(int minIdx, int secondIdx)
{
    int nSameCall = 0;
    float q = 0;
    for (uint32_t i = 0; i < this->bootstrappedLikelihoods[0].size(); ++i)
    {
        // get min and second index within bootstrap sample
        q = this->bootstrappedLikelihoods[secondIdx][i] - this->bootstrappedLikelihoods[minIdx][i];
        this->bootstrappedQualities.push_back(q);
        
        // check wheter the calls are actually the same
        std::vector<float> likelihoods;
        for (uint32_t j = 0; j < this->bootstrappedLikelihoods.size(); ++j)
            likelihoods.push_back(this->bootstrappedLikelihoods[j][i]);
        std::vector<float> tempLikelihoods = likelihoods;
        std::sort(tempLikelihoods.begin(), tempLikelihoods.end());
        int minIndex = std::find(likelihoods.begin(), likelihoods.end(), tempLikelihoods[0]) - likelihoods.begin();
        if (minIndex == minIdx)
            ++nSameCall;
    }
    // determine certainty
    this->callCertainty = (float) nSameCall / (float) this->bootstrappedLikelihoods[minIdx].size();

    //get confidence interval
    std::vector<float> bQs = this->bootstrappedQualities;
    std::sort(bQs.begin(), bQs.end());
    int minPercentileIdx = (int) (0.025 * bQs.size()) - 1;
    int maxPercentileIdx = (int) (0.975 * bQs.size()) - 1;
    this->lowerBoundQuality = bQs[minPercentileIdx];
    this->upperBoundQuality = bQs[maxPercentileIdx];
}

void GenotypeResult::writeBootstrapData(std::string prefix)
{
    std::string bsFilename = prefix + "_bootstrap_data.txt";
    std::ofstream f(bsFilename);
    if (f.is_open())
    {
        for (auto gt : this->genotypeNames)
            f << gt << "\t";
        f << "Quality" << std::endl;

        for (uint32_t i = 0; i < this->bootstrappedLikelihoods[0].size(); ++i)
        {
            std::vector<float> v;
            float minL = this->bootstrappedLikelihoods[0][i];
            for (uint32_t j = 0; j < this->bootstrappedLikelihoods.size(); ++j) {
                float L = this->bootstrappedLikelihoods[j][i];
                v.push_back(L);
                minL =  std::min(minL, L);
            }
            for (uint32_t j = 0; j < v.size(); ++j)
                f << (v[j] - minL) << "\t";
            f << this->bootstrappedQualities[i] << std::endl;
        }
    }
    f.close();
}

float GenotypeResult::mean(std::vector<float> v)
{
    double m = 0;
    for (float x : v)
        m += x;
    m /= v.size();
    return (float) m;
}

float GenotypeResult::sd(std::vector<float> v, float m)
{
    double s = 0;
    for (float x : v)
        s += (x-m)*(x-m);
    s /= v.size();
    s = std::sqrt(s);
    return (float) s;
}

float GenotypeResult::getLikelihood(std::string genotype)
{
    int index = -1;
    for (uint32_t i = 0; i < this->genotypeNames.size(); ++i)
    {
        if (this->genotypeNames[i] == genotype)
        {
            index = (int) i;
            break;
        }
    }
    if (index == -1)
        throw std::runtime_error("Did not find Genotype to adjust!");
    return this->genotypeLikelihoods[index];
}

void GenotypeResult::printAllLikelihoods()
{
    std::cout << std::endl;
    for (uint32_t i = 0; i < this->genotypeNames.size(); ++i)
        std::cout << this->genotypeNames[i] << "\t";
    std::cout << std::endl;
    for (uint32_t i = 0; i < this->genotypeLikelihoods.size(); ++i)
        std::cout << this->genotypeLikelihoods[i] << "\t";
    std::cout << std::endl << std::endl;
}

void GenotypeResult::createOutputString()
{
    determineQualityStats();
    std::string outputString = "";
    outputString.append(this->sampleName);
    outputString.append("\t");
    outputString.append(this->filename);
    outputString.append("\t");
    outputString.append(this->calledGenotype);
    outputString.append("\t");
    outputString.append(std::to_string(this->quality));
    outputString.append("\t");
    outputString.append(std::to_string(this->lowerBoundQuality));
    outputString.append("\t");
    outputString.append(std::to_string(this->upperBoundQuality));
    outputString.append("\t");
    outputString.append(std::to_string(this->callCertainty));
    outputString.append("\t");
    outputString.append(std::to_string(this->observedReads));
    outputString.append("\t");
    outputString.append(std::to_string(this->outlierCount));
    outputString.append("\t");
    outputString.append(std::to_string(this->meanQuality));
    outputString.append("\t");
    outputString.append(std::to_string(this->minQuality));
    outputString.append("\t");
    outputString.append(std::to_string(this->maxQuality));
    outputString.append("\t");
    if (this->meanQuality >= 55 && this->callCertainty >= 0.9999)
	    outputString.append("TRUE");
    else
	    outputString.append("FALSE");
    // outputString.append("\n");
    this->outputString = outputString;
}

std::string GenotypeResult::getOutputString()
{
    return this->outputString;
}

void GenotypeResult::determineQualityStats()
{
    for (auto mapQ : this->mappingQualities)
    {
        if (mapQ < this->minQuality)
            this->minQuality = mapQ;
        if (mapQ > this->maxQuality)
            this->maxQuality = mapQ;
        this->meanQuality += (float) mapQ;
    }
    this->meanQuality /= this->mappingQualities.size();
}

void GenotypeResult::storeEvidence(int64_t insertSize, std::string orientation, std::string junctionString, std::string breakpointString, std::string bridgeString, std::vector<int> mapQs)
{
    ++ this->observedReads;
    for (int mapQ : mapQs)
        this->mappingQualities.push_back(mapQ);

    this->sampleDistribution.addInsertSizeProbability(insertSize, orientation, junctionString, breakpointString, bridgeString, 1.0);
}

void GenotypeResult::writeEvidence(std::string prefix)
{
    this->sampleDistribution.scaleDistribution();
    this->sampleDistribution.writeDistributionBinned(prefix);
}

float GenotypeResult::getQuality()
{
    return this->quality;
}

void GenotypeResult::clearData()
{
    this->sampleDistribution = GenotypeDistribution();
    std::vector<int>().swap(this->mappingQualities);
    std::vector<float>().swap(this->templateWeights);
    std::vector<float>().swap(this->bootstrappedQualities);
    std::vector<std::vector<float>>().swap(this->bootstrappedLikelihoods);
    std::vector<std::vector<float>>().swap(this->templateProbabilities);
    std::vector<std::string>().swap(this->genotypeNames);
    std::vector<float>().swap(genotypeLikelihoods);
    std::vector<float>().swap(genotypeLikelihoods);
}

std::string GenotypeResult::getSampleName()
{
    return this->sampleName;
}

void GenotypeResult::addOutlier(bool outlier)
{
    if (outlier)
        ++this->outlierCount;
    return;
}