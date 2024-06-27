#pragma once
#include "hdbscanResult.hpp"
#include "hdbscanParameters.hpp"
#include <cmath>

class hdbscanRunner
{
public:
    static hdbscanResult run(hdbscanParameters parameters);
};

using namespace std;

inline double ManhattanDistance(std::vector<double> attributesOne, std::vector<double> attributesTwo)
{
    double distance = 0;
    for (size_t i = 0; i < attributesOne.size() && i < attributesTwo.size(); i++)
    {
        distance += std::abs(attributesOne[i] - attributesTwo[i]);
    }

    return distance;
}

inline double EuclideanDistance(std::vector<double> attributesOne, std::vector<double> attributesTwo)
{
    double distance = 0;
    for (size_t i = 0; i < attributesOne.size() && i < attributesTwo.size(); i++)
    {
        distance += std::pow(attributesOne[i] - attributesTwo[i], 2);
    }

    return std::sqrt(distance);
}

inline double MinkowskiDistance(std::vector<double> attributesOne, std::vector<double> attributesTwo, double p)
{
    double distance = 0;
    for (size_t i = 0; i < attributesOne.size() && i < attributesTwo.size(); i++)
    {
        distance += std::pow(std::abs(attributesOne[i] - attributesTwo[i]), p);
    }

    distance = std::pow(distance, 1.0/p);

    return distance;
}
