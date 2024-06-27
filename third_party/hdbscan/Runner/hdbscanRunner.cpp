#include "hdbscanRunner.hpp"
#include "hdbscanResult.hpp"
#include "hdbscanParameters.hpp"

#include "../HdbscanStar/hdbscanAlgorithm.hpp"
#include "../HdbscanStar/undirectedGraph.hpp"
#include "../HdbscanStar/cluster.hpp"
#include "../HdbscanStar/outlierScore.hpp"

using namespace hdbscanStar;
using namespace std;

hdbscanResult hdbscanRunner::run(hdbscanParameters parameters)
{
    int numPoints = parameters.dataset.size() != 0 ? parameters.dataset.size() : parameters.distances.size();

    hdbscanAlgorithm algorithm;
    hdbscanResult result;
    if (parameters.distances.size() == 0)
    {
        std::vector<std::vector<double>> distances(numPoints);
        for (int i = 0; i < numPoints; i++)
        {
            distances[i].resize(numPoints);

            for (int j = 0; j < i; j++)
            {
                if (parameters.distanceFunction.length() == 0)
                {
                    double distance = EuclideanDistance(parameters.dataset[i], parameters.dataset[j]);
                    distances[i][j] = distance;
                    distances[j][i] = distance;

                }
                else if (parameters.distanceFunction == "Euclidean")
                {
                    double distance = EuclideanDistance(parameters.dataset[i], parameters.dataset[j]);
                    distances[i][j] = distance;
                    distances[j][i] = distance;
                }
                else if (parameters.distanceFunction == "Manhattan")
                {
                    double distance = ManhattanDistance(parameters.dataset[i], parameters.dataset[j]);
                    distances[i][j] = distance;
                    distances[j][i] = distance;
                }
                else if (parameters.distanceFunction == "Minkowski")
                {
                    double distance = MinkowskiDistance(parameters.dataset[i], parameters.dataset[j], parameters.MinkowskiP);
                    distances[i][j] = distance;
                    distances[j][i] = distance;
                }
            }
        }

        parameters.distances = distances;
    }

    std::vector <double> coreDistances = algorithm.calculateCoreDistances(
            parameters.distances,
            parameters.minPoints);

    undirectedGraph mst = algorithm.constructMst(
                              parameters.distances,
                              coreDistances,
                              true);
    mst.quicksortByEdgeWeight();

    std::vector<double> pointNoiseLevels(numPoints);
    std::vector<int> pointLastClusters(numPoints);

    std::vector< std::vector <int> > hierarchy;

    std::vector<cluster*> clusters;
    algorithm.computeHierarchyAndClusterTree(
        &mst,
        parameters.minClusterSize,
        parameters.constraints,
        hierarchy,
        pointNoiseLevels,
        pointLastClusters,
        clusters);
    bool infiniteStability = algorithm.propagateTree(clusters);

    std::vector<int> prominentClusters = algorithm.findProminentClusters(clusters, hierarchy, numPoints);
    std::vector<double> membershipProbabilities = algorithm.findMembershipScore(prominentClusters, coreDistances);
    std::vector<outlierScore> scores = algorithm.calculateOutlierScores(
                                           clusters,
                                           pointNoiseLevels,
                                           pointLastClusters,
                                           coreDistances);

    return hdbscanResult(prominentClusters, scores, membershipProbabilities,  infiniteStability);
}
