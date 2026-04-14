#include "LassoRegression.h"
#include "matrix.h"
#include <cmath>
#include <iostream>
#include <fstream>

LassoRegression::LassoRegression(std::vector<std::vector<double> > samples, std::vector<double> target) {
    this->numberOfSamples  = samples.size();
    this->numberOfFeatures = samples[0].size();
    this->features         = featuresMatrix(samples);

    this->features = normalizeFeatures(this->features);
    this->weights  = initialWeights();
    this->target   = targetAsArray(target);
}

// Add destructor
LassoRegression::~LassoRegression() {
    // Delete features matrix
    for (int i = 0; i < numberOfSamples; i++) {
        delete[] features[i];
    }
    delete[] features;

    // Delete weights and target arrays
    delete[] weights;
    delete[] target;
}

double *LassoRegression::predictions() {
    double *result = new double[numberOfSamples];

    for (int sampleIdx = 0; sampleIdx < numberOfSamples; sampleIdx++) {
        double prediction = 0.0;
        for (int featureIdx = 0; featureIdx < numberOfFeatures; featureIdx++) {
            prediction += features[sampleIdx][featureIdx] * weights[featureIdx];
        }
        result[sampleIdx] = prediction;
    }

    return result;
}

double *LassoRegression::ro() {
    double *results = new double[numberOfFeatures];

    for (int idx = 0; idx < numberOfFeatures; idx++) {
        double *feature_val   = feature(idx);
        double *penaltyVector = vectorMultiply(feature_val, numberOfSamples, weights[idx]);
        delete[] feature_val;

        double *predictions_val = predictions();
        double *temp_multiply = vectorMultiply(predictions_val, numberOfSamples, -1); // FIX: Store intermediate result
        double *predictionDiff = vectorAdd(target, temp_multiply, numberOfSamples);
        delete[] predictions_val;
        delete[] temp_multiply; // FIX: Delete intermediate result

        feature_val      = feature(idx);
        double *temp_add = vectorAdd(predictionDiff, penaltyVector, numberOfSamples); // FIX: Store intermediate result
        double *roVector = vectorMultiplyComponentWise(feature_val, temp_add, numberOfSamples);
        delete[] feature_val;
        delete[] temp_add; // FIX: Delete intermediate result

        double roValue = vectorSum(roVector, numberOfSamples);
        results[idx]   = roValue;

        delete[] penaltyVector;
        delete[] predictionDiff;
        delete[] roVector;
    }

    return results;
}

double LassoRegression::coordinateDescentStep(int weightIdx, double alpha) {
    double *roValues = ro();

    double newWeight;
    if (weightIdx == 0) {
        newWeight = roValues[weightIdx];
    } else if (roValues[weightIdx] < (-1.0) * alpha / 2.0) {
        newWeight = roValues[weightIdx] + alpha / 2.0;
    } else if (roValues[weightIdx] > alpha / 2.0) {
        newWeight = roValues[weightIdx] - alpha / 2.0;
    } else {
        newWeight = 0.0;
    }

    delete[] roValues;
    return newWeight;
}

double *LassoRegression::cyclicalCoordinateDescent(double tolerance, double alpha) {
    bool condition = true;
    double maxChange;

    while (condition) {
        maxChange          = 0.0;
        double *newWeights = new double[numberOfFeatures];

        for (int weightIdx = 0; weightIdx < numberOfFeatures; ++weightIdx) {
            double oldWeight        = weights[weightIdx];
            double newWeight        = coordinateDescentStep(weightIdx, alpha);
            newWeights[weightIdx]   = newWeight;
            weights[weightIdx]      = newWeight;
            double coordinateChange = fabs(oldWeight - newWeight);

            if (coordinateChange > maxChange) {
                maxChange = coordinateChange;
            }
        }

        delete[] newWeights; // FIX: Delete allocated array

        if (maxChange < tolerance) {
            condition = false;
        }
    }

    return weights;
}

double **LassoRegression::featuresMatrix(std::vector<std::vector<double> > samples) {
    double **matrix = emptyMatrix();

    for (int sampleIdx = 0; sampleIdx < numberOfSamples; sampleIdx++) {
        for (int featureIdx = 0; featureIdx < numberOfFeatures; featureIdx++) {
            matrix[sampleIdx][featureIdx] = samples[sampleIdx][featureIdx];
        }
    }

    return matrix;
}

double **LassoRegression::normalizeFeatures(double **matrix) {
    for (int featureIdx = 0; featureIdx < numberOfFeatures; ++featureIdx) {
        double *feature_vec = feature(featureIdx); // FIX: Store pointer
        double featureNorm  = norm(feature_vec, numberOfSamples);
        delete[] feature_vec; // FIX: Delete allocated array

        for (int sampleIdx = 0; sampleIdx < numberOfSamples; ++sampleIdx) {
            matrix[sampleIdx][featureIdx] /= featureNorm;
        }
    }

    return matrix;
}

double **LassoRegression::emptyMatrix() {
    double **result = new double *[numberOfSamples];
    for (int sampleIdx = 0; sampleIdx < numberOfSamples; sampleIdx++) {
        result[sampleIdx] = new double[numberOfFeatures];
    }

    return result;
}

double *LassoRegression::initialWeights() {
    double *weights = new double[numberOfFeatures];

    for (int idx = 0; idx < numberOfFeatures; idx++) {
        weights[idx] = 0.5;
    }

    return weights;
}

double *LassoRegression::targetAsArray(std::vector<double> target) {
    double *result = new double[target.size()];

    for (int targetIdx = 0; targetIdx < target.size(); targetIdx++) {
        result[targetIdx] = target[targetIdx];
    }

    return result;
}

double *LassoRegression::feature(int featureIdx) {
    double *result = new double[numberOfSamples];

    for (int idx = 0; idx < numberOfSamples; idx++) {
        result[idx] = features[idx][featureIdx];
    }

    return result;
}

void LassoRegression::dumpWeightsToFile() {
    std::ofstream file;
    file.open("weights.txt");

    for (int weightIdx = 0; weightIdx < numberOfFeatures; ++weightIdx) {
        file << weights[weightIdx] << " ";
    }
    file.close();
}
