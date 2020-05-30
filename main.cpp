#include <stdio.h>
#include <random>
#include <vector>

#define DETERMINISTIC() 1

static const size_t c_Iterative_Tests = 1000;
static const size_t c_Iterative_Iterations = 100000;

static const size_t c_Flat_Tests = 100000;

static const float c_pi = 3.14159265359f;

float Lerp(float A, float B, float t)
{
    return A * (1.0f - t) + B * t;
}

void TestIterative(std::mt19937& rng)
{
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    auto func = [] (float x)
    {
        return sin(x * c_pi);
        //return sin(100.0f / (x + 0.001f)) * 0.5f + 0.5f;
    };

    float monteCarloYAvg = 0.0f;
    float monteCarloYAvgSquared = 0.0f;
    float monteCarloStdDev = 0.0f;

    float russianRouletteYAvg = 0.0f;
    float russianRouletteYAvgSquared = 0.0f;
    float russianRouletteStdDev = 0.0f;
    size_t russianRouletteMinIterations = 1000;
    size_t russianRouletteMaxIterations = 0;
    float russianRouletteAvgIterations = 0.0f;

    // TODO: ditch lerp for averaging!

    for (size_t testIndex = 1; testIndex <= c_Iterative_Tests; ++testIndex)
    {
        // monte carlo
        {
            // take one sample of this recursive integral by doing "a lot" of iterations instead of infinite iterations
            float Y = dist(rng);
            float YAvg = 0.0f;
            for (size_t iterationIndex = 1; iterationIndex <= c_Iterative_Iterations; ++iterationIndex)
            {
                if (iterationIndex == 4028)
                {
                    int ijkl = 0;
                }
                Y = func(Y);
                YAvg += Y / float(c_Iterative_Iterations);
            }

            // integrate the sample
            monteCarloYAvg += YAvg / float(c_Iterative_Tests);
            monteCarloYAvgSquared += (YAvg * YAvg) / float(c_Iterative_Tests);

            // TODO: write out value and variance to file? or maybe store in an array so we can make a better csv?
        }

        // Russian roulette
        {
            // take one sample of this recursive integral by using Russian roulette instead of infinite iterations
            float Y = dist(rng);
            float YAvg = 0.0f;
            float weightTotal = 1.0f;

            float weight = 1.0f;

            // TODO: i'm not convinced this is correct. making p = Y * 2.0 should keep the average value the same but decrease variance. it changes the value though

            size_t iterationIndex = 1;
            while (1)
            {
                Y = func(Y);
                YAvg += Y * weight;

                float p = Y;
                if (dist(rng) > p)
                    break;

                weight *= (1.0f - p);

                ++iterationIndex;
            }

            // TODO: do flattened test, maybe before finishing this and calling it ok

            // TODO: if the iteration count is low, we can change the probability of being killed to be lower. we get more iterations then, so better quality, but still unbiased. I think?

            // integrate the sample
            russianRouletteYAvg = Lerp(russianRouletteYAvg, YAvg, 1.0f / float(testIndex));
            russianRouletteYAvgSquared = Lerp(russianRouletteYAvgSquared, YAvg * YAvg, 1.0f / float(testIndex));

            russianRouletteMinIterations = std::min(russianRouletteMinIterations, iterationIndex);
            russianRouletteMaxIterations = std::max(russianRouletteMaxIterations, iterationIndex);
            russianRouletteAvgIterations = Lerp(russianRouletteAvgIterations, float(iterationIndex), 1.0f / float(testIndex));

            // TODO: write out value and variance to file? or maybe store in an array so we can make a better csv?
        }
    }

    // calculate standard deviations
    monteCarloStdDev = sqrt(abs(monteCarloYAvgSquared) - (monteCarloYAvg * monteCarloYAvg));
    russianRouletteStdDev = sqrt(abs(russianRouletteYAvgSquared) - (russianRouletteYAvg * russianRouletteYAvg));

    // TODO: rename russianRouletteYAvgSquared to russianRouletteYSquaredAvg

    // TODO: print out value and variance of each technique, and loopcounts
    
    printf("Iterative Function (%zu tests):\n", c_Iterative_Tests);
    printf(" Monte Carlo: \n  value = %f\n  stddev = %f\n  iterations = %zu\n", monteCarloYAvg, monteCarloStdDev, c_Iterative_Iterations);
    printf(" Russian Roulette: \n  value = %f\n  stddev = %f\n  iterations min = %zu, max = %zu, avg = %0.2f\n\n\n", russianRouletteYAvg, russianRouletteStdDev, russianRouletteMinIterations, russianRouletteMaxIterations, russianRouletteAvgIterations);

    FILE* file;
    fopen_s(&file, "iterative.csv", "w+t");
    fclose(file);
}

void TestFlat(std::mt19937& rng)
{
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    auto func = [](float x)
    {
        return x * x;
    };

    float monteCarlo_YAvg = 0.0f;
    float monteCarlo_YSquaredAvg = 0.0f;

    float rejectionSampling1_YAvg = 0.0f;
    float rejectionSampling1_YSquaredAvg = 0.0f;
    size_t rejectionSampling1_Attempts = 0;

    float rejectionSampling2_YAvg = 0.0f;
    float rejectionSampling2_YSquaredAvg = 0.0f;
    size_t rejectionSampling2_Attempts = 0;

    float rejectionSampling3_YAvg = 0.0f;
    float rejectionSampling3_YSquaredAvg = 0.0f;
    size_t rejectionSampling3_Attempts = 0;

    struct SamplePoint
    {
        float value;
        float stdDev;
    };

    struct Result
    {
        SamplePoint monteCarlo;
        SamplePoint rejectionSampling1;
        SamplePoint rejectionSampling2;
        SamplePoint rejectionSampling3;
    };

    std::vector<Result> results(c_Flat_Tests);

    for (size_t index = 1; index <= c_Flat_Tests; ++index)
    {
        // monte carlo
        {
            float x = dist(rng);
            float y = func(x);

            monteCarlo_YAvg = Lerp(monteCarlo_YAvg, y, 1.0f / float(index));
            monteCarlo_YSquaredAvg = Lerp(monteCarlo_YSquaredAvg, y*y, 1.0f / float(index));

            results[index - 1].monteCarlo.value = monteCarlo_YAvg;
            results[index - 1].monteCarlo.stdDev = sqrt(abs(monteCarlo_YSquaredAvg) - (monteCarlo_YAvg*monteCarlo_YAvg));
        }

        // Rejection sampling the pdf y = 2x for importance sampling
        {
            float x;
            do
            {
                x = dist(rng);
                rejectionSampling1_Attempts++;
            }
            while(dist(rng) > x);

            float pdf = 2.0f * x;
            float y = func(x) / pdf;

            rejectionSampling1_YAvg = Lerp(rejectionSampling1_YAvg, y, 1.0f / float(index));
            rejectionSampling1_YSquaredAvg = Lerp(rejectionSampling1_YSquaredAvg, y*y, 1.0f / float(index));

            results[index - 1].rejectionSampling1.value = rejectionSampling1_YAvg;
            results[index - 1].rejectionSampling1.stdDev = sqrt(abs(rejectionSampling1_YSquaredAvg) - (rejectionSampling1_YAvg*rejectionSampling1_YAvg));
        }

        // Rejection sampling the pdf y = 2(1-x) for importance sampling - higher variance than the other because the pdf doesn't fit the shape as well
        {
            float x;
            do
            {
                x = dist(rng);
                rejectionSampling2_Attempts++;
            }
            while(dist(rng) > (1-x));

            float pdf = 2.0f * (1.0f - x);
            float y = func(x) / pdf;

            rejectionSampling2_YAvg = Lerp(rejectionSampling2_YAvg, y, 1.0f / float(index));
            rejectionSampling2_YSquaredAvg = Lerp(rejectionSampling2_YSquaredAvg, y*y, 1.0f / float(index));

            results[index - 1].rejectionSampling2.value = rejectionSampling2_YAvg;
            results[index - 1].rejectionSampling2.stdDev = sqrt(abs(rejectionSampling2_YSquaredAvg) - (rejectionSampling2_YAvg*rejectionSampling2_YAvg));
        }

        // Rejection sampling the pdf y = 3(1-x)^2 for importance sampling - higher variance than the other because the pdf doesn't fit the shape as well, but also higher rejection rate because acceptance area is smaller percentage of 2d rectangle!
        {
            float x;
            do
            {
                x = dist(rng);
                rejectionSampling3_Attempts++;
            } while (dist(rng) > (1 - x) * (1 - x));

            float pdf = 3.0f * (1.0f - x) * (1.0f - x);
            float y = func(x) / pdf;

            rejectionSampling3_YAvg = Lerp(rejectionSampling3_YAvg, y, 1.0f / float(index));
            rejectionSampling3_YSquaredAvg = Lerp(rejectionSampling3_YSquaredAvg, y*y, 1.0f / float(index));

            results[index - 1].rejectionSampling3.value = rejectionSampling3_YAvg;
            results[index - 1].rejectionSampling3.stdDev = sqrt(abs(rejectionSampling3_YSquaredAvg) - (rejectionSampling3_YAvg*rejectionSampling3_YAvg));
        }
    }

    // report results
    printf("Flat Function (%zu tests):\n", c_Flat_Tests);
    printf(" Monte Carlo: \n  value = %f\n  stddev = %f\n", monteCarlo_YAvg, results.rbegin()->monteCarlo.stdDev);
    printf(" Rejection Sampling y=2x: \n  value = %f\n  stddev = %f\n  attempts = %zu (%i%%)\n", rejectionSampling1_YAvg, results.rbegin()->rejectionSampling1.stdDev, rejectionSampling1_Attempts, int(float(rejectionSampling1_Attempts) * 100.0f / float(c_Flat_Tests)));
    printf(" Rejection Sampling y=2(1-x): \n  value = %f\n  stddev = %f\n  attempts = %zu (%i%%)\n", rejectionSampling2_YAvg, results.rbegin()->rejectionSampling2.stdDev, rejectionSampling2_Attempts, int(float(rejectionSampling2_Attempts) * 100.0f / float(c_Flat_Tests)));
    printf(" Rejection Sampling y=3(1-x)^2: \n  value = %f\n  stddev = %f\n  attempts = %zu (%i%%)\n", rejectionSampling3_YAvg, results.rbegin()->rejectionSampling3.stdDev, rejectionSampling3_Attempts, int(float(rejectionSampling3_Attempts) * 100.0f / float(c_Flat_Tests)));

    // TODO: write csv of results!
}

int main(int argc, char** argv)
{
#if DETERMINISTIC()
    std::mt19937 rng;
#else
    std::random_device rd;
    std::mt19937 rng(rd());
#endif

    //TestIterative(rng);

    TestFlat(rng);

    return 0;
}

/*

TODO:

* maybe the iterative function is too chaotic. could try something simpler & smoother.

! need to boost values of survivors

- chaotic iterative function. y=sin(1/x)
 * "many iterations" vs "russian roulette" of y value?

- some function z=f(x,y) not iterative.
 * z = (x*x + y) maybe? form 0 to 1.
 * use x as the probability and see how it goes
 * use 1-x as the probability and see how it goes
 * also do with uniform probability.


 * output these things to csvs

 ? maybe graph histogram of accepted values to show how it's like drawing from another PDF & importance sampling?
 * link to monte carlo 1d explanation of yours
 * link to this: http://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/Russian_Roulette_and_Splitting.html
 * link to this for chaotic function: https://computing.dcu.ie/~humphrys/Notes/Neural/chaos.html
 
 ? should i do the test multiple times and get variance besides just error?

 ? should i link to this for how lerp is calculating an averag? https://blog.demofox.org/2020/03/10/how-do-i-calculate-variance-in-1-pass/



 Blog Notes:

 "The connection between russian roulette and rejection sampling / importance sampling"

 ? mention how you can directly draw from the PDF without importance sampling: https://blog.demofox.org/2017/08/05/generating-random-numbers-from-a-specific-distribution-by-inverting-the-cdf/
* link to 1d monte carlo basics: https://blog.demofox.org/2018/06/12/monte-carlo-integration-explanation-in-1d/

 * summarize the PBRT post.
  * russian roulette CAN reduce variance but it might not

 * russian roulette takes a lot fewer steps, that could make it better ultimately if steps were really costly to do.

 * can even do a fixed percentage - but why?

 */