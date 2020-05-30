#include <stdio.h>
#include <random>
#include <vector>

#define DETERMINISTIC() 0

// TODO: up these counts after things are working
static const size_t c_Iterative_Tests = 1000000;
static const size_t c_Iterative_Iterations = 100;

static const size_t c_Flat_Tests = 100000;

static const size_t c_histogramBuckets = 256;

// ===================================================================

static const float c_pi = 3.14159265359f;

template <typename T>
T Clamp(T value, T min, T max)
{
    if (value <= min)
        return min;
    else if (value >= max)
        return max;
    else
        return value;
}

float Lerp(float A, float B, float t)
{
    return A * (1.0f - t) + B * t;
}

void TestIterative(std::mt19937& rng)
{
    std::uniform_real_distribution<float> dist_Neg1_1(-1.0f, 1.0f);
    std::uniform_real_distribution<float> dist_0_1(0.0f, 1.0f);

    auto func = [] (float x)
    {
        return Clamp(sin(x * c_pi / 2.0f), -1.0f, 1.0f);
    };

    float monteCarloYAvg = 0.0f;
    float monteCarloYSquaredAvg = 0.0f;
    float monteCarloStdDev = 0.0f;

    // TODO: clean these up
    float russianRouletteYAvg = 0.0f;
    float russianRouletteYSquaredAvg = 0.0f;
    float russianRouletteStdDev = 0.0f;
    size_t russianRouletteMinIterations = c_Iterative_Iterations;
    size_t russianRouletteMaxIterations = 0;
    float russianRouletteAvgIterations = 0.0f;

    // TODO: ditch lerp for averaging? maybe not. actually no.. go back :P
    // TODO: could do a sum of a function instead of integrating it, for the inner function. might make more sense. could make it like sin(x) + 0.1f so there was an actual non zero value to reach

    for (size_t testIndex = 1; testIndex <= c_Iterative_Tests; ++testIndex)
    {
        // monte carlo
        {
            // take one sample of this recursive integral by doing "a lot" of iterations instead of infinite iterations
            float Y = dist_Neg1_1(rng);
            float YAvg = 0.0f;
            for (size_t iterationIndex = 1; iterationIndex <= c_Iterative_Iterations; ++iterationIndex)
            {
                Y = func(Y);
                YAvg += Y;
                // TODO: rename YAvg to YSum
            }

            // TODO: make this code look like the MC in the flat test

            // integrate the sample
            monteCarloYAvg = Lerp(monteCarloYAvg, YAvg, 1.0f / float(testIndex));
            monteCarloYSquaredAvg = Lerp(monteCarloYSquaredAvg, YAvg*YAvg, 1.0f / float(testIndex));
        }

        // Russian roulette
        {
            // take one sample of this recursive integral by using Russian roulette instead of infinite iterations
            float x = dist_Neg1_1(rng);
            float YAvg = 0.0f;

            float weight = 1.0f;
            size_t russianRoulette_Samples = 0;
            for (size_t iterationIndex = 1; iterationIndex <= c_Iterative_Iterations; ++iterationIndex)
            {
                // probability of killing - kill more towards zero, where the function is smaller, less at the edges where it's larger.
                float p = 1.0f - abs(x);
                if (dist_0_1(rng) < p)
                    break;

                // boost the weighting for this and all future iterations, to account for the ones that got killed
                weight /= (1.0f - p);

                // evaluate the function
                float y = func(x);

                // take the answer as our next input
                x = y;

                // keep track of the YAvg
                YAvg += y * weight;

                russianRoulette_Samples++;
            }

            // integrate the sample
            russianRouletteYAvg = Lerp(russianRouletteYAvg, YAvg, 1.0f / float(testIndex));
            russianRouletteYSquaredAvg = Lerp(russianRouletteYSquaredAvg, YAvg*YAvg, 1.0f / float(testIndex));

            // keep track of min/max/avg samples
            russianRouletteMinIterations = std::min(russianRouletteMinIterations, russianRoulette_Samples);
            russianRouletteMaxIterations = std::max(russianRouletteMaxIterations, russianRoulette_Samples);
            russianRouletteAvgIterations = Lerp(russianRouletteAvgIterations, float(russianRoulette_Samples), 1.0f / float(testIndex));
        }
    }

    // TODO: do same RR test, but make it be 1/10th the probability of killing

    // TODO: use a result structure like the flat test

    // calculate standard deviations
    monteCarloStdDev = sqrt(abs(monteCarloYSquaredAvg) - (monteCarloYAvg * monteCarloYAvg));
    russianRouletteStdDev = sqrt(abs(russianRouletteYSquaredAvg) - (russianRouletteYAvg * russianRouletteYAvg));

    // TODO: print out value and variance of each technique, and loopcounts
    // TODO: write CSV
    
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

    // y = x^2
    auto func = [](float x)
    {
        return x * x;
    };

    // integrated from 0 to 1 is 1/3 or 0.33333...
    const float c_actualValue = 1.0f / 3.0f;

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

    float russianRoulette1_YAvg = 0.0f;
    float russianRoulette1_YSquaredAvg = 0.0f;
    size_t russianRoulette1_Samples = 0;

    float russianRoulette2_YAvg = 0.0f;
    float russianRoulette2_YSquaredAvg = 0.0f;
    size_t russianRoulette2_Samples = 0;

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
        SamplePoint russianRoulette1;
        SamplePoint russianRoulette2;
    };

    std::vector<Result> results(c_Flat_Tests);

    std::vector<size_t> rejectionSampling_histogram(c_histogramBuckets, 0);
    std::vector<size_t> russianRoulette_histogram(c_histogramBuckets, 0);

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

            // keep track of the histogram, to compare the rejection sampling histogram against the russian roulette histogram
            size_t histogramBucket = Clamp<size_t>(size_t(x*float(c_histogramBuckets)), 0, c_histogramBuckets - 1);
            rejectionSampling_histogram[histogramBucket]++;
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

        // Russian Roulette: constant 25% probability of killing a sample
        {
            float x = dist(rng);
            float y = 0.0f;

            // probability of killing = 0.25
            float p = 0.25f;
            if (dist(rng) > p)
            {
                y = func(x) / (1.0f - p);
                russianRoulette1_Samples++;
            }
            else
            {
                y = 0.0f;
            }

            // the averaging happens whether we killed or not
            russianRoulette1_YAvg = Lerp(russianRoulette1_YAvg, y, 1.0f / float(index));
            russianRoulette1_YSquaredAvg = Lerp(russianRoulette1_YSquaredAvg, y*y, 1.0f / float(index));

            results[index - 1].russianRoulette1.value = russianRoulette1_YAvg;
            results[index - 1].russianRoulette1.stdDev = sqrt(abs(russianRoulette1_YSquaredAvg) - (russianRoulette1_YAvg*russianRoulette1_YAvg));
        }

        // Russian Roulette: the probability of killing a sample is 1-x
        // Meaning the probability of keeping a sample is x
        // So, we are rejection sampling the pdf y=x
        {
            float x = dist(rng);
            float y = 0.0f;

            // probability of killing
            float p = (1.0f - x);
            if (dist(rng) > p)
            {
                y = func(x) / (1.0f - p);
                russianRoulette2_Samples++;

                // keep track of the histogram, to compare the rejection sampling histogram against the russian roulette histogram
                size_t histogramBucket = Clamp<size_t>(size_t(x*float(c_histogramBuckets)), 0, c_histogramBuckets - 1);
                russianRoulette_histogram[histogramBucket]++;
            }
            else
            {
                y = 0.0f;
            }

            // the averaging happens whether we killed or not, which is why we have to boost y to compensate for lost samples
            russianRoulette2_YAvg = Lerp(russianRoulette2_YAvg, y, 1.0f / float(index));
            russianRoulette2_YSquaredAvg = Lerp(russianRoulette2_YSquaredAvg, y*y, 1.0f / float(index));

            results[index - 1].russianRoulette2.value = russianRoulette2_YAvg;
            results[index - 1].russianRoulette2.stdDev = sqrt(abs(russianRoulette2_YSquaredAvg) - (russianRoulette2_YAvg*russianRoulette2_YAvg));
        }
    }

    // report results
    printf("Flat Function (%zu tests):\n", c_Flat_Tests);
    printf(" Monte Carlo: \n  value = %f\n  stddev = %f\n", results.rbegin()->monteCarlo.value, results.rbegin()->monteCarlo.stdDev);
    printf(" Rejection Sampling y=2x: \n  value = %f\n  stddev = %f\n  attempts = %zu (%i%%)\n", results.rbegin()->rejectionSampling1.value, results.rbegin()->rejectionSampling1.stdDev, rejectionSampling1_Attempts, int(float(rejectionSampling1_Attempts) * 100.0f / float(c_Flat_Tests)));
    printf(" Rejection Sampling y=2(1-x): \n  value = %f\n  stddev = %f\n  attempts = %zu (%i%%)\n", results.rbegin()->rejectionSampling2.value, results.rbegin()->rejectionSampling2.stdDev, rejectionSampling2_Attempts, int(float(rejectionSampling2_Attempts) * 100.0f / float(c_Flat_Tests)));
    printf(" Rejection Sampling y=3(1-x)^2: \n  value = %f\n  stddev = %f\n  attempts = %zu (%i%%)\n", results.rbegin()->rejectionSampling3.value, results.rbegin()->rejectionSampling3.stdDev, rejectionSampling3_Attempts, int(float(rejectionSampling3_Attempts) * 100.0f / float(c_Flat_Tests)));
    printf(" Russian Roulette p=0.25: \n  value = %f\n  stddev = %f\n  samples = %zu (%i%%)\n", results.rbegin()->russianRoulette1.value, results.rbegin()->russianRoulette1.stdDev, russianRoulette1_Samples, int(float(russianRoulette1_Samples) * 100.0f / float(c_Flat_Tests)));
    printf(" Russian Roulette p=1-x: \n  value = %f\n  stddev = %f\n  samples = %zu (%i%%)\n", results.rbegin()->russianRoulette2.value, results.rbegin()->russianRoulette2.stdDev, russianRoulette2_Samples, int(float(russianRoulette2_Samples) * 100.0f / float(c_Flat_Tests)));

    // write results to csv
    {
        FILE* file;
        fopen_s(&file, "flat.csv", "w+t");
        fprintf(file, "\"Sample\"");
        fprintf(file, ",\"Monte Carlo Value\",\"Monte Carlo Abs Error\",\"Monte Carlo StdDev\"");
        fprintf(file, ",\"RS 2x Value\",\"RS 2x Abs Error\",\"RS 2x StdDev\"");
        fprintf(file, ",\"RS 2(1-x) Value\",\"RS 2(1-x) Abs Error\",\"RS 2(1-x) StdDev\"");
        fprintf(file, ",\"RS 3(1-x)^2 Value\",\"RS 3(1-x)^2 Abs Error\",\"RS 3(1-x)^2 StdDev\"");
        fprintf(file, ",\"RR 0.25 Value\",\"RR 0.25 Abs Error\",\"RR 0.25 StdDev\"");
        fprintf(file, ",\"RR 1-x Value\",\"RR 1-x Abs Error\",\"RR 1-x StdDev\"");
        fprintf(file, "\n");
        size_t index = 1;
        for (const Result& result : results)
        {
            fprintf(file, "\"%zu\"", index);
            fprintf(file, ",\"%f\",\"%f\",\"%f\"", result.monteCarlo.value, abs(result.monteCarlo.value - c_actualValue), result.monteCarlo.stdDev);
            fprintf(file, ",\"%f\",\"%f\",\"%f\"", result.rejectionSampling1.value, abs(result.rejectionSampling1.value - c_actualValue), result.rejectionSampling1.stdDev);
            fprintf(file, ",\"%f\",\"%f\",\"%f\"", result.rejectionSampling2.value, abs(result.rejectionSampling2.value - c_actualValue), result.rejectionSampling2.stdDev);
            fprintf(file, ",\"%f\",\"%f\",\"%f\"", result.rejectionSampling3.value, abs(result.rejectionSampling3.value - c_actualValue), result.rejectionSampling3.stdDev);
            fprintf(file, ",\"%f\",\"%f\",\"%f\"", result.russianRoulette1.value, abs(result.russianRoulette1.value - c_actualValue), result.russianRoulette1.stdDev);
            fprintf(file, ",\"%f\",\"%f\",\"%f\"", result.russianRoulette2.value, abs(result.russianRoulette2.value - c_actualValue), result.russianRoulette2.stdDev);
            fprintf(file, "\n");
            index++;
        }
        fclose(file);
    }

    // Write the histograms to show that rejection sampling y=2x is the same as russian roulette with a probability of keeping a sample as y=x 
    // Note that a PDF is non negative everywere and integrates to 1, but it is relative probabilities, not absolute, so some values can be greater than 1
    // but that doesn't mean it's > 100% probability.
    {
        FILE* file;
        fopen_s(&file, "histogram.csv", "w+t");
        fprintf(file, "\"y=2x\",\"Rejection Sampling y=2x\",\"Russian Roulette p=(1-x)\"\n");
        for (size_t index = 0; index < c_histogramBuckets; ++index)
        {
            fprintf(file, "\"%f\",", float(index * 2) / float(c_histogramBuckets));
            fprintf(file, "\"%f\",", float(c_histogramBuckets) * float(rejectionSampling_histogram[index]) / float(c_Flat_Tests));
            fprintf(file, "\"%f\"\n", float(c_histogramBuckets) * float(russianRoulette_histogram[index]) / float(russianRoulette2_Samples));
        }
        fclose(file);
    }
}

int main(int argc, char** argv)
{
#if DETERMINISTIC()
    std::mt19937 rng;
#else
    std::random_device rd;
    std::mt19937 rng(rd());
#endif

    // TODO: uncomment before checkin
    //TestFlat(rng);
    TestIterative(rng);

    return 0;
}

/*

TODO:

* maybe the iterative function is too chaotic. could try something simpler & smoother.

! need to boost values of survivors

- chaotic iterative function. y=sin(1/x)
 * "many iterations" vs "russian roulette" of y value?





 Blog Notes:

 "The connection between russian roulette and rejection sampling / importance sampling"

 Maybe start the explanation with the non recursive function.
 * If you do 4 samples of the function y=1, and reject 25% of the time, you get an average of 3/4 instead of 1.
 * You need to divide value by (1-0.25) aka 0.75 or 3/4.  Aka you need to multiply by 4/3.
 * 4/3 times the 3 samples that survived is 12/3 or 4.
 * Averaging as if there are 4 samples. You divide 4 by 4 and get 1.
 * The 3 samples each made themselves count for a little bit more, to make up for the lost sample.


 ? mention how you can directly draw from the PDF without importance sampling: https://blog.demofox.org/2017/08/05/generating-random-numbers-from-a-specific-distribution-by-inverting-the-cdf/
* link to 1d monte carlo basics: https://blog.demofox.org/2018/06/12/monte-carlo-integration-explanation-in-1d/
 * link to this: http://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/Russian_Roulette_and_Splitting.html
 * link to this for chaotic function: https://computing.dcu.ie/~humphrys/Notes/Neural/chaos.html

  ? should i link to this for how lerp is calculating an averag? https://blog.demofox.org/2020/03/10/how-do-i-calculate-variance-in-1-pass/


 * summarize the PBRT post.
  * russian roulette CAN reduce variance but it might not

 * russian roulette takes a lot fewer steps, that could make it better ultimately if steps were really costly to do.

 * can even do a fixed percentage for russian roulette- but why? i did it to show it

 */