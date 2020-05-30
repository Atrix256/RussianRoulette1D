#include <stdio.h>
#include <random>
#include <vector>

static const size_t c_Iterative_Tests = 1000;
static const size_t c_Iterative_Iterations = 100000;

float Lerp(float A, float B, float t)
{
    return A * (1.0f - t) + B * t;
}

void TestIterative(std::mt19937& rng)
{
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    auto func = [] (float x)
    {
        return sin(100.0f / (x + 0.001f)) * 0.5f + 0.5f;
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
        /*
        // monte carlo
        {
            // take one sample of this recursive integral by doing "a lot" of iterations instead of infinite iterations
            float Y = dist(rng);
            float YAvg = 0.0f;
            for (size_t iterationIndex = 1; iterationIndex <= c_Iterative_Iterations; ++iterationIndex)
            {
                Y = func(Y);
                YAvg += Y / float(c_Iterative_Iterations);
            }

            // integrate the sample
            monteCarloYAvg += YAvg / float(c_Iterative_Tests);
            monteCarloYAvgSquared += (YAvg * YAvg) / float(c_Iterative_Tests);

            // TODO: write out value and variance to file? or maybe store in an array so we can make a better csv?
        }
        */

        // Russian roulette
        {
            // take one sample of this recursive integral by using Russian roulette instead of infinite iterations
            float Y = dist(rng);
            float YAvg = 0.0f;
            float weightTotal = 1.0f;

            // TODO: i'm not convinced this is correct. making p = Y * 2.0 should keep the average value the same but decrease variance. it changes the value though

            size_t iterationIndex = 1;
            while (1)
            {
                Y = func(Y);
                YAvg += Y;

                float p = Y;
                if (dist(rng) > p)
                    break;

                weightTotal += p;

                ++iterationIndex;
            }
            YAvg /= weightTotal;

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

    // TODO: print out value and variance of each technique, and loopcounts
    
    printf("Iterative Function (%zu tests):\n", c_Iterative_Tests);
    printf(" Monte Carlo: \n  value = %f\n  variance = %f\n  iterations = %zu\n", monteCarloYAvg, monteCarloStdDev, c_Iterative_Iterations);

    printf(" Russian Roulette: \n  value = %f\n  variance = %f\n  iterations min = %zu, max = %zu, avg = %0.2f\n\n\n", russianRouletteYAvg, russianRouletteStdDev, russianRouletteMinIterations, russianRouletteMaxIterations, russianRouletteAvgIterations);

    FILE* file;
    fopen_s(&file, "iterative.csv", "w+t");
    fclose(file);
}


int main(int argc, char** argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());

    TestIterative(rng);

    return 0;
}

/*

TODO:

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

 * summarize the PBRT post.
  * russian roulette CAN reduce variance but it might not

 * russian roulette takes a lot fewer steps, that could make it better ultimately if steps were really costly to do.

 * can even do a fixed percentage - but why?

 */