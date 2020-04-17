// Copyright (c) ppy Pty Ltd <contact@ppy.sh>. Licensed under the MIT Licence.
// See the LICENCE file in the repository root for full licence text.

using System;
using System.Linq;
using System.Collections.Generic;
using MathNet.Numerics.Distributions;

namespace osu.Game.Rulesets.Osu.Difficulty.MathUtil
{
    /// <summary>
    /// approximate poisson binomial CDF defined by miss probabilities
    /// see "Refined Normal Approximation (RNA)" from
    /// https://www.researchgate.net/publication/257017356_On_computing_the_distribution_function_for_the_Poisson_binomial_distribution
    /// </summary>
    public class PoissonBinomial
    {
        private readonly double mu, sigma, v;

        public PoissonBinomial(IList<double> probabilities)
        {
            mu = probabilities.Sum();

            sigma = 0;
            double gamma = 0;

            foreach (double p in probabilities)
            {
                sigma += p * (1 - p);
                gamma += p * (1 - p) * (1 - 2 * p);
            }

            sigma = Math.Sqrt(sigma);

            v = gamma / (6 * Math.Pow(sigma, 3));
        }

        public double Cdf(double count)
        {
            double k = (count + 0.5 - mu) / sigma;

            double result = Normal.CDF(0, 1, k) + v * (1 - k * k) * Normal.PDF(0, 1, k);

            if (result < 0) return 0;
            if (result > 1) return 1;

            return result;
        }
    }

    public class IterativePoissonBinomial
    {
        private double mu=0, var=0, gamma=0;

        public void AddProbability(double p)
        {
            mu += p;
            var += p * (1 - p);
            gamma += p * (1 - p) * (1 - 2 * p);
        }

        public double Cdf(double count)
        {
            if (var == 0)
                return mu <= count ? 1 : 0;

            double sigma = Math.Sqrt(var);
            double v = gamma / (6 * Math.Pow(sigma, 3));
            double k = (count + 0.5 - mu) / sigma;

            double result = Normal.CDF(0, 1, k) + v * (1 - k * k) * Normal.PDF(0, 1, k);

            if (result < 0) return 0;
            if (result > 1) return 1;

            return result;
        }


    }
}
