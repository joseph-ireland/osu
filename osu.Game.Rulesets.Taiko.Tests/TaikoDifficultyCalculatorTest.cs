// Copyright (c) ppy Pty Ltd <contact@ppy.sh>. Licensed under the MIT Licence.
// See the LICENCE file in the repository root for full licence text.

using NUnit.Framework;
using osu.Game.Beatmaps;
using osu.Game.Rulesets.Difficulty;
using osu.Game.Rulesets.Taiko.Difficulty;
using osu.Game.Tests.Beatmaps;

namespace osu.Game.Rulesets.Taiko.Tests
{
    public class TaikoDifficultyCalculatorTest : DifficultyCalculatorTest
    {
        protected override string ResourceAssembly => "osu.Game.Rulesets.Taiko";

        [TestCase(2.9811336589467095, "diffcalc-test")]
        [TestCase(2.9811336589467095, "diffcalc-test-strong")]
        public void Test(double expected, string name)
            => base.Test(expected, name);

        protected override LegacyDifficultyCalculator CreateDifficultyCalculator(WorkingBeatmap beatmap) => new TaikoLegacyDifficultyCalculator(new TaikoRuleset(), beatmap);

        protected override Ruleset CreateRuleset() => new TaikoRuleset();
    }
}
