// Copyright (c) ppy Pty Ltd <contact@ppy.sh>. Licensed under the MIT Licence.
// See the LICENCE file in the repository root for full licence text.

using System;
using System.Collections.Generic;
using System.Linq;
using osu.Game.Beatmaps;
using osu.Game.Rulesets.Difficulty;
using osu.Game.Rulesets.Difficulty.Preprocessing;
using osu.Game.Rulesets.Difficulty.Skills;
using osu.Game.Rulesets.Mods;
using osu.Game.Rulesets.Osu.Difficulty.Preprocessing;
using osu.Game.Rulesets.Osu.Difficulty.Skills;
using osu.Game.Rulesets.Osu.Mods;
using osu.Game.Rulesets.Osu.Objects;

namespace osu.Game.Rulesets.Osu.Difficulty
{
    public class OsuDifficultyCalculator : DifficultyCalculator
    {

        public OsuDifficultyCalculator(Ruleset ruleset, WorkingBeatmap beatmap)
            : base(ruleset, beatmap)
        {
        }


        private IEnumerable<OsuHitObjectDifficulty> hitObjectDifficulties(Skill[] skills)
        {
            var timesIt = skills[0].Timestamps.GetEnumerator();
            var aimStarsIt = skills[0].HitObjectStars().GetEnumerator();
            var aimCumStarsIt = skills[0].CumulativeHitObjectStars().GetEnumerator();

            var speedStarsIt = skills[1].HitObjectStars().GetEnumerator();
            var speedCumStarsIt = skills[1].CumulativeHitObjectStars().GetEnumerator();

            while(timesIt.MoveNext() && aimStarsIt.MoveNext() && aimCumStarsIt.MoveNext() && speedStarsIt.MoveNext() && speedCumStarsIt.MoveNext())
            {
                yield return new OsuHitObjectDifficulty
                {
                    Time = timesIt.Current,
                    AimStars = aimStarsIt.Current,
                    AimCumulativeStars = aimCumStarsIt.Current,
                    SpeedStars = speedStarsIt.Current,
                    SpeedCumulativeStars = speedCumStarsIt.Current
                };
            }


        }

        protected override DifficultyAttributes CreateDifficultyAttributes(IBeatmap beatmap, Mod[] mods, Skill[] skills, double clockRate)
        {
            if (beatmap.HitObjects.Count == 0)
                return new OsuDifficultyAttributes { Mods = mods };

            IList<double> aimComboSR = skills[0].ComboSR;
            IList<double> aimMissCounts = skills[0].MissCounts;

            IList<double> speedComboSR = skills[1].ComboSR;
            IList<double> speedMissCounts = skills[1].MissCounts;

            double missSrIncrement = Skill.MissSRIncrement;


            double aimRating = aimComboSR.Last();
            double speedRating = speedComboSR.Last();
            double starRating = aimRating + speedRating + Math.Abs(aimRating - speedRating) / 2;

            // Todo: These int casts are temporary to achieve 1:1 results with osu!stable, and should be removed in the future
            double hitWindowGreat = (int)(beatmap.HitObjects.First().HitWindows.Great / 2) / clockRate;
            double preempt = (int)BeatmapDifficulty.DifficultyRange(beatmap.BeatmapInfo.BaseDifficulty.ApproachRate, 1800, 1200, 450) / clockRate;

            int maxCombo = beatmap.HitObjects.Count;
            // Add the ticks + tail of the slider. 1 is subtracted because the head circle would be counted twice (once for the slider itself in the line above)
            maxCombo += beatmap.HitObjects.OfType<Slider>().Sum(s => s.NestedHitObjects.Count - 1);

            return new OsuDifficultyAttributes
            {
                StarRating = starRating,
                Mods = mods,
                MissSRIncrement = missSrIncrement,
                AimStrain = aimRating,
                AimComboSR = aimComboSR,
                AimMissCounts = aimMissCounts,
                SpeedStrain = speedRating,
                SpeedComboSR = speedComboSR,
                SpeedMissCounts = speedMissCounts,
                ApproachRate = preempt > 1200 ? (1800 - preempt) / 120 : (1200 - preempt) / 150 + 5,
                OverallDifficulty = (80 - hitWindowGreat) / 6,
                MaxCombo = maxCombo,
                HitObjectDifficulties = hitObjectDifficulties(skills).Where(x => x.AimStars != 0).ToList(),
            };
        }

        protected override IEnumerable<DifficultyHitObject> CreateDifficultyHitObjects(IBeatmap beatmap, double clockRate)
        {
            // The first jump is formed by the first two hitobjects of the map.
            // If the map has less than two OsuHitObjects, the enumerator will not return anything.
            for (int i = 1; i < beatmap.HitObjects.Count; i++)
            {
                var lastLast = i > 1 ? beatmap.HitObjects[i - 2] : null;
                var last = beatmap.HitObjects[i - 1];
                var current = beatmap.HitObjects[i];

                yield return new OsuDifficultyHitObject(current, lastLast, last, clockRate);
            }
        }

        protected override Skill[] CreateSkills(IBeatmap beatmap) => new Skill[]
        {
            new Aim(),
            new Speed()
        };

        protected override Mod[] DifficultyAdjustmentMods => new Mod[]
        {
            new OsuModDoubleTime(),
            new OsuModHalfTime(),
            new OsuModEasy(),
            new OsuModHardRock(),
        };
    }
}
