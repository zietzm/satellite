import Data.Complex (Complex ((:+)))
import qualified Data.Complex as C
import Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as V
import GHC.Float (double2Float, float2Double, int2Float)
import qualified Norm
import qualified Signals
import qualified Sync
import Test.Tasty
import Test.Tasty.HUnit

main :: IO ()
main = defaultMain tests

tests :: TestTree
tests = testGroup "Tests" [dspTests, syncTests, normTests]

dspTests :: TestTree
dspTests =
  testGroup
    "Signal processing"
    [ testCase "FFT expected length" $
        V.length (Signals.fft wave) @?= V.length wave,
      testCase "Envelope expected length" $
        V.length (Signals.getEnvelope wave) @?= V.length wave,
      testCase "FFT expected result" $
        let maxDiffFft = maxDiffD (Signals.fft wave) waveFft
         in assertBool ("FFT differed from Python: " ++ show maxDiffFft) (maxDiffFft < 1e-6),
      testCase "Wave envelope" $
        let maxDiffEnv = maxDiffF (Signals.getEnvelope wave) waveEnvelope
         in assertBool ("Envelope differed from Python: " ++ show maxDiffEnv) (maxDiffEnv < 1e-6),
      -- TODO: We could make this more accurate in the future. For now, it works.
      testCase "Low pass filter" $
        let t = V.enumFromN 0 2000 :: Vector Float
            xLow = V.map (sin . (* (2 * 5 * pi))) t
            xHigh = V.map (sin . (* (2 * 250 * pi))) t
            x = V.zipWith (+) xLow xHigh
            xL = map float2Double $ V.toList x
            filtered = Signals.lowpass 2000.0 0.125 xL
            filtered' = V.fromList $ map double2Float filtered
            maxDiffFilt = maxDiffF (V.drop 10 filtered') (V.drop 10 xLow)
         in assertBool ("Filter error" ++ show maxDiffFilt) (maxDiffFilt < 5e-3),
      testCase "Downsampling length" $
        V.length rDown @?= V.length expectDown,
      testCase "Upsampling length" $
        V.length rUp @?= V.length expectUp,
      testCase "Downsampling result" $
        assertBool
          ("Downsampling error" ++ show maxDiffDown)
          (maxDiffDown < 1e-5),
      testCase
        "Upampling result"
        $ assertBool ("Upsampling error" ++ show maxDiffUp) (maxDiffUp < 1e-5),
      testCase "Interpolation" $
        let xOld = V.fromList [1.0, 3.0, 5.0] :: V.Vector Float
            yOld = V.fromList [1.0, 9.0, 25.0]
            xNew = V.fromList [2.0, 4.0]
            result = Signals.interp xOld yOld xNew
            expected = V.fromList [5.0, 17.0]
            maxDiff = V.maximum $ V.map abs $ V.zipWith (-) result expected
         in assertBool ("Interpolation error: " ++ show maxDiff) (maxDiff < 1e-6),
      testCase "Interpolation outside range" $
        let xOld = V.fromList [1.0, 3.0, 5.0] :: V.Vector Float
            yOld = V.fromList [1.0, 9.0, 25.0]
            xNew = V.fromList [0.0, 6.0]
            result = Signals.interp xOld yOld xNew
            expected = V.fromList [1.0, 25.0]
            maxDiff = V.maximum $ V.map abs $ V.zipWith (-) result expected
         in assertBool ("Interpolation error: " ++ show maxDiff) (maxDiff < 1e-6)
    ]
  where
    -- Test Fourier resampling
    t2 = V.generate 100 ((* 0.01) . int2Float) :: Vector Float
    x2 = V.map (sin . (* (2 * pi * 10))) t2
    rDown = Signals.fourierResample 50 x2
    maxDiffDown = maxDiffF rDown expectDown
    rUp = Signals.fourierResample 200 x2
    maxDiffUp = maxDiffF rUp expectUp

syncTests :: TestTree
syncTests =
  testGroup
    "Synchronization"
    [ testCase "Self-convolve = length for syncA" $
        Sync.crossCorr Sync.syncA Sync.syncA @?= V.fromList [40.0],
      testCase "Self-convolve = length for syncB" $
        Sync.crossCorr Sync.syncB Sync.syncB @?= V.fromList [40.0],
      testCase "Adjust pattern 2 [1,2] -> [1,1,2,2]" $
        let raw = V.fromList [1.0, 2.0]
            expected = V.fromList [1.0, 1.0, 2.0, 2.0]
         in Sync.upsamplePattern 2 raw @?= expected,
      testCase "Find basic peaks" $ basicPeaks @?= V.fromList [1, 3, 5],
      testCase "Peak heights" $ peakHeights @?= V.fromList [2, 5, 9],
      testCase "Prioritize highest" $ filteredPeaks @?= V.fromList [1, 3, 5],
      testCase "Find peaks" $ peaks @?= V.fromList [1, 3, 5]
    ]
  where
    -- Test peak finding
    heights = V.fromList [1, 2, 1, 5, 1, 9, 1] :: Vector Int
    basicPeaks = Sync.findBasicPeaks 0 heights
    peakHeights = V.map (heights V.!) basicPeaks
    filteredPeaks = Sync.prioritizeHighest 2 peakHeights basicPeaks
    peaks = Sync.findPeaks 0 2 heights

normTests :: TestTree
normTests =
  testGroup
    "Normalization"
    [ testCase "sortStorableVector sorts correctly" $
        let unsorted = V.fromList [3.0, 1.0, 4.0, 1.5, 9.0] :: Vector Float
            expected = V.fromList [1.0, 1.5, 3.0, 4.0, 9.0]
         in Norm.sortStorableVector unsorted @?= expected,
      testCase "percentile computes 25th and 75th percentiles" $
        let input = V.fromList [1.0, 2.0, 3.0, 4.0, 5.0]
            result = Norm.percentile [25.0, 75.0] input
            expected = [2.0, 4.0] -- Approximate, assuming linear interpolation
         in assertApproxEqual "Percentiles" 1e-5 result expected,
      testCase "small range percentiles" $
        let input = V.fromList [1.0, 1.1, 1.2, 1.3, 1.4]
            result = Norm.percentile [20.0, 80.0] input
            expected = [1.08, 1.32]
         in assertApproxEqual "Percentiles" 1e-5 result expected,
      testCase "rangeNorm clamps and normalizes between 0 and 1" $
        let input = V.fromList [1.0, 2.0, 3.0, 4.0, 5.0]
            result = Norm.rangeNorm 25.0 75.0 input
            expected = V.fromList [0.0, 0.0, 0.5, 1.0, 1.0] -- Based on 25th=2, 75th=4
         in assertApproxEqual "Normalized values" 1e-5 (V.toList result) (V.toList expected),
      testCase "rangeNorm handles edge case with small range" $
        let input = V.fromList [1.0, 1.1, 1.2, 1.3, 1.4]
            result = Norm.rangeNorm 20.0 80.0 input
            expected = V.fromList [0.0, 1 / 12, 1 / 2, 11 / 12, 1] -- Approx, small range
         in assertApproxEqual "Small range normalization" 1e-5 (V.toList result) (V.toList expected)
    ]

-- Helper function for approximate equality due to floating-point imprecision
assertApproxEqual :: String -> Double -> [Double] -> [Double] -> Assertion
assertApproxEqual msg eps actual expected =
  let pairs = zip actual expected
      diffs = map (\(a, e) -> abs (a - e) <= eps) pairs
   in assertBool (msg ++ ": " ++ show actual ++ " not approx equal to " ++ show expected) (and diffs)

maxDiffF :: Vector Float -> Vector Float -> Float
maxDiffF x y = V.maximum $ V.zipWith (\a b -> abs (a - b)) x y

maxDiffD :: Vector (Complex Double) -> Vector (Complex Double) -> Double
maxDiffD x y = V.maximum $ V.map C.magnitude $ V.zipWith (-) x y

-- Built in Python using:
-- t = np.linspace(0, 1, 32)
-- signal = np.sin(2 * np.pi * 5 * t) + 0.5 * np.sin(2 * np.pi * 20 * t)
{- ORMOLU_DISABLE -}
wave :: Vector Float
wave = V.fromList [
  0.00000000e+00,  4.53256389e-01,  1.38184310e+00, -9.60096056e-02,
  -1.03342672e+00, -4.43517970e-01, -5.63694914e-01,  6.74208626e-01,
   1.39239925e+00, -1.69512943e-01, -5.01690921e-01, -7.02834217e-01,
  -8.93714109e-01,  8.96954456e-01,  1.09936577e+00,  3.63996927e-02,
  -3.63996927e-02, -1.09936577e+00, -8.96954456e-01,  8.93714109e-01,
   7.02834217e-01,  5.01690921e-01,  1.69512943e-01, -1.39239925e+00,
  -6.74208626e-01,  5.63694914e-01,  4.43517970e-01,  1.03342672e+00,
   9.60096056e-02, -1.38184310e+00, -4.53256389e-01, -3.67394040e-15]
{- ORMOLU_ENABLE -}

-- Built in Python using:
-- np.fft.fft(signal)
waveFft :: Vector (Complex Double)
waveFft =
  V.fromList
    [ negate 2.88657986e-15 :+ 0.0,
      1.56605153e-02 :+ negate 0.15900388,
      7.19912997e-02 :+ negate 0.3619247,
      2.13316394e-01 :+ negate 0.70320991,
      6.53756462e-01 :+ negate 1.57830772,
      7.01262326e+00 :+ negate 13.11969534,
      negate 1.91335470e+00 :+ 2.86353766,
      negate 1.22009858e+00 :+ 1.48669442,
      negate 1.12514917e+00 :+ 1.12514917,
      negate 1.26692084e+00 :+ 1.03973506,
      negate 1.83502514e+00 :+ 1.2261246,
      negate 6.04812457e+00 :+ 3.23278993,
      2.96641940e+00 :+ negate 1.22873115,
      1.07930981e+00 :+ negate 0.32740505,
      6.58027278e-01 :+ negate 0.13088976,
      5.05431552e-01 :+ negate 0.04978066,
      4.64274046e-01 :+ 0.0,
      5.05431552e-01 :+ 0.04978066,
      6.58027278e-01 :+ 0.13088976,
      1.07930981e+00 :+ 0.32740505,
      2.96641940e+00 :+ 1.22873115,
      negate 6.04812457e+00 :+ negate 3.23278993,
      negate 1.83502514e+00 :+ negate 1.2261246,
      negate 1.26692084e+00 :+ negate 1.03973506,
      negate 1.12514917e+00 :+ negate 1.12514917,
      negate 1.22009858e+00 :+ negate 1.48669442,
      negate 1.91335470e+00 :+ negate 2.86353766,
      7.01262326e+00 :+ 13.11969534,
      6.53756462e-01 :+ 1.57830772,
      2.13316394e-01 :+ 0.70320991,
      7.19912997e-02 :+ 0.3619247,
      1.56605153e-02 :+ 0.15900388
    ]

-- Built in Python using:
-- np.abs(scipy.signal.hilbert(signal))
{- ORMOLU_DISABLE -}
waveEnvelope :: Vector Float
waveEnvelope =
  V.fromList
    [ 0.41780733, 0.86659862, 1.41418247, 1.50467192, 1.04570146, 0.50550353,
      0.85786936, 1.34758746, 1.48181892, 1.15394755, 0.57287563, 0.75230654,
      1.30184956, 1.49352365, 1.22515443, 0.65818438, 0.65818438, 1.22515443,
      1.49352365, 1.30184956, 0.75230654, 0.57287563, 1.15394755, 1.48181892,
      1.34758746, 0.85786936, 0.50550353, 1.04570146, 1.50467192, 1.41418247,
      0.86659862, 0.41780733
    ]
{- ORMOLU_ENABLE -}

-- From Python:
-- fs = 100  # Original sampling rate (Hz)
-- t = np.linspace(0, 1, fs, endpoint=False)  # 100 samples over 1 second
-- freq = 10  # Frequency of the sine wave (Hz)
-- signal_orig = np.sin(2 * np.pi * freq * t)
-- signal_down = signal.resample(signal_orig, 50)
{- ORMOLU_DISABLE -}
expectDown :: Vector Float
expectDown =
  V.fromList
    [ -1.91395591e-15,  9.51056516e-01,  5.87785252e-01, -5.87785252e-01,
      -9.51056516e-01, -9.32860229e-17,  9.51056516e-01,  5.87785252e-01,
      -5.87785252e-01, -9.51056516e-01, -5.51744758e-16,  9.51056516e-01,
       5.87785252e-01, -5.87785252e-01, -9.51056516e-01, -1.06695564e-15,
       9.51056516e-01,  5.87785252e-01, -5.87785252e-01, -9.51056516e-01,
      -4.19778198e-16,  9.51056516e-01,  5.87785252e-01, -5.87785252e-01,
      -9.51056516e-01, -2.06574631e-16,  9.51056516e-01,  5.87785252e-01,
      -5.87785252e-01, -9.51056516e-01, -3.35356809e-15,  9.51056516e-01,
       5.87785252e-01, -5.87785252e-01, -9.51056516e-01,  1.80989786e-15,
       9.51056516e-01,  5.87785252e-01, -5.87785252e-01, -9.51056516e-01,
       3.12889646e-16,  9.51056516e-01,  5.87785252e-01, -5.87785252e-01,
      -9.51056516e-01, -6.20717256e-16,  9.51056516e-01,  5.87785252e-01,
      -5.87785252e-01, -9.51056516e-01
    ]
{- ORMOLU_ENABLE -}

-- From Python:
-- fs = 100  # Original sampling rate (Hz)
-- t = np.linspace(0, 1, fs, endpoint=False)  # 100 samples over 1 second
-- freq = 10  # Frequency of the sine wave (Hz)
-- signal_orig = np.sin(2 * np.pi * freq * t)
-- signal_up = signal.resample(signal_orig, 200)
{- ORMOLU_DISABLE -}
expectUp :: Vector Float
expectUp =
  V.fromList
    [ 4.41762107e-31,  3.09016994e-01,  5.87785252e-01,  8.09016994e-01,
        9.51056516e-01,  1.00000000e+00,  9.51056516e-01,  8.09016994e-01,
        5.87785252e-01,  3.09016994e-01,  8.60966643e-17, -3.09016994e-01,
       -5.87785252e-01, -8.09016994e-01, -9.51056516e-01, -1.00000000e+00,
       -9.51056516e-01, -8.09016994e-01, -5.87785252e-01, -3.09016994e-01,
       -1.92678075e-16,  3.09016994e-01,  5.87785252e-01,  8.09016994e-01,
        9.51056516e-01,  1.00000000e+00,  9.51056516e-01,  8.09016994e-01,
        5.87785252e-01,  3.09016994e-01,  4.41368032e-16, -3.09016994e-01,
       -5.87785252e-01, -8.09016994e-01, -9.51056516e-01, -1.00000000e+00,
       -9.51056516e-01, -8.09016994e-01, -5.87785252e-01, -3.09016994e-01,
       -4.89858720e-16,  3.09016994e-01,  5.87785252e-01,  8.09016994e-01,
        9.51056516e-01,  1.00000000e+00,  9.51056516e-01,  8.09016994e-01,
        5.87785252e-01,  3.09016994e-01,  5.75955384e-16, -3.09016994e-01,
       -5.87785252e-01, -8.09016994e-01, -9.51056516e-01, -1.00000000e+00,
       -9.51056516e-01, -8.09016994e-01, -5.87785252e-01, -3.09016994e-01,
       -6.82536794e-16,  3.09016994e-01,  5.87785252e-01,  8.09016994e-01,
        9.51056516e-01,  1.00000000e+00,  9.51056516e-01,  8.09016994e-01,
        5.87785252e-01,  3.09016994e-01, -2.62148693e-15, -3.09016994e-01,
       -5.87785252e-01, -8.09016994e-01, -9.51056516e-01, -1.00000000e+00,
       -9.51056516e-01, -8.09016994e-01, -5.87785252e-01, -3.09016994e-01,
       -9.79717439e-16,  3.09016994e-01,  5.87785252e-01,  8.09016994e-01,
        9.51056516e-01,  1.00000000e+00,  9.51056516e-01,  8.09016994e-01,
        5.87785252e-01,  3.09016994e-01,  1.06581410e-15, -3.09016994e-01,
       -5.87785252e-01, -8.09016994e-01, -9.51056516e-01, -1.00000000e+00,
       -9.51056516e-01, -8.09016994e-01, -5.87785252e-01, -3.09016994e-01,
       -1.17239551e-15,  3.09016994e-01,  5.87785252e-01,  8.09016994e-01,
        9.51056516e-01,  1.00000000e+00,  9.51056516e-01,  8.09016994e-01,
        5.87785252e-01,  3.09016994e-01, -2.13162821e-15, -3.09016994e-01,
       -5.87785252e-01, -8.09016994e-01, -9.51056516e-01, -1.00000000e+00,
       -9.51056516e-01, -8.09016994e-01, -5.87785252e-01, -3.09016994e-01,
       -1.46957616e-15,  3.09016994e-01,  5.87785252e-01,  8.09016994e-01,
        9.51056516e-01,  1.00000000e+00,  9.51056516e-01,  8.09016994e-01,
        5.87785252e-01,  3.09016994e-01, -1.99704086e-15, -3.09016994e-01,
       -5.87785252e-01, -8.09016994e-01, -9.51056516e-01, -1.00000000e+00,
       -9.51056516e-01, -8.09016994e-01, -5.87785252e-01, -3.09016994e-01,
        5.44317312e-15,  3.09016994e-01,  5.87785252e-01,  8.09016994e-01,
        9.51056516e-01,  1.00000000e+00,  9.51056516e-01,  8.09016994e-01,
        5.87785252e-01,  3.09016994e-01,  5.46365787e-15, -3.09016994e-01,
       -5.87785252e-01, -8.09016994e-01, -9.51056516e-01, -1.00000000e+00,
       -9.51056516e-01, -8.09016994e-01, -5.87785252e-01, -3.09016994e-01,
       -1.95943488e-15,  3.09016994e-01,  5.87785252e-01,  8.09016994e-01,
        9.51056516e-01,  1.00000000e+00,  9.51056516e-01,  8.09016994e-01,
        5.87785252e-01,  3.09016994e-01,  5.59824522e-15, -3.09016994e-01,
       -5.87785252e-01, -8.09016994e-01, -9.51056516e-01, -1.00000000e+00,
       -9.51056516e-01, -8.09016994e-01, -5.87785252e-01, -3.09016994e-01,
       -2.15211295e-15,  3.09016994e-01,  5.87785252e-01,  8.09016994e-01,
        9.51056516e-01,  1.00000000e+00,  9.51056516e-01,  8.09016994e-01,
        5.87785252e-01,  3.09016994e-01, -1.15191077e-15, -3.09016994e-01,
       -5.87785252e-01, -8.09016994e-01, -9.51056516e-01, -1.00000000e+00,
       -9.51056516e-01, -8.09016994e-01, -5.87785252e-01, -3.09016994e-01
    ]
{- ORMOLU_ENABLE -}
