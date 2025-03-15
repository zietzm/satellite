import Data.Complex (Complex ((:+)))
import qualified Data.Complex as C
import qualified Data.Vector.Storable as V
import qualified Data.Vector.Unboxed as U
import GHC.Float (double2Float, float2Double)
import qualified Signals
import qualified Sync
import Test.Tasty
import Test.Tasty.HUnit

main :: IO ()
main = defaultMain tests

tests :: TestTree
tests = testGroup "Tests" [dspTests, syncTests]

dspTests :: TestTree
dspTests =
  testGroup
    "Signal processing"
    [ testCase "FFT expected length" $
        U.length (Signals.fft wave) @?= V.length wave,
      testCase "Envelope expected length" $
        V.length (Signals.getEnvelope wave) @?= V.length wave,
      testCase "FFT expected result" $
        assertBool ("FFT differed from Python: " ++ show maxDiffFft) (maxDiffFft < 1e-6),
      testCase "Wave envelope" $
        assertBool ("Envelope differed from Python: " ++ show maxDiffEnv) (maxDiffEnv < 1e-6),
      -- TODO: We could make this more accurate in the future. For now, it works.
      testCase "Low pass filter" $
        assertBool ("Filter didn't work right" ++ show maxDiffFilt) (maxDiffFilt < 2e-3)
    ]
  where
    maxDiffFft = maxDiffD (Signals.fft wave) waveFft
    maxDiffEnv = maxDiffF (Signals.getEnvelope wave) waveEnvelope

    -- Test the low pass filter
    t = V.enumFromN 0 2000 :: V.Vector Float
    xLow = V.map (sin . (* (2 * 5 * pi))) t
    xHigh = V.map (sin . (* (2 * 250 * pi))) t
    x = V.zipWith (+) xLow xHigh
    xL = map float2Double $ V.toList x
    filtered = Signals.lowpass 2000.0 0.125 xL
    filtered' = V.fromList $ map double2Float filtered
    maxDiffFilt = maxDiffF (V.drop 10 filtered') (V.drop 10 xLow)

syncTests :: TestTree
syncTests =
  testGroup
    "Synchronization"
    [ testCase "Self-convolve = length for syncA" $
        Sync.crossCorr Sync.syncA Sync.syncA @?= V.fromList [40.0],
      testCase "Self-convolve = length for syncB" $
        Sync.crossCorr Sync.syncB Sync.syncB @?= V.fromList [40.0],
      testCase "Adjust pattern 2 [1,2] -> [1,1,2,2]" $
        Sync.adjustPattern 2 (V.fromList [1.0, 2.0]) @?= V.fromList [1.0, 1.0, 2.0, 2.0]
    ]

maxDiffF :: V.Vector Float -> V.Vector Float -> Float
maxDiffF x y = V.maximum $ V.zipWith (\a b -> abs a - b) x y

maxDiffD :: U.Vector (Complex Double) -> U.Vector (Complex Double) -> Double
maxDiffD x y = U.maximum $ U.map C.magnitude $ U.zipWith (-) x y

-- Built in Python using:
-- t = np.linspace(0, 1, 32)
-- signal = np.sin(2 * np.pi * 5 * t) + 0.5 * np.sin(2 * np.pi * 20 * t)
{- ORMOLU_DISABLE -}
wave :: V.Vector Float
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
waveFft :: U.Vector (Complex Double)
waveFft =
  U.fromList
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
waveEnvelope :: V.Vector Float
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
