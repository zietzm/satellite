module Main (main) where

import Data.Function ((&))
import qualified Data.Vector.Storable as V
import qualified Norm
import qualified Signals
import qualified Sync
import qualified Wav

actualDecode :: IO ()
actualDecode = do
  (info, rawSamples) <- Wav.readWAV "/Users/zietzm/projects/sdr/noaa/noaa19_short.wav"
  putStrLn $ "Info: " ++ show info

  let fs = fromIntegral $ Wav.sampleRate info :: Int
  let sampleRate = fromIntegral fs
  putStrLn $ "Sample rate: " ++ show sampleRate

  let factor = 4 * 2 * 2080 / sampleRate
  let upsampleFactor = fromIntegral (V.length rawSamples) * factor
  putStrLn $ "Upsample factor: " ++ show upsampleFactor

  let processed =
        rawSamples
          & Norm.maxNorm
          & Signals.getEnvelope
          -- This step is messed up
          & V.toList
          & Signals.lowpass sampleRate 2080
          & V.fromList
          -- Back on track
          & Signals.fourierResample (round upsampleFactor)
          & Norm.meanNorm

  putStrLn $ "Processed: " ++ show (V.take 5 processed)
  putStrLn $ "N processed: " ++ show (V.length processed)

  let aPulses = Sync.crossCorr Sync.syncA processed
  let bPulses = Sync.crossCorr Sync.syncB processed

  putStrLn $ "Sync A: " ++ show (V.take 5 aPulses)
  putStrLn $ "Sync B: " ++ show (V.take 5 bPulses)

  let newFs = sampleRate * factor / 2
  -- let samplesPerLine = newFs - newFs `mod` 2
  let samplesPerLine = newFs

  let aMinH = 1.5 * V.sum aPulses / fromIntegral (V.length aPulses)
  let bMinH = 1.5 * V.sum bPulses / fromIntegral (V.length bPulses)
  let distance = round (0.8 * samplesPerLine)
  let aPeaks = Sync.findPeaks aMinH distance aPulses
  let bPeaks = Sync.findPeaks bMinH distance bPulses
  putStrLn $ "Peaks A: " ++ show (V.take 5 aPeaks)
  putStrLn $ "Peaks B: " ++ show (V.take 5 bPeaks)

  let (aPeaks', bPeaks') = Sync.alignPeaks aPeaks bPeaks
  putStrLn $ "Aligned A: " ++ show (V.take 5 aPeaks')
  putStrLn $ "Aligned B: " ++ show (V.take 5 bPeaks')

  putStrLn $ "Inputs: " ++ show (round samplesPerLine, V.head aPeaks', V.head bPeaks')
  let sampleCoords = Sync.getSampleCoords (round samplesPerLine) 2080 aPeaks' bPeaks' :: V.Vector Float
  putStrLn $ "Sample coordinates: " ++ show (V.take 5 sampleCoords)

-- Interpolate the image at some given points (downsample)

main :: IO ()
main = actualDecode
