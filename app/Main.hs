module Main (main) where

import Data.Function ((&))
import qualified Data.Vector.Storable as V
import qualified Export
import qualified Norm
import qualified Signals
import qualified Sync
import qualified Wav

actualDecode :: IO ()
actualDecode = do
  let upFactor = 4
  (info, rawSamples) <- Wav.readWAV "/Users/zietzm/projects/sdr/noaa/noaa19_short.wav"
  putStrLn $ "Info: " ++ show info

  let fs = fromIntegral $ Wav.sampleRate info :: Int
  let sampleRate = fromIntegral fs :: Float
  putStrLn $ "Sample rate: " ++ show sampleRate

  let factor = upFactor * 2 * Sync.wordsPerLine / sampleRate
  let nResampled = round $ fromIntegral (V.length rawSamples) * factor
  putStrLn $ "Number of resampled points: " ++ show nResampled

  let processed =
        rawSamples
          & Norm.maxNorm
          & Signals.getEnvelope
          & Signals.lowpassVec 1.5 sampleRate Sync.wordsPerLine
          & Signals.fourierResample nResampled
          & Norm.meanNorm

  putStrLn $ "Processed: " ++ show (V.take 5 processed)
  putStrLn $ "N processed: " ++ show (V.length processed)

  let filterA = Sync.upsamplePattern 4 Sync.syncA
  let filterB = Sync.upsamplePattern 4 Sync.syncB
  let aPulses = Sync.crossCorr filterA processed
  let bPulses = Sync.crossCorr filterB processed

  putStrLn $ "Sync A: " ++ show (V.take 5 aPulses)
  putStrLn $ "Sync B: " ++ show (V.take 5 bPulses)
  putStrLn $ "Sync A length: " ++ show (V.length aPulses)
  putStrLn $ "Sync B length: " ++ show (V.length bPulses)

  let newFs = sampleRate * factor / Sync.linesPerSecond
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
  putStrLn $ "Sample coordinates: " ++ show (V.take 5 sampleCoords) ++ show (V.last sampleCoords)
  putStrLn $ "Sample coordinates length: " ++ show (V.length sampleCoords)

  -- Interpolate the image at some given points (downsample)
  let xOld = V.enumFromN 0 (V.length processed)
  let interpolated = Signals.interp xOld processed sampleCoords
  putStrLn $ "Interpolated: " ++ show (V.take 5 interpolated) ++ show (V.last interpolated)
  putStrLn $ "Interpolated length: " ++ show (V.length interpolated)

  let normed = Norm.rangeNorm 1 99 interpolated
  let final = V.map (\x -> fromInteger $ floor (x * 255)) normed
  putStrLn $ "Final: " ++ show (V.take 5 final) ++ " " ++ show (V.last final) ++ " " ++ show (V.maximum final)

  let img = Export.resolveImage final 2080 (V.length final `div` 2080)
  result <- Export.saveImg "/Users/zietzm/projects/satellite/TEST_HS.png" img
  case result of
    Left err -> putStrLn $ "Error writing PNG: " ++ err
    Right True -> putStrLn "PNG written successfully!"
    Right False -> putStrLn "Failed to write PNG, but no specific error."

main :: IO ()
main = actualDecode
