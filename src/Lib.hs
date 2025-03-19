module Lib (decodeApt) where

import Codec.Picture.Types as Img
import Data.Function ((&))
import qualified Data.Vector.Storable as V
import qualified Export
import qualified Norm
import qualified Signals
import qualified Sync
import qualified Wav

data DecodeError
  = InvalidUpsample Int
  | NoPeaks
  | InvalidSampleRate Float
  | InvalidSampleCount Int
  | ProcessingError String
  | InvalidInterpolated String
  deriving (Show)

decodeApt :: Int -> Wav.WAVInfo -> V.Vector Wav.Sample -> Either DecodeError Img.DynamicImage
decodeApt upsample info samples = do
  if upsample <= 0
    then Left $ InvalidUpsample upsample
    else Right ()

  -- Compute some size factors used for resampling the signal
  let fs = fromIntegral $ Wav.sampleRate info :: Float
      nyquist = fs / 2
      factor = fromIntegral upsample * Sync.wordsPerLine / nyquist
      nSamples = V.length samples
      nResampled = round $ fromIntegral nSamples * factor

      -- Initial signal processing: normalization, Hilbert transform,
      -- lowpass filter, resample, then normalize again.
      processed =
        samples
          & Norm.maxNorm
          & Signals.getEnvelope
          & Signals.lowpassVec 1.5 fs Sync.wordsPerLine
          & Signals.fourierResample nResampled
          & Norm.meanNorm

      newFs = fs * factor / Sync.linesPerSecond
      (aPeaks, bPeaks) = findSyncPeaks upsample newFs processed

  if V.length aPeaks <= 0 || V.length bPeaks <= 0
    then Left NoPeaks
    else Right ()

  -- Interpolate between peaks to downsample without aliasing
  let xNew = Sync.getSampleCoords (round newFs) 2080 aPeaks bPeaks
      interpolated = Signals.interpCoords processed xNew

  if V.length interpolated <= 0
    then Left $ InvalidInterpolated "Interpolated vector is empty"
    else Right ()

  -- Clip brightness values and scale to increase contrast
  let normed = Norm.rangeNorm 1 99 interpolated
      -- Convert to grayscale image
      final = V.map (\x -> fromInteger $ floor (x * 255)) normed
      img = Export.resolveImage final 2080 (V.length final `div` 2080)

  Right img

findSyncPeaks :: (RealFrac a, V.Storable a) => Int -> a -> V.Vector a -> (V.Vector Int, V.Vector Int)
findSyncPeaks upfactor fs xs = (aPeaks', bPeaks')
  where
    -- Convolve with sync patterns to find line starts/stops
    filterA = Sync.upsamplePattern upfactor Sync.syncA
    filterB = Sync.upsamplePattern upfactor Sync.syncB
    aPulses = Sync.crossCorr filterA xs
    bPulses = Sync.crossCorr filterB xs

    -- Find peaks in the convolution outputs
    aMinH = 1.5 * V.sum aPulses / fromIntegral (V.length aPulses)
    bMinH = 1.5 * V.sum bPulses / fromIntegral (V.length bPulses)
    distance = round (0.8 * fs)
    aPeaks = Sync.findPeaks aMinH distance aPulses
    bPeaks = Sync.findPeaks bMinH distance bPulses
    (aPeaks', bPeaks') = Sync.alignPeaks aPeaks bPeaks
