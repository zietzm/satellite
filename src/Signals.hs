module Signals
  ( getEnvelope,
    fft,
    ifft,
    toFftType,
    fromFftType,
    lowpass,
    fourierResample,
  )
where

import DSP.Filter.IIR.Cookbook (lpf)
import Data.Complex (Complex ((:+)))
import qualified Data.Complex as C
import qualified Data.Vector.FFT as FFT
import qualified Data.Vector.Storable as V
import qualified Data.Vector.Unboxed as U
import GHC.Float (double2Float, float2Double)
import qualified Numeric.FFT.Vector.Invertible as FFTW

fft :: V.Vector Float -> U.Vector (Complex Double)
fft = FFT.fft . toFftType

ifft :: U.Vector (Complex Double) -> V.Vector Float
ifft = fromFftType . FFT.ifft

toFftType :: V.Vector Float -> U.Vector (Complex Double)
toFftType = U.convert . V.map ((:+ 0) . float2Double)

fromFftType :: U.Vector (Complex Double) -> V.Vector Float
fromFftType = V.convert . U.map (double2Float . C.magnitude)

-- | Compute the envelope of a signal using FFT-based Hilbert transform
getEnvelope :: V.Vector Float -> V.Vector Float
getEnvelope signal = envelope
  where
    -- Compute FFT
    fftResult = fft signal
    len = U.length fftResult -- Will be next power of 2 from original length

    -- Create mask for the Hilbert transform
    multiplier :: U.Vector (Complex Double)
    multiplier = U.generate len value
      where
        value i
          | i == 0 = 1 -- DC
          | i < len `div` 2 = 2 -- Positive frequencies
          | i == len `div` 2 = 1 -- Nyquist
          | otherwise = 0 -- Negative frequencies

    -- Apply multiplier to get analytic signal in frequency domain
    analyticFreq = U.zipWith (*) fftResult multiplier

    -- Inverse FFT to get analytic signal in time domain
    analyticSignal = FFT.ifft analyticFreq

    -- Slice to get the same number of outputs as inputs (not smallest 2^n above)
    sliced = U.take (V.length signal) analyticSignal
    envelope = fromFftType sliced

-- | Apply a low pass filter to data with a given sampling and cutoff frequency
lowpass :: Double -> Double -> [Double] -> [Double]
lowpass fs fc signal = reverse $ lpf bw w $ reverse $ lpf bw w signal
  where
    w = 2 * pi * fc / fs -- Normalized angular frequency
    bw = 0.707

fourierResample :: Int -> V.Vector Float -> V.Vector Float
fourierResample nNew xs = result
  where
    nOrig = V.length xs
    mOrig = nOrig `div` 2 + 1
    mNew = nNew `div` 2 + 1
    xs' = V.map float2Double xs
    fftX = FFTW.run FFTW.dftR2C xs'
    changeLen vec
      | mNew <= mOrig = V.take mNew vec
      | otherwise = vec V.++ V.replicate (mNew - mOrig) 0.0
    fftX' = changeLen fftX
    resampled = FFTW.run FFTW.dftC2R fftX'
    -- result = V.map double2Float resampled
    -- resampled' = V.map (* (fromIntegral nOrig / fromIntegral nNew)) resampled
    resampled' = V.map (* (fromIntegral nNew / fromIntegral nOrig)) resampled
    result = V.map double2Float resampled'
