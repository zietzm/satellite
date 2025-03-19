module Signals
  ( getEnvelope,
    fft,
    ifftMagnitude,
    filtFilt,
    lowpassVec,
    fourierResample,
    interp,
  )
where

import DSP.Filter.IIR.Cookbook (lpf)
import Data.Complex (Complex ((:+)))
import qualified Data.Complex as C
import Data.Vector.Storable (Storable, Vector)
import qualified Data.Vector.Storable as V
import GHC.Float (double2Float, float2Double)
import qualified Numeric.FFT.Vector.Invertible as FFTW

fft :: Vector Float -> Vector (Complex Double)
fft = FFTW.run FFTW.dft . toFftType

ifftMagnitude :: Vector (Complex Double) -> Vector Float
ifftMagnitude = getMagnitudes . FFTW.run FFTW.idft

toFftType :: Vector Float -> Vector (Complex Double)
toFftType = V.map ((:+ 0) . float2Double)

getMagnitudes :: Vector (Complex Double) -> Vector Float
getMagnitudes = V.map (double2Float . C.magnitude)

-- | Compute the envelope of a signal using FFT-based Hilbert transform
getEnvelope :: Vector Float -> Vector Float
getEnvelope signal = envelope
  where
    -- Compute FFT
    fftResult = fft signal
    n = V.length fftResult -- Should be same as input length per docs

    -- Create mask for the Hilbert transform
    hilbertMask :: Int -> Complex Double -> Complex Double
    hilbertMask i x
      | i == 0 = x -- DC
      | i < n `div` 2 = 2.0 * x -- Positive frequencies
      | i == n `div` 2 = x -- Nyquist
      | otherwise = 0.0 -- Negative frequencies

    -- Apply multiplier to get analytic signal in frequency domain
    analyticFreq = V.imap hilbertMask fftResult

    -- Inverse FFT to get analytic signal in time domain, get magnitude
    envelope = ifftMagnitude analyticFreq

-- | Apply a low pass filter to data with a given sampling and cutoff frequency
filtFilt :: (Floating a) => a -> a -> a -> [a] -> [a]
filtFilt bw fs fc = reverse . lpf bw w . reverse . lpf bw w
  where
    w = 2 * pi * fc / fs -- Normalized angular frequency

lowpassVec :: (Floating a, Storable a) => a -> a -> a -> Vector a -> Vector a
lowpassVec bw fs fc = V.fromList . filtFilt bw fs fc . V.toList

fourierResample :: Int -> Vector Float -> Vector Float
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

-- | @interp@ linearly interpolate a 1D signal at new coordinates
interp :: (Ord a, Fractional a, Storable a) => Vector a -> Vector a -> Vector a -> Vector a
interp xOld yOld xNew
  | V.null xOld || V.null yOld = V.fromList [] -- Handle empty inputs
  | V.length xOld /= V.length yOld = V.fromList [] -- Handle mismatched lengths
  | otherwise = V.fromList $ go 0 0 []
  where
    nNew = V.length xNew

    go newP oldP ys
      | newP >= nNew = reverse ys
      | x <= V.head xOld = go (newP + 1) oldP (V.head yOld : ys)
      | x >= V.last xOld = go (newP + 1) oldP (V.last yOld : ys)
      | xOld V.! oldP <= x && x <= xOld V.! (oldP + 1) = go (newP + 1) oldP (yNew : ys)
      | otherwise = go newP (oldP + 1) ys
      where
        x = xNew V.! newP
        slope =
          ((yOld V.! (oldP + 1)) - (yOld V.! oldP))
            / ((xOld V.! (oldP + 1)) - (xOld V.! oldP))
        yNew = (yOld V.! oldP) + (x - (xOld V.! oldP)) * slope
