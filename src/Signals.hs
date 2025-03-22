module Signals
  ( getEnvelope,
    fft,
    fftReal,
    ifft,
    ifftMagnitude,
    fftN,
    fftNReal,
    ifftReal,
    filtFilt,
    lowpassVec,
    fourierResample,
    interpCoords,
  )
where

import DSP.Filter.IIR.Cookbook (lpf)
import Data.Complex (Complex ((:+)))
import qualified Data.Complex as C
import Data.Vector.Storable (Storable, Vector)
import qualified Data.Vector.Storable as V
import qualified Numeric.FFT.Vector.Invertible as FFTW

fft :: (Real a, Storable a) => Vector a -> Vector (Complex Double)
fft = FFTW.run FFTW.dft . toFftType

fftReal :: (Real a, Storable a) => Vector a -> Vector (Complex Double)
fftReal = FFTW.run FFTW.dftR2C . V.map realToFrac

ifft :: (RealFloat a, Storable a) => Vector (Complex Double) -> Vector a
ifft = V.map (realToFrac . C.realPart) . FFTW.run FFTW.idft

ifftMagnitude :: (RealFloat a, Storable a) => Vector (Complex Double) -> Vector a
ifftMagnitude = V.map (realToFrac . C.magnitude) . FFTW.run FFTW.idft

toFftType :: (Real a, Storable a) => V.Vector a -> V.Vector (Complex Double)
toFftType = V.map toComplexDouble

toComplexDouble :: (Real a) => a -> Complex Double
toComplexDouble x = realToFrac x :+ 0

fftN :: (Real a, Storable a) => Int -> Vector a -> Vector (Complex Double)
fftN n xs = fft padded
  where
    nPad = n - V.length xs
    padded = xs V.++ V.replicate nPad 0

fftNReal :: (RealFrac a, Storable a) => Int -> Vector a -> Vector (Complex Double)
fftNReal n xs = FFTW.run FFTW.dftR2C padded
  where
    nPad = n - V.length xs
    padded = V.map realToFrac xs V.++ V.replicate nPad 0

ifftReal :: (RealFrac a, Storable a) => Vector (Complex Double) -> Vector a
ifftReal = V.map realToFrac . FFTW.run FFTW.dftC2R

-- | Compute the envelope of a signal using FFT-based Hilbert transform
getEnvelope :: (RealFloat a, Storable a) => Vector a -> Vector a
getEnvelope signal = envelope
  where
    n = V.length signal
    fftResult = fftReal signal
    nFft = n `div` 2 + 1
    nToPad = n - nFft

    -- Create mask for the Hilbert transform
    hilbertMask :: Int -> Complex Double -> Complex Double
    hilbertMask i x
      | i == 0 = x -- DC
      | i < n `div` 2 = 2.0 * x -- Positive frequencies
      | i == n `div` 2 = x -- Nyquist
      | otherwise = 0.0 -- Negative frequencies

    -- Apply multiplier to get analytic signal in frequency domain
    analyticFreq = V.imap hilbertMask fftResult V.++ V.replicate nToPad (0 :+ 0)

    -- Inverse FFT to get analytic signal in time domain, get magnitude
    envelope = ifftMagnitude analyticFreq

-- | Apply a low pass filter to data with a given sampling and cutoff frequency
filtFilt :: (Floating a) => a -> a -> a -> [a] -> [a]
filtFilt bw fs fc = reverse . lpf bw w . reverse . lpf bw w
  where
    w = 2 * pi * fc / fs -- Normalized angular frequency

-- | Apply a low pass filter to data. Params are bandwidth and sampling and cutoff
-- frequencies.
lowpassVec :: (Floating a, Storable a) => a -> a -> a -> Vector a -> Vector a
lowpassVec bw fs fc = V.fromList . filtFilt bw fs fc . V.toList

fourierResample :: (RealFrac a, Storable a) => Int -> Vector a -> Vector a
fourierResample nNew xs = result
  where
    nOrig = V.length xs
    mOrig = nOrig `div` 2 + 1
    mNew = nNew `div` 2 + 1
    xs' = V.map realToFrac xs
    fftX = FFTW.run FFTW.dftR2C xs'
    changeLen vec
      | mNew <= mOrig = V.take mNew vec
      | otherwise = vec V.++ V.replicate (mNew - mOrig) 0.0
    fftX' = changeLen fftX
    resampled = FFTW.run FFTW.dftC2R fftX'
    resampled' = V.map (* (fromIntegral nNew / fromIntegral nOrig)) resampled
    result = V.map realToFrac resampled'

-- | @interpCoords@ interpolate a 1D signal at new coordinates, assuming that the index of the
-- input vector is the x-coordinate, and the new x-coordinates are on this same
-- scale.
interpCoords :: (RealFrac a, Storable a) => Vector a -> Vector a -> Vector a
interpCoords yOld xNew
  | V.null yOld || V.null xNew = V.empty
  | otherwise = V.map getNearest xNew
  where
    nOld = fromIntegral $ V.length yOld
    getNearest x
      | x <= 0.0 = V.head yOld
      | x >= nOld - 1 = V.last yOld
      | otherwise = y0 + (x - x0) * slope
      where
        x0 = fromIntegral (floor x :: Int)
        y0 = yOld V.! floor x
        y1 = yOld V.! (floor x + 1)
        slope = y1 - y0
