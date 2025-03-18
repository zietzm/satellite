module Norm
  ( maxNorm,
    meanNorm,
    rangeNorm,
    percentile,
    sortStorableVector,
  )
where

import Control.Monad.ST (runST)
import qualified Data.Vector.Algorithms.Intro as ALG
import Data.Vector.Storable (Storable, Vector)
import qualified Data.Vector.Storable as V
import qualified Signals

-- | Normalize to have a max absolute value of one
maxNorm :: (Ord a, Fractional a, Storable a) => Vector a -> Vector a
maxNorm xs = V.map (/ maxAbs) xs
  where
    maxAbs = V.maximum $ V.map abs xs

-- | Normalize to have zero mean
meanNorm :: (Fractional a, Storable a) => Vector a -> Vector a
meanNorm xs = V.map (subtract meanVal) xs
  where
    n = fromIntegral $ V.length xs
    meanVal = V.sum xs / n

-- | Clamp values between percentiles and set the range to 0 - 1
rangeNorm :: (Fractional a, Ord a, Storable a) => a -> a -> Vector a -> Vector a
rangeNorm lowPct highPct xs = V.map adjustVal xs
  where
    percentiles = percentile [lowPct, highPct] xs
    lowVal = head percentiles
    highVal = last percentiles
    valRange = highVal - lowVal + 1e-8
    adjustVal x
      | x < lowVal = 0
      | x > highVal = 1
      | otherwise = (x - lowVal) / valRange

-- | Compute percentiles of a vector
percentile :: (Fractional a, Ord a, Storable a) => [a] -> Vector a -> [a]
percentile ps xs = V.toList newY
  where
    n = V.length xs
    nF = (fromRational . toRational) n
    psIx = map (\x -> (x / 100.0) * (nF - 1)) ps
    xNew = V.fromList psIx
    xOld = V.enumFromN 0 n
    yOld = sortStorableVector xs
    newY = Signals.interp xOld yOld xNew

-- | Sort a vector
sortStorableVector :: (Ord a, Storable a) => Vector a -> Vector a
sortStorableVector vec = runST $ do
  mvec <- V.thaw vec
  ALG.sort mvec
  V.freeze mvec
