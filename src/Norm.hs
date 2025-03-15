module Norm
  ( maxNorm,
    meanNorm,
  )
where

import Data.Vector.Storable (Storable, Vector)
import qualified Data.Vector.Storable as V

-- | Normalize to have a max absolute value of one
maxNorm :: (Ord a, Fractional a, Storable a) => Vector a -> Vector a
maxNorm xs = V.map (/ maxAbs) xs
  where
    maxAbs = V.maximum $ V.map abs xs

-- | Normalize to have zero mean
meanNorm :: (Num a, Storable a) => Vector a -> Vector a
meanNorm xs = V.map (subtract meanVal) xs
  where
    meanVal = V.sum xs
