module Sync
  ( syncA,
    syncB,
    adjustPattern,
    crossCorr,
  )
where

import Data.Vector.Storable (Storable, Vector)
import qualified Data.Vector.Storable as V

-- | @wordsPerLine@ the number of words per line in the NOAA APT format
wordsPerLine :: Int
wordsPerLine = 2080

-- | @linesPerSecond@ the number of lines per second in the NOAA APT format
linesPerSecond :: Int
linesPerSecond = 2

{- ORMOLU_DISABLE -}
-- | @syncA@ the A synchronization pattern for NOAA APT images
syncA :: Vector Float
syncA = V.fromList [
       -1, -1, -1, -1,
       1, 1, -1, -1,  -- Pulse 1
       1, 1, -1, -1,  -- Pulse 2
       1, 1, -1, -1,  -- Pulse 3
       1, 1, -1, -1,  -- Pulse 4
       1, 1, -1, -1,  -- Pulse 5
       1, 1, -1, -1,  -- Pulse 6
       1, 1, -1, -1,  -- Pulse 7
       -1, -1, -1, -1,
       -1, -1, -1, -1
   ]

-- | @syncB@ the B synchronization pattern for NOAA APT images
syncB :: Vector Float
syncB = V.fromList [
       -1, -1, -1, -1,
       1, 1, 1, -1, -1,  -- Pulse 1
       1, 1, 1, -1, -1,  -- Pulse 2
       1, 1, 1, -1, -1,  -- Pulse 3
       1, 1, 1, -1, -1,  -- Pulse 4
       1, 1, 1, -1, -1,  -- Pulse 5
       1, 1, 1, -1, -1,  -- Pulse 6
       1, 1, 1, -1, -1,  -- Pulse 7
       -1
   ]
{- ORMOLU_ENABLE -}

-- | @adjustPattern@ increase the size of a pattern by repeating elements
-- e.g. adjustPattern 2 (a,b) -> (a,a,b,b)
adjustPattern :: Int -> Vector Float -> Vector Float
adjustPattern upFactor = V.concatMap (V.replicate upFactor)

-- | @crossCorr@ computes the cross-correlation between two vectors
crossCorr :: (Num a, Storable a) => Vector a -> Vector a -> Vector a
crossCorr fs xs = V.generate n getElem
  where
    n1 = V.length fs
    n2 = V.length xs
    n = max n1 n2 - min n1 n2 + 1
    getElem i = V.sum $ V.zipWith (*) fs (V.drop i xs)
