module Sync
  ( syncA,
    syncB,
    upsamplePattern,
    crossCorr,
    findPeaks,
    findBasicPeaks,
    prioritizeHighest,
  )
where

import Data.List (sortBy)
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
upsamplePattern :: Int -> Vector Float -> Vector Float
upsamplePattern upFactor = V.concatMap (V.replicate upFactor)

-- | @crossCorr@ computes the cross-correlation between two vectors
crossCorr :: (Num a, Storable a) => Vector a -> Vector a -> Vector a
crossCorr fs xs = V.generate n getElem
  where
    n1 = V.length fs
    n2 = V.length xs
    n = max n1 n2 - min n1 n2 + 1
    getElem i = V.sum $ V.zipWith (*) fs (V.drop i xs)

-- | @findPeaks@ finds peaks in a 1d signal that are at least height tall and
-- distance from another peak. Mimics scipy.signal.find_peaks (with fewer options)
findPeaks :: (Ord a, Storable a) => a -> Int -> Vector a -> Vector Int
findPeaks height distance xs = peaks
  where
    basicPeaks = findBasicPeaks height xs
    peakHeights = V.map (xs V.!) basicPeaks
    peaks = prioritizeHighest distance peakHeights basicPeaks

-- | Find indices that are higher than height and both adjacent neighbors
findBasicPeaks :: (Ord a, Storable a) => a -> Vector a -> Vector Int
findBasicPeaks height xs = V.filter (> 0) $ V.izipWith3 isPeak xs (V.drop 1 xs) (V.drop 2 xs)
  where
    isPeak i p c n
      | c >= height && p < c && n < c = i + 1
      | otherwise = 0 -- Zero indicates no peak (above is never zero)

-- | Filter peaks down to those which are a minimum distance apart. Prioritize
-- keeping higher peaks.
prioritizeHighest :: (Ord a, Storable a) => Int -> Vector a -> Vector Int -> Vector Int
prioritizeHighest distance heights peaks = V.fromList filtered
  where
    peakHeights = zip (V.toList peaks) (V.toList heights)
    sortedPairs = sortBy (\(_, h1) (_, h2) -> compare h2 h1) peakHeights
    sortedPeaks = map fst sortedPairs

    filterPeaks :: [Int] -> [Int] -> [Int]
    filterPeaks [] _ = []
    filterPeaks (p : ps) kept =
      let tooClose = [i | i <- kept, abs (p - i) < distance]
       in if null tooClose
            then p : filterPeaks ps (p : kept)
            else filterPeaks ps kept

    filtered = reverse $ filterPeaks sortedPeaks []
