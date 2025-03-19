module Sync
  ( syncA,
    syncB,
    upsamplePattern,
    crossCorr,
    findPeaks,
    findBasicPeaks,
    prioritizeHighest,
    alignPeaks,
    getSampleCoords,
    linesPerSecond,
    wordsPerLine,
  )
where

import Data.List (sort, sortBy)
import Data.Vector.Storable (Storable, Vector)
import qualified Data.Vector.Storable as V

-- | @wordsPerLine@ the number of words per line in the NOAA APT format
wordsPerLine :: Float
wordsPerLine = 2080

-- | @linesPerSecond@ the number of lines per second in the NOAA APT format
linesPerSecond :: Float
linesPerSecond = 2

{- ORMOLU_DISABLE -}
-- | @syncA@ the A synchronization pattern for NOAA APT images
syncA :: (Num a, Storable a) => Vector a
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
syncB :: (Num a, Storable a) => Vector a
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
upsamplePattern :: (Storable a) => Int -> Vector a -> Vector a
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

    filtered = sort $ filterPeaks sortedPeaks []

-- | @alignPeaks@ ensure that A peaks precede B peaks
alignPeaks :: (Ord a, Storable a) => Vector a -> Vector a -> (Vector a, Vector a)
alignPeaks a b
  | V.null a || V.null b = (a, b)
  | V.head a >= V.head b = alignPeaks a (V.drop 1 b)
  | otherwise = (a, b)

getSampleCoords :: Int -> Int -> Vector Int -> Vector Int -> Vector Float
getSampleCoords rawNPerLine nPerLine a b = go a b V.empty
  where
    nPerSide = nPerLine `div` 2
    step = fromIntegral rawNPerLine / fromIntegral nPerLine :: Float

    go x y acc
      | V.null x || V.null y = acc
      | otherwise = go x' y' (V.concat [acc, newValsX, newValsY])
      where
        xPeak = fromIntegral $ V.head x :: Float
        yPeak = fromIntegral $ V.head y :: Float
        newValsX = V.generate nPerSide (\i -> xPeak + fromIntegral i * step)
        newValsY = V.generate nPerSide (\i -> yPeak + fromIntegral i * step)
        x' = V.drop 1 x
        y' = V.drop 1 y
