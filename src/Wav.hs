{-# LANGUAGE DeriveGeneric #-}

module Wav
  ( -- * Types
    WAVInfo (..),
    Sample,

    -- * Loading
    readWAV,

    -- * Utilities
    duration,
    channels,
    nSamples,
  )
where

import qualified Codec.Audio.Wave as W
import Control.Monad.IO.Class (MonadIO, liftIO)
import qualified Data.ByteString as BS
import qualified Data.Serialize as S
import Data.Vector.Storable (Storable)
import qualified Data.Vector.Storable as V
import Data.Word (Word16, Word32)
import GHC.Generics (Generic)
import System.IO

-- | WAV file metadata
data WAVInfo = WAVInfo
  { -- | Samples per second (Hz)
    sampleRate :: Word32,
    -- | Number of channels
    channelCount :: Word16,
    -- | Bytes per second
    byteRate :: Word32,
    -- | Bits per second
    bitRate :: Double
  }
  deriving (Show, Generic)

-- | Sample type - only 32-bit float
type Sample = Float

-- | Read a WAV file, error if not 32-bit float
readWAV :: (MonadIO m) => FilePath -> m (WAVInfo, V.Vector Sample)
readWAV path = liftIO $ do
  wave <- W.readWaveFile path
  let offset = fromIntegral $ W.waveDataOffset wave
      size = fromIntegral $ W.waveDataSize wave
  audioBytes <- withFile path ReadMode $ \handle -> do
    hSeek handle AbsoluteSeek offset
    BS.hGet handle size
  case W.waveSampleFormat wave of
    W.SampleFormatIeeeFloat32Bit -> do
      case decodeFloat32Samples audioBytes of
        Left err -> error $ "Failed to decode samples: " ++ err
        Right samples -> pure (toWAVInfo wave, samples)
    (W.SampleFormatPcmInt _) -> do
      case decodePcmInt16Samples audioBytes of
        Left err -> error $ "Failed to decode samples: " ++ err
        Right samples -> pure (toWAVInfo wave, samples)
    fmt -> error $ "WAV file must be 32-bit float format, but found: " ++ show fmt

-- | Convert Wave metadata to WAVInfo
toWAVInfo :: W.Wave -> WAVInfo
toWAVInfo wave =
  WAVInfo
    { sampleRate = W.waveSampleRate wave,
      channelCount = W.waveChannels wave,
      byteRate = W.waveByteRate wave,
      bitRate = W.waveBitRate wave
    }

-- | Decode ByteString into a Vector of Float samples
decodeFloat32Samples :: BS.ByteString -> Either String (V.Vector Sample)
decodeFloat32Samples bs = S.runGet (V.replicateM nSamp S.getFloat32le) bs
  where
    nSamp = BS.length bs `div` 4 -- 4 bytes per Float32

-- | Decode 16-bit PCM samples and convert to Float
decodePcmInt16Samples :: BS.ByteString -> Either String (V.Vector Sample)
decodePcmInt16Samples bs = word16Vec
  where
    nSamp = BS.length bs `div` 2 -- 2 bytes per Int16
    word16Vec = S.runGet (V.replicateM nSamp getInt16AsFloat) bs

getInt16AsFloat :: S.Get Float
getInt16AsFloat = do
  w16 <- S.getInt16le -- Read a 16-bit unsigned integer (little-endian)
  return (fromIntegral w16 :: Float)

-- | Calculate duration in seconds
duration :: WAVInfo -> V.Vector Sample -> Float
duration info samples = fromIntegral (V.length samples) / fromIntegral (sampleRate info)

-- | Get the number of channels
channels :: WAVInfo -> Word16
channels = channelCount

nSamples :: (Storable a) => V.Vector a -> Int
nSamples = V.length
