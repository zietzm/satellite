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
  case W.waveSampleFormat wave of
    W.SampleFormatIeeeFloat32Bit -> do
      let offset = fromIntegral $ W.waveDataOffset wave
          size = fromIntegral $ W.waveDataSize wave
      audioBytes <- withFile path ReadMode $ \handle -> do
        hSeek handle AbsoluteSeek offset
        BS.hGet handle size
      case decodeFloat32Samples audioBytes of
        Left err -> error $ "Failed to decode samples: " ++ err
        Right samples -> pure (toWAVInfo wave, samples)
    _ -> error "WAV file must be 32-bit float format"

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
decodeFloat32Samples :: BS.ByteString -> Either String (V.Vector Float)
decodeFloat32Samples bs = S.runGet (V.replicateM nSamp S.getFloat32le) bs
  where
    nSamp = BS.length bs `div` 4 -- 4 bytes per Float32

-- | Calculate duration in seconds
duration :: WAVInfo -> V.Vector Sample -> Float
duration info samples = fromIntegral (V.length samples) / fromIntegral (sampleRate info)

-- | Get the number of channels
channels :: WAVInfo -> Word16
channels = channelCount

nSamples :: (Storable a) => V.Vector a -> Int
nSamples = V.length
