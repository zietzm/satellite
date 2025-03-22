{-# LANGUAGE DeriveDataTypeable #-}

module Main (main) where

import Control.Monad (unless)
import Control.Monad.Trans.Class (lift)
import Control.Monad.Trans.Except (ExceptT, except, runExceptT, throwE)
import Data.Bifunctor (first)
import qualified Export
import qualified Lib
import System.Console.CmdArgs ((&=))
import qualified System.Console.CmdArgs as Cmd
import qualified Wav

data Decoder = Decoder
  { input :: FilePath,
    output :: FilePath,
    verbose :: Bool,
    upsampleFactor :: Int
  }
  deriving (Show, Cmd.Data, Cmd.Typeable)

decoder :: Decoder
decoder =
  Decoder
    { input = Cmd.def &= Cmd.argPos 0 &= Cmd.typ "INPUT",
      output = Cmd.def &= Cmd.argPos 1 &= Cmd.typ "OUTPUT",
      verbose = Cmd.def &= Cmd.help "Print verbose output",
      upsampleFactor = 4 &= Cmd.opt "4" &= Cmd.help "Upsampling level at which to perform computations"
    }
    &= Cmd.summary "SatelliteDecoder v1.0 - Decode NOAA APT files from WAV to PNG"
    &= Cmd.program "satellite-decoder"

main :: IO ()
main = runExceptT happyPath >>= either (putStrLn . ("Error: " ++)) (const $ putStrLn "PNG written successfully!")
  where
    happyPath :: ExceptT String IO ()
    happyPath = do
      args <- lift $ Cmd.cmdArgs decoder
      let infile = input args
      (info, rawSamples) <- lift $ Wav.readWAV infile
      let upFactor = upsampleFactor args
      img <- except $ first show $ Lib.decodeApt upFactor info rawSamples
      let outfile = output args
      success <- except =<< lift (Export.saveImg outfile img)
      unless success $ throwE "Failed to write PNG"
