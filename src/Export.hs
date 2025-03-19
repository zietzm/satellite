module Export (resolveImage, saveImg) where

import Codec.Picture.Png as Png
import Codec.Picture.Types as Img
import qualified Data.Vector.Storable as V
import Data.Word (Word8)

saveImg :: FilePath -> Img.DynamicImage -> IO (Either String Bool)
saveImg = Png.writeDynamicPng

type Color = Word8

resolveImage :: V.Vector Color -> Int -> Int -> Img.DynamicImage
resolveImage input w h =
  Img.ImageY8 $
    Img.Image
      { imageWidth = w,
        imageHeight = h,
        imageData = V.take (w * h) input
      }
