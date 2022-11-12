#pragma once

#include "Image/ImageAlgorithms.h"

std::shared_ptr<baseImage> DoBlur(std::shared_ptr<baseImage> curImg, int filtWid, float imageStDev);
std::shared_ptr<baseImage> DoChanConvert(std::shared_ptr<baseImage> curImg, int nchan);
std::shared_ptr<baseImage> DoDespeckle(std::shared_ptr<baseImage> curImg);
std::shared_ptr<baseImage> DoGrow(std::shared_ptr<baseImage> curImg);
std::shared_ptr<baseImage> DoNoise(std::shared_ptr<baseImage> curImg, float mean, float stdev);
std::shared_ptr<baseImage> DoResize(std::shared_ptr<baseImage> curImg, int newWid, int newHgt);
std::shared_ptr<baseImage> DoResizeCA(std::shared_ptr<baseImage> curImg, int newWid, int newHgt);
std::shared_ptr<baseImage> DoRotate90(std::shared_ptr<baseImage> curImg);
std::shared_ptr<baseImage> DoRingDespeckleFilter(std::shared_ptr<baseImage> curImg, int filtWid, int clearWid, float threshold, int iterations);
std::shared_ptr<baseImage> DoVCD(std::shared_ptr<baseImage> curImg, int filtWid, float imageStDev, float colorStDev, int iterations);
template <class Image_T> void DoDiff(Image_T& imgA, const Image_T& imgB, const float scale);
template <class Image_T> void DoFlatten(Image_T& dstImg, const Image_T& srcImg, const float shrinkFac, const float biasFac);
template <class Image_T> void DoFlattenHoriz(Image_T& dstImg, const Image_T& srcImg, const float a, const float b);
template <class Image_T> void DoLinearMap(Image_T& Img, const f3Pixel p0, const f3Pixel p1);
template <class Image_T>
void DoMatMul(Image_T& Img, float m0, float m1, float m2, float m3, float m4, float m5, float m6, float m7, float m8, float m9, float ma, float mb, float mc,
              float md, float me, float mf);
template <class Image_T> void DoOptimizeLinearMap(const Image_T& Img);
template <class Image_T> void DoPrintTransect(Image_T& Img, const int y);
template <class Image_T> void DoThreshold(Image_T& Img, const typename Image_T::PixType threshold);
void DoCopyBlock(std::shared_ptr<baseImage> curImg, const std::shared_ptr<baseImage> srcImg, int dx, int dy, int sx, int sy, int sw, int sh, int mode,
                 uc4Pixel key, float alpha);
void DoFill(std::shared_ptr<baseImage> curImg, float r, float g, float b, float a);
void DoFlip(std::shared_ptr<baseImage> curImg, bool isVert);
void DoGradient(std::shared_ptr<baseImage> curImg, int c, bool isVert, int minx, int maxx);
