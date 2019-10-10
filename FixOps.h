
#pragma once

#include "Image/ImageAlgorithms.h"

template <class old_T, class new_T>
new_T replaceImage(old_T oldImg, new_T newImg);

void DoCopyBlock(baseImage* curImg, std::string &blitFName, int dx, int dy, int sx, int sy, int sw, int sh, int mode, uc4Pixel key, float alpha);
void DoFill(baseImage* curImg, float r, float g, float b, float a);
void DoFlip(baseImage* curImg, bool isVert);
void DoGradient(baseImage* curImg, int c, bool isVert, int minx, int maxx);
baseImage* DoBlur(baseImage* curImg, int filtWid, float imageStDev);
baseImage* DoChanConvert(baseImage* curImg, int nchan);
baseImage* DoTemplateMatchFilter(baseImage* curImg, int iterations, int filtWid, int clearWid, float threshold);
baseImage* DoDespeckle(baseImage* curImg);
baseImage* DoNoise(baseImage* curImg, float mean, float stdev);
baseImage* DoResize(baseImage* curImg, int newWid, int newHgt);
baseImage* DoCAResize(baseImage* curImg, int newWid, int newHgt);
baseImage* DoVCD(baseImage* curImg, int filtWid, float imageStDev, float colorStDev, int iterations);
template <class Image_T> void DoFlatten(Image_T& dstImg, const Image_T& srcImg, const float shrinkFac, const float biasFac);
template <class Image_T> void DoHorizFlatten(Image_T& dstImg, const Image_T& srcImg, const float a, const float b);
template <class Image_T> void DoMatMul(Image_T& Img, float m0, float m1, float m2, float m3, float m4, float m5, float m6, float m7, float m8, float m9, float ma, float mb, float mc, float md, float me, float mf);
template <class Image_T> void DoPrintTransect(Image_T& Img, const int y);
template <class Image_T> void DoThreshold(Image_T& Img, const typename Image_T::PixType threshold);
template <class Image_T> void DoLinearMap(Image_T& Img, const f3Pixel p0, const f3Pixel p1);
template <class Image_T> void DoOptimizeLinearMap(const Image_T& Img);
