#include "FixOps.h"

#include "Math/DownSimplex.h"
#include "Math/Matrix44.h"
#include "Math/Random.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

template <class Image_T> void Noise(Image_T& Img, float mean, float stdev)
{
    for (int i = 0; i < Img.size(); i++) {
        float v = NRandf(mean, stdev);
        Img[i] = typename Image_T::PixType(v);
    }
}

// Copy pixels from srcImg into dstImg
// key is used for mode 2 and alpha is used for mode 4.
template <class Image_T>
void CopyBlock(Image_T* dstImg, const Image_T* srcImg, int dx, int dy, int sx, int sy, int sw, int sh, int mode, typename Image_T::PixType key, float alpha)
{
    ASSERT_R(dstImg);
    ASSERT_R(srcImg);

    if (sw <= 0) sw = 100000; // If sw and sh are unspecified it indicates maximum size
    if (sh <= 0) sh = 100000;

    // Clip source box to source image
    int nsw = srcImg->w_virtual();
    int nsh = srcImg->h_virtual();
    if (sx + sw > nsw) sw = nsw - sx;
    if (sy + sh > nsh) sh = nsh - sy;

    // Clip source box to dest image
    int ndw = dstImg->w_virtual();
    int ndh = dstImg->h_virtual();
    if (dx + sw > ndw) sw = ndw - dx;
    if (dy + sh > ndh) sh = ndh - dy;

    ASSERT_R(dx + sw <= ndw && dy + sh <= ndh);
    ASSERT_R(sw > 0 && sh > 0);

    // Copy the source image into the dest image
    switch (mode) {
    case 0: std::cerr << "Copying"; break;
    case 1: std::cerr << "Over blending with pixel alpha"; break;
    case 2: std::cerr << "Keyed blitting (key=" << key << ')'; break;
    case 3: std::cerr << "Adobe Multiply lerping "; break;
    case 4: std::cerr << "Over blending (alpha=" << alpha << ')'; break;
    case 5: std::cerr << "Adding "; break;
    case 6: std::cerr << "Copying channel " << key[0]; break;
    default: std::cerr << "ERROR"; break;
    }

    std::cerr << " to " << dx << ',' << dy << " from " << sx << ',' << sy << " size " << sw << 'x' << sh << " mode " << mode << '\n';
    CopyRect(*dstImg, *srcImg, sx, sy, dx, dy, sw, sh, mode, key, alpha);
}

// True if all channels are below threshold
template <class Pixel_T> inline bool isDark(const typename Pixel_T p, const typename Pixel_T::ElType threshold)
{
    for (int c = 0; c < Pixel_T::Chan; c++)
        if (p[c] > threshold) return false;

    return true;
}

// Search a ring around each pixel; if the ring is all light then paint the box inside the ring light
// Good for removing isolated dark speckles within a light background
template <class Image_T>
void RingDespeckleFilter(Image_T& dstImg, const Image_T& srcImg, const int filtWid, const int ringThick, const typename Image_T::PixType::ElType threshold)
{
    const int NF = filtWid / 2;    // NF is outer radius
    const int NC = NF - ringThick; // NC is inner radius

    dstImg.SetSize(srcImg.w(), srcImg.h());
    typename Image_T::PixType fillColor(0);
    dstImg.fill(fillColor);

    for (int y = 0; y < srcImg.h(); y++) {
        for (int x = 0; x < srcImg.w(); x++) {
            if (dstImg(x, y) != fillColor) continue; // Already been set as part of a previous speckle
            dstImg(x, y) = srcImg(x, y);
            if (!isDark(srcImg(x, y), threshold)) continue; // If not dark then no need to check the ring

            // See whether this pixel is a speckle or an object
            bool ringIsLight = true;
            typename Image_T::PixType::MathPixType ringColorAccum(0);
            typename Image_T::PixType::MathType ringPixCount = 0;

            for (int by = y - NF; by <= y + NF && ringIsLight; by++) {
                if (by < 0 || by >= srcImg.h()) continue;
                for (int bx = x - NF; bx <= x + NF && ringIsLight; bx++) {
                    if (bx < 0 || bx >= srcImg.w()) continue;

                    if (abs(y - by) >= NC || abs(x - bx) >= NC) { // If template pixel is in outer ring
                        ringIsLight = ringIsLight && !isDark(srcImg(bx, by), threshold);
                        ringColorAccum += static_cast<typename Image_T::PixType::MathPixType>(srcImg(bx, by));
                        ringPixCount++;
                    }
                }
            }

            if (ringIsLight) {
                // Paint the whole center box inside the ring
                typename Image_T::PixType ringColor = ringColorAccum / ringPixCount;
                for (int by = y - NC; by <= y + NC; by++) {
                    if (by < 0 || by >= srcImg.h()) continue;
                    for (int bx = x - NC; bx <= x + NC; bx++) {
                        if (bx < 0 || bx >= srcImg.w()) continue;
                        dstImg(bx, by) = ringColor;
                    }
                }
            }
        }
    }
}

template <class Image_T> void TransposeImage(Image_T& Out, const Image_T& In)
{
    int iw = In.w(), ih = In.h();

    Out.SetSize(ih, iw);
    for (int y = 0; y < iw; y++) {
        for (int x = 0; x < ih; x++) { Out(x, y) = In(y, x); }
    }
}

// Seam carving / content-aware resize
// This function removes one column
template <class Image_T> void ContentAwareResize(Image_T& Out, const Image_T& In, const int dw, const int dh)
{
    const int chan = Image_T::PixType::Chan;
    typedef tPixel<float, chan> FPT;

    int iw = In.w(), ih = In.h();
    Out.SetSize(iw - 1, ih);

    // Energy is sum of absolute differences of channels vs. pixel to the right of this.
    f1Image PixEnergy(iw, ih);
    for (int y = 0; y < ih; y++) {
        for (int x = 0; x < iw - 1; x++) { PixEnergy(x, y) = (Abs<float, chan>(static_cast<FPT>(In(x, y)) - static_cast<FPT>(In(x + 1, y))).sum_chan()); }
        PixEnergy(iw - 1, y) = 255.0f;
    }
    // PixEnergy.Save(std::string("Pix" + std::to_string(PixEnergy.w()) + "x" + std::to_string(PixEnergy.h()) + ".png").c_str());

    f1Image CumPixEnergy(iw, ih);
    for (int x = 0; x < iw; x++) CumPixEnergy(x, 0) = PixEnergy(x, 0);
    for (int y = 1; y < ih; y++) {
        CumPixEnergy(0, y) = PixEnergy(0, y) + std::min(CumPixEnergy(0, y - 1), CumPixEnergy(1, y - 1));
        for (int x = 1; x < iw - 1; x++) {
            CumPixEnergy(x, y) = PixEnergy(x, y) + std::min(CumPixEnergy(x - 1, y - 1), std::min(CumPixEnergy(x, y - 1), CumPixEnergy(x + 1, y - 1)));
        }
        CumPixEnergy(iw - 1, y) = PixEnergy(iw - 1, y) + std::min(CumPixEnergy(iw - 2, y - 1), CumPixEnergy(iw - 1, y - 1));
    }
    // CumPixEnergy.Save(std::string("Cum" + std::to_string(CumPixEnergy.w()) + "x" + std::to_string(CumPixEnergy.h()) + ".png").c_str());

    // Count number of minimum paths
    int cntmin = 0;
    float e = 255.0f;
    for (int x = 0; x < iw; x++) {
        if (CumPixEnergy(x, ih - 1) < e) {
            e = CumPixEnergy(x, ih - 1);
            cntmin = 1;
        } else if (CumPixEnergy(x, ih - 1) == e)
            cntmin++;
    }

    // Randomly choose a minimum path
    int chc = LRand(cntmin);
    int xi = 0;
    cntmin = 0;
    for (int x = 0; x < iw; x++) {
        if (CumPixEnergy(x, ih - 1) == e) {
            xi = x;
            if (chc == cntmin) break;
            cntmin++;
        }
    }

    // Trace the seam up from the bottom copying the pixels as it goes
    for (int y = ih - 1; y >= 0; y--) {
        if (xi > 0) memcpy(&Out(0, y), &In(0, y), sizeof(Image_T::PixType) * xi);
        if (xi < iw - 1) memcpy(&Out(xi, y), &In(xi + 1, y), sizeof(Image_T::PixType) * (iw - xi - 1));

        int xn = xi;
        if (y > 0 && xi > 0 && CumPixEnergy(xi - 1, y - 1) < CumPixEnergy(xn, y - 1)) xn = xi - 1;
        if (y > 0 && xi < iw - 1 && CumPixEnergy(xi + 1, y - 1) < CumPixEnergy(xn, y - 1)) xn = xi + 1;
        xi = xn;
    }
    // Out.Save(std::string("Out" + std::to_string(Out.w()) + "x" + std::to_string(Out.h()) + ".png").c_str());
}
}; // namespace

// Copy a block into curImg from blitImage. If sw and sh are 0 it means maximum copy size
void DoCopyBlock(std::shared_ptr<baseImage> curImg, std::shared_ptr<baseImage> blitImg, int dx, int dy, int sx, int sy, int sw, int sh, int mode, uc4Pixel key,
                 float alpha)
{
    std::cerr << "Copy: " << '\n';
    ASSERT_R(typeid(*curImg) == typeid(*blitImg));

    if (dynamic_cast<const uc1Image*>(curImg.get())) {
        CopyBlock(dynamic_cast<uc1Image*>(curImg.get()), dynamic_cast<const uc1Image*>(blitImg.get()), dx, dy, sx, sy, sw, sh, mode, key, alpha);
    } else if (dynamic_cast<const uc3Image*>(curImg.get())) {
        CopyBlock(dynamic_cast<uc3Image*>(curImg.get()), dynamic_cast<const uc3Image*>(blitImg.get()), dx, dy, sx, sy, sw, sh, mode, key, alpha);
    } else if (dynamic_cast<const uc4Image*>(curImg.get())) {
        CopyBlock(dynamic_cast<uc4Image*>(curImg.get()), dynamic_cast<const uc4Image*>(blitImg.get()), dx, dy, sx, sy, sw, sh, mode, key, alpha);
    } else if (dynamic_cast<const f1Image*>(curImg.get())) {
        CopyBlock(dynamic_cast<f1Image*>(curImg.get()), dynamic_cast<const f1Image*>(blitImg.get()), dx, dy, sx, sy, sw, sh, mode, key, alpha);
    } else if (dynamic_cast<const f3Image*>(curImg.get())) {
        CopyBlock(dynamic_cast<f3Image*>(curImg.get()), dynamic_cast<const f3Image*>(blitImg.get()), dx, dy, sx, sy, sw, sh, mode, key, alpha);
    } else if (dynamic_cast<const f4Image*>(curImg.get())) {
        CopyBlock(dynamic_cast<f4Image*>(curImg.get()), dynamic_cast<const f4Image*>(blitImg.get()), dx, dy, sx, sy, sw, sh, mode, key, alpha);
    } else {
        throw DMcError("Unsupported image type.\n");
    }
}

void DoFill(std::shared_ptr<baseImage> curImg, float r, float g, float b, float a)
{
    std::cerr << "Fill: " << r << ',' << g << ',' << b << ',' << a << '\n';

    if (f1Image* f1I = dynamic_cast<f1Image*>(curImg.get())) {
        f1I->fill(f1Pixel(r));
    } else if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
        f3I->fill(f3Pixel(r, g, b));
    } else if (f4Image* f4I = dynamic_cast<f4Image*>(curImg.get())) {
        f4I->fill(f4Pixel(r, g, b, a));
    } else if (uc1Image* uc1I = dynamic_cast<uc1Image*>(curImg.get())) {
        uc1I->fill(uc1Pixel(r));
    } else if (uc3Image* uc3I = dynamic_cast<uc3Image*>(curImg.get())) {
        uc3I->fill(uc3Pixel(r, g, b));
    } else if (uc4Image* uc4I = dynamic_cast<uc4Image*>(curImg.get())) {
        uc4I->fill(uc4Pixel(r, g, b, a));
    } else {
        throw DMcError("Unsupported image type.\n");
    }
}

void DoFlip(std::shared_ptr<baseImage> curImg, bool isVert)
{
    std::cerr << "Flip: " << (isVert ? "vertical" : "horizontal") << std::endl;

    if (isVert) {
        // Modifies image in place
        size_t linelen = curImg->w_virtual() * curImg->size_pixel_virtual();
        void* line = new unsigned char[linelen];
        for (int y = 0; y < curImg->h_virtual() / 2; y++) {
            void* pf = const_cast<void*>(curImg->pv_virtual(0, y));
            void* pt = const_cast<void*>(curImg->pv_virtual(0, curImg->h_virtual() - y - 1));
            memcpy(line, pf, linelen);
            memcpy(pf, pt, linelen);
            memcpy(pt, line, linelen);
        }
        delete[] line;
    } else {
        for (int y = 0; y < curImg->h_virtual(); y++) {
            for (int x = 0; x < curImg->w_virtual() / 2; x++) {
                void* pf = const_cast<void*>(curImg->pv_virtual(x, y));
                void* pt = const_cast<void*>(curImg->pv_virtual(curImg->w_virtual() - x - 1, y));
                for (int i = 0; i < curImg->size_pixel_virtual(); i++) {
                    unsigned char c = *(((unsigned char*)(pf)) + i);
                    *(((unsigned char*)(pf)) + i) = *(((unsigned char*)(pt)) + i);
                    *(((unsigned char*)(pt)) + i) = c;
                }
            }
        }
    }
}

void DoGradient(std::shared_ptr<baseImage> curImg, int c, bool isVert, int minx, int maxx)
{
    std::cerr << "Gradient: channel " << c << (isVert ? " vertical" : " horizontal") << " minx " << minx << " maxx " << maxx << std::endl;

    // Modifies image in place
    if (curImg->size_element_virtual() > 1 || !curImg->is_integer_virtual()) throw DMcError("Unsupported image type.\n");

    bool up = maxx > minx;
    float delta = 255.0f / float(maxx - minx);
    float valy = up ? 0.0f : 255.0f;
    if (!up) std::swap(minx, maxx);
    if (minx < 0) valy += delta * -minx;

    for (int y = 0; y < curImg->h_virtual(); y++) {
        float valx = up ? 0.0f : 255.0f;
        if (minx < 0) valx += delta * -minx;

        for (int x = 0; x < curImg->w_virtual(); x++) {
            unsigned char* pv = (unsigned char*)(curImg->pv_virtual(x, y));
            pv[c] = (unsigned char)(isVert ? valy : valx);

            if (x >= minx && x < maxx) valx += delta;
        }
        if (y >= minx && y < maxx) valy += delta;
    }
}

std::shared_ptr<baseImage> DoBlur(std::shared_ptr<baseImage> curImg, int filtWid, float imageStDev)
{
    std::cerr << "Blur: " << filtWid << " " << imageStDev << std::endl;

    if (std::shared_ptr<f1Image> f1I = dynamic_pointer_cast<f1Image>(curImg)) {
        std::shared_ptr<f1Image> dstImg(new f1Image);
        GaussianBlur(*dstImg, *f1I, filtWid, imageStDev);
        curImg = dstImg;
    } else if (std::shared_ptr<f3Image> f3I = dynamic_pointer_cast<f3Image>(curImg)) {
        std::shared_ptr<f3Image> dstImg(new f3Image);
        GaussianBlur(*dstImg, *f3I, filtWid, imageStDev);
        curImg = dstImg;
    } else if (std::shared_ptr<f4Image> f4I = dynamic_pointer_cast<f4Image>(curImg)) {
        std::shared_ptr<f4Image> dstImg(new f4Image);
        GaussianBlur(*dstImg, *f4I, filtWid, imageStDev);
        curImg = dstImg;
    } else {
        throw DMcError("Unsupported image type.\n");
    }

    return curImg;
}

std::shared_ptr<baseImage> DoGrow(std::shared_ptr<baseImage> curImg)
{
    std::cerr << "Grow black region" << std::endl;

    if (std::shared_ptr<f1Image> srcIm = dynamic_pointer_cast<f1Image>(curImg)) {
        f1Image* dstIm(new f1Image(srcIm->w(), srcIm->h()));
        for (int y = 0; y < srcIm->h(); y++)
            for (int x = 0; x < srcIm->w(); x++) {
                bool nearBlack = false;
                for (int sy = std::max(0, y - 1); sy <= y + 1 && sy < srcIm->h(); sy++)
                    for (int sx = std::max(0, x - 1); sx <= x + 1 && sx < srcIm->w(); sx++) // If I have a black neighbor set me to black
                        nearBlack = nearBlack || ((*srcIm)(sx, sy) < 0.5f);
                (*dstIm)(x, y) = nearBlack ? f1Pixel(0.f) : (*srcIm)(x, y);
            }
        return std::shared_ptr<baseImage>(dstIm);
    } else {
        throw DMcError("Unsupported image type.\n");
    }

    return curImg;
}

namespace {
template <class Elem_T, int srcChan> std::shared_ptr<baseImage> SrcChanConv(const tImage<tPixel<Elem_T, srcChan>>* src, int dstChan)
{
    std::cerr << "ChanConv " << srcChan << " to " << dstChan << " channel " << typeid(Elem_T).name() << " image\n";

    switch (dstChan) {
    case 1: return std::shared_ptr<baseImage>(new tImage<tPixel<Elem_T, 1>>(*src)); break;
    case 2: return std::shared_ptr<baseImage>(new tImage<tPixel<Elem_T, 2>>(*src)); break;
    case 3: return std::shared_ptr<baseImage>(new tImage<tPixel<Elem_T, 3>>(*src)); break;
    case 4: return std::shared_ptr<baseImage>(new tImage<tPixel<Elem_T, 4>>(*src)); break;
    default: return std::shared_ptr<baseImage>(new tImage<tPixel<Elem_T, 1>>(*src)); break;
    }
}

template <class Elem_T> std::shared_ptr<baseImage> ElTypeConv(const std::shared_ptr<baseImage> curImg, int dstChan)
{
    switch (curImg->chan_virtual()) {
    case 1: return SrcChanConv<Elem_T, 1>(dynamic_cast<tImage<tPixel<Elem_T, 1>>*>(curImg.get()), dstChan); break;
    case 2: return SrcChanConv<Elem_T, 2>(dynamic_cast<tImage<tPixel<Elem_T, 2>>*>(curImg.get()), dstChan); break;
    case 3: return SrcChanConv<Elem_T, 3>(dynamic_cast<tImage<tPixel<Elem_T, 3>>*>(curImg.get()), dstChan); break;
    case 4: return SrcChanConv<Elem_T, 4>(dynamic_cast<tImage<tPixel<Elem_T, 4>>*>(curImg.get()), dstChan); break;
    default: return SrcChanConv<Elem_T, 1>(dynamic_cast<tImage<tPixel<Elem_T, 1>>*>(curImg.get()), dstChan); break;
    }
}
} // namespace

std::shared_ptr<baseImage> DoChanConvert(std::shared_ptr<baseImage> curImg, int dstChan)
{
    baseImage* cI = curImg.get();
    if (cI->size_element_virtual() == 4 && !cI->is_integer_virtual()) return ElTypeConv<float>(curImg, dstChan);
    if (cI->size_element_virtual() == 1 && cI->is_integer_virtual()) return ElTypeConv<unsigned char>(curImg, dstChan);

    throw DMcError("Unsupported image element type.\n");

    return curImg;
}

std::shared_ptr<baseImage> DoRingDespeckleFilter(std::shared_ptr<baseImage> curImg, int filtWid, int ringThick, float threshold, int iterations)
{
    std::cerr << "RingDespeckleFilter: " << filtWid << " " << ringThick << " " << threshold << " " << iterations << std::endl;

    for (int l = 0; l < iterations; l++) {
        int curFiltWid = (l == 0 ? filtWid : (rand() % filtWid)) | 1; // First iter is precise; rest are random

        if (uc1Image* uc1I = dynamic_cast<uc1Image*>(curImg.get())) {
            uc1Image* dstImg = new uc1Image;
            RingDespeckleFilter(*dstImg, *uc1I, curFiltWid, ringThick, static_cast<unsigned char>(threshold));
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else if (uc3Image* uc3I = dynamic_cast<uc3Image*>(curImg.get())) {
            uc3Image* dstImg = new uc3Image;
            RingDespeckleFilter(*dstImg, *uc3I, curFiltWid, ringThick, static_cast<unsigned char>(threshold));
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else {
            throw DMcError("Unsupported image type.\n");
        }
    }
    return curImg;
}

std::shared_ptr<baseImage> DoDespeckle(std::shared_ptr<baseImage> curImg)
{
    std::cerr << "Despeckle\n";

    if (f1Image* f1I = dynamic_cast<f1Image*>(curImg.get())) {
        f1Image* dstImg = new f1Image;
        Despeckle(*dstImg, *f1I);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
        f3Image* dstImg = new f3Image;
        Despeckle(*dstImg, *f3I);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (f4Image* f4I = dynamic_cast<f4Image*>(curImg.get())) {
        f4Image* dstImg = new f4Image;
        Despeckle(*dstImg, *f4I);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (uc1Image* uc1I = dynamic_cast<uc1Image*>(curImg.get())) {
        uc1Image* dstImg = new uc1Image;
        Despeckle(*dstImg, *uc1I);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (uc3Image* uc3I = dynamic_cast<uc3Image*>(curImg.get())) {
        uc3Image* dstImg = new uc3Image;
        Despeckle(*dstImg, *uc3I);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (uc4Image* uc4I = dynamic_cast<uc4Image*>(curImg.get())) {
        uc4Image* dstImg = new uc4Image;
        Despeckle(*dstImg, *uc4I);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else {
        throw DMcError("Unsupported image type.\n");
    }

    return curImg;
}

template <class Image_T>
void DoMatMul(Image_T& Img, float m0, float m1, float m2, float m3, float m4, float m5, float m6, float m7, float m8, float m9, float ma, float mb, float mc,
              float md, float me, float mf)
{
    std::cerr << "Matrix multiply on pixel values\n";

    float im[16] = {m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, ma, mb, mc, md, me, mf};
    Matrix44<f3vec> Mat(im);

    // Modifies image in place
    for (int p = 0; p < Img.size(); p++) {
        f4Pixel pf = Img[p]; // Converts from Image_T's pixel type to f4Pixel
        Mat.Project(pf.r(), pf.g(), pf.b(), pf.a());
        Img[p] = pf; // Converts from f4Pixel to Image_T's pixel type
    }
}
template void DoMatMul(f3Image& Img, float m0, float m1, float m2, float m3, float m4, float m5, float m6, float m7, float m8, float m9, float ma, float mb,
                       float mc, float md, float me, float mf);
template void DoMatMul(f4Image& Img, float m0, float m1, float m2, float m3, float m4, float m5, float m6, float m7, float m8, float m9, float ma, float mb,
                       float mc, float md, float me, float mf);

std::shared_ptr<baseImage> DoNoise(std::shared_ptr<baseImage> curImg, float mean, float stdev)
{
    std::cerr << "Noise: " << mean << " " << stdev << '\n';

    if (uc1Image* uc1I = dynamic_cast<uc1Image*>(curImg.get())) {
        Noise(*uc1I, mean, stdev);
    } else if (uc3Image* uc3I = dynamic_cast<uc3Image*>(curImg.get())) {
        Noise(*uc3I, mean, stdev);
    } else if (uc4Image* uc4I = dynamic_cast<uc4Image*>(curImg.get())) {
        Noise(*uc4I, mean, stdev);
    } else if (f1Image* f1I = dynamic_cast<f1Image*>(curImg.get())) {
        Noise(*f1I, mean, stdev);
    } else if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
        Noise(*f3I, mean, stdev);
    } else if (f4Image* f4I = dynamic_cast<f4Image*>(curImg.get())) {
        Noise(*f4I, mean, stdev);
    } else {
        throw DMcError("Unsupported image type.\n");
    }

    return curImg;
}

std::shared_ptr<baseImage> DoResize(std::shared_ptr<baseImage> curImg, int newWid, int newHgt)
{
    std::cerr << "Resize: " << newWid << 'x' << newHgt << std::endl;

    if (f1Image* f1I = dynamic_cast<f1Image*>(curImg.get())) {
        f1Image* dstImg = new f1Image();
        Resample(*dstImg, *f1I, newWid, newHgt);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
        f3Image* dstImg = new f3Image();
        Resample(*dstImg, *f3I, newWid, newHgt);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (f4Image* f4I = dynamic_cast<f4Image*>(curImg.get())) {
        f4Image* dstImg = new f4Image();
        Resample(*dstImg, *f4I, newWid, newHgt);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (uc1Image* uc1I = dynamic_cast<uc1Image*>(curImg.get())) {
        uc1Image* dstImg = new uc1Image();
        Resample(*dstImg, *uc1I, newWid, newHgt);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (uc3Image* uc3I = dynamic_cast<uc3Image*>(curImg.get())) {
        uc3Image* dstImg = new uc3Image();
        Resample(*dstImg, *uc3I, newWid, newHgt);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (uc4Image* uc4I = dynamic_cast<uc4Image*>(curImg.get())) {
        uc4Image* dstImg = new uc4Image();
        Resample(*dstImg, *uc4I, newWid, newHgt);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else {
        throw DMcError("Unsupported image type.\n");
    }

    return curImg;
}

// Shrink by one pixel in one dimension
std::shared_ptr<baseImage> CAResizeBranch(std::shared_ptr<baseImage> curImg, int newWid, int newHgt)
{
    ASSERT_R(newWid == curImg->w_virtual() || newHgt == curImg->h_virtual());

    if (f1Image* f1I = dynamic_cast<f1Image*>(curImg.get())) {
        f1Image* dstImg = new f1Image();
        ContentAwareResize(*dstImg, *f1I, newWid, newHgt);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
        f3Image* dstImg = new f3Image();
        ContentAwareResize(*dstImg, *f3I, newWid, newHgt);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (f4Image* f4I = dynamic_cast<f4Image*>(curImg.get())) {
        f4Image* dstImg = new f4Image();
        ContentAwareResize(*dstImg, *f4I, newWid, newHgt);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (uc1Image* uc1I = dynamic_cast<uc1Image*>(curImg.get())) {
        uc1Image* dstImg = new uc1Image();
        ContentAwareResize(*dstImg, *uc1I, newWid, newHgt);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (uc3Image* uc3I = dynamic_cast<uc3Image*>(curImg.get())) {
        uc3Image* dstImg = new uc3Image();
        ContentAwareResize(*dstImg, *uc3I, newWid, newHgt);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (uc4Image* uc4I = dynamic_cast<uc4Image*>(curImg.get())) {
        uc4Image* dstImg = new uc4Image();
        ContentAwareResize(*dstImg, *uc4I, newWid, newHgt);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else {
        throw DMcError("Unsupported image type.\n");
    }

    return curImg;
}

std::shared_ptr<baseImage> DoResizeCA(std::shared_ptr<baseImage> curImg, int newWid, int newHgt)
{
    std::cerr << "Content-Aware Resize: " << newWid << 'x' << newHgt << std::endl;
    const int ITERS_PER_DIM = 16;

    while (curImg->w_virtual() > newWid || curImg->h_virtual() > newHgt) {
        for (int i = 0; i < ITERS_PER_DIM && curImg->w_virtual() > newWid; i++) { curImg = CAResizeBranch(curImg, newWid, curImg->h_virtual()); }

        // Transpose it to do vertical shrinks
        if (f1Image* f1I = dynamic_cast<f1Image*>(curImg.get())) {
            f1Image* dstImg = new f1Image();
            TransposeImage(*dstImg, *f1I);
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
            f3Image* dstImg = new f3Image();
            TransposeImage(*dstImg, *f3I);
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else if (f4Image* f4I = dynamic_cast<f4Image*>(curImg.get())) {
            f4Image* dstImg = new f4Image();
            TransposeImage(*dstImg, *f4I);
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else if (uc1Image* uc1I = dynamic_cast<uc1Image*>(curImg.get())) {
            uc1Image* dstImg = new uc1Image();
            TransposeImage(*dstImg, *uc1I);
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else if (uc3Image* uc3I = dynamic_cast<uc3Image*>(curImg.get())) {
            uc3Image* dstImg = new uc3Image();
            TransposeImage(*dstImg, *uc3I);
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else if (uc4Image* uc4I = dynamic_cast<uc4Image*>(curImg.get())) {
            uc4Image* dstImg = new uc4Image();
            TransposeImage(*dstImg, *uc4I);
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else {
            throw DMcError("Unsupported image type.\n");
        }

        // It's transposed so tread carefully.
        for (int i = 0; i < ITERS_PER_DIM && curImg->w_virtual() > newHgt; i++) { curImg = CAResizeBranch(curImg, newHgt, curImg->h_virtual()); }

        // Transpose it back
        if (f1Image* f1I = dynamic_cast<f1Image*>(curImg.get())) {
            f1Image* dstImg = new f1Image();
            TransposeImage(*dstImg, *f1I);
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
            f3Image* dstImg = new f3Image();
            TransposeImage(*dstImg, *f3I);
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else if (f4Image* f4I = dynamic_cast<f4Image*>(curImg.get())) {
            f4Image* dstImg = new f4Image();
            TransposeImage(*dstImg, *f4I);
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else if (uc1Image* uc1I = dynamic_cast<uc1Image*>(curImg.get())) {
            uc1Image* dstImg = new uc1Image();
            TransposeImage(*dstImg, *uc1I);
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else if (uc3Image* uc3I = dynamic_cast<uc3Image*>(curImg.get())) {
            uc3Image* dstImg = new uc3Image();
            TransposeImage(*dstImg, *uc3I);
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else if (uc4Image* uc4I = dynamic_cast<uc4Image*>(curImg.get())) {
            uc4Image* dstImg = new uc4Image();
            TransposeImage(*dstImg, *uc4I);
            curImg = std::shared_ptr<baseImage>(dstImg);
        } else {
            throw DMcError("Unsupported image type.\n");
        }
    }

    return curImg;
}

std::shared_ptr<baseImage> DoVCD(std::shared_ptr<baseImage> curImg, int filtWid, float imageStDev, float colorStDev, int iterations)
{
    std::cerr << "VCD: " << filtWid << " " << imageStDev << " " << colorStDev << " " << iterations << std::endl;

    if (f1Image* f1I = dynamic_cast<f1Image*>(curImg.get())) {
        f1Image* dstImg = new f1Image;
        VCD(*dstImg, *f1I, filtWid, imageStDev, colorStDev, iterations);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
        f3Image* dstImg = new f3Image;
        VCD(*dstImg, *f3I, filtWid, imageStDev, colorStDev, iterations);
        curImg = std::shared_ptr<baseImage>(dstImg);
    } else {
        throw DMcError("Unsupported image type.\n");
    }
    return curImg;
}

template <class Image_T> void DoFlattenHoriz(Image_T& dstImg, const Image_T& srcImg, const float a, const float b)
{
    std::cerr << "HorizFlatten: " << a << " " << b << '\n';

    float* avg = new float[srcImg.h()];
    float* savg = new float[srcImg.h()];

    for (int y = 0; y < srcImg.h(); y++) {
        float s = 0;
        for (int x = 0; x < srcImg.w(); x++) {
            if (Image_T::PixType::Chan == 1)
                s += srcImg(x, y)[0];
            else
                s += srcImg(x, y).luminance();
        }
        avg[y] = s / float(srcImg.w());
    }

    int kd2 = a;

    for (int y = kd2; y < srcImg.h() - kd2; y++) {
        float s = 0;
        for (int x = -kd2; x <= kd2; x++) { s += avg[y + x]; }
        savg[y] = s / float(kd2 + kd2 + 1);
    }

    dstImg.SetSize(srcImg.w(), srcImg.h());

    for (int y = 0; y < srcImg.h(); y++) {
        float scale = savg[y] / avg[y];
        // if(y % 100 == 0)
        // scale = 0.0f;
        // std::cerr << avg[y] << " " << savg[y] << " " << scale << '\n';
        for (int x = 0; x < srcImg.w(); x++) {
            if (y % 100 == 0)
                dstImg(x, y) = Image_T::PixType(0);
            else
                dstImg(x, y) = srcImg(x, y) * scale;
        }
    }

    delete[] avg;
}
template void DoFlattenHoriz(f1Image& dstImg, const f1Image& srcImg, const float a, const float b);
template void DoFlattenHoriz(f3Image& dstImg, const f3Image& srcImg, const float a, const float b);

template <class Image_T> void DoFlatten(Image_T& dstImg, const Image_T& srcImg, const float shrinkFac, const float biasFac)
{
    std::cerr << "Flatten: shrink " << shrinkFac << " bias " << biasFac << '\n';

    // Get the blurry image by downsampling to 512x512-ish, then blurring, then bicubic upsampling.
    std::cerr << "Shrinking.\n";
    Image_T* shrunkImg = new Image_T(srcImg);

    while (shrunkImg->w() >= 512 || shrunkImg->h() >= 512) {
        // Downsample by two and replace shrunkImg.
        Image_T* smallerImg = new Image_T;
        Downsample2x2(*smallerImg, *shrunkImg);
        delete shrunkImg;
        shrunkImg = smallerImg;

        std::cerr << "Shrunk to " << shrunkImg->w() << "x" << shrunkImg->h() << std::endl;
    }

    int KernelSize = int(sqrtf(srcImg.w() * srcImg.h()) * shrinkFac);
    KernelSize = (KernelSize * shrunkImg->w()) / srcImg.w();
    if (!(KernelSize & 1)) KernelSize++;
    if (KernelSize < 3) KernelSize = 3;

    Image_T blurredImg, shrunkBlurredImg;

    std::cerr << "Blurring. shrinkFac = " << shrinkFac << " KernelSize = " << KernelSize << "\n";
    GaussianBlur(shrunkBlurredImg, *shrunkImg, KernelSize, KernelSize * 0.333f);
    std::cerr << "Upsampling.\n";
    // BlurImg will be full-size blurred.
    Resample(blurredImg, shrunkBlurredImg, srcImg.w(), srcImg.h());

    // shrunkBlurredImg.Save("shrunkBlurredImg.jpg");
    // blurredImg.Save("blurredImg.png");
    // shrunkImg->Save("fOpImg.jpg");

    // Do image average in double precision, since float fails about 16MP.
    tPixel<double, Image_T::PixType::Chan> csum(0);
    for (int i = 0; i < srcImg.size(); i++) csum += tPixel<double, Image_T::PixType::Chan>(srcImg[i]);

    typename Image_T::PixType Avg = csum / double(srcImg.size());

    std::cerr << "Applying flatten. Avg = " << Avg << "\n";

    dstImg.SetSize(blurredImg.w(), blurredImg.h());

    for (int i = 0; i < dstImg.size(); i++) {
        typename Image_T::PixType Biased(srcImg[i] + (Avg - blurredImg[i])); // Flatten by biasing
        typename Image_T::PixType Scaled((srcImg[i] / blurredImg[i]) * Avg); // Flatten by scaling
        dstImg[i] = Biased * biasFac + Scaled * (1.0f - biasFac);
    }

    delete shrunkImg;
}
template void DoFlatten(f1Image& dstImg, const f1Image& srcImg, const float shrinkFac, const float biasFac);
template void DoFlatten(f3Image& dstImg, const f3Image& srcImg, const float shrinkFac, const float biasFac);

// -transect 128 > fname.csv
template <class Image_T> void DoPrintTransect(Image_T& Img, const int y)
{
    std::cerr << "Transect: y " << y << '\n';

    for (int x = 0; x < Img.w(); x++) {
        std::cout << x << ',' << static_cast<float>(Img(x, y).r()) << ',' << static_cast<float>(Img(x, y).g()) << ',' << static_cast<float>(Img(x, y).b())
                  << std::endl;
    }
}
template void DoPrintTransect(uc1Image& Img, const int y);
template void DoPrintTransect(f1Image& Img, const int y);
template void DoPrintTransect(uc3Image& Img, const int y);
template void DoPrintTransect(f3Image& Img, const int y);

template <class Image_T> void DoThreshold(Image_T& Img, const typename Image_T::PixType threshold)
{
    std::cerr << "Threshold: " << threshold << '\n';

    for (int i = 0; i < Img.size(); i++) {
        for (int c = 0; c < Image_T::PixType::Chan; c++) {
            Img[i][c] = Img[i][c] >= threshold[c] ? element_traits<typename Image_T::PixType::ElType>::one() : 0;
        }
    }
}
template void DoThreshold(uc1Image& Img, const uc1Image::PixType threshold);
template void DoThreshold(uc3Image& Img, const uc3Image::PixType threshold);

template <class Image_T> float LinearMapError(const HVector<float>& OVec, void* imgp)
{
    for (const float f : OVec)
        if (f < 0.f || f > 1.0f) return std::numeric_limits<float>::max();

    Image_T& Img = *reinterpret_cast<Image_T*>(imgp);

    f3vec v0(OVec[0], OVec[1], OVec[2]), v1(OVec[3], OVec[4], OVec[5]);
    float err = 0;

    f3vec dir = v1 - v0;
    dir.normalize();

    for (int i = 0; i < Img.size(); i++) {
        f3vec v(Img[i].r(), Img[i].g(), Img[i].b());
        f3vec vv0 = v - v0;
        f3vec vp = dir * dot(dir, vv0) + v0;
        err += (vp - v).lenSqr();
    }

    // std::cerr << OVec << " Error = " << err << '\n';
    return err;
}
template float LinearMapError<f3Image>(const HVector<float>& OVec, void* imgp);
template float LinearMapError<f4Image>(const HVector<float>& OVec, void* imgp);

template <class Image_T> void DoOptimizeLinearMap(const Image_T& Img)
{
    std::vector<HVector<float>> OVecs(7, HVector<float>(6));
    for (auto& i : OVecs) i.rand();

    int nfunk = 0;
    float err = DownSimplex(&OVecs[0], 6, 0.01f, LinearMapError<Image_T>, (void*)&Img, nfunk);
    std::cerr << OVecs[0] << "Error = " << err << '\n';
}
template void DoOptimizeLinearMap(const f3Image& Img);
template void DoOptimizeLinearMap(const f4Image& Img);

// Compute the two endpoints in RGB color space of a gradient that's the best fit color map for the pixels in Img
template <class Image_T> void DoLinearMap(Image_T& Img, const f3Pixel p0, const f3Pixel p1)
{
    f3vec v0(p0.r(), p0.g(), p0.b()), v1(p1.r(), p1.g(), p1.b());
    f3vec dir = v1 - v0;
    dir.normalize();

    for (int i = 0; i < Img.size(); i++) {
        f3vec v(Img[i].r(), Img[i].g(), Img[i].b());
        f3vec vv0 = v - v0;
        f3vec vp = dir * dot(dir, vv0) + v0;

        Img[i][0] = vp[0];
        Img[i][1] = vp[1];
        Img[i][2] = vp[2];
    }
}
template void DoLinearMap(f3Image& Img, const f3Pixel p0, const f3Pixel p1);
template void DoLinearMap(f4Image& Img, const f3Pixel p0, const f3Pixel p1);

template <class Image_T> void DoDiff(Image_T& imgA, const Image_T& imgB, const float scale)
{
    typedef typename Image_T::PixType::ElType ChType;
    typedef element_traits<ChType>::MathType MType;
    MType maxChan = 0;
    typename Image_T::PixType difPix;

    for (int y = 0; y < imgA.h(); y++) {
        for (int x = 0; x < imgA.w(); x++) {
            bool printit = false;
            for (int c = 0; c < imgA.chan(); c++) {
                MType cA = imgA(x, y)[c];
                MType cB = imgB(x, y)[c];
                MType dif = std::abs(cA - cB);
                difPix[c] = dif;
                if (dif > maxChan) { // Print each pixel if its diff is the largest so far
                    maxChan = dif;
                    printit = true;
                }
            }

            if (printit) {
                std::cerr << "(" << x << "," << y << ") " << difPix << " = " << imgA(x, y) << " - " << imgB(x, y) << '\n';
                printit = false;
            }

            imgA(x, y) = static_cast<typename Image_T::PixType::FloatMathPixType>(difPix) * scale;
        }
    }
}
template void DoDiff(uc1Image& imgA, const uc1Image& imgB, const float scale);
template void DoDiff(uc3Image& imgA, const uc3Image& imgB, const float scale);
template void DoDiff(uc4Image& imgA, const uc4Image& imgB, const float scale);
template void DoDiff(f1Image& imgA, const f1Image& imgB, const float scale);
template void DoDiff(f3Image& imgA, const f3Image& imgB, const float scale);
template void DoDiff(f4Image& imgA, const f4Image& imgB, const float scale);
