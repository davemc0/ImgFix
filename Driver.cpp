// This program loads images, applies algorithms to them, and saves them.
// Use -ld and -st to store images in numbered registers
// Copyright David K. McAllister, 2009-2022.

// ImgFix -in in.jpg -despecklering 13 1 190 1 -threshold 189 -blur 3 1 -out out.png
// For making alpha cut-out of signature written with black ink on white paper
// ImgFix -in orig.png -threshold 95 -chan 1 -float -grow -st 1 -in orig.png -float -chan 1 -comp 1 0 0 0 0 0 0 5 -chan 4 -uchar -float -st 3 -fill 0 0 0 1
// -comp 3 0 0 0 0 0 0 6 3 -matmul 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 1 -out sig.jpg

// TODO:
// Add rotation

#include "FixOps.h"
#include "Helpers.h"
#include "Image/ImageAlgorithms.h"
#include "Math/Random.h"

#include <fstream>
#include <iostream>
#include <string>

static void Usage(const char* message = NULL, const bool Exit = true)
{
    if (message) std::cerr << "\nERROR: " << message << std::endl;

    std::cerr << // One big string
        "Usage: <options>\n\n"
        "-in fname.jpg               Load the current image\n"
        "-out fname.jpg              Save the current image\n"
        "-ld I                       Load image from slot I to current\n"
        "-st I                       Save a copy of curImage to slot I\n"
        "-chan N                     Convert image to N channel\n"
        "-float                      Convert image to float\n"
        "-uchar                      Convert image to uchar\n"
        "-sp mcp mif mi te map       Set LoadSaveParams: maxColorPalette, maxItersFast, maxIters, targetErr, makeArtisticPalette\n\n"

        "-blur f s                   Blur filt_width img_stdev; try 5 2.0\n"
        "-comp I dx dy sx sy w h M   Composite a subrect of the given image to an arbitrary place in the cur image\n"
        "                                M:0=copy 1=over blend 2=key (int) r g b a 3=adobe multiply 4=over blend const alpha 5=add 6=copy channel c\n"
        "-despeckle                  Clamp each pixel to be within the bounding box of its eight neighbors\n"
        "-despecklering f r t i      Ring Despeckle filt_width ring_thickness threshold iters; random smaller rings if iter>1; try 15 2 180 1\n"
        "-diff I scale               Cur = (I - Cur) * scale\n"
        "-fill r g b a               Fill image with specified color (0..255 or 0..1.0)\n"
        "-findlinmap                 Find best fit line to project f3Image pixel values\n"
        "-flatten shrink bias        Flatten; shrink is blurred image rel. size; bias blends btw. bias and scale methods; try 0.1 0.5\n"
        "-flattenh a b               Flatten horizontal rows to remove scanner artifacts; ignores params\n"
        "-flip h|v                   Flip image horizontally or vertically\n"
        "-gradient C h|v minx maxx   Set channel C to a horiz or vert gradient with 0 and 255 at given scanlines\n"
        "-grow                       Set pixel to black if any neighboring pixel is black, otherwise no change\n"
        "-linmap r0 g0 b0 r1 g1 b1   Project f3Image pixel values onto a line in RGB space\n"
        "-makef  N w h               Make an fImage of N channels, width and height; fill with black\n"
        "-makeuc N w h               Make a ucImage of N channels, width and height; fill with black\n"
        "-matmul 0123456789ab[cdef]  Multiply each pixel value by this 4x3 (for f3Image) or 4x4 (for f4Image) matrix\n"
        "-noise mean stdev           Fill image with noise\n"
        "-quant                      Quantize the image to pal colors (set params with -sp)\n"
        "-resize w|'x' h|'y'|'X'     Resize the given image to w x h. Bicubic upsampling; box filter then bicubic downsampling\n"
        "                                'X' => scale factor; 'x'|'y' => compute this dim\n"
        "-resizeca w|'x' h|'y'|'X'   Content-Aware Resize the given image to w x h. Seam carving; currently just shrinking\n"
        "                                'X' => scale factor; 'x'|'y' => compute this dim\n"
        "-threshold N                Convert image to two-tone at given threshold value\n"
        "-transect y                 Print channel 0 values of scanline y to stdout in CSV format\n"
        "-vcd f s c i                VCD filt_width img_stdev col_stdev iterations; try 5 3.0 0.25 1\n";

    if (Exit) exit(1);
}

static void Driver(int& argc, char** argv)
{
    ImgSlots Slots;
    std::shared_ptr<baseImage> curImg = NULL;
    LoadSaveParams SP;

    if (argc < 2) Usage();

    for (int i = 1; i < argc; i++) {
        std::string starg(argv[i]);

        if (curImg) std::cerr << "curImg: " << (*curImg) << '\n';

        if (starg == "-h" || starg == "-help") {
            Usage();
        } else if (starg == "-in") {
            if (argc <= i + 1) Usage();

            std::string srcFName(argv[i + 1]);
            curImg = std::shared_ptr<baseImage>(LoadtImage(srcFName.c_str()));
            std::cerr << "Load: " << srcFName << " size=" << curImg->w_virtual() << "x" << curImg->h_virtual() << " chan=" << curImg->chan_virtual() << '\n';

            RemoveArgs(argc, argv, i, 2);
        } else if (starg == "-makef") {
            if (argc <= i + 3) Usage();

            int chan = atoi(argv[i + 1]), w = atoi(argv[i + 2]), h = atoi(argv[i + 3]);

            std::cerr << "Make: f" << chan << "Image\n";

            std::shared_ptr<baseImage> newImg = NULL;
            if (chan == 1) newImg = std::shared_ptr<baseImage>(new f1Image(w, h, f1Pixel(0)));
            if (chan == 3) newImg = std::shared_ptr<baseImage>(new f3Image(w, h, f3Pixel(0)));
            if (chan == 4) newImg = std::shared_ptr<baseImage>(new f4Image(w, h, f4Pixel(0)));
            curImg = newImg;

            RemoveArgs(argc, argv, i, 4);
        } else if (starg == "-makeuc") {
            if (argc <= i + 3) Usage();

            int chan = atoi(argv[i + 1]), w = atoi(argv[i + 2]), h = atoi(argv[i + 3]);

            std::cerr << "Make: uc" << chan << "Image\n";

            std::shared_ptr<baseImage> newImg = NULL;
            if (chan == 1) newImg = std::shared_ptr<baseImage>(new uc1Image(w, h, uc1Pixel(0)));
            if (chan == 3) newImg = std::shared_ptr<baseImage>(new uc3Image(w, h, uc3Pixel(0)));
            if (chan == 4) newImg = std::shared_ptr<baseImage>(new uc4Image(w, h, uc4Pixel(0)));
            curImg = newImg;

            RemoveArgs(argc, argv, i, 4);
        } else if (starg == "-noise") {
            if (argc <= i + 2) Usage();

            float mean = atof(argv[i + 1]), stdev = atof(argv[i + 2]);
            curImg = DoNoise(curImg, mean, stdev);

            RemoveArgs(argc, argv, i, 3);
        } else if (starg == "-fill") {
            if (argc <= i + 4) Usage();
            float r = atof(argv[i + 1]);
            float g = atof(argv[i + 2]);
            float b = atof(argv[i + 3]);
            float a = atof(argv[i + 4]);
            DoFill(curImg, r, g, b, a);

            RemoveArgs(argc, argv, i, 5);
        } else if (starg == "-comp") {
            if (argc <= i + 8) Usage();
            int dx = atoi(argv[i + 2]), dy = atoi(argv[i + 3]), sx = atoi(argv[i + 4]), sy = atoi(argv[i + 5]);
            int sw = atoi(argv[i + 6]), sh = atoi(argv[i + 7]), mode = atoi(argv[i + 8]);
            uc4Pixel key(element_traits<uc4Pixel::ElType>::one());
            float alpha = 1.0f;

            int subArgc = 0;
            if (mode == 2) {
                subArgc = curImg->chan_virtual();
                if (argc <= i + 8 + subArgc) Usage(("Need " + std::to_string(subArgc) + " ints for key color channels.").c_str());
                for (int c = 0; c < subArgc; c++) { key[c] = atoi(argv[i + 9 + c]); }
            } else if (mode == 4) {
                subArgc = 1;
                if (argc <= i + 8 + subArgc) Usage("Need a float for alpha");
                alpha = atof(argv[i + 8 + subArgc]);
            } else if (mode == 6) {
                subArgc = 1;
                if (argc <= i + 8 + subArgc) Usage("Need an int for channel number");
                key[0] = atoi(argv[i + 8 + subArgc]); // Store channel number to copy in key[0]
            }

            DoCopyBlock(curImg, Slots.get(atoi(argv[i + 1])), dx, dy, sx, sy, sw, sh, mode, key, alpha);

            RemoveArgs(argc, argv, i, 9 + subArgc);
        } else if (starg == "-diff") {
            if (argc <= i + 2) Usage();
            float scale = atof(argv[i + 2]);

            std::shared_ptr<baseImage> imgB = Slots.get(atoi(argv[i + 1]));
            ASSERT_R(typeid(*curImg) == typeid(*imgB));
            ASSERT_R(curImg->w_virtual() == imgB->w_virtual() && curImg->h_virtual() == imgB->h_virtual() && curImg->chan_virtual() == imgB->chan_virtual() &&
                     curImg->size_element_virtual() == imgB->size_element_virtual());
            // Modifies curImg in place

            if (uc1Image* imA = dynamic_cast<uc1Image*>(curImg.get())) {
                uc1Image* imB = dynamic_cast<uc1Image*>(imgB.get());
                DoDiff(*imA, *imB, scale);
            } else if (uc3Image* imA = dynamic_cast<uc3Image*>(curImg.get())) {
                uc3Image* imB = dynamic_cast<uc3Image*>(imgB.get());
                DoDiff(*imA, *imB, scale);
            } else if (uc4Image* imA = dynamic_cast<uc4Image*>(curImg.get())) {
                uc4Image* imB = dynamic_cast<uc4Image*>(imgB.get());
                DoDiff(*imA, *imB, scale);
            } else if (f1Image* imA = dynamic_cast<f1Image*>(curImg.get())) {
                f1Image* imB = dynamic_cast<f1Image*>(imgB.get());
                DoDiff(*imA, *imB, scale);
            } else if (f3Image* imA = dynamic_cast<f3Image*>(curImg.get())) {
                f3Image* imB = dynamic_cast<f3Image*>(imgB.get());
                DoDiff(*imA, *imB, scale);
            } else if (f4Image* imA = dynamic_cast<f4Image*>(curImg.get())) {
                f4Image* imB = dynamic_cast<f4Image*>(imgB.get());
                DoDiff(*imA, *imB, scale);
            } else
                throw DMcError("Unsupported image type.\n");

            RemoveArgs(argc, argv, i, 3);
        } else if (starg == "-matmul") {
            int subArgc = 0;
            if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
                subArgc = 12;
                if (argc <= i + subArgc) Usage();
                DoMatMul(*f3I, atof(argv[i + 1]), atof(argv[i + 2]), atof(argv[i + 3]), atof(argv[i + 4]), atof(argv[i + 5]), atof(argv[i + 6]),
                         atof(argv[i + 7]), atof(argv[i + 8]), atof(argv[i + 9]), atof(argv[i + 10]), atof(argv[i + 11]), atof(argv[i + 12]), 0, 0, 0, 1.0f);
            } else if (f4Image* f4I = dynamic_cast<f4Image*>(curImg.get())) {
                subArgc = 16;
                if (argc <= i + subArgc) Usage();
                DoMatMul(*f4I, atof(argv[i + 1]), atof(argv[i + 2]), atof(argv[i + 3]), atof(argv[i + 4]), atof(argv[i + 5]), atof(argv[i + 6]),
                         atof(argv[i + 7]), atof(argv[i + 8]), atof(argv[i + 9]), atof(argv[i + 10]), atof(argv[i + 11]), atof(argv[i + 12]),
                         atof(argv[i + 13]), atof(argv[i + 14]), atof(argv[i + 15]), atof(argv[i + 16]));
            } else {
                throw DMcError("Unsupported image type.\n");
            }

            RemoveArgs(argc, argv, i, subArgc + 1);
        } else if (starg == "-gradient") {
            // Modifies image in place
            if (argc <= i + 4) Usage();
            int c = atoi(argv[i + 1]), minx = atoi(argv[i + 3]), maxx = atoi(argv[i + 4]);
            bool isVert = argv[i + 2][0] == 'v';
            if (minx == maxx) Usage("Gradient must have different min and max scanline numbers");
            if (c < 0 || c >= curImg->chan_virtual()) Usage("Nonexistent color channel for gradient");
            DoGradient(curImg, c, isVert, minx, maxx);

            RemoveArgs(argc, argv, i, 5);
        } else if (starg == "-grow") {
            curImg = DoGrow(curImg);

            RemoveArgs(argc, argv, i);
        } else if (starg == "-flip") {
            if (argc <= i + 1) Usage();
            bool isVert = (argv[i + 1][0] == 'v' || argv[i + 1][0] == 'V');
            DoFlip(curImg, isVert);

            RemoveArgs(argc, argv, i, 2);
        } else if (starg == "-transect") {
            if (argc <= i + 1) Usage();
            int yval = atoi(argv[i + 1]);

            if (uc1Image* uc1I = dynamic_cast<uc1Image*>(curImg.get())) {
                DoPrintTransect(*uc1I, yval);
            } else if (uc3Image* uc3I = dynamic_cast<uc3Image*>(curImg.get())) {
                DoPrintTransect(*uc3I, yval);
            } else if (f1Image* f1I = dynamic_cast<f1Image*>(curImg.get())) {
                DoPrintTransect(*f1I, yval);
            } else if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
                DoPrintTransect(*f3I, yval);
            } else {
                throw DMcError("Unsupported image type.\n");
            }

            RemoveArgs(argc, argv, i, 2);
        } else if (starg == "-flatten") {
            if (argc <= i + 2) Usage();
            float shrinkFac = atof(argv[i + 1]);
            float biasFac = atof(argv[i + 2]);

            std::shared_ptr<baseImage> dstImg = NULL;
            if (f1Image* f1I = dynamic_cast<f1Image*>(curImg.get())) {
                f1Image* dstImg = new f1Image;
                DoFlatten(*dstImg, *f1I, shrinkFac, biasFac);
                curImg = std::shared_ptr<baseImage>(dstImg);
            } else if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
                f3Image* dstImg = new f3Image;
                DoFlatten(*dstImg, *f3I, shrinkFac, biasFac);
                curImg = std::shared_ptr<baseImage>(dstImg);
            } else {
                throw DMcError("Unsupported image type.\n");
            }

            RemoveArgs(argc, argv, i, 3);
        } else if (starg == "-flattenh") {
            if (argc <= i + 2) Usage();
            float a = atof(argv[i + 1]);
            float b = atof(argv[i + 2]);

            std::shared_ptr<baseImage> dstImg = NULL;
            if (f1Image* f1I = dynamic_cast<f1Image*>(curImg.get())) {
                f1Image* dstImg = new f1Image;
                DoFlattenHoriz(*dstImg, *f1I, a, b);
                curImg = std::shared_ptr<baseImage>(dstImg);
            } else if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
                f3Image* dstImg = new f3Image;
                DoFlattenHoriz(*dstImg, *f3I, a, b);
                curImg = std::shared_ptr<baseImage>(dstImg);
            } else {
                throw DMcError("Unsupported image type.\n");
            }

            RemoveArgs(argc, argv, i, 3);
        } else if (starg == "-chan") {
            if (argc <= i + 1) Usage();
            curImg = DoChanConvert(curImg, atoi(argv[i + 1]));

            RemoveArgs(argc, argv, i, 2);
        } else if (starg == "-float") {
            std::cerr << "Convert: to float\n";
            if (uc1Image* uc = dynamic_cast<uc1Image*>(curImg.get()))
                curImg = std::shared_ptr<baseImage>(new f1Image(*uc));
            else if (uc3Image* uc = dynamic_cast<uc3Image*>(curImg.get()))
                curImg = std::shared_ptr<baseImage>(new f3Image(*uc));
            else if (uc4Image* uc = dynamic_cast<uc4Image*>(curImg.get()))
                curImg = std::shared_ptr<baseImage>(new f4Image(*uc));
            else if (!curImg->is_integer_virtual())
                std::cerr << "Already float\n";
            else
                throw DMcError("Unsupported image type.\n");

            RemoveArgs(argc, argv, i);
        } else if (starg == "-uchar") {
            std::cerr << "Convert: to uchar\n";
            if (f1Image* f = dynamic_cast<f1Image*>(curImg.get()))
                curImg = std::shared_ptr<baseImage>(new uc1Image(*f));
            else if (f3Image* f = dynamic_cast<f3Image*>(curImg.get()))
                curImg = std::shared_ptr<baseImage>(new uc3Image(*f));
            else if (f4Image* f = dynamic_cast<f4Image*>(curImg.get()))
                curImg = std::shared_ptr<baseImage>(new uc4Image(*f));
            else if (!curImg->is_integer_virtual())
                std::cerr << "Already uchar\n";
            else
                throw DMcError("Unsupported image type.\n");

            RemoveArgs(argc, argv, i);
        } else if (starg == "-sp") {
            std::cerr << "LoadSaveParams:"
                      << " maxColorPalette=" << SP.QP.maxColorPalette << " maxItersFast=" << SP.QP.maxItersFast << " maxIters=" << SP.QP.maxIters
                      << " targetErr=" << SP.QP.targetErr << " makeArtisticPalette=" << SP.QP.makeArtisticPalette << '\n';

            if (argc <= i + 5) Usage();
            SP.QP.maxColorPalette = atoi(argv[i + 1]);
            SP.QP.maxItersFast = atoi(argv[i + 2]);
            SP.QP.maxIters = atoi(argv[i + 3]);
            SP.QP.targetErr = atof(argv[i + 4]);
            SP.QP.makeArtisticPalette = atoi(argv[i + 5]);

            RemoveArgs(argc, argv, i, 6);
        } else if (starg == "-threshold") {
            if (argc <= i + 1) Usage();
            float threshold = atof(argv[i + 1]);

            if (uc1Image* uc1I = dynamic_cast<uc1Image*>(curImg.get())) {
                uc1Pixel threshPix(static_cast<unsigned char>(threshold));
                DoThreshold(*uc1I, threshPix);
            } else if (uc3Image* uc3I = dynamic_cast<uc3Image*>(curImg.get())) {
                uc3Pixel threshPix(static_cast<unsigned char>(threshold));
                DoThreshold(*uc3I, threshPix);
            } else {
                throw DMcError("Unsupported image type.\n");
            }

            RemoveArgs(argc, argv, i, 2);
        } else if (starg == "-resize") {
            if (argc <= i + 2) Usage();

            int newWid = atoi(argv[i + 1]);
            int newHgt = atoi(argv[i + 2]);

            if (newWid == 0) newWid = int(float(newHgt) * float(curImg->w_virtual()) / float(curImg->h_virtual()));
            if (newHgt == 0) {
                if (std::string(argv[i + 2]) == "X") {
                    float s = atof(argv[i + 1]);
                    newWid = int(s * float(curImg->w_virtual()));
                    newHgt = int(s * float(curImg->h_virtual()));
                } else
                    newHgt = int(float(newWid) * float(curImg->h_virtual()) / float(curImg->w_virtual()));
            }

            curImg = DoResize(curImg, newWid, newHgt);

            RemoveArgs(argc, argv, i, 3);
        } else if (starg == "-resizeca") {
            if (argc <= i + 2) Usage();

            int newWid = atoi(argv[i + 1]);
            int newHgt = atoi(argv[i + 2]);

            if (newWid == 0) newWid = int(float(newHgt) * float(curImg->w_virtual()) / float(curImg->h_virtual()));
            if (newHgt == 0) {
                if (std::string(argv[i + 2]) == "X") {
                    float s = atof(argv[i + 1]);
                    newWid = int(s * float(curImg->w_virtual()));
                    newHgt = int(s * float(curImg->h_virtual()));
                } else
                    newHgt = int(float(newWid) * float(curImg->h_virtual()) / float(curImg->w_virtual()));
            }

            curImg = DoResizeCA(curImg, newWid, newHgt);

            RemoveArgs(argc, argv, i, 3);
        } else if (starg == "-blur") {
            if (argc <= i + 2) Usage();

            int filtWid = atoi(argv[i + 1]);
            float imageStDev = atof(argv[i + 2]);
            curImg = DoBlur(curImg, filtWid, imageStDev);

            RemoveArgs(argc, argv, i, 3);
        } else if (starg == "-vcd") {
            if (argc <= i + 4) Usage();

            int filtWid = atoi(argv[i + 1]);
            float imageStDev = atof(argv[i + 2]);
            float colorStDev = atof(argv[i + 3]);
            int iterations = atoi(argv[i + 4]);
            curImg = DoVCD(curImg, filtWid, imageStDev, colorStDev, iterations);

            RemoveArgs(argc, argv, i, 5);
        } else if (starg == "-despecklering") {
            if (argc <= i + 4) Usage();

            int filtWid = atoi(argv[i + 1]);
            int ringThick = atoi(argv[i + 2]);
            float threshold = atof(argv[i + 3]);
            int iterations = atoi(argv[i + 4]);
            curImg = DoRingDespeckleFilter(curImg, filtWid, ringThick, threshold, iterations);

            RemoveArgs(argc, argv, i, 5);
        } else if (starg == "-findlinmap") {
            if (argc <= i) Usage();

            if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
                DoOptimizeLinearMap(*f3I);
            } else if (f4Image* f4I = dynamic_cast<f4Image*>(curImg.get())) {
                DoOptimizeLinearMap(*f4I);
            } else {
                throw DMcError("Unsupported image type.\n");
            }

            RemoveArgs(argc, argv, i, 7);
        } else if (starg == "-linmap") {
            if (argc <= i + 6) Usage();

            float r0 = atof(argv[i + 1]);
            float g0 = atof(argv[i + 2]);
            float b0 = atof(argv[i + 3]);
            float r1 = atof(argv[i + 4]);
            float g1 = atof(argv[i + 5]);
            float b1 = atof(argv[i + 6]);
            if (f3Image* f3I = dynamic_cast<f3Image*>(curImg.get())) {
                DoLinearMap(*f3I, f3Pixel(r0, g0, b0), f3Pixel(r1, g1, b1));
            } else if (f4Image* f4I = dynamic_cast<f4Image*>(curImg.get())) {
                DoLinearMap(*f4I, f3Pixel(r0, g0, b0), f3Pixel(r1, g1, b1));
            } else {
                throw DMcError("Unsupported image type.\n");
            }

            RemoveArgs(argc, argv, i, 7);
        } else if (starg == "-despeckle") {
            if (argc <= i + 0) Usage();

            curImg = DoDespeckle(curImg);

            RemoveArgs(argc, argv, i, 1);
        } else if (starg == "-out") {
            if (argc <= i + 1) Usage();

            std::string dstFName(argv[i + 1]);
            std::cerr << "Save: " << dstFName << '\n';
            curImg->Save(dstFName.c_str(), SP);

            RemoveArgs(argc, argv, i, 2);
        } else if (starg == "-ld") {
            if (argc <= i + 1) Usage();
            int slot = atoi(argv[i + 1]);
            std::cerr << "Getting image from slot " << slot << ".\n";
            curImg = Slots.get(slot, true);

            RemoveArgs(argc, argv, i, 2);
        } else if (starg == "-st") {
            if (argc <= i + 1) Usage();
            int slot = atoi(argv[i + 1]);
            std::cerr << "Storing copy of curImage in slot " << slot << ".\n";
            Slots.set(slot, std::shared_ptr<baseImage>(curImg->Copy()));

            RemoveArgs(argc, argv, i, 2);
        } else {
            Usage(("Unknown option:" + starg).c_str());
        }
    }
}

int main(int argc, char** argv)
{
    std::cerr << "Starting " << argv[0] << "\n";
    try {
        Driver(argc, argv);
    }
    catch (DMcError& Er) {
        std::cerr << "DMcError caught: " << Er.Er << std::endl;
    }
    catch (std::string& Er) {
        std::cerr << "String caught: " << Er << std::endl;
    }
    catch (...) {
        std::cerr << "Caught exception.\n";
    }

    return 0;
}
