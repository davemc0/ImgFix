// This program does one or more chained edit commands on srcFName and saves the result to dstFName
// Copyright David K. McAllister, 2009-2012.

// Release/ImgFix.exe -in vcd2.png -templatefilter 35 1 189 -threshold 189 -blur 3 1 -out out.png

#include "FixOps.h"

#include "Image/ImageAlgorithms.h"
#include "Math/Random.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

static void Usage(const char *message = NULL, const bool Exit = true)
{
    if(message)
        cerr << "\nERROR: " << message << endl;

    cerr << "Usage: <options>\n\n";
    cerr << "-in fname.jpg                      Load the current image\n";
    cerr << "-blur f s                          Blur filt_width img_stdev; try 5 2.0\n";
    cerr << "-chan N                            Convert image to N channel\n";
    cerr << "-copy fname.jpg dx dy sx sy w h M  Copy a subrect of the given image to an arbitrary place in the cur image;\n";
    cerr << "                                        M:0=copy 1=over blend 2=key (int) r g b a 3=multiply 4=blend a\n";                                        // Add a copy mode to copy just a given channel
    cerr << "-templatefilter f c t i            Template despeckle box_width ring_thickness threshold iters; random smaller boxes if iter>1; try 15 2 180 1\n";
    cerr << "-despeckle                         Clamp each pixel to be within the bounding box of its eight neighbors\n";
    cerr << "-fill r g b a                      Fill image with specified color (0..255 or 0..1.0)\n";
    cerr << "-flatten shrink bias               Flatten; shrink is blurred image rel. size; bias blends btw. bias and scale methods; try 0.1 0.5\n";
    cerr << "-flip h|v                          Flip image horizontally or vertically\n";
    cerr << "-float                             Convert image to float\n";
    cerr << "-gradient N h|v minx maxx          Set channel N to a horiz or vert gradient with 0 and 255 at given scanlines\n";
    cerr << "-hflatten a b                      Flatten horizontal rows to remove scanner artifacts; ignores params\n";
    cerr << "-makeuc N w h                      Make a current ucImage of N channels, width and height; fill with black\n";
    cerr << "-matmul 0123456789ab[cdef]         Multiply each pixel value by this 4x3 (for f3Image) or 4x4 (for f4Image) matrix\n";
    cerr << "-noise mean stdev                  Fill image with noise\n";
    cerr << "-resize w|'x' h|'y'|'X'            Resize the given image to w x h. Bicubic upsampling; box filter then bicubic downsampling\n";
    cerr << "                                       'X' => scale factor; 'x'|'y' => compute this dim\n";
    cerr << "-caresize w|'x' h|'y'|'X'          Content-Aware Resize the given image to w x h. Seam carving; currently just shrinking\n";
    cerr << "                                       'X' => scale factor; 'x'|'y' => compute this dim\n";
	cerr << "-linmap r0 g0 b0 r1 g1 b1          Project f3Image pixel values onto a line in RGB space\n";
	cerr << "-findlinmap                        Find best fit line to project f3Image pixel values\n";
	cerr << "-threshold N                       Convert image to two-tone at given threshold value\n";
	cerr << "-transect y                        Print channel 0 values of scanline y to stdout in CSV format\n";
    cerr << "-vcd f s c i                       VCD filt_width img_stdev col_stdev iterations; try 5 3.0 0.25 1\n";
    cerr << "-uchar                             Convert image to uchar\n";
    cerr << "-out fname.jpg                     Save the current image\n";
    cerr << "-outbin fname.bin                  Dump the pixel bytes to a raw file\n";

    if(Exit)
        exit(1);
}

static void Driver(int& argc, char** argv)
{
    baseImage* curImg = NULL;

    if(argc < 2)
        Usage();

    for(int i = 1; i < argc; i++) {
        if(curImg) {
            std::cerr << "curImg: " << " size=" << curImg->w_virtual() << "x" << curImg->h_virtual()
                << (curImg->is_integer_virtual() ? " INT" : " FLOAT") << " chan=" << curImg->chan_virtual()
                << " bytes_per_element=" << curImg->size_element_virtual() << '\n';
        }

        if(string(argv[i]) == "-h" || string(argv[i]) == "-help") {
            Usage();
        } else if(string(argv[i]) == "-in") {
            if(argc <= i + 1) Usage();

            std::string srcFName(argv[i + 1]);
            curImg = LoadtImage(srcFName.c_str());
            std::cerr << "Load: " << srcFName << " size=" << curImg->w_virtual() << "x" << curImg->h_virtual() << " chan=" << curImg->chan_virtual() << '\n';

            RemoveArgs(argc, argv, i, 2);
        } else if(string(argv[i]) == "-makeuc") {
            if(argc <= i + 3) Usage();

            int chan = atoi(argv[i + 1]), w = atoi(argv[i + 2]), h = atoi(argv[i + 3]);

            cerr << "Make: uc" << chan << "Image\n";

            baseImage* newImg = NULL;
            if(chan == 1) newImg = new uc1Image(w, h, uc1Pixel(0));
            if(chan == 3) newImg = new uc3Image(w, h, uc3Pixel(0));
            if(chan == 4) newImg = new uc4Image(w, h, uc4Pixel(0));
            curImg = replaceImage(curImg, newImg);
            RemoveArgs(argc, argv, i, 4);
        } else if(string(argv[i]) == "-noise") {
            if(argc <= i + 2) Usage();

            float mean = atof(argv[i + 1]), stdev = atof(argv[i + 2]);
            curImg = DoNoise(curImg, mean, stdev);

            RemoveArgs(argc, argv, i, 3);
        } else if(string(argv[i]) == "-fill") {
            if(argc <= i + 4) Usage();
            float r = atof(argv[i + 1]);
            float g = atof(argv[i + 2]);
            float b = atof(argv[i + 3]);
            float a = atof(argv[i + 4]);
            DoFill(curImg, r, g, b, a);

            RemoveArgs(argc, argv, i, 5);
        } else if(string(argv[i]) == "-copy") {
            if(argc <= i + 8) Usage();
            int dx = atoi(argv[i + 2]), dy = atoi(argv[i + 3]),
                sx = atoi(argv[i + 4]), sy = atoi(argv[i + 5]),
                sw = atoi(argv[i + 6]), sh = atoi(argv[i + 7]),
                md = atoi(argv[i + 8]);
            uc4Pixel key(element_traits<uc4Pixel::ElType>::one());
            float alpha = 1.0f;

            int NArg = 0;
            if(md == 2) {
                NArg = curImg->chan_virtual();
                if(argc <= i + 8 + NArg) Usage("Need ints for each channel for key color");
                for(int c = 0; c < NArg; c++) {
                    key[c] = atoi(argv[i + 9 + c]);
                }
            } else if(md == 4) {
                NArg = 1;
                if(argc <= i + 8 + NArg) Usage("Need a float for alpha");
                alpha = atof(argv[i + 8 + NArg]);
            }

            DoCopyBlock(curImg, std::string(argv[i + 1]), dx, dy, sx, sy, sw, sh, md, key, alpha);

            RemoveArgs(argc, argv, i, 9 + NArg);
        } else if(string(argv[i]) == "-matmul") {
            int NArg = 0;
            if(f3Image* f3I = dynamic_cast<f3Image*>(curImg)) {
                NArg = 12;
                if(argc <= i + NArg) Usage();
                DoMatMul(*f3I, atof(argv[i + 1]), atof(argv[i + 2]), atof(argv[i + 3]), atof(argv[i + 4]), atof(argv[i + 5]), atof(argv[i + 6]), atof(argv[i + 7]), atof(argv[i + 8]), atof(argv[i + 9]), atof(argv[i + 10]), atof(argv[i + 11]), atof(argv[i + 12]), 0, 0, 0, 1.0f);
            } else if(f4Image* f4I = dynamic_cast<f4Image*>(curImg)) {
                NArg = 16;
                if(argc <= i + NArg) Usage();
                DoMatMul(*f4I, atof(argv[i + 1]), atof(argv[i + 2]), atof(argv[i + 3]), atof(argv[i + 4]), atof(argv[i + 5]), atof(argv[i + 6]), atof(argv[i + 7]), atof(argv[i + 8]), atof(argv[i + 9]), atof(argv[i + 10]), atof(argv[i + 11]), atof(argv[i + 12]), atof(argv[i + 13]), atof(argv[i + 14]), atof(argv[i + 15]), atof(argv[i + 16]));
            } else {
                throw DMcError("Unsupported image type.\n");
            }

            RemoveArgs(argc, argv, i, NArg + 1);
        } else if(string(argv[i]) == "-gradient") {
            // Modifies image in place
            if(argc <= i + 4) Usage();
            int c = atoi(argv[i + 1]), minx = atoi(argv[i + 3]), maxx = atoi(argv[i + 4]);
            bool isVert = argv[i + 2][0] == 'v';
            if(minx == maxx)
                Usage("Gradient must have different min and max scanline numbers");
            DoGradient(curImg, c, isVert, minx, maxx);

            RemoveArgs(argc, argv, i, 5);
        } else if(string(argv[i]) == "-flip") {
            if(argc <= i + 1) Usage();
            bool isVert = (argv[i + 1][0] == 'v' || argv[i + 1][0] == 'V');
            DoFlip(curImg, isVert);

            RemoveArgs(argc, argv, i, 2);
        } else if(string(argv[i]) == "-transect") {
            if(argc <= i + 1) Usage();
            int yval = atoi(argv[i + 1]);

            if(uc1Image* uc1I = dynamic_cast<uc1Image*>(curImg)) {
                DoPrintTransect(*uc1I, yval);
            } else if(uc3Image* uc3I = dynamic_cast<uc3Image*>(curImg)) {
                DoPrintTransect(*uc3I, yval);
            } else {
                throw DMcError("Unsupported image type.\n");
            }

            RemoveArgs(argc, argv, i, 2);
        } else if(string(argv[i]) == "-flatten") {
            if(argc <= i + 2) Usage();
            float shrinkFac = atof(argv[i + 1]);
            float biasFac = atof(argv[i + 2]);

            baseImage* dstImg = NULL;
            if(f1Image* f1I = dynamic_cast<f1Image*>(curImg)) {
                f1Image* dstImg = new f1Image;
                DoFlatten(*dstImg, *f1I, shrinkFac, biasFac);
                curImg = replaceImage(curImg, dstImg);
            } else if(f3Image* f3I = dynamic_cast<f3Image*>(curImg)) {
                f3Image* dstImg = new f3Image;
                DoFlatten(*dstImg, *f3I, shrinkFac, biasFac);
                curImg = replaceImage(curImg, dstImg);
            } else {
                throw DMcError("Unsupported image type.\n");
            }

            RemoveArgs(argc, argv, i, 3);
        } else if(string(argv[i]) == "-hflatten") {
            if(argc <= i + 2) Usage();
            float a = atof(argv[i + 1]);
            float b = atof(argv[i + 2]);

            baseImage* dstImg = NULL;
            if(f1Image* f1I = dynamic_cast<f1Image*>(curImg)) {
                f1Image* dstImg = new f1Image;
                DoHorizFlatten(*dstImg, *f1I, a, b);
                curImg = replaceImage(curImg, dstImg);
            } else if(f3Image* f3I = dynamic_cast<f3Image*>(curImg)) {
                f3Image* dstImg = new f3Image;
                DoHorizFlatten(*dstImg, *f3I, a, b);
                curImg = replaceImage(curImg, dstImg);
            } else {
                throw DMcError("Unsupported image type.\n");
            }

            RemoveArgs(argc, argv, i, 3);
        } else if(string(argv[i]) == "-chan") {
            if(argc <= i + 1) Usage();
            curImg = DoChanConvert(curImg, atoi(argv[i + 1]));

            RemoveArgs(argc, argv, i, 2);
        } else if(string(argv[i]) == "-float") {
            cerr << "Convert: to float\n";
            if(uc1Image* uc = dynamic_cast<uc1Image*>(curImg))
                curImg = replaceImage(curImg, new f1Image(*uc));
            else if(uc3Image* uc = dynamic_cast<uc3Image*>(curImg))
                curImg = replaceImage(curImg, new f3Image(*uc));
            else if(uc4Image* uc = dynamic_cast<uc4Image*>(curImg))
                curImg = replaceImage(curImg, new f4Image(*uc));
            else if(!curImg->is_integer_virtual())
                std::cerr << "Already float\n";
            else
                throw DMcError("Unsupported image type.\n");

            RemoveArgs(argc, argv, i);
        } else if(string(argv[i]) == "-uchar") {
            cerr << "Convert: to uchar\n";
            if(f1Image* f = dynamic_cast<f1Image*>(curImg))
                curImg = replaceImage(curImg, new uc1Image(*f));
            else if(f3Image* f = dynamic_cast<f3Image*>(curImg))
                curImg = replaceImage(curImg, new uc3Image(*f));
            else if(f4Image* f = dynamic_cast<f4Image*>(curImg))
                curImg = replaceImage(curImg, new uc4Image(*f));
            else if(!curImg->is_integer_virtual())
                std::cerr << "Already float\n";
            else
                throw DMcError("Unsupported image type.\n");

            RemoveArgs(argc, argv, i);
        } else if(string(argv[i]) == "-threshold") {
            if(argc <= i + 1) Usage();
            float threshold = atof(argv[i + 1]);

            if(uc1Image* uc1I = dynamic_cast<uc1Image*>(curImg)) {
                uc1Pixel threshPix(static_cast<unsigned char>(threshold));
                DoThreshold(*uc1I, threshPix);
            } else if(uc3Image* uc3I = dynamic_cast<uc3Image*>(curImg)) {
                uc3Pixel threshPix(static_cast<unsigned char>(threshold));
                DoThreshold(*uc3I, threshPix);
            } else {
                throw DMcError("Unsupported image type.\n");
            }

            RemoveArgs(argc, argv, i, 2);
        } else if(string(argv[i]) == "-resize") {
            if(argc <= i + 2) Usage();

            int newWid = atoi(argv[i + 1]);
            int newHgt = atoi(argv[i + 2]);

            if(newWid == 0)
                newWid = int(float(newHgt) * float(curImg->w_virtual()) / float(curImg->h_virtual()));
            if(newHgt == 0) {
                if(string(argv[i + 2]) == "X") {
                    float s = atof(argv[i + 1]);
                    newWid = int(s * float(curImg->w_virtual()));
                    newHgt = int(s * float(curImg->h_virtual()));
                } else
                    newHgt = int(float(newWid) * float(curImg->h_virtual()) / float(curImg->w_virtual()));
            }

            curImg = DoResize(curImg, newWid, newHgt);

            RemoveArgs(argc, argv, i, 3);
        } else if(string(argv[i]) == "-caresize") {
            if(argc <= i + 2) Usage();

            int newWid = atoi(argv[i + 1]);
            int newHgt = atoi(argv[i + 2]);

            if(newWid == 0)
                newWid = int(float(newHgt) * float(curImg->w_virtual()) / float(curImg->h_virtual()));
            if(newHgt == 0) {
                if(string(argv[i + 2]) == "X") {
                    float s = atof(argv[i + 1]);
                    newWid = int(s * float(curImg->w_virtual()));
                    newHgt = int(s * float(curImg->h_virtual()));
                } else
                    newHgt = int(float(newWid) * float(curImg->h_virtual()) / float(curImg->w_virtual()));
            }

            curImg = DoCAResize(curImg, newWid, newHgt);

            RemoveArgs(argc, argv, i, 3);
        } else if(string(argv[i]) == "-blur") {
            if(argc <= i + 2) Usage();

            int filtWid = atoi(argv[i + 1]);
            float imageStDev = atof(argv[i + 2]);
            curImg = DoBlur(curImg, filtWid, imageStDev);

            RemoveArgs(argc, argv, i, 3);
        } else if(string(argv[i]) == "-vcd") {
            if(argc <= i + 4) Usage();

            int filtWid = atoi(argv[i + 1]);
            float imageStDev = atof(argv[i + 2]);
            float colorStDev = atof(argv[i + 3]);
            int iterations = atoi(argv[i + 4]);
            curImg = DoVCD(curImg, filtWid, imageStDev, colorStDev, iterations);

            RemoveArgs(argc, argv, i, 5);
		}
		else if (string(argv[i]) == "-templatematch") {
			if (argc <= i + 4) Usage();

			int filtWid = atoi(argv[i + 1]);
			int clearWid = atof(argv[i + 2]); // XXX Why is atof assigned to int?
			float threshold = atof(argv[i + 3]);
			int iterations = atoi(argv[i + 4]);
			curImg = DoTemplateMatchFilter(curImg, iterations, filtWid, clearWid, threshold);

			RemoveArgs(argc, argv, i, 5);
		}
		else if (string(argv[i]) == "-findlinmap") {
			if (argc <= i) Usage();

			if (f3Image* f3I = dynamic_cast<f3Image*>(curImg)) {
				DoOptimizeLinearMap(*f3I);
			}
			else if (f4Image* f4I = dynamic_cast<f4Image*>(curImg)) {
				DoOptimizeLinearMap(*f4I);
			}
			else {
				throw DMcError("Unsupported image type.\n");
			}

			RemoveArgs(argc, argv, i, 7);
		}
		else if (string(argv[i]) == "-linmap") {
			if (argc <= i + 6) Usage();

			float r0 = atof(argv[i + 1]);
			float g0 = atof(argv[i + 2]);
			float b0 = atof(argv[i + 3]);
			float r1 = atof(argv[i + 4]);
			float g1 = atof(argv[i + 5]);
			float b1 = atof(argv[i + 6]);
			if (f3Image* f3I = dynamic_cast<f3Image*>(curImg)) {
				DoLinearMap(*f3I, f3Pixel(r0, g0, b0), f3Pixel(r1, g1, b1));
			}
			else if (f4Image* f4I = dynamic_cast<f4Image*>(curImg)) {
				DoLinearMap(*f4I, f3Pixel(r0, g0, b0), f3Pixel(r1, g1, b1));
			}
			else {
				throw DMcError("Unsupported image type.\n");
			}

			RemoveArgs(argc, argv, i, 7);
		}
		else if (string(argv[i]) == "-despeckle") {
            if(argc <= i + 0) Usage();

            curImg = DoDespeckle(curImg);

            RemoveArgs(argc, argv, i, 1);
        } else if(string(argv[i]) == "-out") {
            if(argc <= i + 1) Usage();

            std::string dstFName(argv[i + 1]);
            std::cerr << "Save: " << dstFName << '\n';
            curImg->Save(dstFName.c_str());

            RemoveArgs(argc, argv, i, 2);
        } else if(string(argv[i]) == "-outbin") {
            if(argc <= i + 1) Usage();

            std::string dstFName(argv[i + 1]);
            std::cerr << "Save: " << dstFName << " " << curImg->size_bytes_virtual() << " bytes\n";
            std::ofstream OutFile(dstFName.c_str(), ios::out | ios::binary);
            if(!OutFile.is_open()) throw DMcError("Could not open file for writing: " + dstFName);
            OutFile.write((char *)curImg->pv_virtual(), curImg->size_bytes_virtual());
            OutFile.close();

            RemoveArgs(argc, argv, i, 2);
        } else {
            Usage(("Unknown option:" + string(argv[i])).c_str());
        }
    }
}

int main(int argc, char** argv)
{
    cerr << "Starting " << argv[0] << "\n";
    try {
        Driver(argc, argv);
    } catch(DMcError &Er) {
        cerr << "DMcError caught: " << Er.Er << endl;
    } catch(std::string &Er) {
        cerr << "String caught: " << Er << endl;
    } catch(...) {
        cerr << "Caught exception.\n";
    }

    return 0;
}