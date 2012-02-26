//------------------------------------------------------------------------------
//
// Contributors: 
//             1) Tommy Hinks
//
//------------------------------------------------------------------------------

#include "salsa_wimage.hpp"
#include <D:/las.git/include/las.hpp>
#include <il.h>
#include <omp.h>
#include <limits>
#include <algorithm>
#include <iostream>
#include <iomanip>

//------------------------------------------------------------------------------

typedef salsa::image<unsigned char>  image;
typedef salsa::image<float>          imagef;
typedef salsa::wimage<unsigned char> wimage;
typedef salsa::wimage<float>         wimagef;

using std::vector;
using std::size_t;
using std::numeric_limits;
using std::min;
using std::max;
using std::string;
using std::stringstream;
using std::cout;
using std::cerr;
using std::fixed;
using std::setprecision;
using std::exception;
using std::abort;

//------------------------------------------------------------------------------

const string ftFileNames[] = {
//"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856085600-856085657_corridor.las",
//"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856085992-856086070_corridor.las",
//"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856086261-856086321_corridor.las",
//"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856086674-856086754_corridor.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856086934-856086971.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856087073-856087074.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856087076-856087142.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856087221-856087287.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856087387-856087470.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856087577-856087642.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856087786-856087874.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856087989-856088055.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856088166-856088257.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856088379-856088442.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856088527-856088621.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856088624-856088628.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856088727-856088787.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856088866-856088949.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856089019-856089071.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856089167-856089243.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856089352-856089404.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856089512-856089591.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856089678-856089731.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856089823-856089900.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856089990-856090026.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856090058-856090097.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856094688-856094717.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856094803-856094857.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856094925-856094977.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856095097-856095169.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856095232-856095317.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856095407-856095489.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856095595-856095684.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856095784-856095862.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856095980-856096056.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856096144-856096214.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856096328-856096399.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856096508-856096571.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856096712-856096790.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856096909-856096986.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856097081-856097144.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856097262-856097268.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856097270-856097320.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856097413-856097447.las"
};

//------------------------------------------------------------------------------

BEGIN_STD_NAMESPACE
bool
operator<(const las::point &lhs, const las::point &rhs)
{ return lhs.gps_time < rhs.gps_time; }
END_STD_NAMESPACE

//------------------------------------------------------------------------------

//! TODO: check errors!
void
saveImage(const image &img, const string &fileName)
{
    ILuint imgHandle = 0;
    ilGenImages(1, &imgHandle);
    ilBindImage(imgHandle);
    ilTexImage(static_cast<ILuint>(img.width()), 
               static_cast<ILuint>(img.height()), 
               1, // z  = 1.
               sizeof(unsigned char), 
               IL_LUMINANCE, 
               IL_UNSIGNED_BYTE, 
               reinterpret_cast<void*>(const_cast<unsigned char*>(&img[0])));
    ilEnable(IL_FILE_OVERWRITE);
    ilSaveImage(fileName.c_str());
    ilDeleteImages(1, &imgHandle); 
    // TODO: check errors!!
}

//! Find min/max values in image.
template<typename P>
void 
imgMinMax(const salsa::image<P> &img, P &imgMin, P &imgMax)
{
    const int n = static_cast<int>(img.size());
    P shMin;
    P shMax;
#pragma omp parallel
    {
    P thMin =  (numeric_limits<P>::max)();
    P thMax = -(numeric_limits<P>::max)();
#pragma omp for nowait
    for(int i = 0; i < n; ++i) {
        thMin = min<P>(img[i], thMin);
        thMax = max<P>(img[i], thMax);
    }
#pragma omp critical
    {
        shMin = min<P>(shMin, thMin);
        shMax = max<P>(shMax, thMax);
    }
    }
    imgMin = shMin;
    imgMax = shMax;
}

void
saveImage(const imagef &imgf, const string &fileName)
{
    float imgMin =  (numeric_limits<float>::max)();
    float imgMax = -(numeric_limits<float>::max)();
    imgMinMax(imgf, imgMin, imgMax);

    image img(imgf.width(), imgf.height(), 0);
    const float div = 1.f/(imgMax - imgMin);
    const int size = static_cast<int>(imgf.size());
    int i;
#pragma omp parallel for
    for (i = 0; i < size; ++i) {
        img[i] = static_cast<unsigned char>(
            (numeric_limits<unsigned char>::max)()*((imgf[i] - imgMin)*div));
    }
    saveImage(img, fileName);
}

void
saveImage(const image  &img, 
          const string &fileName,
          const string &suffix, 
          const int     frame,
          const string &ext = "png")
{
    stringstream ss;
    ss << fileName << (suffix.empty() ? "" : "_") << suffix 
       << frame << "." << ext;
    saveImage(img, ss.str());
}


void
saveImage(const imagef &imgf, 
          const string &fileName,
          const string &suffix, 
          const int     frame,
          const string &ext = "png")
{
    stringstream ss;
    ss << fileName << (suffix.empty() ? "" : "_") << suffix 
       << frame << "." << ext;
    saveImage(imgf, ss.str());
}

//------------------------------------------------------------------------------

void
addElevation(const vector<las::point> &p, 
             const double wzmin, const double wzmax, 
             wimage &wimg)
{
    const double f = 1./(wzmax - wzmin);
    const size_t size = p.size();
    for (size_t i = 0; i < size; ++i) {
        if (p[i].return_num == 1) { // First return only.
            const unsigned char pval1 = 
                static_cast<unsigned char>(
                   (numeric_limits<unsigned char>::max)()*((p[i].z - wzmin)*f));
            if (pval1 > wimg(p[i].x, p[i].y)) {
                wimg(p[i].x, p[i].y) = pval1;
            }
        }
    }
}

//------------------------------------------------------------------------------

//! Transfer visibility to main image.
void
accumVisibility(const wimagef &vbuf, wimagef &vis)
{
    if (vbuf.image().size() == vis.image().size()) {
        int i;
        const int size = static_cast<int>(vbuf.image().size());
#pragma omp parallel for
        for (i = 0; i < size; ++i) {
            vis[i] += min<float>(vbuf[i], 1.f);
        }
    }
}

//! Add visibility for a chunk of points.
void
addVisibility(const vector<las::point> &p, 
              const size_t              i0, 
              const size_t              i1, 
              wimagef                  &vbuf,
              wimagef                  &vis)
{
    if (i0 <= i1 && i1 <= p.size()) {
        vbuf.image().clear(0.f);
        for (size_t i = i0; i < i1; ++i) {
            if (p[i].return_num <= p[i].num_returns) {
                vbuf(p[i].x, p[i].y) += p[i].return_num/p[i].num_returns;
            }
        }
        accumVisibility(vbuf, vis);
    }
}

//------------------------------------------------------------------------------

void
addScanLineVisibility(vector<las::point> &ft, 
                      const size_t        fts0,
                      const size_t        fts1,
                      const size_t        fts,
                      wimagef            &vbuf,
                      wimagef            &vaccum)
{
    if (fts0 <= fts1 && fts1 <= ft.size()) {
        size_t i = fts0;
        double t0 = ft[i].gps_time;
        int i0 = i;
        int i1 = i0;
        for (; i < fts1; ++i) {
            const double t = ft[i].gps_time;
            if (t == t0) {
                ++i1;
            }
            else {
                addVisibility(ft, i0, i1, vbuf, vaccum);
                t0 = t;
                i0 = i + 1;
                i1 = i0;
            }
        }
        saveImage(vaccum.image(), "occlusion", "sl", fts);
    }
}


void
addSegmentVisibility(vector<las::point> &ft, 
                     size_t             &fts,
                     wimagef            &vbuf,
                     wimagef            &ftsAccum,
                     wimagef            &slAccum)
{
    if (!ft.empty()) {
        sort(ft.begin(), ft.end()); // TODO: parallel!

        size_t i = 0;
        double t0 = floor(ft[i].gps_time);
        int i0 = i;
        int i1 = i0;
        const size_t n = ft.size();
        for (; i < n; ++i) {
            const double t = floor(ft[i].gps_time);
            if (t == t0) {
                ++i1;
            }
            else {
                addVisibility(ft, i0, i1, vbuf, ftsAccum);
                saveImage(ftsAccum.image(), "occlusion", "fts", fts);

                addScanLineVisibility(ft, i0, i1, fts, vbuf, slAccum);

                t0 = t;
                i0 = i + 1;
                i1 = i0;
                ++fts;
            }
        }
    }
}





void
trackToSegments(vector<las::point>          &ft, 
                vector<vector<las::point> > &fts,
                const size_t                 scanRate,
                const size_t                 pulseRate,
                const size_t                 numEchoes)
{
    if (!ft.empty()) {
        sort(ft.begin(), ft.end());
        fts.push_back(vector<las::point>(1, ft[0]));
        fts.back().reserve(scanRate*pulseRate*numEchoes);
        const size_t size = ft.size();
        for (size_t i = 0; i < size; ++i) {
            if (floor(ft[i].gps_time) == floor(fts.back().back().gps_time)) {
                fts.back().push_back(ft[i]);
            }
            else {
                fts.push_back(vector<las::point>(1, ft[i]));
                fts.back().reserve(scanRate*pulseRate*numEchoes);
            }
        }
    }
}

void
trackToScanLines(vector<las::point>          &ft, 
                 vector<vector<las::point> > &sl,
                 const size_t                 pulseRate,
                 const size_t                 numEchoes)
{
    if (!ft.empty()) {
        sort(ft.begin(), ft.end());
        sl.push_back(vector<las::point>(1, ft[0]));
        sl.back().reserve(pulseRate*numEchoes);
        const size_t ftSize = ft.size();
        for (size_t i = 0; i < ftSize; ++i) {
            if (ft[i].gps_time == sl.back().back().gps_time) {
                sl.back().push_back(ft[i]);
            }
            else {
                sl.push_back(vector<las::point>(1, ft[i]));
                sl.back().reserve(pulseRate*numEchoes);
            }
        }
    }
}

void
segmentToScanLines(vector<las::point>          &fts, 
                   vector<vector<las::point> > &sl,
                   const size_t                 pulseRate,
                   const size_t                 numEchoes)
{
    if (!fts.empty()) {
        sort(fts.begin(), fts.end());
        sl.push_back(vector<las::point>(1, fts[0]));
        sl.back().reserve(pulseRate*numEchoes);
        const size_t ftsSize = fts.size();
        for (size_t i = 0; i < ftsSize; ++i) {
            if (fts[i].gps_time == sl.back().back().gps_time) {
                sl.back().push_back(fts[i]);
            }
            else {
                sl.push_back(vector<las::point>(1, fts[i]));
                sl.back().reserve(pulseRate*numEchoes);
            }
        }
    }
}

//------------------------------------------------------------------------------

// Input must be a set of LAS files, where each files is exactly one 
// flight track.

int
main(int argc, char *argv[]) 
{
    try {
        // Settings. TODO: should be command line arguments.

        const size_t scanRate = 150;
        const size_t pulseRate = 1000;
        const size_t numEchoes = 5;
        const double wdx = 0.5;
        const double wdy = 0.5;

        const double programStart = omp_get_wtime();
        cout << fixed << setprecision(3);

        ilInit();

        double wxmin =  (numeric_limits<double>::max)();
        double wxmax = -(numeric_limits<double>::max)();
        double wymin =  (numeric_limits<double>::max)();
        double wymax = -(numeric_limits<double>::max)();
        double wzmin =  (numeric_limits<double>::max)();
        double wzmax = -(numeric_limits<double>::max)();
        for (int ft = 0; ft < 40; ++ft) {
            las10::public_header_block phb;
            las10::read::header(ftFileNames[ft], phb);
            wxmin = min<double>(wxmin, phb.min_x);
            wxmax = max<double>(wxmax, phb.max_x);
            wymin = min<double>(wymin, phb.min_y);
            wymax = max<double>(wymax, phb.max_y);
            wzmin = min<double>(wzmin, phb.min_z);
            wzmax = max<double>(wzmax, phb.max_z);
        }

        wimage ftElevation(wxmin, wxmax, wymin, wymax, wdx, wdy, 0);
        wimage ftsElevation(wxmin, wxmax, wymin, wymax, wdx, wdy, 0);
        cout << ftElevation << "\n"; // Just to show image dims.

        wimagef ftVisAccum(wxmin, wxmax, wymin, wymax, wdx, wdy, 0.f);
        wimagef ftsVisAccum(wxmin, wxmax, wymin, wymax, wdx, wdy, 0.f);
        wimagef slVisAccum(wxmin, wxmax, wymin, wymax, wdx, wdy, 0.f);
        wimagef vbuf(wxmin, wxmax, wymin, wymax, wdx, wdy, 0.f);


        size_t ftsIndex = 0;
        int slIndex = 0;
        for (int ft = 0; ft < 40; ++ft) {
            const double ftStart = omp_get_wtime();
            cout << "[" << ft + 1 << "|40]: '" << ftFileNames[ft] << "'\n";
            vector<las::point> ftPoints;
            las10::read::points(ftFileNames[ft], ftPoints);

            addElevation(ftPoints, wzmin, wzmax, ftElevation);
            saveImage(ftElevation.image(), "elevation", "ft", ft);

            addVisibility(ftPoints, 0, ftPoints.size(), vbuf, ftVisAccum);
            saveImage(ftVisAccum.image(), "occlusion", "ft", ft);

            addSegmentVisibility(ftPoints, ftsIndex, vbuf, ftsVisAccum, slVisAccum);
            //saveImage(ftsVisAccum.image(), "occlusion", "fts", ftsIndex);

            const double ftEnd = omp_get_wtime();
            cout << ftPoints.size() << " points in " 
                 << ftEnd - ftStart << " [s]\n";

            /*
            vector<vector<las::point> > ftsPoints;
            trackToSegments(ftPoints, ftsPoints, scanRate, pulseRate, numEchoes);
            const size_t ftsSize = ftsPoints.size();
            for (size_t fts = 0; fts < ftsSize; ++fts) {
                addElevation(ftsPoints[fts], wzmin, wzmax, ftsElevation);
                saveImage(ftsElevation.image(), "elevation", "fts", ftsIndex);

                addVisibility(ftsPoints[fts], ftsVisAccum);
                saveImage(ftsVisAccum.image(), "occlusion", "fts", ftsIndex);

                vector<vector<las::point> > slPoints;
                segmentToScanLines(ftPoints, slPoints, pulseRate, numEchoes);
                const size_t slSize = slPoints.size();
                for (size_t sl = 0; sl < slSize; ++sl) {
                    addVisibility(slPoints[sl], slVisAccum);
                    //saveImage(slVisAccum.image(), "occlusion", "sl", slIndex);
                    ++slIndex;
                }
                saveImage(slVisAccum.image(), "occlusion", "sl", ftsIndex);
                ++ftsIndex;
            }

            */
        }

        const double programEnd = omp_get_wtime();
        cout << "Program execution time: " 
             << programEnd - programStart << " [s]\n";

        return 0;
    }   
    catch (exception& ex) {
        cerr << "Exception in main(): " << ex.what() << "\n";
        abort();
    }
    catch (...) {
        cerr << "Unknown exception. Aborting...\n";
        abort();
    }
}
