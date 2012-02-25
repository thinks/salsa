//------------------------------------------------------------------------------
//
// Contributors: 
//             1) Tommy Hinks
//
//------------------------------------------------------------------------------

#include "salsa_wimage.hpp"
#include <D:/las.git/include/las.hpp>
#include <il.h>
#include <limits>
#include <algorithm>
#include <iostream>

//------------------------------------------------------------------------------

const std::string ft[] = {
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856085600-856085657_corridor.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856085992-856086070_corridor.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856086261-856086321_corridor.las",
"D:/Dropbox/lidar/dublin-2007/tracks/las/ft/FT_856086674-856086754_corridor.las",
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

//! TODO: check errors!
void
saveImage(const salsa::image<unsigned char> &img, const std::string &fileName)
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

//------------------------------------------------------------------------------

void
elevation(const std::vector<las::point> &points, 
          const double                   zmin, 
          const double                   zmax, 
          salsa::wimage<unsigned char>  &wimg)
{
    const std::size_t size = points.size();
    for (std::size_t i = 0; i < size; ++i) {
        if (points[i].return_num == 1) {

            const double nz = (points[i].z - zmin)/(zmax - zmin);
            const unsigned char pval1 = std::numeric_limits<unsigned char>::max()*nz;
            const unsigned char pval0 = wimg(points[i].x, points[i].y);
            if (pval1 > wimg(points[i].x, points[i].y)) {
                wimg(points[i].x, points[i].y) = pval1;
            }
        }
    }
}


void
addVisibility(const std::vector<las::point> &p, 
              salsa::wimage<float>          &wimg)
{
    salsa::wimage<float> vis(wimg.wxmin(), wimg.wxmax(), 
                             wimg.wymin(), wimg.wymax(), 
                             wimg.wdx(), wimg.wdy(), 0.f);

    const std::size_t size = p.size();
    for (std::size_t i = 0; i < size; ++i) {
        if (p[i].return_num <= p[i].num_returns) {
            vis(p[i].x, p[i].y) += p[i].return_num/p[i].num_returns;
        }
    }

    // Transfer visibility to main image.

    for (std::size_t i = 0; i < wimg.image().size(); ++i) {
        wimg[i] += std::min<float>(vis[i], 1.f);
    }
}

void
occlusionPass1(const std::vector<las::point> &p, 
               salsa::wimage<float>          &wimg)
{
    salsa::wimage<float> tmp(wimg.wxmin(), wimg.wxmax(), 
                             wimg.wymin(), wimg.wymax(), 
                             wimg.wdx(), wimg.wdy(), 0.f);

    const std::size_t size = p.size();
    for (std::size_t i = 0; i < size; ++i) {
        if (p[i].return_num <= p[i].num_returns) {
            tmp(p[i].x, p[i].y) += p[i].return_num/p[i].num_returns;
        }
    }

    for (std::size_t x = 0; x < wimg.image().width(); ++x) {
        for (std::size_t y = 0; y < wimg.image().height(); ++y) {
            wimg(x, y) += std::min<float>(tmp(x, y), 1.f);
        }
    }
}

//! Normalize pixel values and convert image to unsigned byte.
void
occlusionPass2(const salsa::image<float>   &img0,
               salsa::image<unsigned char> &img1)
{
    float vmin =  (std::numeric_limits<float>::max)();
    float vmax = -(std::numeric_limits<float>::max)();
    for (std::size_t i = 0; i < img0.size(); ++i) {
        vmin = std::min<float>(vmin, img0[i]);
        vmax = std::max<float>(vmax, img0[i]);
    }

    for (std::size_t i = 0; i < img0.size(); ++i) {
        img1[i] = (std::numeric_limits<unsigned char>::max)()*((img0[i] - vmin)/(vmax - vmin));
    }
}


namespace std {
    bool
    operator<(const las::point &lhs, const las::point &rhs)
    { return lhs.gps_time < rhs.gps_time; }
}

void
trackToSegments(std::vector<las::point>       &ft, 
        std::vector<std::vector<las::point> > &fts,
        const std::size_t                      scanRate  = 150,
        const std::size_t                      pulseRate = 1000)
{
    if (!ft.empty()) {
        std::sort(ft.begin(), ft.end());

        fts.push_back(std::vector<las::point>(1, ft[0]));

        const std::size_t size = ft.size();
        for (std::size_t i = 0; i < size; ++i) {
            if (std::floor(ft[i].gps_time) == std::floor(fts.back().back().gps_time)) {
                fts.back().push_back(ft[i]);
            }
            else {
                fts.push_back(std::vector<las::point>(1, ft[i]));
                fts.back().reserve(scanRate*pulseRate);
            }
        }
    }
}

//ftsToSl(const std::vector<las::point>         &ft, 
//        std::vector<std::vector<las::point> > &fts)
//{
//}

//------------------------------------------------------------------------------

// Input must be a set of LAS files, where each files is exactly one 
// flight track.

int
main(int argc, char *argv[]) 
{
    try {
        ilInit();

        double wxmin =  (std::numeric_limits<double>::max)();
        double wxmax = -(std::numeric_limits<double>::max)();
        double wymin =  (std::numeric_limits<double>::max)();
        double wymax = -(std::numeric_limits<double>::max)();
        double wzmin =  (std::numeric_limits<double>::max)();
        double wzmax = -(std::numeric_limits<double>::max)();
        for (int i = 0; i < 44; ++i) {
            las10::public_header_block phb;
            las10::read::header(ft[i], phb);
            wxmin = std::min<double>(wxmin, phb.min_x);
            wxmax = std::max<double>(wxmax, phb.max_x);
            wymin = std::min<double>(wymin, phb.min_y);
            wymax = std::max<double>(wymax, phb.max_y);
            wzmin = std::min<double>(wzmin, phb.min_z);
            wzmax = std::max<double>(wzmax, phb.max_z);
        }

        const double wdx = 0.25;
        const double wdy = 0.25;
    
        salsa::wimage<unsigned char> wimg(wxmin, wxmax, wymin, wymax, wdx, wdy, 0);
        std::cerr << wimg << "\n";

        salsa::wimage<float> ftVisAccum(wxmin, wxmax, wymin, wymax, wdx, wdy, 0);
        salsa::wimage<float> ftsVisAccum(wxmin, wxmax, wymin, wymax, wdx, wdy, 0);
        salsa::wimage<float> slVisAccum(wxmin, wxmax, wymin, wymax, wdx, wdy, 0);

        for (int i = 0; i < 44; ++i) {
            std::vector<las::point> points;
            las10::read::points(ft[i], points);
            //elevation(points, wzmin, wzmax, wimg);

            std::vector<std::vector<las::point> > fts;
            trackToSegments(points, fts);
            for (std::size_t i = 0; i < fts.size(); ++i) {
                addVisibility(fts[i], ftsVisAccum);
            }

            addVisibility(points, ftVisAccum);
        }

        salsa::image<unsigned char> ftImg(ftVisAccum.image().width(), ftVisAccum.image().height(), 0);
        occlusionPass2(ftVisAccum.image(), ftImg);
        saveImage(ftImg, "occlusion_ft.png");

        salsa::image<unsigned char> ftsImg(ftsVisAccum.image().width(), ftsVisAccum.image().height(), 0);
        occlusionPass2(ftsVisAccum.image(), ftsImg);
        saveImage(ftsImg, "occlusion_fts.png");


        saveImage(wimg.image(), "elevation.png");
        // return success;

        return 0;
    }   
    catch (std::exception& ex) {
        std::cerr << ex.what() << "\n";
        std::abort();
    }
    catch (...) {
        std::cerr << "Unknown exception. Aborting...\n";
        std::abort();
    }
}
