#ifndef LEVEL_H
#define LEVEL_H

#include <stdexcept>
#include <vector>
#include <unordered_set>
#include <array>
#include <unordered_map>
#include <algorithm>

#include "PolarCoordinates.h"

class BoxTree;

class Level {

    friend class BoxTree;

    private:

        enum ConeRefinementStrategy {

            Laplace = 0,
            Helmholtz = 1

        };

        const int level_;
        const int nlevels_;
        const double min_x_;
        const double min_y_;
        const double boxsize_;
        const long long nboxesperside_;

        std::vector<long long> relconesegment2cocenteredmortonboxid_;
        std::vector<long long> relconesegment2nonrelconesegment_;

        double scaling_r2s_;
        long long nconeselevation_;
        std::vector<double> radial_intervals_cone_;

        std::unordered_set<long long> mortonidofrelboxes_;
        std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_;

        std::unordered_set<long long> propagationrequiredrelconesegments_;
        std::unordered_set<long long> interpolationrequiredrelconesegments_;

        std::unordered_map<long long, std::unordered_map<long long, long long>> mortonboxnonrelcone2relcone_;

        std::unordered_map<long long, std::array<std::unordered_map<long long, std::vector<int>>, 4>> protobox_; 

    public:

        Level(const int level, const int nlevels, 
              const double min_x, const double min_y,
              const double boxsize, const double wavenumber) :
              level_(level), nlevels_(nlevels), min_x_(min_x), min_y_(min_y),
              boxsize_(boxsize), nboxesperside_(1 << level_) {

            InitializeConeSegments(wavenumber);            

        }

        ~Level() {

        }

        template <ConeRefinementStrategy S = Helmholtz>
        void InitializeConeSegments(const double wavenumber) {

            const double minsize = (3.0/2.0 * boxsize_) * (1.0-1e-10);
            const double maxsize = (std::sqrt(2.0) * boxsize_ * nboxesperside_ - std::sqrt(2.0)/2.0 * boxsize_) * (1.0+1e-10);
            constexpr long long nconesradial = 1;
            constexpr long long nconeselevation = 2; 
            long long nconeslevelscaling = 1;

            if (S == Laplace) {

            } else if (S == Helmholtz) {

                if (level_ > 1)
                    nconeslevelscaling = std::max<long long>(1, std::ceil(wavenumber * boxsize_ * M_1_PI * 1.2)); //?

            } else
                throw std::invalid_argument("Unknown cone refinement strategy.");

            nconeselevation_ = nconeslevelscaling * nconeselevation;
            const long long nconesradialonlevel = nconeslevelscaling * nconesradial;

            scaling_r2s_ = std::sqrt(2.0)/2.0 * boxsize_;
            const double eta = std::sqrt(2.0)/3.0-1e-5;

            const double tmpsizeinr = (maxsize - minsize + 1e-10);
            const double tmpsizeins = std::min(eta/nconesradialonlevel, scaling_r2s_*(1.0/minsize - 1.0/(minsize+tmpsizeinr)));
            
            const long long nactualcones = std::ceil(scaling_r2s_*(1.0/minsize - 1.0/maxsize)/tmpsizeins);
            radial_intervals_cone_.resize(nactualcones+1);
            radial_intervals_cone_[0] = minsize;
    
            const double finalsizeins = (scaling_r2s_*(1.0/minsize - 1.0/maxsize)/nactualcones);
            for (long long i = nactualcones-1; i >= 0; i--) {
                const long long iter = nactualcones - i;       
                radial_intervals_cone_[iter] = scaling_r2s_*maxsize/(scaling_r2s_ + maxsize*i*finalsizeins);
            }

        }

        static long long Point2Morton (const double x, const double y,
                                       const double min_x, const double min_y,
                                       const double boxsize, const int level) {

            const long long posx = static_cast<long long>((x - min_x)/boxsize);
            const long long posy = static_cast<long long>((y - min_y)/boxsize);

            if (posx < 0 || posy < 0 || posx >= 1 << level || posy >= 1 << level)
                throw std::logic_error("The 3D box index cannot be < 0 or larger than the number of boxes");

            return Box2Morton(posx, posy, level);

        }

        long long  Point2Morton(const double x, const double y) const {

            return Point2Morton(x, y, min_x_, min_y_, boxsize_, level_);

        }

        static long long Box2Morton(const long long x, const long long y, const int level) {

            uint64_t answer = 0;
            answer |= splitBy2(x, level) | splitBy2(y, level) << 1;

            return answer;

        }

        long long Box2Morton(const long long x, const long long y) const {

            return Box2Morton(x, y, level_);

        }

        static uint64_t splitBy2(const unsigned int a, const int level) {

            uint64_t x = a;
            x &= (static_cast<uint64_t>(1) << (level+1)) -1;
            x = (x | (x << 16)) & 0xFFFF0000FFFF;
            x = (x | (x << 8)) & 0xFF00FF00FF00FF;
            x = (x | (x << 4)) & 0xF0F0F0F0F0F0F0F;
            x = (x | (x << 2)) & 0x3333333333333333;
            x = (x | (x << 1)) & 0x5555555555555555;
            
            return x;

        }

        inline uint64_t splitBy2(const unsigned int a) const {

            return splitBy2(a, level_);

        } 

        static unsigned int mergeFrom2(uint64_t x) {

            x &= 0x5555555555555555;
            x = (x ^ (x >> 1)) & 0x3333333333333333;
            x = (x ^ (x >> 2)) & 0xF0F0F0F0F0F0F0F;
            x = (x ^ (x >> 4)) & 0xFF00FF00FF00FF;
            x = (x ^ (x >> 8)) & 0xFFFF0000FFFF;
            x = (x ^ (x >> 16)) & 0xffffffff;

            return x;

        }

        static void Morton2Box(long long morton_box, long long& x, long long& y, const int level) {

            if (morton_box >= 1 << (2*level))
                throw std::invalid_argument("Current level cannot have this morton order.");

            x = mergeFrom2(morton_box);
            y = mergeFrom2(morton_box >> 1);

        }

        void Morton2Box(long long morton_box, long long& x, long long& y) const {

            Morton2Box(morton_box, x, y, level_);

        }

        static long long even (const long long x) {

            return static_cast<long long>(x/2)*2;

        }

        bool IsRelevant(long long morton_box) const {

            return mortonidofrelboxes_.count(morton_box) > 0;

        }

        std::vector<long long> GetCousins(const long long morton_box) const {

            std::vector<long long> cousins;
            cousins.reserve(189);

            long long x, y;
            Morton2Box(morton_box, x, y);
            long long xlo = std::max<long long>(even(x)-2, 0);
            long long xup = std::min<long long>(nboxesperside_-1, even(x)+3);
            long long ylo = std::max<long long>(even(y)-2, 0);
            long long yup = std::min<long long>(nboxesperside_-1, even(y)+3);
            
            for (int j = ylo; j <= yup; j++) {
                for (int i = xlo; i <= xup; i++) {
                    if (std::abs(i - x) > 1 || 
                        std::abs(j - y) > 1) {

                            const long long cousin_morton = Box2Morton(i, j);
                            if (IsRelevant(cousin_morton))
                                cousins.push_back(cousin_morton);

                    }
                }
            }

            std::sort(cousins.begin(), cousins.end());
            return cousins;

        }

        void GetBoxCenter(long long morton_box, double& x, double& y) const {

            long long posx, posy;
            Morton2Box(morton_box, posx, posy);
            x = min_x_ + posx * boxsize_ + 0.5 * boxsize_;
            y = min_y_ + posy * boxsize_ + 0.5 * boxsize_;

        }

        long long GetNRadialIntervals() const  {return radial_intervals_cone_.size() - 1;}

        long long Twodimensionalconesegment2nonrelconesegment(long long r, long long theta) const {

            return theta * GetNRadialIntervals() + r;

        }

        long long Cart2Nonrelcone(const long long mortonbox, double x, double y) const {

            double centerx, centery;
            GetBoxCenter(mortonbox, centerx, centery);
            x -= centerx;
            y -= centery;
            Functions::CartToPol(x, y);
            const long long posr = static_cast<long long>(std::upper_bound(radial_intervals_cone_.begin(), radial_intervals_cone_.end(), x)
            - radial_intervals_cone_.begin()) - 1;
            const long long postheta = static_cast<long long>(y * nconeselevation_ / M_PI);
            return Twodimensionalconesegment2nonrelconesegment(posr, postheta);

        }

        long long MortonboxNonrelconesegment2Relconesegment(long long morton_box, long long nonrelconesegment) const {

            return mortonboxnonrelcone2relcone_.at(morton_box).at(nonrelconesegment);
            
        }

        long long Nonrelconesegment2radialpos(const long long nonrelconesegment) const {

            return nonrelconesegment % GetNRadialIntervals();

        }

        void Nonrelconesegment2twodimensionalconesegment(const long long nonrelconesegment, long long & r, long long & theta) const {

            r = Nonrelconesegment2radialpos(nonrelconesegment);

            theta = nonrelconesegment / GetNRadialIntervals();

        }

        double ChebPosr2Radius(long long nonrelposinrdirection, const double r) const {

            double min = radial_intervals_cone_[nonrelposinrdirection];

            double reciprocal_min = scaling_r2s_ / radial_intervals_cone_[nonrelposinrdirection+1];
            double dr = scaling_r2s_/min-reciprocal_min;

            return scaling_r2s_/((r + 1.0)/2.0*dr + reciprocal_min);

        }

        void Cheb2Cart(const long long mortonbox, long long nonrelconesegment, double & x, double & y) const {

            long long posr, postheta;
            
            Nonrelconesegment2twodimensionalconesegment(nonrelconesegment, posr, postheta);

            const double dangle = M_PI/nconeselevation_;

            x = ChebPosr2Radius(posr, x);
            y = (y + 1.0)/2.0*dangle + dangle * postheta;

            double centerx, centery;
            GetBoxCenter(mortonbox, centerx, centery);
            Functions::PolToCart(x, y);
            x += centerx;
            y += centery;

        }

        std::vector<long long> GetNeighbours(long long morton_box) {

            std::vector<long long> neighbours;
            neighbours.reserve(9);
            long long x, y;
            Morton2Box(morton_box, x, y);
            long long xlo = std::max<long long>(0, x-1);
            long long xup = std::min<long long>(nboxesperside_-1, x+1);
            long long ylo = std::max<long long>(0, y-1);
            long long yup = std::min<long long>(nboxesperside_-1, y+1);

            for (int j = ylo; j <= yup; j++) {
                for (int i = xlo; i <= xup; i++) {

                    const long long neighbour_morton = Box2Morton(i, j);

                    if (IsRelevant(neighbour_morton))
                        neighbours.push_back(neighbour_morton);

                }
            }

            std::sort(neighbours.begin(), neighbours.end());
            return neighbours;

        }

        double Cheb2Radius(long long nonrelconesegment, const double r) const {

            const long long nonrelposindirection = Nonrelconesegment2radialpos(nonrelconesegment);
            return ChebPosr2Radius(nonrelposindirection, r);

        }

        double PolRadius2ChebRadius(long long nonrelconesegment_rpos, const double r) const {

            double min = radial_intervals_cone_[nonrelconesegment_rpos];

            double reciprocal_min = scaling_r2s_/radial_intervals_cone_[nonrelconesegment_rpos+1];
            double dr = scaling_r2s_/min-reciprocal_min;

            return (scaling_r2s_/r - reciprocal_min)/dr*2.0 - 1.0;
            
        }

        void Cart2Cheb(const long long mortonbox, const long long nonrelconesegment, double& x, double& y) const {

            double centerx, centery;
            GetBoxCenter(mortonbox, centerx, centery);
            x -= centerx;
            y -= centery;
            Functions::CartToPol(x, y);
            long long posr, postheta;
            Nonrelconesegment2twodimensionalconesegment(nonrelconesegment, posr, postheta);
            const double dangle = M_PI/nconeselevation_;
            x = PolRadius2ChebRadius(posr, x);
            y = (y - dangle*postheta)/dangle*2.0 - 1.0;

        }

        void Cart2Cheb(const long long mortonbox, double& x, double& y, long long& nonrelconesegment_new) const {

            double centerx, centery;
            GetBoxCenter(mortonbox, centerx, centery);
            x -= centerx;
            y -= centery;
            Functions::CartToPol(x, y);
            const long long posr = static_cast<long long>(std::upper_bound(radial_intervals_cone_.begin(), radial_intervals_cone_.end(), x) - radial_intervals_cone_.begin()) - 1;
            const long long postheta = static_cast<long long>(y * nconeselevation_ / M_PI);
            nonrelconesegment_new = Twodimensionalconesegment2nonrelconesegment(posr, postheta);
            const double dangle = M_PI/nconeselevation_;
            x = PolRadius2ChebRadius(posr, x);
            y = (y - dangle*postheta)/dangle*2.0 - 1.0;

        }

};

#endif