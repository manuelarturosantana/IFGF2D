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

    // This means that BoxTree can access the private and protected members of Level
    friend class BoxTree;

    private:

        // An enum user defined type whose possible values are a set of named constants.
        enum ConeRefinementStrategy {

            Laplace = 0,
            Helmholtz = 1

        };

        const int level_;
        const int nlevels_;
        const double min_x_; // The x value of the lower left corner of the level d = 1 box.
        const double min_y_; // The y value of the lower left corner of the level d = 1 box.
        const double boxsize_;
        const long long nboxesperside_;

        std::vector<long long> relconesegment2cocenteredmortonboxid_;
        // I believe that this also tracks the number relevant cones.
        std::vector<long long> relconesegment2nonrelconesegment_;

        // little h in the paper. s = h / r
        double scaling_r2s_;
        // Number of cones in the theta direction.
        long long nconeselevation_;
        // Interval endpoints in s of the radial interval.
        std::vector<double> radial_intervals_cone_;

        // morton ID of relevant boxes
        /*
        * unordered_set does not store values, only if an element exits
        * unordered_map is a hash table, stores keys and values.
        */
        std::unordered_set<long long> mortonidofrelboxes_;
        // There is an assumption here that the points are orderd in such a way that
        // all the points in a box go in ascending order (done in teh sort box function). 
        // Each element holds the start index of the points in the morton box, as well as
        // the number of points in that box.
        std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_;

        std::unordered_set<long long> propagationrequiredrelconesegments_;
        std::unordered_set<long long> interpolationrequiredrelconesegments_;

        // This tracks given the morton box, and non-relvant cone the ID of the relevant cone.
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

            // The closest radial point in the boxes in the paper
            const double minsize = (3.0/2.0 * boxsize_) * (1.0-1e-10); 
            // The distance of anypoint or box size in the paper.
            // Level comes in from nboxesperside
            const double maxsize = (std::sqrt(2.0) * boxsize_ * nboxesperside_ - std::sqrt(2.0)/2.0 * boxsize_) * (1.0+1e-10);
            constexpr long long nconesradial = 1;
            constexpr long long nconeselevation = 2; // Twice as many in the theta direction
            long long nconeslevelscaling = 1;

            if (S == Laplace) {

            } else if (S == Helmholtz) {

                if (level_ > 1)
                    // TODO: This looks like a hueristic to determine how many initial cones, but I'm not sure where it cones from.
                    // Note M_1_PI is 1 / PI
                    nconeslevelscaling = std::max<long long>(1, std::ceil(wavenumber * boxsize_ * M_1_PI * 1.2)); //?

            } else {
                throw std::invalid_argument("Unknown cone refinement strategy.");
            }

            nconeselevation_ = nconeslevelscaling * nconeselevation;
            const long long nconesradialonlevel = nconeslevelscaling * nconesradial;

            // Scaling for s change of variables.
            scaling_r2s_ = std::sqrt(2.0)/2.0 * boxsize_;
            // For all source points x' and target points x we have |x'| / |x| < eta
            const double eta = std::sqrt(2.0)/3.0-1e-5;

            const double tmpsizeinr = (maxsize - minsize + 1e-10);
            // S is interpolated from [0, eta].
            // The second term is 1/minsize - 1/maxsize, the distance in the s interval.
            // So this says use the smaller of the two, the whole size, or an interval in eta
            const double tmpsizeins = std::min(eta/nconesradialonlevel, scaling_r2s_*(1.0/minsize - 1.0/(minsize+tmpsizeinr)));
            // the s distance divided by the size of each s interval.
            const long long nactualcones = std::ceil(scaling_r2s_*(1.0/minsize - 1.0/maxsize)/tmpsizeins);
            
            radial_intervals_cone_.resize(nactualcones+1);
            radial_intervals_cone_[0] = minsize;
    
            const double finalsizeins = (scaling_r2s_*(1.0/minsize - 1.0/maxsize)/nactualcones);
            for (long long i = nactualcones-1; i >= 0; i--) {
                const long long iter = nactualcones - i;
                // Like 27 in the IFGF paper but I don't quite understand this part :(
                // Appears like this is stored in the r variable.
                radial_intervals_cone_[iter] = scaling_r2s_*maxsize/(scaling_r2s_ + maxsize*i*finalsizeins);
            }

        }

        // Returns the morton index of the box which the point (x,y) is in
        static long long Point2Morton (const double x, const double y,
                                       const double min_x, const double min_y,
                                       const double boxsize, const int level) {

            const long long posx = static_cast<long long>((x - min_x)/boxsize);
            const long long posy = static_cast<long long>((y - min_y)/boxsize);

            if (posx < 0 || posy < 0 || posx >= 1 << level || posy >= 1 << level)
                throw std::logic_error("The 2D box index cannot be < 0 or larger than the number of boxes");

            return Box2Morton(posx, posy, level);

        }

        // Returns the morton index of the box which the point (x,y) is in
        long long  Point2Morton(const double x, const double y) const {

            return Point2Morton(x, y, min_x_, min_y_, boxsize_, level_);

        }

        // Returns the index of the morton order of the box given the x and y cartesian
        // coordinate of the box in the grid, starting from the lower left corner.
        static long long Box2Morton(const long long x, const long long y, const int level) {

            uint64_t answer = 0;
            // x and y get split into bits with zeros between them, and then
            // get interleaved.
            answer |= splitBy2(x, level) | splitBy2(y, level) << 1;

            return answer;

        }

        long long Box2Morton(const long long x, const long long y) const {

            return Box2Morton(x, y, level_);

        }

        // 
        static uint64_t splitBy2(const unsigned int a, const int level) {

            uint64_t x = a;
            // The rightside in 2^(level+1) - 1, which means that in binary it if level was 3 we would have
            // the right side evaluating to (in binary) 01111. Then &= grabs means that
            // x now lowest level + 1 bits of a.
            x &= (static_cast<uint64_t>(1) << (level+1)) -1;
            // The rest of this code spreads zeros between the bits of x
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

        // Grab every other bit from x
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
                throw std::invalid_argument(" Morton2Box: Current level cannot have this morton order.");

            x = mergeFrom2(morton_box);
            y = mergeFrom2(morton_box >> 1);

        }

        // Given a value in the morton_ordering, return the x and y coordinates of the 
        // relevant box
        void Morton2Box(long long morton_box, long long& x, long long& y) const {

            Morton2Box(morton_box, x, y, level_);

        }

        // Return closest even number <= x
        static long long even (const long long x) {

            return static_cast<long long>(x/2)*2;

        }

        bool IsRelevant(long long morton_box) const {

            return mortonidofrelboxes_.count(morton_box) > 0;

        }

        // Should be called get Relevant cousins. Doesn't return any other.
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

        // The non-relevant cones repsent the relative location of the cone to the center of the box
        // in a certain numbering
        long long Cart2Nonrelcone(const long long mortonbox, double x, double y) const {

            double centerx, centery;
            GetBoxCenter(mortonbox, centerx, centery);
            x -= centerx;
            y -= centery;
            Functions::CartToPol(x, y);
              // upper_bound returns an iterator pointing to the first element greater than the value being searched for.
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
            // Max and min make sure we stay in bounds
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