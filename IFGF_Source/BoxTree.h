#ifndef BOXTREE_H
#define BOXTREE_H

#include <vector>
#include <limits>
#include <map>
#include <set>
#include <iostream>
#include <numeric>

#include "Level.h"
#include "Vector.h"
#include "Interpolator.h"

#include "../utils/DebugUtils.hpp"

class Level;
class BoxTree;

const Interpolator<Chebyshev> IPSCHEME_;
// TODO: This is a terrible hack :(..
// Should use a template
// const int P_ = 3*5;
const int P_ = 3*5;



class BoxTree 
{

    private:

        std::vector<double> x_; // x coordinates of all values to be evaluated in the IFGF, sorted in Morton order
        std::vector<double> y_; // y coordinates of all values to be evaluated in the IFGF, sorted in Morton order

        std::vector<long long> sorting_; // List as map to go from oringinal ordering to morton box ordering
        int nlevels_;
        double wavenumber_;
        long long N_; // Total number of points in IFGF.

        // int P_;

        std::vector<Level*> levels_;

        std::vector<double> solution_real_;
        std::vector<double> solution_imag_;

        std::vector<double> conesegments_current_real_;
        std::vector<double> conesegments_current_imag_;
        std::vector<double> conesegments_prev_real_;
        std::vector<double> conesegments_prev_imag_;

    private:

        void Initialize() 
        {

            InitializeBoxesAndLevels();

            InitializeRelevantConeSegments();

            InitializeSolutionAndConeSegmentCoefficients();

        }

        // Sort points according to the Morton ordering
        void SortBox(const std::array<double, 2>& min, const double boxsize) 
        {

            // std::iota(sorting_.begin(), sorting_.end(), 0);
            sorting_.reserve(N_);
            for (long long i = 0; i < N_; i++) {

                sorting_.push_back(i);

            }


            std::vector<long long> morton_code(N_);

            for (long long i = 0; i < N_; i++) {
                // Point2Morton Returns the morton index of the box which the point (x,y) is in
                // So points are sorted into their boxes.
                const long long morton = Level::Point2Morton(x_[i], y_[i], min[0], min[1], boxsize, nlevels_-1);
                morton_code[i] = morton;

            }

            new_sort(sorting_, morton_code);

            
            std::vector<double> tmp_xy(N_);

            for (long long i = 0; i < N_; i++) {

                tmp_xy[i] = x_[sorting_[i]];

            }

            x_ = tmp_xy;

            for (long long i = 0; i < N_; i++) {

                tmp_xy[i] = y_[sorting_[i]];

            }

            y_ = tmp_xy;

        }

        // Compute the start index of the points in each box and the number of points in 
        // each box. Assumes that the points are sorted into their Morton boxes
        // which is done in the sort box function
        void InitializeLevelDBoxesData() 
        {

            long long old_morton = -1;
            long long npoints;

            for (long long i = 0; i < N_; i++) {

                const long long morton = levels_.back()->Point2Morton(x_[i], y_[i]);

                if (morton != old_morton) {

                    // Records the index of the previous points
                    // TODO: Does this work if all the points in the box are not orderd 
                    // right next to each other. Could test by shuffling the points.
                    if (i != 0) {

                        levels_.back()->mortonbox2discretizationpoints_[old_morton] = {i - npoints, npoints};

                    }

                    old_morton = morton;
                    npoints = 1;
                    levels_.back()->mortonidofrelboxes_.insert(morton);
                    
                } else {

                    npoints++;

                }

            }

            levels_.back()->mortonbox2discretizationpoints_[old_morton] = {N_ - npoints, npoints};

        }

        // In each level flags which boxes are relevant
        void InitializeRelBoxesAllLevels() 
        {

            std::unordered_set<long long> relevantparentmorton;
            std::unordered_set<long long> relevantmorton = levels_.back()->mortonidofrelboxes_;

            for (int i = nlevels_-1; i >= 2; i--) {

                // Store the relevant boxes for this level
                if (i != nlevels_-1)
                    levels_[i]->mortonidofrelboxes_.insert(relevantmorton.begin(), relevantmorton.end());

                for (auto j = relevantmorton.begin(); j != relevantmorton.end(); j++) {
                    // j is an iterator. *j give the value that j points. Dividing by 4 = 2^2
                    // moves you one level up the tree since you are moving 2 bits over, one level
                    // of the tree
                    relevantparentmorton.insert(static_cast<long long>(*j/4));

                }

                std::swap(relevantmorton, relevantparentmorton);
                relevantparentmorton.clear();

            }

        }

        // I believe that this function is compute bounding box
        // Computes the values min and max such that
        // min (resp max) contains values below (resp. above) the smallest x and y value.
        void ComputeBB(const std::vector<double>& pointsx, const std::vector<double>& pointsy,
                    std::array<double, 2>& min, std::array<double ,2>& max) const
        {

            min = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
            max = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};

            for (long long i = 0; i < N_; i++) {

                if (pointsx[i] < min[0]) min[0] = pointsx[i];
                if (pointsy[i] < min[1]) min[1] = pointsy[i];
                
                if (pointsx[i] > max[0]) max[0] = pointsx[i];
                if (pointsy[i] > max[1]) max[1] = pointsy[i];
                
            }

            max[0] += 1e-10;
            max[1] += 1e-10;

            min[0] -= 1e-10;
            min[1] -= 1e-10;
            
        }

        // Flags relevant boxes on all levels, creates the first bounding box
        void InitializeBoxesAndLevels() 
        {

            std::array<double, 2> min, max;

            ComputeBB(x_, y_, min, max);          

            // Choose the initial box to have side lengths slightly larger than the distance
            // of the x and y coordinates  or the points furthest from the origin
            double boxsize = std::max(max[0] - min[0], max[1] - min[1]);

            // Sort points according to morton ordering at the lowest level.
            SortBox(min, boxsize / (1 << (nlevels_ - 1)));

            levels_.resize(nlevels_, nullptr);

            for (int i = 1; i < nlevels_; i++) {

                levels_[i] = new Level(i, nlevels_, min[0], min[1], boxsize/(1 << i), wavenumber_);

            }

            InitializeLevelDBoxesData();

            InitializeRelBoxesAllLevels();

        }

        // Returns an ordered map where if it is given the cousin box as a key, it returns a vector containing
        // all cone segements relevant to that box
        void GetRelevantConeSegmentsDueToCousinSurface(const int level, std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment) const {

            long long old_morton_box = -1;
            std::vector<long long> cousins;
            std::unordered_map<long long, std::unordered_set<long long>> tmpmortonboxnonrelconesegment; 

            for (long long i = 0; i < N_; i++) {

                const double x = x_[i];
                const double y = y_[i];

                const long long morton_box = levels_[level]->Point2Morton(x, y);

                if (morton_box != old_morton_box) {
                    // Get Cousins only returns the relevant cousins
                    cousins = levels_[level]->GetCousins(morton_box);
                    old_morton_box = morton_box;

                }

                // For each cousin, compute the cone segments that matter to it.
                for (size_t j = 0; j < cousins.size(); j++) {
                    const long long morton_cousin = cousins[j];
                    const long long nonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_cousin, x, y);
                    // store the cone segement in the unordered set relating to the cousin.
                    // Note if the cone segement is already flagged nothing happens when
                    // inserting it again here.
                    tmpmortonboxnonrelconesegment[morton_cousin].insert(nonrel_conesegment);

                }

            }

            for (auto i = tmpmortonboxnonrelconesegment.begin(); i != tmpmortonboxnonrelconesegment.end(); i++) {
                // second is the value of the iterator i, which is itself and iterator
                std::vector<long long> tmp(i->second.begin(), i->second.end());
                std::sort(tmp.begin(), tmp.end());
                // first is the key of the iterator
                mortonboxnonrelconesegment[i->first] = tmp;

            }

        }

        std::vector<long long> GetMortonRelChildren(const int level, const long long morton_box) const {

            std::vector<long long> children;

            if (level >= nlevels_ - 1)
                return children;

            children.reserve(4);

            for (long long morton_children = morton_box*4; morton_children < (morton_box+1)*4; morton_children++) {

                if (levels_[level+1]->IsRelevant(morton_children))
                    children.push_back(morton_children);

            }

            return children;

        }

        void GetRelevantConeSegmentsDueToParents(const int level, std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment) {

            std::unordered_map<long long, std::unordered_set<long long>> mortonboxnonrelconesegment_all;

            double x, y;
            long long old_cocentered_morton_box = -1;
            std::vector<long long> morton_children;
            std::unordered_map<long long, std::array<std::unordered_map<long long, std::vector<int>>, 4>> protobox;

            for (size_t iter = 0; iter < levels_[level-1]->relconesegment2nonrelconesegment_.size(); iter++) {

                long long rel_conesegment = iter; 

                const long long cocentered_morton_box = levels_[level-1]->relconesegment2cocenteredmortonboxid_[rel_conesegment];
                const long long nonrel_conesegment = levels_[level-1]->relconesegment2nonrelconesegment_[rel_conesegment];

                if (cocentered_morton_box != old_cocentered_morton_box) {

                    morton_children = GetMortonRelChildren(level-1, cocentered_morton_box);

                    if (morton_children.size() == 0) {

                        throw std::logic_error("Every cocentered morton box must have children in get relevant cone segments function.");

                    }

                    old_cocentered_morton_box = cocentered_morton_box;

                }

                for (size_t childiter = 0; childiter < morton_children.size(); childiter++) {

                    const long long morton_child = morton_children[childiter];
                    const long long child_pos = morton_child % 4;

                    if (protobox.count(nonrel_conesegment) == 0 || protobox.at(nonrel_conesegment)[child_pos].size() == 0) {

                        std::unordered_map<long long, std::vector<int>> & tmpchildproto = protobox[nonrel_conesegment][child_pos];

                        for (int interpiter = 0; interpiter < IPSCHEME_.GetNInterpPoints(); interpiter++) {

                            IPSCHEME_.GetInterpolationPoint(interpiter, x, y);

                            levels_[level-1]->Cheb2Cart(cocentered_morton_box, nonrel_conesegment, x, y);

                            long long locnonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_child, x, y);

                            tmpchildproto[locnonrel_conesegment].push_back(interpiter);

                        }

                    }

                }

                for (auto iter = protobox.begin(); iter != protobox.end(); iter++) {

                    if (levels_[level]->protobox_.count(iter->first) == 0) {

                        levels_[level]->protobox_[iter->first] = iter->second;

                    } else {

                        for (int childiter = 0; childiter < 4; childiter++) {

                            if (iter->second[childiter].size() != 0 && levels_[level]->protobox_[iter->first][childiter].size() == 0) {

                                levels_[level]->protobox_[iter->first][childiter] = iter->second[childiter];

                            }

                        }

                    }

                }

                std::unordered_map<long long, std::array<std::unordered_map<long long, std::vector<int>>, 4>>().swap(protobox);

            }

            old_cocentered_morton_box = -1;

            for (size_t iter = 0; iter < levels_[level-1]->relconesegment2nonrelconesegment_.size(); iter++) {

                long long rel_conesegment = iter; 

                const long long cocentered_morton_box = levels_[level-1]->relconesegment2cocenteredmortonboxid_[rel_conesegment];
                const long long nonrel_conesegment = levels_[level-1]->relconesegment2nonrelconesegment_[rel_conesegment];

                if (cocentered_morton_box != old_cocentered_morton_box) {

                    morton_children = GetMortonRelChildren(level-1, cocentered_morton_box);

                    if (morton_children.size() == 0) {

                        throw std::logic_error("Every cocentered morton box must have children in get relevant cone segments function.");

                    }

                    old_cocentered_morton_box = cocentered_morton_box;

                }

                for (size_t childiter = 0; childiter < morton_children.size(); childiter++) {

                    const long long morton_child = morton_children[childiter];
                    const long long child_pos = morton_child % 4;

                    const std::unordered_map<long long, std::vector<int>> & allnonrelconesegments = levels_[level]->protobox_.at(nonrel_conesegment)[child_pos];

                    for (auto i = allnonrelconesegments.begin(); i != allnonrelconesegments.end(); i++) {

                        mortonboxnonrelconesegment_all[morton_child].insert(i->first);

                    }

                }

            }

            for (auto iter = mortonboxnonrelconesegment_all.begin(); iter != mortonboxnonrelconesegment_all.end(); iter++) {

                std::vector<long long> tmp(iter->second.begin(), iter->second.end());
                std::sort(tmp.begin(), tmp.end());
                mortonboxnonrelconesegment[iter->first] = std::move(tmp);

            }

        }

        std::unordered_map<long long, std::vector<long long>> MergeRelevantConeSegments(
            const std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment1,
            const std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment2) const {

            std::unordered_map<long long, std::vector<long long>> merged_map = mortonboxnonrelconesegment1;

            for (auto iter = mortonboxnonrelconesegment2.begin(); iter != mortonboxnonrelconesegment2.end(); iter++) {

                if (merged_map.count(iter->first) == 0) 
                    merged_map.insert({iter->first, std::vector<long long>{}});

            }

            for (auto iter = mortonboxnonrelconesegment2.begin(); iter != mortonboxnonrelconesegment2.end(); iter++) {

                std::vector<long long> & vec1 = merged_map.at(iter->first);

                std::vector<long long> vec;
                vec.reserve(iter->second.size() + vec1.size());
                std::merge(vec1.begin(), vec1.end(), iter->second.begin(), iter->second.end(), std::back_inserter(vec));

                auto last = std::unique(std::begin(vec), std::end(vec));
                vec.erase(last, std::end(vec));
                vec1 = std::move(vec);

            }

            return merged_map;

        }

        void LoadBalanceRelevantConeSegments(const int level, const std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment) {

            long long nallconesegments = 0;

            std::vector<long long> relcone2nonrelcone;
            std::vector<long long> relcone2cocenteredmorton;

            for (auto iter = mortonboxnonrelconesegment.begin(); iter != mortonboxnonrelconesegment.end(); iter++) {

                nallconesegments += iter->second.size();

            }

            relcone2nonrelcone.resize(nallconesegments);
            relcone2cocenteredmorton.resize(nallconesegments);

            long long counter = 0;

            for (auto iter = mortonboxnonrelconesegment.begin(); iter != mortonboxnonrelconesegment.end(); iter++) {
                for (auto coneiter = iter->second.begin(); coneiter != iter->second.end(); coneiter++) {

                    relcone2nonrelcone[counter] = *coneiter;
                    relcone2cocenteredmorton[counter] = iter->first;
                    counter++;

                }
            }

            levels_[level]->relconesegment2nonrelconesegment_ = relcone2nonrelcone;
            levels_[level]->relconesegment2cocenteredmortonboxid_ = relcone2cocenteredmorton;

            std::vector<long long>().swap(relcone2nonrelcone);
            std::vector<long long>().swap(relcone2cocenteredmorton);

            for (size_t iter = 0; iter < levels_[level]->relconesegment2nonrelconesegment_.size(); iter++) {

                const long long cocenteredmortonbox = levels_[level]->relconesegment2cocenteredmortonboxid_[iter];
                const long long nonrelconeid = levels_[level]->relconesegment2nonrelconesegment_[iter];
                levels_[level]->mortonboxnonrelcone2relcone_[cocenteredmortonbox][nonrelconeid] = iter;

            }

            for (auto iter = mortonboxnonrelconesegment.begin(); iter != mortonboxnonrelconesegment.end(); iter++) {
                for (auto coneiter = iter->second.begin(); coneiter != iter->second.end(); coneiter++) {

                    if (!(levels_[level]->mortonboxnonrelcone2relcone_.count(iter->first) > 0 && 
                        levels_[level]->mortonboxnonrelcone2relcone_.at(iter->first).count(*coneiter) > 0)) {
                            std::cout << "WARNING\n";
                        }

                }
            }

        }

        void SetPropagationAndInterpolationRequiredConeSegments(const int level, 
            const std::unordered_map<long long, std::vector<long long>> & interpolation_conesegments,
            const std::unordered_map<long long, std::vector<long long>> & propagation_conesegments) {

            levels_[level]->interpolationrequiredrelconesegments_.clear();
            levels_[level]->propagationrequiredrelconesegments_.clear();

            for (auto iter = interpolation_conesegments.begin(); iter != interpolation_conesegments.end(); iter++) {

                const std::vector<long long> nonrelconesegments = iter->second;

                for (size_t csiter = 0; csiter < nonrelconesegments.size(); csiter++) {

                    const long long relconesegment = levels_[level]->MortonboxNonrelconesegment2Relconesegment(iter->first, nonrelconesegments[csiter]);

                    levels_[level]->interpolationrequiredrelconesegments_.insert(relconesegment);

                }

            }

            for (auto iter = propagation_conesegments.begin(); iter != propagation_conesegments.end(); iter++) {

                const std::vector<long long> nonrelconesegments = iter->second;

                for (size_t csiter = 0; csiter < nonrelconesegments.size(); csiter++) {

                    const long long relconesegment = levels_[level]->MortonboxNonrelconesegment2Relconesegment(iter->first, nonrelconesegments[csiter]);

                    levels_[level]->propagationrequiredrelconesegments_.insert(relconesegment);

                }

            }

        }

        void InitializeRelevantConeSegments() {

            std::unordered_map<long long, std::vector<long long>> interpolation_conesegments;
            std::unordered_map<long long, std::vector<long long>> propagation_conesegments;

            for (int level = 2; level < nlevels_; level++) {

                GetRelevantConeSegmentsDueToCousinSurface(level, interpolation_conesegments);

            
                if (level > 2) {

                    GetRelevantConeSegmentsDueToParents(level, propagation_conesegments);

                }

                if (interpolation_conesegments.size() == 0 && propagation_conesegments.size() == 0) {
                    std::cout << "Level: " << level << std::endl;
                    throw std::logic_error("The number of cone segments on any level cannot be 0");
                }

                LoadBalanceRelevantConeSegments(level, MergeRelevantConeSegments(interpolation_conesegments, propagation_conesegments));

                SetPropagationAndInterpolationRequiredConeSegments(level, interpolation_conesegments, propagation_conesegments);

                interpolation_conesegments.clear();
                propagation_conesegments.clear();  

            }

        }

        void InitializeSolutionAndConeSegmentCoefficients() {

            solution_real_ = std::vector<double>(N_, 0.0);
            solution_imag_ = std::vector<double>(N_, 0.0);

            long long maxncoeffs = 0;
            for (int level = 2; level < nlevels_; level++) {

                maxncoeffs = std::max<long long>(maxncoeffs, levels_[level]->relconesegment2nonrelconesegment_.size());

            }
            maxncoeffs *= IPSCHEME_.GetNInterpPoints();

            conesegments_current_real_ = std::vector<double>(maxncoeffs, 0.0);
            conesegments_current_imag_ = std::vector<double>(maxncoeffs, 0.0);
            conesegments_prev_real_ = std::vector<double>(maxncoeffs, 0.0);
            conesegments_prev_imag_ = std::vector<double>(maxncoeffs, 0.0);
            
        }

        void ZeroSolution() 
        {

            for (long long i = 0; i < N_; i++) {

                solution_real_[i] = 0.0;
                solution_imag_[i] = 0.0;

            }

        }

        template <void _kernel(const double, const double, const double, const double, const double, const double, double&, double&)>
        void SingularInteractions(std::vector<double> & density_real, std::vector<double> & density_imag) {

            long long old_mortonbox = -1;
            std::vector<long long> neighbours;
            long long npoints = 0;
            long long points_begin = 0;

            for (long long i = 0; i < N_; i++) {

                const long long mortonbox = levels_.back()->Point2Morton(x_[i], y_[i]);
                double result_real = solution_real_[i];
                double result_imag = solution_imag_[i];
                double tmp_result_real, tmp_result_imag;
                const double x = x_[i];
                const double y = y_[i];

                if (mortonbox != old_mortonbox) {

                    neighbours = levels_.back()->GetNeighbours(mortonbox);

                    if (neighbours.size() > 9)
                        throw std::logic_error("Cannot have more than 9 neighbours.");
                    
                    old_mortonbox = mortonbox;

                }

                for (size_t neighbouriter = 0; neighbouriter < neighbours.size(); neighbouriter++) {

                    const long long neighbourmorton = neighbours[neighbouriter];

                    const std::array<long long, 2> & points_data = levels_.back()->mortonbox2discretizationpoints_.at(neighbourmorton);

                    points_begin = points_data[0];
                    npoints = points_data[1];

                    for (long long sourceiter = 0; sourceiter < npoints; sourceiter++) {

                        // TODO: ADD HERE IF NEAR SINGULAR
                        if (points_begin + sourceiter == i)
                            continue;

                        double locx = x_[points_begin + sourceiter];
                        double locy = y_[points_begin + sourceiter];

                        double locdensity_real = density_real[points_begin + sourceiter];
                        double locdensity_imag = density_imag[points_begin + sourceiter];

                        _kernel(locx, locy, x, y, locdensity_real, locdensity_imag, tmp_result_real, tmp_result_imag);

                        result_real += tmp_result_real;
                        result_imag += tmp_result_imag;

                    }

                }

                solution_real_[i] = result_real;
                solution_imag_[i] = result_imag;

            }

        }

        template <void _kernel(const double, const double, const double, const double, const double, const double, double&, double&),
                  void _factorization(const double, double&, double&)>
        void LevelDEvaluations(std::vector<double> & density_real, std::vector<double> & density_imag) {

            long long old_cocentered_morton_box = -1;
            double locdensity_real;
            double locdensity_imag;
            double locx;
            double locy;

            long long npoints;
            long long points_begin;

            for (size_t i = 0; i < levels_.back()->relconesegment2nonrelconesegment_.size(); i++) {

                const long long cocentered_morton_box = levels_.back()->relconesegment2cocenteredmortonboxid_[i];
                const long long nonrel_conesegment = levels_.back()->relconesegment2nonrelconesegment_[i];
                double x, y;
                const long long coefficients_begin = i * IPSCHEME_.GetNInterpPoints(); 
    
                if (old_cocentered_morton_box != cocentered_morton_box) {

                    const std::array<long long, 2>& points_data = levels_.back()->mortonbox2discretizationpoints_.at(cocentered_morton_box);
                    points_begin = points_data[0];
                    npoints = points_data[1];

                    old_cocentered_morton_box = cocentered_morton_box;

                }

                for (int j = 0; j < IPSCHEME_.GetNInterpPoints(); j++) {

                    IPSCHEME_.GetInterpolationPoint(j, x, y);
                    
                    const double radius = levels_.back()->Cheb2Radius(nonrel_conesegment, x);
                    
                    levels_.back()->Cheb2Cart(cocentered_morton_box, nonrel_conesegment, x, y);
                    
                    double result_real = 0;
                    double result_imag = 0;
                    double kernel_real, kernel_imag;
                    
                    // This computes what is known in the paper as F_S
                    for (long long k = 0; k < npoints; k++) {

                        locx = x_[points_begin + k];
                        locy = y_[points_begin + k];
                        locdensity_real = density_real[points_begin + k];
                        locdensity_imag = density_imag[points_begin + k];

                        _kernel(locx, locy, x, y, locdensity_real, locdensity_imag, kernel_real, kernel_imag);
                        result_real += kernel_real;
                        result_imag += kernel_imag;

                    }

                    double fac_real, fac_imag;
                    _factorization(radius, fac_real, fac_imag);
                    const double dabsfac = 1.0/(fac_real*fac_real + fac_imag*fac_imag);                    

                    const long long coefficients_pos = coefficients_begin + j;
                    // Divide by the factorization to get the slowy varying function
                    conesegments_prev_real_[coefficients_pos] = (result_real * fac_real + result_imag * fac_imag) * dabsfac;
                    conesegments_prev_imag_[coefficients_pos] = (result_imag * fac_real - result_real * fac_imag) * dabsfac;
                
                }

                // Now do the interpolation 
                IPSCHEME_.GenerateInterpolant(&conesegments_prev_real_[coefficients_begin]);
                IPSCHEME_.GenerateInterpolant(&conesegments_prev_imag_[coefficients_begin]);
            
            }            
        
            for (size_t i = 0; i < conesegments_current_real_.size(); i++) {

                conesegments_current_real_[i] = conesegments_prev_real_[i];
                conesegments_current_imag_[i] = conesegments_prev_imag_[i];

            }
        
        }

        // Having already computed the the relevant cones compute the interpolants
        // in the cones for the next level up.
        template<void _factorization(const double, double&, double&)>
        void Propagation(const int level) {

            std::array<double, P_> x, y, radii, value_real, value_imag;
            double locx, locy;
            std::array<double const *, 4> childcoeffs_real, childcoeffs_imag;
            long long old_cocentered_morton_box = -1;
            std::vector<long long> morton_children;
            std::array<long long, 4> old_locrelconesegment = {-1, -1, -1, -1};

            for (size_t i = 0; i < levels_[level-1]->relconesegment2nonrelconesegment_.size(); i++) {

                long long rel_conesegment = i; 

                const long long cocentered_morton_box = levels_[level-1]->relconesegment2cocenteredmortonboxid_[rel_conesegment];
                const long long nonrel_conesegment = levels_[level-1]->relconesegment2nonrelconesegment_[rel_conesegment];

                const long long new_coeffs_begin_id = rel_conesegment * IPSCHEME_.GetNInterpPoints(); 
                
                // TODO: Seems like a silly if statement. Won't this always be true by the previous line
                if (new_coeffs_begin_id % IPSCHEME_.GetNInterpPoints() != 0)
                    throw std::logic_error("Coeffs need to start at a multiple of interpolation points.");

                if (cocentered_morton_box != old_cocentered_morton_box) {

                    morton_children = GetMortonRelChildren(level-1, cocentered_morton_box);

                    old_cocentered_morton_box = cocentered_morton_box;

                }

                const std::array<std::unordered_map<long long, std::vector<int>>, 4>& allchildnonrelcones = levels_[level]->protobox_.at(nonrel_conesegment);

                for (int interpid = 0; interpid < IPSCHEME_.GetNInterpPoints(); interpid++) {

                    IPSCHEME_.GetInterpolationPoint(interpid, x[interpid], y[interpid]);
                    radii[interpid] = levels_[level-1]->Cheb2Radius(nonrel_conesegment, x[interpid]);

                    levels_[level-1]->Cheb2Cart(cocentered_morton_box, nonrel_conesegment, x[interpid], y[interpid]);

                    value_real[interpid] = 0.0;
                    value_imag[interpid] = 0.0;

                }

                // For each child of a parent evaluate at new interpolation points,
                for (size_t childiter = 0; childiter < morton_children.size(); childiter++) {

                    const long long morton_child = morton_children[childiter];
                    // Order the four children
                    const long long child_pos = morton_child % 4;
                    const std::unordered_map<long long, std::vector<int>>& childnonrelcones = allchildnonrelcones.at(child_pos);
                    
                    for (auto childconeiter = childnonrelcones.begin(); childconeiter != childnonrelcones.end(); childconeiter++) {
                        // TODO: Does the vector in child cone iter holds points to interpolate at
                        for (size_t interpiter = 0; interpiter < childconeiter->second.size(); interpiter++) {

                            const int interpid = childconeiter->second[interpiter];
                            locx = x[interpid];
                            locy = y[interpid];
                            levels_[level]->Cart2Cheb(morton_child, childconeiter->first, locx, locy);
                            
                            const long long locrelconesegment = levels_[level]->mortonboxnonrelcone2relcone_.at(morton_child).at(childconeiter->first);
                            const double locradius = levels_[level]->Cheb2Radius(childconeiter->first, locx);
                            double locfac_real, locfac_imag;
                            _factorization(locradius, locfac_real, locfac_imag);

                            if (locrelconesegment != old_locrelconesegment[childiter]) {

                                const long long coeffs_begin_id = locrelconesegment * IPSCHEME_.GetNInterpPoints();
                                childcoeffs_real[childiter] = &conesegments_prev_real_[coeffs_begin_id];
                                childcoeffs_imag[childiter] = &conesegments_prev_imag_[coeffs_begin_id];

                                old_locrelconesegment[childiter] = locrelconesegment;

                            }

                            const double tmpreal = IPSCHEME_.Interpolate(locx, locy, childcoeffs_real[childiter]);
                            const double tmpimag = IPSCHEME_.Interpolate(locx, locy, childcoeffs_imag[childiter]);
                            
                            value_real[interpid] += tmpreal * locfac_real - tmpimag * locfac_imag;
                            value_imag[interpid] += tmpimag * locfac_real + tmpreal * locfac_imag;
                            
                        }

                    }

                }

                for (int interpid = 0; interpid < IPSCHEME_.GetNInterpPoints(); interpid++) {

                    double fac_real, fac_imag;
                    _factorization(radii[interpid], fac_real, fac_imag);
                    const double dabsfac = 1.0/(fac_real * fac_real + fac_imag * fac_imag);

                    conesegments_current_real_[new_coeffs_begin_id+interpid] = (value_real[interpid] * fac_real + value_imag[interpid] * fac_imag) * dabsfac;
                    conesegments_current_imag_[new_coeffs_begin_id+interpid] = (value_imag[interpid] * fac_real - value_real[interpid] * fac_imag) * dabsfac;
                    
                }

                IPSCHEME_.GenerateInterpolant(&conesegments_current_real_[new_coeffs_begin_id]);
                IPSCHEME_.GenerateInterpolant(&conesegments_current_imag_[new_coeffs_begin_id]);

            }

        }


        // Evaluate the interpolants. For every point aggregate all the
        // contributions from all boxes for which x lives in a cousin box from.
        // We loop over points, and then boxes for which the box x lives in is a 
        // cousin box to.
        template<void _factorization(const double, double&, double&)>
        void Interpolation(const int level) {

            double locx, locy;
            long long old_morton_box = -1;
            std::array<long long, 189> old_relconesegment; old_relconesegment.fill(-1);
            std::array<double const *, 189> cousincoeffs_real, cousincoeffs_imag;
            std::vector<long long> morton_cousinboxes;

            for (long long pointiter = 0; pointiter < N_; pointiter++) {

                const double x = x_[pointiter];
                const double y = y_[pointiter];

                double value_real = solution_real_[pointiter];
                double value_imag = solution_imag_[pointiter];

                const long long morton_box = levels_[level]->Point2Morton(x, y);

                if (morton_box != old_morton_box) {

                    morton_cousinboxes = levels_[level]->GetCousins(morton_box);
                    old_morton_box = morton_box;

                }

                for (size_t cousiniter = 0; cousiniter < morton_cousinboxes.size(); cousiniter++) {

                    locx = x;
                    locy = y;

                    const long long morton_cousinbox = morton_cousinboxes[cousiniter];
                    
                    long long relconesegment = -1;
                    long long nonrelconesegment;
                    levels_[level]->Cart2Cheb(morton_cousinbox, locx, locy, nonrelconesegment);
                    relconesegment = levels_[level]->mortonboxnonrelcone2relcone_.at(morton_cousinbox).at(nonrelconesegment);
                    const double radius = levels_[level]->Cheb2Radius(nonrelconesegment, locx);
                    double fac_real, fac_imag;
                    _factorization(radius, fac_real, fac_imag);
                    
                    if (relconesegment != old_relconesegment[cousiniter]) {

                        const long long coeffs_begin_id = relconesegment * IPSCHEME_.GetNInterpPoints(); 
                        cousincoeffs_real[cousiniter] = &conesegments_prev_real_[coeffs_begin_id];
                        cousincoeffs_imag[cousiniter] = &conesegments_prev_imag_[coeffs_begin_id];

                        old_relconesegment[cousiniter] = relconesegment;

                    }
                    
                    const double tmpreal = IPSCHEME_.Interpolate(locx, locy, cousincoeffs_real[cousiniter]);
                    const double tmpimag = IPSCHEME_.Interpolate(locx, locy, cousincoeffs_imag[cousiniter]);

                    // Complex multiplication
                    value_real += tmpreal * fac_real - tmpimag * fac_imag;
                    value_imag += tmpreal * fac_imag + tmpimag * fac_real;

                }

                solution_real_[pointiter] = value_real;
                solution_imag_[pointiter] = value_imag;

            }

        } 

        void SwapCoefficient() {

            for (size_t i = 0; i < conesegments_current_real_.size(); i++) {

                double tmp = conesegments_current_real_[i];
                conesegments_current_real_[i] = conesegments_prev_real_[i];
                conesegments_prev_real_[i] = tmp;

                tmp = conesegments_current_imag_[i];
                conesegments_current_imag_[i] = conesegments_prev_imag_[i];
                conesegments_prev_imag_[i] = tmp;

            }

        }

        void SortDensity(std::vector<double>& density_real, std::vector<double>& density_imag) {

            std::vector<double> tmp_real(density_real);
            std::vector<double> tmp_imag(density_imag);

            for (long long i = 0; i < N_; i++) {
                density_real[i] = tmp_real[sorting_[i]];
                density_imag[i] = tmp_imag[sorting_[i]];

            }

        }

    public:

    // The initializer list x_(x) syntax means x is copied into x_(x), unless x_ is a reference.
        BoxTree(const std::vector<double>& x, const std::vector<double>& y, int nlevels, double wavenumber):
        x_(x), y_(y), nlevels_(nlevels), wavenumber_(wavenumber), N_(x_.size())
        {

            Initialize();

        }

        ~BoxTree() 
        {

        }

        template <void _kernel(const double, const double, const double, const double, const double, const double, double&, double&),
                  void _factorization(const double, double&, double&)>
        void Solve(const bool with_singular_interactions, std::vector<double> & density_real, std::vector<double> & density_imag) {

            if (nlevels_ < 3)
                throw std::invalid_argument("IFGF accelerator requires at least 3 levels.");

            const int D = nlevels_ - 1;

            // P_ = IPSCHEME_.GetNInterpPoints();

            ZeroSolution();

            // TODO SPEEDUP: If we sort the points in the RP method according to the morton order
            // then we don't have to do all this sorting here. However, it is easier to track
            // and debug this way, especially when assigning the quadrature and jacobian weights.
            // This leads to O(N) extra work every iteration, although there is of a necessity
            // non-contiguous memory access pattern in this operation.

            // Put the density in the correct order of the points.
            SortDensity(density_real, density_imag);


            if (with_singular_interactions) 
                SingularInteractions<_kernel>(density_real, density_imag);

            LevelDEvaluations<_kernel, _factorization>(density_real, density_imag);

            for (int level = D; level >= 2; level--) {

                if (level > 2) {

                    Propagation<_factorization>(level);

                }

                Interpolation<_factorization>(level);

                SwapCoefficient();

            }

            for (long long i = 0; i < N_; i++) {

                // Save the density in the original ordering
                // Draw out a picture with three squences to see why this works.
                density_real[sorting_[i]] = solution_real_[i];
                density_imag[sorting_[i]] = solution_imag_[i];

            }

        }

};

#endif