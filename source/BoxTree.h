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

class Level;
class BoxTree;

const Interpolator<Chebyshev> IPSCHEME_;
const int P_ = 3*5;

class BoxTree 
{

    private:

        std::vector<double> x_; // x coordinates of all values to be evaluated in the IFGF
        std::vector<double> y_; // y coordinates of all values to be evaluated in the IFGF

        int nlevels_;
        double wavenumber_;
        long long N_; // Total number of points in IFGF.

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

            // TODO: This allocation could be better by using
            // #include <numeric>
            // std::vector<long long> sorting(N_);
            // std::iota(sorting.begin(), sorting.end(), 0);
            // About twice as fast but won't make a noticeable different in the code.
            std::vector<long long> sorting;

            for (long long i = 0; i < N_; i++) {

                sorting.push_back(i);

            }

            std::vector<long long> morton_code(N_);

            for (long long i = 0; i < N_; i++) {

                const long long morton = Level::Point2Morton(x_[i], y_[i], min[0], min[1], boxsize, nlevels_-1);
                morton_code[i] = morton;

            }

            new_sort(sorting, morton_code);

            std::vector<double> tmp_xy(N_);

            for (long long i = 0; i < N_; i++) {

                tmp_xy[i] = x_[sorting[i]];

            }

            x_ = tmp_xy;

            for (long long i = 0; i < N_; i++) {

                tmp_xy[i] = y_[sorting[i]];

            }

            y_ = tmp_xy;

        }

        void InitializeLevelDBoxesData() 
        {

            long long old_morton = -1;
            long long npoints;

            for (long long i = 0; i < N_; i++) {

                const long long morton = levels_.back()->Point2Morton(x_[i], y_[i]);

                if (morton != old_morton) {

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

        void InitializeRelBoxesAllLevels() 
        {

            std::unordered_set<long long> relevantparentmorton;
            std::unordered_set<long long> relevantmorton = levels_.back()->mortonidofrelboxes_;

            for (int i = nlevels_-1; i >= 2; i--) {

                if (i != nlevels_-1)
                    levels_[i]->mortonidofrelboxes_.insert(relevantmorton.begin(), relevantmorton.end());

                for (auto j = relevantmorton.begin(); j != relevantmorton.end(); j++) {
                    // j is an iterator. *j give the value that j points to
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

        void InitializeBoxesAndLevels() 
        {

            std::array<double, 2> min, max;

            ComputeBB(x_, y_, min, max);          

            // Choose the initial box to have side lengths slightly larger than the distance
            // of the x and y coordinates  or the points furthest from the origin
            double boxsize = std::max(max[0] - min[0], max[1] - min[1]);

            // Sort according to morton ordering. TODO: MORTON ORDERING
            SortBox(min, boxsize / (1 << (nlevels_ - 1)));

            levels_.resize(nlevels_, nullptr);

            for (int i = 1; i < nlevels_; i++) {

                levels_[i] = new Level(i, nlevels_, min[0], min[1], boxsize/(1 << i), wavenumber_);

            }

            //
            InitializeLevelDBoxesData();

            InitializeRelBoxesAllLevels();

        }

        void GetRelevantConeSegmentsDueToCousinSurface(const int level, std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment) const {

            long long old_morton_box = -1;
            std::vector<long long> cousins;
            std::unordered_map<long long, std::unordered_set<long long>> tmpmortonboxnonrelconesegment; 

            for (long long i = 0; i < N_; i++) {

                const double x = x_[i];
                const double y = y_[i];

                const long long morton_box = levels_[level]->Point2Morton(x, y);

                if (morton_box != old_morton_box) {

                    cousins = levels_[level]->GetCousins(morton_box);
                    old_morton_box = morton_box;

                }

                for (int j = 0; j < cousins.size(); j++) {

                    const long long morton_cousin = cousins[j];
                    const long long nonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_cousin, x, y);
                    tmpmortonboxnonrelconesegment[morton_cousin].insert(nonrel_conesegment);
    
                }

            }

            for (auto i = tmpmortonboxnonrelconesegment.begin(); i != tmpmortonboxnonrelconesegment.end(); i++) {

                std::vector<long long> tmp(i->second.begin(), i->second.end());
                std::sort(tmp.begin(), tmp.end());
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

            for (long long iter = 0; iter < levels_[level-1]->relconesegment2nonrelconesegment_.size(); iter++) {

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

                for (int childiter = 0; childiter < morton_children.size(); childiter++) {

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

            for (long long iter = 0; iter < levels_[level-1]->relconesegment2nonrelconesegment_.size(); iter++) {

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

                for (int childiter = 0; childiter < morton_children.size(); childiter++) {

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

            for (long long iter = 0; iter < levels_[level]->relconesegment2nonrelconesegment_.size(); iter++) {

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

                for (long long csiter = 0; csiter < nonrelconesegments.size(); csiter++) {

                    const long long relconesegment = levels_[level]->MortonboxNonrelconesegment2Relconesegment(iter->first, nonrelconesegments[csiter]);

                    levels_[level]->interpolationrequiredrelconesegments_.insert(relconesegment);

                }

            }

            for (auto iter = propagation_conesegments.begin(); iter != propagation_conesegments.end(); iter++) {

                const std::vector<long long> nonrelconesegments = iter->second;

                for (long long csiter = 0; csiter < nonrelconesegments.size(); csiter++) {

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

                if (interpolation_conesegments.size() == 0 && propagation_conesegments.size() == 0)
                    throw std::logic_error("The number of cone segments on any level cannot be 0");

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

                for (int neighbouriter = 0; neighbouriter < neighbours.size(); neighbouriter++) {

                    const long long neighbourmorton = neighbours[neighbouriter];

                    const std::array<long long, 2> & points_data = levels_.back()->mortonbox2discretizationpoints_.at(neighbourmorton);

                    points_begin = points_data[0];
                    npoints = points_data[1];

                    for (long long sourceiter = 0; sourceiter < npoints; sourceiter++) {

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

            for (long long i = 0; i < levels_.back()->relconesegment2nonrelconesegment_.size(); i++) {

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
                    conesegments_prev_real_[coefficients_pos] = (result_real * fac_real + result_imag * fac_imag) * dabsfac;
                    conesegments_prev_imag_[coefficients_pos] = (result_imag * fac_real - result_real * fac_imag) * dabsfac;
                
                }

                IPSCHEME_.GenerateInterpolant(&conesegments_prev_real_[coefficients_begin]);
                IPSCHEME_.GenerateInterpolant(&conesegments_prev_imag_[coefficients_begin]);
            
            }            
        
            for (long long i = 0; i < conesegments_current_real_.size(); i++) {

                conesegments_current_real_[i] = conesegments_prev_real_[i];
                conesegments_current_imag_[i] = conesegments_prev_imag_[i];

            }
        
        }

        template<void _factorization(const double, double&, double&)>
        void Propagation(const int level) {

            std::array<double, P_> x, y, radii, value_real, value_imag;
            double locx, locy;
            std::array<double const *, 4> childcoeffs_real, childcoeffs_imag;
            long long old_cocentered_morton_box = -1;
            std::vector<long long> morton_children;
            std::array<long long, 4> old_locrelconesegment = {-1, -1, -1, -1};

            for (long long i = 0; i < levels_[level-1]->relconesegment2nonrelconesegment_.size(); i++) {

                long long rel_conesegment = i; 

                const long long cocentered_morton_box = levels_[level-1]->relconesegment2cocenteredmortonboxid_[rel_conesegment];
                const long long nonrel_conesegment = levels_[level-1]->relconesegment2nonrelconesegment_[rel_conesegment];

                const long long new_coeffs_begin_id = rel_conesegment * IPSCHEME_.GetNInterpPoints(); 
                
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

                for (int childiter = 0; childiter < morton_children.size(); childiter++) {

                    const long long morton_child = morton_children[childiter];
                    const long long child_pos = morton_child % 4;
                    const std::unordered_map<long long, std::vector<int>>& childnonrelcones = allchildnonrelcones.at(child_pos);
                    
                    for (auto childconeiter = childnonrelcones.begin(); childconeiter != childnonrelcones.end(); childconeiter++) {
                        
                        for (int interpiter = 0; interpiter < childconeiter->second.size(); interpiter++) {

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

                for (int cousiniter = 0; cousiniter < morton_cousinboxes.size(); cousiniter++) {

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

                    value_real += tmpreal * fac_real - tmpimag * fac_imag;
                    value_imag += tmpreal * fac_imag + tmpimag * fac_real;

                }

                solution_real_[pointiter] = value_real;
                solution_imag_[pointiter] = value_imag;

            }

        } 

        void SwapCoefficient() {

            for (long long i = 0; i < conesegments_current_real_.size(); i++) {

                double tmp = conesegments_current_real_[i];
                conesegments_current_real_[i] = conesegments_prev_real_[i];
                conesegments_prev_real_[i] = tmp;

                tmp = conesegments_current_imag_[i];
                conesegments_current_imag_[i] = conesegments_prev_imag_[i];
                conesegments_prev_imag_[i] = tmp;

            }

        }

    public:

        BoxTree(std::vector<double>& x, std::vector<double>& y, int nlevels, double wavenumber):
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

                density_real[i] = solution_real_[i];
                density_imag[i] = solution_imag_[i];

            }

        }

};

#endif