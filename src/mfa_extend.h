#pragma once
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>
#include <mfa/mfa.hpp>
#include "opts.h"
#include "block.hpp"
#include "../utility/utility_function.h"


namespace mfa_extend
{
    //prepare the decoder ct(), cs()
    template<typename T>
    void decoder_params(MatrixXi&  ct, VectorXi& cs, const mfa::MFA_Data<T>&  mfa_data)
    {
        int tot_iters = (mfa_data.p + VectorXi::Ones(mfa_data.dom_dim)).prod();
        cs = VectorXi::Ones(mfa_data.dom_dim);
        ct.resize(tot_iters, mfa_data.p.size());
        for (size_t i = 0; i < mfa_data.p.size(); i++)   // for all dims
        {
            if (i > 0)
            {
                cs[i] = cs[i - 1] * mfa_data.tmesh.tensor_prods[0].nctrl_pts[i - 1];
            }
        }



        for (int i = 0; i < tot_iters; i++)      // 1-d flattening all n-d nested loop computations
        {
            int div = tot_iters;
            int i_temp = i;
            for (int j = mfa_data.p.size() - 1; j >= 0; j--)
            {
                div      /= (mfa_data.p(j) + 1);
                ct(i, j) =  i_temp / div;
                i_temp   -= (ct(i, j) * div);
            }
        }
    }

    // compute a point from a NURBS n-d volume at a given parameter value
    // slower version for single points
    // based on a seleted span
    template<typename T>
    void VolPt(
            const VectorX<T>&       param,      // parameter value in each dim. of desired point
            VectorX<T>&             out_pt,     // (output) point, allocated by caller
            const TensorProduct<T>& tensor,     // tensor product to use for decoding
            vector<int>&            span,       //span index should be in [p, n]
            const mfa::MFA_Data<T>& mfa_data,
            MatrixXi&  ct, VectorXi& cs,
            const VectorXi&         derivs)     // derivative to take in each domain dim. (0 = value, 1 = 1st deriv, 2 = 2nd deriv, ...)
                                                // pass size-0 vector if unused
    {
        int last = mfa_data.tmesh.tensor_prods[0].ctrl_pts.cols() - 1;      // last coordinate of control point
        if (derivs.size())                                                  // extra check for derivatives, won't slow down normal point evaluation
        {
            if (derivs.size() != mfa_data.p.size())
            {
                fprintf(stderr, "Error: size of derivatives vector is not the same as the number of domain dimensions\n");
                exit(0);
            }
            for (auto i = 0; i < mfa_data.p.size(); i++)
                if (derivs(i) > mfa_data.p(i))
                    fprintf(stderr, "Warning: In dimension %d, trying to take derivative %d of an MFA with degree %d will result in 0. This may not be what you want",
                            i, derivs(i), mfa_data.p(i));
        }

        // init
        vector <MatrixX<T>> N(mfa_data.p.size());                           // basis functions in each dim.
        vector<VectorX<T>>  temp(mfa_data.p.size());                        // temporary point in each dim.
        // vector<int>         span(mfa_data.p.size());                        // span in each dim.
        VectorX<T>          ctrl_pt(last + 1);                              // one control point
        int                 ctrl_idx;                                       // control point linear ordering index
        VectorX<T>          temp_denom = VectorX<T>::Zero(mfa_data.p.size());// temporary rational NURBS denominator in each dim

        // set up the volume iterator
        VectorXi npts = mfa_data.p + VectorXi::Ones(mfa_data.dom_dim);      // local support is p + 1 in each dim.
        mfa::VolIterator vol_iter(npts);                                         // for iterating in a flat loop over n dimensions

        // basis funs
        for (size_t i = 0; i < mfa_data.dom_dim; i++)                       // for all dims
        {
            temp[i]    = VectorX<T>::Zero(last + 1);
            // span[i]    = mfa_data.tmesh.FindSpan(i, param(i), tensor);
            N[i]       = MatrixX<T>::Zero(1, tensor.nctrl_pts(i));
            if (derivs.size() && derivs(i))
            {
#ifndef MFA_TMESH   // original version for one tensor product
                MatrixX<T> Ders = MatrixX<T>::Zero(derivs(i) + 1, tensor.nctrl_pts(i));
                mfa_data.DerBasisFuns(i, param(i), span[i], derivs(i), Ders);
                N[i].row(0) = Ders.row(derivs(i));
#endif
            }
            else
            {
#ifndef MFA_TMESH   // original version for one tensor product
                mfa_data.OrigBasisFuns(i, param(i), span[i], N[i], 0);
#else               // tmesh version
                mfa_data.BasisFuns(i, param(i), span[i], N[i], 0);
#endif
            }
        }

        // linear index of first control point
        ctrl_idx = 0;
        for (int j = 0; j < mfa_data.p.size(); j++)
            ctrl_idx += (span[j] - mfa_data.p(j) + ct(0, j)) * cs[j];
        size_t start_ctrl_idx = ctrl_idx;

        while (!vol_iter.done())
        {
            // always compute the point in the first dimension
            ctrl_pt = tensor.ctrl_pts.row(ctrl_idx);
            T w     = tensor.weights(ctrl_idx);

#ifdef WEIGH_ALL_DIMS                                                           // weigh all dimensions
            temp[0] += (N[0])(0, vol_iter.idx_dim(0) + span[0] - mfa_data.p(0)) * ctrl_pt * w;
#else                                                                           // weigh only range dimension
            for (auto j = 0; j < last; j++)
                (temp[0])(j) += (N[0])(0, vol_iter.idx_dim(0) + span[0] - mfa_data.p(0)) * ctrl_pt(j);
            (temp[0])(last) += (N[0])(0, vol_iter.idx_dim(0) + span[0] - mfa_data.p(0)) * ctrl_pt(last) * w;
#endif

            temp_denom(0) += w * N[0](0, vol_iter.idx_dim(0) + span[0] - mfa_data.p(0));

            vol_iter.incr_iter();                                           // must call near bottom of loop, but before checking for done span below

            // for all dimensions except last, check if span is finished
            ctrl_idx = start_ctrl_idx;
            for (size_t k = 0; k < mfa_data.p.size(); k++)
            {
                if (vol_iter.cur_iter() < vol_iter.tot_iters())
                    ctrl_idx += ct(vol_iter.cur_iter(), k) * cs[k];         // ctrl_idx for the next iteration
                if (k < mfa_data.dom_dim - 1 && vol_iter.done(k))
                {
                    // compute point in next higher dimension and reset computation for current dim
                    // use prev_idx_dim because iterator was already incremented above
                    temp[k + 1]        += (N[k + 1])(0, vol_iter.prev_idx_dim(k + 1) + span[k + 1] - mfa_data.p(k + 1)) * temp[k];
                    temp_denom(k + 1)  += temp_denom(k) * N[k + 1](0, vol_iter.prev_idx_dim(k + 1) + span[k + 1] - mfa_data.p(k + 1));
                    temp_denom(k)       = 0.0;
                    temp[k]             = VectorX<T>::Zero(last + 1);
                }
            }
        }

        T denom;                                                            // rational denominator
        if (derivs.size() && derivs.sum())
            denom = 1.0;                                                    // TODO: weights for derivatives not implemented yet
        else
            denom = temp_denom(mfa_data.p.size() - 1);

#ifdef WEIGH_ALL_DIMS                                                           // weigh all dimensions
        out_pt = temp[mfa_data.p.size() - 1] / denom;
#else                                                                           // weigh only range dimension
        out_pt   = temp[mfa_data.p.size() - 1];
        out_pt(last) /= denom;
#endif


    }

            // Decode variable model at single point
    template<typename T>
    void DecodeVar(
            const VectorX<T>&   param,
            VectorX<T>&         out_point,
            const mfa::MFA_Data<T>& mfa_data,
            vector<int>&            span,       //span index should be in [p, n]
            const VectorXi&     derivs = VectorXi())
    {
        MatrixXi ct;
        VectorXi cs;
        decoder_params(ct, cs, mfa_data);
        VolPt(param, out_point, mfa_data.tmesh.tensor_prods[0], span, mfa_data,ct,cs, derivs);
    }

    template<typename T>
    void move_point_to_span(VectorX<T>&   param, VectorX<T>& correct_param)
    {
        correct_param = param;
        for(int i = 0; i < param.size(); ++i)
        {
            if (correct_param(i) < 0)
            {
                correct_param(i) = 0;
            }
            else if (correct_param(i) > 1)
            {
                correct_param(i) = 1;
            }
        }
    }

    // we extend the query to entire space. If a point is outside the domain, use the info of the closest span.
    template<typename T>
    void recover_mfa(
            const Block<T>* b,               // mfa data model
            const VectorX<T>&   point,                  // parameters of point to decode
            VectorX<T>&         cpt,const VectorXi&     derivs = VectorXi())
    {
        VectorX<T> domain_range = b->core_maxs - b->core_mins;
        VectorX<T> param = (point-b->core_mins).cwiseQuotient(domain_range);
        cpt.resize(b->mfa->nvars());

        // if(utility::InDomain(param))
        // {
        b->mfa->DecodeVar(0,param,cpt,derivs);
        // }
        // else
        // {
        //     VectorX<T> correct_param;
        //     move_point_to_span(param, correct_param);
        //     const mfa::MFA_Data<T>& mfa_data = b->mfa->var(0);
        //     vector<int>         span(mfa_data.p.size());  
        //     for(int i = 0; i < param.size(); ++i)
        //     {
        //         span[i]    = mfa_data.tmesh.FindSpan(i, correct_param(i), mfa_data.tmesh.tensor_prods[0]);
        //     }

        //     DecodeVar(param, cpt, mfa_data, span, derivs);
        // }

        if (derivs.size() > 0)
        {
            for(int i=0;i<derivs.size();++i)
            {
                if(derivs[i]>0)
                {
                    cpt = cpt/pow(domain_range[i],derivs[i]);
                }
            }
        }
    }

}