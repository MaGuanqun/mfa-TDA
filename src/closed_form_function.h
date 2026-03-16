#pragma once

#include <Eigen/Dense>


namespace closed_form_function
{
    int initial_func_type(string function_name)
    {
        if(function_name=="quartic_potential")
        {
            return 1;
        }
        if(function_name=="quartic_potential_2")
        {
            return 2;
        }
        if(function_name=="rotating_quartic_multiwell")
        {
            return 3;
        }
        if(function_name=="quartic_potential_3d")
        {
            return 4;
        }
        std::cout<<"error: invalid function name"<<std::endl;
        exit(0);
        return -1;

    }


    VectorXd domain_min(int func_type)
    {
        VectorXd result(3);
        switch (func_type)
        {
        case 1:
            result << -2.0,-2.0,0.0;
            break;
        case 2:
            result << -2.0,-2.0,0.0;
            break;
        case 3:
            result << -2.0,-2.0,0.0;
            break;
        case 4:
            result.resize(4);
            result << -2.0,-2.0,-2.0,0.0;
            break;
        default:
            break;
        }
        return result;
    }

    VectorXd domain_max(int func_type)
    {
        VectorXd result(3);
        switch (func_type)
        {
        case 1:
            result << 2.0,2.0,4.0;
            break;
        case 2:
            result << 2.0,2.0,4.0;
            break;
        case 3:
            result << 2.0,2.0,4.0;
            break;
        case 4:
            result.resize(4);
            result << 2.0,2.0,2.0,4.0;
            break;
        default:
            break;
        }
        return result;      
    }

    VectorXi block_num(int func_type) //number of blocks that splits the domain in each dimension
    {
        VectorXi result(3);
        switch (func_type)
        {
        case 1:
            result << 10,10,10;
            break;
        case 2:
            result << 10,10,10;
            break;
        case 3:
            result << 10,10,10;
            break;
        case 4:
            result.resize(4);
            result << 10,10,10,10;
            break;
        default:
            break;
        }
        return result;
    }


    VectorXi point_num_in_block(int func_type) //number of initial points in each block in each dimension
    {
        VectorXi result(3);
        switch (func_type)
        {
        case 1:
            result << 5,5,5;
            break;
        case 2:
            result << 5,5,5;
            break;
        case 3:
            result << 5,5,5;
            break;
        case 4:
            result.resize(4);
            result << 4,4,4,4;
            break;
        default:
            break;
        }
        return result;
    }

    //$$f(x,y,t)=\frac{x^4}{4}+\frac{1-t}{2}\,x^2+\frac{1}{2}\,(y-\sin t)^4$$
    //the third derivatives only contain /partial_t = 1 or 0
    template<typename T>
    void quartic_potential(const VectorX<T>&   point,VectorX<T>&         result,const VectorXi&     derivs = VectorXi())
    {
        result.resize(1);
        if(derivs.size()==0)
        {
            result(0)=0.25*pow(point(0),4)+0.5*(1-point(2))*point(0)*point(0)+0.5*pow(point(1)-sin(point(2)),4);
            return;
        }
        if(derivs.sum()==1)
        {
            if(derivs(0)==1)
            {
                result(0)=point(0)*point(0)*point(0)+point(0)*(1-point(2));
                return;
            }
            if(derivs(1)==1)
            {
                result(0)=2.0*pow(point(1)-sin(point(2)),3);
                return;
            }
            // derivs(2)==1
            result(0)=-0.5*point(0)*point(0)-2.0*pow(point(1)-sin(point(2)),3)*cos(point(2));
            return;
            
        }
        if(derivs.sum()==2)
        {
            if(derivs(2)==2)
            {
                T sin_t=sin(point(2));
                result(0)=2.0*pow(point(1)-sin_t,3)*sin_t+6.0*pow(point(1)-sin_t,2)*pow(cos(point(2)),2);
                return;
            }
            if(derivs(2)==1)
            {
                if(derivs(1)==1) //(0,1,1)
                {
                    result(0)=-6.0*pow(point(1)-sin(point(2)),2)*cos(point(2));
                    return;
                }
                // derivs(0)==1 (1,0,1)
                result(0)=-point(0);
                return;
            }
            // derivs(2)==0
            if(derivs(1)==2) // (0,2,0)
            {
                result(0)=6.0*pow(point(1)-sin(point(2)),2);
                return;
            }

            if(derivs(0)==2) //(2,0,0)
            {
                result(0)=3.0*point(0)*point(0)+1.0-point(2);
                return;
            }

            //(1,1,0)
            result(0)=0.0;
            return;
            
        }
        if(derivs.sum()==3)
        {
            if(derivs(2)==0)
            {
                if(derivs(1)==0) //(3,0,0)
                {
                    result(0)=6.0*point(0);
                    return;
                }
                if(derivs(0)==0) //(0,3,0)
                {
                    result(0)=12.0*(point(1)-sin(point(2)));
                    return;
                }
                //(2,1,0) // (1,2,0)
                result(0)=0.0;
                return;
            }

            // the following must be derivs(2)==1 
            //(2,0,1)
            if(derivs(1)==0)
            {
               result(0)=-1.0;
                return;
            }
            if(derivs(1)==1)  // (1,1,1)
            {
                result(0)=0.0;
                return;
            }
            // derivs(1)==2  // (0,2,1)
                result(0)=-12.0*(point(1)-sin(point(2)))*cos(point(2));
                return;
        }

        std::cout<<"Error: the order of derivative is larger than 3, which is not supported!"<<std::endl;
        exit(0);
    }

//$$f(x,y,t)=\frac{x^4}{4}+\frac{1-t}{2}\,x^2+\frac{y^4}{4}+\frac{\cos t}{2}y^2$$
    //the third derivatives only contain /partial_t = 1 or 0
    template<typename T>
    void quartic_potential_2(const VectorX<T>&   point,VectorX<T>&         result,const VectorXi&     derivs = VectorXi())
    {
        result.resize(1);
        if(derivs.size()==0)
        {
            result(0)=0.25*pow(point(0),4)+0.5*(1-point(2))*point(0)*point(0)+0.25*pow(point(1),4)+0.5*cos(point(2))*point(1)*point(1);
            return;
        }
        if(derivs.sum()==1)
        {
            if(derivs(0)==1)
            {
                result(0)=point(0)*point(0)*point(0)+point(0)*(1-point(2));
                return;
            }
            if(derivs(1)==1)
            {
                result(0)=point(1)*point(1)*point(1)+point(1)*cos(point(2));
                return;
            }
            // derivs(2)==1
            result(0)=-0.5*point(0)*point(0)-0.5*point(1)*point(1)*sin(point(2));
            return;
            
        }
        if(derivs.sum()==2)
        {
            if(derivs(2)==2)
            {
                result(0)=-0.5*point(1)*point(1)*cos(point(2));
                return;
            }
            if(derivs(2)==1)
            {
                if(derivs(1)==1) //(0,1,1)
                {
                    result(0)=-point(1)*sin(point(2));
                    return;
                }
                // derivs(0)==1 (1,0,1)
                result(0)=-point(0);
                return;
            }
            // derivs(2)==0
            if(derivs(1)==2) // (0,2,0)
            {
                result(0)=3.0*point(1)*point(1)+cos(point(2));
                return;
            }

            if(derivs(0)==2) //(2,0,0)
            {
                result(0)=3.0*point(0)*point(0)+1.0-point(2);
                return;
            }

            //(1,1,0)
            result(0)=0.0;
            return;
            
        }
        if(derivs.sum()==3)
        {
            if(derivs(2)==0)
            {
                if(derivs(1)==0) //(3,0,0)
                {
                    result(0)=6.0*point(0);
                    return;
                }
                if(derivs(0)==0) //(0,3,0)
                {
                    result(0)=6.0*point(1);
                    return;
                }
                //(2,1,0) // (1,2,0)
                result(0)=0.0;
                return;
            }

            // the following must be derivs(2)==1 
            //(2,0,1)
            if(derivs(1)==0)
            {
               result(0)=-1.0;
                return;
            }
            if(derivs(1)==1)  // (1,1,1)
            {
                result(0)=0.0;
                return;
            }
            // derivs(1)==2  // (0,2,1)
                result(0)=-sin(point(2));
                return;
        }

        std::cout<<"Error: the order of derivative is larger than 3, which is not supported!"<<std::endl;
        exit(0);
    }




    //$$f(x,y,t)=\frac{x^4}{4}+\frac{1-t}{2}x^2+\frac{y^4}{4}+\frac{\cos t}{2}y^2 +\frac{z^4}{4}+\frac{\cos t}{2}z^2$$
    //the third derivatives only contain /partial_t = 1 or 0
    template<typename T>
    void quartic_potential_3d(const VectorX<T>&   point,VectorX<T>&         result,const VectorXi&     derivs = VectorXi())
    {
        result.resize(1);
        if(derivs.size()==0)
        {
            result(0)=0.25*pow(point(0),4)+0.5*(1-point(3))*point(0)*point(0)+0.25*pow(point(1),4)+0.5*cos(point(3))*point(1)*point(1)+0.25*pow(point(2),4)+0.5*cos(point(3))*point(2)*point(2);
            return;
        }
        if(derivs.sum()==1)
        {
            if(derivs(0)==1)
            {
                result(0)=point(0)*point(0)*point(0)+point(0)*(1-point(3));
                return;
            }
            if(derivs(1)==1)
            {
                result(0)=point(1)*point(1)*point(1)+point(1)*cos(point(3));
                return;
            }
            if(derivs(2)==1)
            {
                result(0)=point(2)*point(2)*point(2)+point(2)*cos(point(3));
                return;
            }
            // derivs(3)==1
            result(0)=-0.5*point(0)*point(0)-0.5*point(1)*point(1)*sin(point(3))-0.5*point(2)*point(2)*sin(point(3));
            return;
            
        }
        if(derivs.sum()==2)
        {
            if(derivs(3)==2)
            {
                result(0)=-0.5*point(1)*point(1)*cos(point(3))-0.5*point(2)*point(2)*cos(point(3));
                return;
            }
            if(derivs(3)==1)
            {
                if(derivs(1)==1 && derivs(2)==0) //(0,1,0,1)
                {
                    result(0)=-point(1)*sin(point(3));
                    return;
                }
                if(derivs(1)==0 && derivs(2)==1) //(0,0,1,1)
                {
                    result(0)=-point(2)*sin(point(3));
                    return;
                }
                // derivs(0)==1 (1,0,0,1)
                result(0)=-point(0);
                return;
            }
            // derivs(3)==0
            if(derivs(1)==2) // (0,2,0,0)
            {
                result(0)=3.0*point(1)*point(1)+cos(point(3));
                return;
            }

            if(derivs(0)==2) //(2,0,0,0)
            {
                result(0)=3.0*point(0)*point(0)+1.0-point(3);
                return;
            }
            if(derivs(2)==2) //(0,0,2,0)
            {
                result(0)=3.0*point(2)*point(2)+cos(point(3));
                return;
            }
            //(1,1,0,0) // (0,1,1,0) // (1,0,1,0)
            result(0)=0.0;
            return;
            
        }
        if(derivs.sum()==3)
        {
            if(derivs(3)==0)
            {
                if(derivs(0)==3) //(3,0,0,0)
                {
                    result(0)=6.0*point(0);
                    return;
                }
                if(derivs(1)==3) //(0,3,0,0)
                {
                    result(0)=6.0*point(1);
                    return;
                }
                if(derivs(2)==3) //(0,0,3,0)
                {
                    result(0)=6.0*point(2);
                    return;
                }
                //(2,1,0,0) // (1,2,0,0) // (0,1,2,0) (1,1,1,0)
                result(0)=0.0;
                return;
            }

            // the following must be derivs(3)==1 
            //(2,0,0,1)
            if(derivs(0)==2)
            {
               result(0)=-1.0;
                return;
            }
            if(derivs(1)==2)
            {
                result(0)=-sin(point(3));
                return;
            }
            if(derivs(2)==2)
            {
                result(0)=-sin(point(3));
                return;
            }
            // (1,1,0,1) (1,0,1,1) (0,1,1,1)
                result(0)=0.0;
                
                return;
        }

        std::cout<<"Error: the order of derivative is larger than 3, which is not supported!"<<std::endl;
        exit(0);
    }

// Rotating quartic multi-well potential (lazy derivative evaluation)
// f(x,y,t) = 1/4 * [ (x cos t + y sin t)^2 - 1 ]^2
//          + 1/4 * [ (-x sin t + y cos t)^2 - 1 ]^2
//
// point(0) = x, point(1) = y, point(2) = t
// derivs = (dx, dy, dt), sum(derivs) <= 3
// For order 3, only dt = 0 or 1 is supported.

    template<typename T>
    void rotating_quartic_multiwell(const VectorX<T>& point,
                                    VectorX<T>&       result,
                                    const VectorXi&   derivs = VectorXi())
    {
        result.resize(1);

        const T x = point(0);
        const T y = point(1);
        const T t = point(2);

        const T c = std::cos(t);
        const T s = std::sin(t);

        // rotated coordinates (always needed)
        const T u =  x * c + y * s;
        const T v = -x * s + y * c;

        const int order = (derivs.size() == 0 ? 0 : derivs.sum());
        const int dx = (derivs.size() > 0 ? derivs(0) : 0);
        const int dy = (derivs.size() > 1 ? derivs(1) : 0);
        const int dt = (derivs.size() > 2 ? derivs(2) : 0);

        // ---------------------
        // 0th order: value only
        // ---------------------
        if(order == 0) {
            const T u2 = u * u;
            const T v2 = v * v;
            result(0) = T(0.25) * ((u2 - T(1)) * (u2 - T(1))
                                + (v2 - T(1)) * (v2 - T(1)));
            return;
        }

        // convenience: basic spatial Jacobian (cheap, always used for derivatives)
        const T ux = c;
        const T uy = s;
        const T vx = -s;
        const T vy = c;

        // ---------------------
        // 1st order derivatives
        // ---------------------
        if(order == 1) {
            // only compute fu,fv when needed
            const T u2 = u * u;
            const T v2 = v * v;
            const T fu = (u2 - T(1)) * u;
            const T fv = (v2 - T(1)) * v;

            if(dx == 1 && dy == 0 && dt == 0) {
                // f_x = f_u u_x + f_v v_x
                result(0) = fu * ux + fv * vx;
                return;
            }
            if(dx == 0 && dy == 1 && dt == 0) {
                // f_y = f_u u_y + f_v v_y
                result(0) = fu * uy + fv * vy;
                return;
            }
            if(dx == 0 && dy == 0 && dt == 1) {
                // f_t = f_u u_t + f_v v_t
                // u_t = v, v_t = -u  (omega = 1)
                const T ut = v;
                const T vt = -u;
                result(0) = fu * ut + fv * vt;
                return;
            }

            std::cerr << "Error: unsupported first-order multi-index!\n";
            std::exit(0);
        }

        // ---------------------
        // 2nd order derivatives
        // ---------------------
        if(order == 2) {
            const T u2 = u * u;
            const T v2 = v * v;
            const T fu  = (u2 - T(1)) * u;
            const T fv  = (v2 - T(1)) * v;
            const T fuu = T(3) * u2 - T(1);
            const T fvv = T(3) * v2 - T(1);

            const T ut  = v;   // u_t
            const T vt  = -u;  // v_t
            const T uxt = vx;  // u_xt
            const T uyt = vy;  // u_yt
            const T vxt = -ux; // v_xt
            const T vyt = -uy; // v_yt
            const T utt = -u;  // u_tt
            const T vtt = -v;  // v_tt

            // // spatial second derivatives
            // const T f_xx = fuu * ux * ux + fvv * vx * vx;
            // const T f_xy = fuu * ux * uy + fvv * vx * vy;
            // const T f_yy = fuu * uy * uy + fvv * vy * vy;

            // // mixed with t
            // const T f_xt = fuu * ut * ux + fu * uxt
            //             + fvv * vt * vx + fv * vxt;

            // const T f_yt = fuu * ut * uy + fu * uyt
            //             + fvv * vt * vy + fv * vyt;

            // const T f_tt = fuu * ut * ut + fvv * vt * vt
            //             + fu * utt + fv * vtt;

            if(dt == 0) {
                if(dx == 2 && dy == 0) { result(0) = fuu * ux * ux + fvv * vx * vx; return; }
                if(dx == 1 && dy == 1) { result(0) = fuu * ux * uy + fvv * vx * vy; return; }
                if(dx == 0 && dy == 2) { result(0) = fuu * uy * uy + fvv * vy * vy; return; }
            } else if(dt == 1) {
                if(dx == 1 && dy == 0) { result(0) = fuu * ut * ux + fu * uxt + fvv * vt * vx + fv * vxt; return; }
                if(dx == 0 && dy == 1) { result(0) = fuu * ut * uy + fu * uyt + fvv * vt * vy + fv * vyt; return; }
            } else if(dt == 2) {
                if(dx == 0 && dy == 0) { result(0) = fuu * ut * ut + fvv * vt * vt + fu * utt + fv * vtt; return; }
            }

            std::cerr << "Error: unsupported second-order multi-index!\n";
            std::exit(0);
        }

        // ---------------------
        // 3rd order derivatives (dt = 0 or 1)
        // ---------------------
        if(order == 3) {
            if(dt > 1) {
                std::cerr << "Error: third derivatives with dt > 1 are not supported!\n";
                std::exit(0);
            }

            const T u2 = u * u;
            const T v2 = v * v;
            const T fu   = (u2 - T(1)) * u;
            const T fv   = (v2 - T(1)) * v;
            const T fuu  = T(3) * u2 - T(1);
            const T fvv  = T(3) * v2 - T(1);
            const T fuuu = T(6) * u;
            const T fvvv = T(6) * v;

            const T ut  = v;
            const T vt  = -u;
            const T uxt = vx;
            const T uyt = vy;
            const T vxt = -ux;
            const T vyt = -uy;

            // // purely spatial third derivatives
            // const T f_xxx = fuuu * ux * ux * ux + fvvv * vx * vx * vx;
            // const T f_xxy = fuuu * ux * ux * uy + fvvv * vx * vx * vy;
            // const T f_xyy = fuuu * ux * uy * uy + fvvv * vx * vy * vy;
            // const T f_yyy = fuuu * uy * uy * uy + fvvv * vy * vy * vy;

            // // mixed (one t)
            // const T f_xxt = fuuu * ut * ux * ux
            //             + T(2) * fuu * ux * uxt
            //             + fvvv * vt * vx * vx
            //             + T(2) * fvv * vx * vxt;

            // const T f_xyt = fuuu * ut * ux * uy
            //             + fuu * (uxt * uy + ux * uyt)
            //             + fvvv * vt * vx * vy
            //             + fvv * (vxt * vy + vx * vyt);

            // const T f_yyt = fuuu * ut * uy * uy
            //             + T(2) * fuu * uy * uyt
            //             + fvvv * vt * vy * vy
            //             + T(2) * fvv * vy * vyt;

            if(dt == 0) {
                if(dx == 3 && dy == 0) { result(0) = fuuu * ux * ux * ux + fvvv * vx * vx * vx; return; }
                if(dx == 2 && dy == 1) { result(0) = fuuu * ux * ux * uy + fvvv * vx * vx * vy; return; }
                if(dx == 1 && dy == 2) { result(0) = fuuu * ux * uy * uy + fvvv * vx * vy * vy; return; }
                if(dx == 0 && dy == 3) { result(0) = fuuu * uy * uy * uy + fvvv * vy * vy * vy; return; }

                std::cerr << "Error: unsupported spatial third-order multi-index!\n";
                std::exit(0);
            } else { // dt == 1
                if(dx == 2 && dy == 0) { result(0) = fuuu * ut * ux * ux
                        + T(2) * fuu * ux * uxt
                        + fvvv * vt * vx * vx
                        + T(2) * fvv * vx * vxt; return; }
                if(dx == 1 && dy == 1) { result(0) = fuuu * ut * ux * uy
                        + fuu * (uxt * uy + ux * uyt)
                        + fvvv * vt * vx * vy
                        + fvv * (vxt * vy + vx * vyt); return; }
                if(dx == 0 && dy == 2) { result(0) = fuuu * ut * uy * uy
                        + T(2) * fuu * uy * uyt
                        + fvvv * vt * vy * vy
                        + T(2) * fvv * vy * vyt; return; }

                std::cerr << "Error: unsupported third-order multi-index with dt=1!\n";
                std::exit(0);
            }
        }

        std::cerr << "Error: derivative order > 3 is not supported!\n";
        std::exit(0);
    }
    // template<typename T>
    // void closed_form_function(const VectorX<T>& point,VectorX<T>& result, const int function_type=0, const VectorXi& derivs = VectorXi())
    // {
    //     switch (function_type)
    //     {
    //     case 1:
    //         quartic_potential(point,result,derivs);
    //         break;
    //     case 2:
    //         quartic_potential_2(point,result,derivs);
    //         break;
    //     default:
    //         std::cout<<"error: invalid function type"<<std::endl;
    //         exit(0);
    //         break;
    //     }
    // }

}