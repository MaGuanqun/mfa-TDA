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

    template<typename T>
    void closed_form_function(const VectorX<T>& point,VectorX<T>& result, const int function_type=0, const VectorXi& derivs = VectorXi())
    {
        switch (function_type)
        {
        case 1:
            quartic_potential(point,result,derivs);
            break;
        default:
            std::cout<<"error: invalid function type"<<std::endl;
            exit(0);
            break;
        }
    }

}