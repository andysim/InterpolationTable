#include <stdio.h>
#include <math.h>
#include <algorithm>

class ERFCInterpolator
{
  private:
    /// The number of elements in the lookup table
    int tabledim_;
    /// Distance (in argument space) between table entries
    double spacing_;
    /// Inverse of the table spacing, to avoid many divisions
    double invspace_;
    /// Sets up the object
    void init();
    /// The table containing values to be interpolated
    double *table_ = NULL;
    /// The maximum value to interpolate
    double maxinterp_;
    /// The maximum value in the table
    double max_;
  public:
    ERFCInterpolator();
    ~ERFCInterpolator();
    /// Uses cubic interpolation to approximate erfc(kR)
    double compute(double kR) const; 
};




void
ERFCInterpolator::init()
{
    // This can be hardwired to some reasonable value.  This choice yields a
    // maximum error of 
    tabledim_ = 800;
    // erfc(6) = 2.15E-17, so we're done after this.
    max_ = 6.0;
    table_ = new double[tabledim_ + 1]; // Don't forget to add 1, because we loop from 0 <= x <= npts
    spacing_ = max_ / (double)tabledim_;
    maxinterp_ = max_ - spacing_;
    invspace_ = 1.0/spacing_;
    for(int i = 0; i <= tabledim_; ++i)
        table_[i] = erfc((double)i*spacing_);
}


ERFCInterpolator::~ERFCInterpolator()
{
    if(table_)
        delete[] table_;
}

ERFCInterpolator::ERFCInterpolator()
{
    init();
}


double
ERFCInterpolator::compute(const double kR) const
{
    // Handle the case where we're sampling around the last element, or beyond
    if(kR >= maxinterp_)
        return 0.0;
    double x = kR*invspace_;
    int intx = int(x+0.5);
    // Make sure we have some wiggle room to get three points to sample
    if(intx == 0)
        ++intx;
    double remainder = x - (double)intx;
    // Read three contiguous values.  Write it using incremented pointers, to help the
    // compiler identify contiguous memory access.
    double *ptr = &table_[intx-1];
    double valm = *ptr; ++ptr;
    double val0 = *ptr; ++ptr;
    double valp = *ptr;
    double d1 = 0.5*(valm - valp);
    double d2 = (val0+val0-valm-valp)*remainder;
    // Linear interpolation
    double erfckR = val0 - d1*remainder;
    // Add cubic terms
    erfckR -= 0.5*d2*remainder;
    return erfckR;
}


int main()
{
    ERFCInterpolator erfcapprox;
    // Maximum value to sample
    double tablemax = 6.2;
    // Number of (evenly-spaced) points to sample
    int npts = 20000;
    double del = tablemax / (double)npts;
    double maxerror = 0.0;
    for(int i = 0; i <= npts; ++i){
       double val = (double)i * del;
       double me = erfcapprox.compute(val);
       double exact = erfc(val);
       double error = me - exact;
       printf("%8.4f: ME: %20.10f, EXACT: %20.10f ERROR: %20.16f\n", val, me, exact, error);
       maxerror = fabs(error) > maxerror ? fabs(error) : maxerror;
    }
    printf("\nMax error: %10.2E\n", maxerror);
    return 0;
}
