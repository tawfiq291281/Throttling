#ifndef VECTOR_TENSOR_OPS_H
#define VECTOR_TENSOR_OPS_H

#include "Vector.H"
#include "Tensor.H"

namespace Foam
{
inline Vector<double> operator-(const Vector<double>& v, const Tensor<double>& t)
{
    return Vector<double>(v.x() - t.xx(), v.y() - t.yy(), v.z() - t.zz());
}
}

#endif
