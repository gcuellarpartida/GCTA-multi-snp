
/*
 * Interface to the extended function for EIGEN lib
 *
 * 2014 by Jian Yang <jian.yang@uq.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#ifndef _EIGENFUNC_H
#define _EIGENFUNC_H

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#include "CommFunc.h"
#include "StatFunc.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/IterativeSolvers>

#include <map>
#include <vector>
#include <algorithm>

using namespace Eigen;
using namespace std;

namespace eigen_func
{
    void rank(VectorXf &x, VectorXf &rank);
    void inverse_norm_rank_transform(VectorXf &x);
}

#endif

