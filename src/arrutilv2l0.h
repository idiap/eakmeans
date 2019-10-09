/*
Copyright (c) 2015-2018 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <james.newling@gmail.com>
All rights reserved.

eakmeans is a library for exact and approximate k-means written in C++ and
Python. This file is part of eakmeans. See file COPYING for more details.

This file is part of eakmeans.

eakmeans is free software: you can redistribute it and/or modify
it under the terms of the 3-Clause BSD Licence. See
https://opensource.org/licenses/BSD-3-Clause for more details.

eakmeans is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See file
COPYING for more details.
*/

#ifndef ARRUTILV2L0_H
#define ARRUTILV2L0_H

#include <memory>
#include <exception>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <set>
#include <chrono>
#include <map>
#include <limits>

/* Herein all things doable directly by blas
 * Herein all distances, distances squared, mins, maxs, combos therof (which if using blas would be different)
 * 
 * rules of functions to make easier to use / remember:
 * (1) dimensions of arrays must appear before arrays, but as late as possible
 * (2) functions which return must be getxxx
 * (3) functions which set_ must be set_xxx (or subtractfrom , addto , update , something obvious)
 * (4) thing(s) being set_ should come as late as possible (excluding flag like parameters, background increment parameters etc.) without violating above rules 
 * (5) if array being set_ is dimension d, there should be d trailing 's' to function name
 * (6) if operation is on 1-D and 2-D array, should have r/c somewhere telling whether row or column
 * (7) if operation on 2-D and 2-D array should have rr/rc/cr/cc as above (unless a flag like bool asrow)
 * (8) nrows before ncols in parameter list
 * for [TFloat = double, TInt = unsigned] autogeneration of functions to arrutilv2.cpp is done by python function 
 * */

#ifdef WITHBLAS
#include "arrutilv2l0withblas.h"
#else
#include "arrutilv2l0blasless.h"
#endif

#ifdef _MSC_VER
// add setenv under Windows
int setenv(const char *name, const char *value, int overwrite)
{
	int errcode = 0;
	if (!overwrite) {
		size_t envsize = 0;
		errcode = getenv_s(&envsize, NULL, 0, name);
		if (errcode || envsize) return errcode;
	}
	return _putenv_s(name, value);
}
#endif // _MSC_VER

#endif

