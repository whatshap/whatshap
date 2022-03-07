#ifndef MATMUL_H
#define MATMUL_H

#include "vector2d.h"
#include <iostream>
#include <vector>
#include <cassert>
using namespace std;

vector<long double> MatMul(vector<long double> &A, Vector2D<long double> &B, bool transpose) {
    
    vector<long double> res;
    if (!transpose) {
        assert (A.size() == B.get_size0());    
        for (int i = 0; i < B.get_size1(); i++) {
            long double value = 0.0L;
            for (int j = 0; j < A.size(); j++) {
                value = value + A[j]*B.at(j,i);
            }
            res.push_back(value);
        }
        return res;
    }
    else {
        assert (A.size() == B.get_size1());    
        for (int i = 0; i < B.get_size0(); i++) {
            long double value = 0.0L;
            for (int j = 0; j < A.size(); j++) {
                value = value + A[j]*B.at(i,j);
            }
            res.push_back(value);
        }
        return res;
    }
    
}

#endif