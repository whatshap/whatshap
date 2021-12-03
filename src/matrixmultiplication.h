#ifndef MATMUL_H
#define MATMUL_H

#include "vector2d.h"
#include <iostream>
#include <vector>
#include <cassert>
using namespace std;

vector<long double> MatMul(vector<long double> &A, Vector2D<double> &B, bool transpose = false) {
    
    vector<long double> res;
    double value;
    if (!transpose) {
        assert (A.size() == B.get_size0());    
        for (int i = 0; i < B.get_size1(); i++) {
            value = 0;
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
            value = 0;
            for (int j = 0; j < A.size(); j++) {
                value = value + A[j]*B.at(i,j);
            }
            res.push_back(value);
        }
        return res;
    }
    
}

#endif