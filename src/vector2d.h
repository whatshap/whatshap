#ifndef VECTOR_2D_H
#define VECTOR_2D_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>

template <typename T>
class Vector2D {
public:
	Vector2D(size_t size0 = 0, size_t size1 = 0, const T& initial_value = T() ) : size0(size0), size1(size1), v(size0*size1, initial_value) {}

	const T& at(size_t index0, size_t index1) const {
		return v[index0*size1 + index1];
	}

	T& at(size_t index0, size_t index1) {
		return v[index0*size1 + index1];
	}

	void set(size_t index0, size_t index1, const T& value) {
		v[index0*size1 + index1] = value;
	}

	void set_all(const T& value) {
		v.assign(size0*size1, value);
	}

	size_t get_size0(){
	  return size0;
	}

	size_t get_size1(){
	  return size1;
	}

	void divide_entries_by(T val) {
	  std::transform(v.begin(), v.end(), v.begin(),
               std::bind2nd(std::divides<T>(),val));
	}


	friend std::ostream& operator<<(std::ostream& out, const Vector2D<T>& v) {
		out << "       ";
		for (size_t j=0; j<v.size1; ++j) {
			out << ' ' << std::setw(11) << j;
		}
		out << std::endl;
		for (size_t i=0; i<v.size0; ++i) {
			if (i==0) out << "[[";
			else out << " [";
			out << std::setw(4) << i << ':';
			for (size_t j=0; j<v.size1; ++j) {
				out << ' ' << std::setw(11) << v.at(i,j);
			}
			if (i==v.size0-1) out << "]]";
			else out << "]" << std::endl;
		}
		return out;
	}

private:
	size_t size0;
	size_t size1;
	std::vector<T> v;
};

#endif
