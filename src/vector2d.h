#ifndef VECTOR_2D_H
#define VECTOR_2D_H

#include <vector>

template <typename T> 
class Vector2D {
public:
	Vector2D(size_t size0 = 0, size_t size1 = 0, const T& initial_value = T() ) : size0(size0), size1(size1), v(size0*size1, initial_value) {}

	const T& at(size_t index0, size_t index1) const {
		return v[index0*size1 + index1];
	}

	void set(size_t index0, size_t index1, const T& value) {
		v[index0*size1 + index1] = value;
	}

	void set_all(const T& value) {
		v.assign(size0*size1, value);
	}

private:
	size_t size0;
	size_t size1;
	std::vector<T> v;
};

#endif
