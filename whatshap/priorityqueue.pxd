from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.pair cimport pair
from libcpp cimport bool

ctypedef int item_type
ctypedef vector[int] priority_type
ctypedef priority_type* priority_type_ptr
ctypedef pair[priority_type_ptr,item_type] queue_entry_type


cdef class PriorityQueue:
	cdef vector[queue_entry_type] heap
	cdef unordered_map[item_type,int] positions

	cdef void c_push(self, priority_type_ptr score, int item)
	cdef void _swap(self, int index1, int index2)
	cdef bool _score_lower(self, int index1, int index2)
	cdef void _sift_up(self, int index)
	cdef void _sift_down(self, int index)
	cdef queue_entry_type c_pop(self)
	cdef void c_change_score(self, item_type item, priority_type_ptr c_new_score)
	cdef priority_type_ptr c_get_score_by_item(self, item_type item)
	cdef int size(self)
	cdef bool c_is_empty(self)
