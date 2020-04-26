from libc.stdlib cimport rand, RAND_MAX

cdef double p_attention = 0.1

cpdef void set_p_attention(double new_p_attention):
    global p_attention
    p_attention = new_p_attention

cpdef int attention():
    return rand() < RAND_MAX * p_attention

cdef extern from "stdlib.h":
    double drand48()
    void srand48(long int seedval)

cdef extern from "time.h":
    long int time(int)

# srand48(time(0))
srand48(100)
# TODO: this is a seed to reproduce bugs, put to line of code above for
# production
drand48() #This gives a float in range [0,1)
