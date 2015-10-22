// Provide functions to allocate and free memory that is alligned to the
// cache line.

// The following code was taken from:
//   http://stackoverflow.com/questions/1919183/how-to-allocate-and-free-aligned-memory-in-c
//   09. April. 2014

#define ALIGN 64
#define CACHE_LINE_SIZE ALIGN

void *aligned_malloc(size_t size) {
    void *mem = malloc(size+ALIGN+sizeof(void*));
    void **ptr = (void**)((size_t)(((char*)mem)+ALIGN+sizeof(void*)) & ~(ALIGN-1));
    ptr[-1] = mem;
    return ptr;
}

void aligned_free(void *ptr) {
    free(((void**)ptr)[-1]);
}
