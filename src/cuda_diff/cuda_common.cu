// cuda_common.cu — host I/O helpers used by every CLI driver.
#include "cuda_common.h"

float *read_floats_or_die(const char *path, size_t n_floats) {
    FILE *f = std::fopen(path, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); std::exit(1); }
    float *buf = (float *)std::malloc(n_floats * sizeof(float));
    if (std::fread(buf, sizeof(float), n_floats, f) != n_floats) {
        fprintf(stderr, "short read on %s\n", path); std::exit(1);
    }
    std::fclose(f);
    return buf;
}

void write_floats_or_die(const char *path, const float *buf,
                         size_t n_floats) {
    FILE *f = std::fopen(path, "wb");
    if (!f) { fprintf(stderr, "cannot open %s for write\n", path);
              std::exit(1); }
    std::fwrite(buf, sizeof(float), n_floats, f);
    std::fclose(f);
}
