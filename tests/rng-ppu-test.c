#include <stdio.h>
#include <libspe2.h>

int main(void) {
    spe_program_handle_t *image = spe_image_open("./rng-test-spu");
    spe_context_ptr_t context;
    spe_stop_info_t stop = {0};
    unsigned int entry = SPE_DEFAULT_ENTRY;
    int result;
    if (!image) return 1;
    context = spe_context_create(0, NULL);
    if (!context) return 1;
    if (spe_program_load(context, image)) return 1;
    spe_image_close(image);
    result = spe_context_run(context, &entry, 0, NULL, NULL, &stop);
    if (result < 0 || stop.stop_reason != SPE_EXIT ||
        stop.result.spe_exit_code != 0) {
        fprintf(stderr, "SPU RNG test failed\n");
        spe_context_destroy(context);
        return 1;
    }
    if (spe_context_destroy(context)) return 1;
    puts("RNG launcher passed");
    return 0;
}
