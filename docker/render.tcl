# Host paths /project and /build are supplied by docker/run-render.sh.
source $env(LIB_DIR)/cell/systemsim.tcl
fast_wait_for_linux_prompt
mysim console create UART0 in string "callthru source /build/forwardFloatSPU > forwardFloatSPU"
mysim console create UART0 in string "callthru source /build/smallptPPU > smallptPPU"
mysim console create UART0 in string "chmod +x forwardFloatSPU smallptPPU"
mysim console create UART0 in string "./smallptPPU $env(CELL_SAMPLES) $env(CELL_SPUS) $env(CELL_SEED)"
mysim console create UART0 in string "if test -s image.ppm; then callthru sink /build/image.ppm < image.ppm; fi"
mysim console create UART0 in string "callthru exit"
mysim go
exit
