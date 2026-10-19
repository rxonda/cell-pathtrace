# Run with: systemsim -cell -n -f /project/docker/smoke.tcl
source $env(LIB_DIR)/cell/systemsim.tcl
fast_wait_for_linux_prompt
mysim console create UART0 in string "uname -m"
mysim console create UART0 in string "mount | grep spufs"
mysim console create UART0 in string "command -v gcc; command -v spu-gcc; command -v ppu-gcc"
mysim console create UART0 in string "callthru exit"
mysim go
exit
