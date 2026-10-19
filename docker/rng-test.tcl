source $env(LIB_DIR)/cell/systemsim.tcl
fast_wait_for_linux_prompt
mysim console create UART0 in string "callthru source /build/rng-test-spu > rng-test-spu"
mysim console create UART0 in string "callthru source /build/rng-test-ppu > rng-test-ppu"
mysim console create UART0 in string "chmod +x rng-test-spu rng-test-ppu"
mysim console create UART0 in string "./rng-test-ppu"
mysim console create UART0 in string "callthru exit"
mysim go
exit
