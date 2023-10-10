---
Title: Startup
Weight: 1
---


General Startup procedure of the Octopus Code
=============================================

The flow chart for most calculations has many steps in common, most of them related to setting up the basic data structures of the code.


The 'main' program only performs tasks which are independent of the
actual calculation and the physical system:

[`main/main.F90:`](https://gitlab.com/octopus-code/octopus/blob/develop/src/main/main.F90)

In a strongly simplified way, we have:
```Fortran
program main

  [...] 

  ! "constructors":           ! start code components

  call global_init()               ! initialize the mpi, clocks, etc.
  call parser_init()               ! initialize the input parser
  call messages_init()             ! initialize the message system
  call walltimer_init()            ! initialize the timer module
  call io_init()                   ! initialize the I/O subsystem
  call calc_mode_par_init()        ! initialize parallelization strategy
  call profiling_init()            ! initialize and start the profiling system

  call run(inp_calc_mode)          ! pass control to the 'actual code' running the calculation

  ! "destructors":            ! stop code components in reverse order

  call profiling_end()
  call calc_mode_par_end()
  call io_end()
  call walltimer_end() 
  call messages_end()
  call parser_end() 
  call global_end()

end programme
```


The actual calculation is started from the routine `run()`, which initialized the data structures representing the actual system and calculation mode, before starting the corresponding calculation:

[`main/run.F90`](https://gitlab.com/octopus-code/octopus/blob/develop/src/main/run.F90) 
defines the module `run_oct_m`:

This module does not contain own data (apart from constants).

{{% expand "Definition of run()" %}}
```Fortran
#include_subroutine run
```
{{% /expand %}}


`scf/ground_state.F90`:

```Fortran
#include_subroutine ground_state_run
```

```Fortran
#include_subroutine ground_state_run_legacy
```

```Fortran
#include_subroutine electrons_ground_state_run
```


