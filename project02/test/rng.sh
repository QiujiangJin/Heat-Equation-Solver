#!./libs/bin/bats

@test "verify solution of 1d 2nd order with Gauss-Seidel method to reference" {
  run ../src/solver ./input/input_12g.dat
    [ "$status" -eq 0 ]
  run h5diff -d 0.01 output.h5 ./ref/ref_12g.h5 "numerical_solution" "analytical_solution"
    [ "$status" -eq 0 ]
  run rm output.h5
}

@test "verify solution of 1d 2nd order with Jacobi  method to reference" {
  run ../src/solver ./input/input_12j.dat
    [ "$status" -eq 0 ]
  run h5diff -d 0.01 output.h5 ./ref/ref_12j.h5 "numerical_solution" "analytical_solution"
    [ "$status" -eq 0 ]
  run rm output.h5
}

@test "verify solution of 1d 4th order with Gauss-Seidel method to reference" {
  run ../src/solver ./input/input_14g.dat
    [ "$status" -eq 0 ]
  run h5diff -d 0.01 output.h5 ./ref/ref_14g.h5 "numerical_solution" "analytical_solution"
    [ "$status" -eq 0 ]
  run rm output.h5
}

@test "verify solution of 2d 2nd order with Gauss-Seidel method to reference" {
  run ../src/solver ./input/input_22g.dat
    [ "$status" -eq 0 ]
  run h5diff -d 0.01 output.h5 ./ref/ref_22g.h5 "numerical_solution" "analytical_solution"
    [ "$status" -eq 0 ]
  run rm output.h5
}

@test "verify solution of 2d 2nd order with Jacobi method to reference" {
  run ../src/solver ./input/input_22j.dat
    [ "$status" -eq 0 ]
  run h5diff -d 0.01 output.h5 ./ref/ref_22j.h5 "numerical_solution" "analytical_solution"
    [ "$status" -eq 0 ]
  run rm output.h5
}

@test "verify solution of 2d 4th order with Gauss-Seidel method to reference" {
  run ../src/solver ./input/input_24g.dat
    [ "$status" -eq 0 ]
  run h5diff -d 0.01 output.h5 ./ref/ref_24g.h5 "numerical_solution" "analytical_solution"
    [ "$status" -eq 0 ]
  run rm output.h5
}
