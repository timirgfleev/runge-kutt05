# Runge-Kutta Method Implementation in C++

This project is a C++ implementation of the Runge-Kutta method, a numerical technique used to solve ordinary differential equations.


## Project Structure

The project is organized as follows:
- 'data.txt': contains integration settings for main: X_START, X_END, X_0, Y(X_0), H_START, ERROR_TARGET
- `main.cpp`: It reads input data from `data.txt`, tries to integrate with target accuracy, and writes the output to `rez.txt`.
- `include/`: This directory contains header files for the project.
  - `file_io.h`: Contains functions for reading from and writing to files.
  - `math_function.h`: Contains definition of mathematical function f(x, y) for integration
  - `runge.h`: Contains the implementation of the Runge-Kutta method.
- `test/`: This directory contains test files for the project.
  - `test_1.cpp`: Contains tests for the Runge-Kutta method.
  - `additional_test.cpp`: Contains additional tests for the project.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

You need to have CMake installed on your system to build this project.

### Installing

1. Clone the repository
2. Navigate to the project directory
3. Run the following commands:

```sh
cmake -B build
cmake --build build
```

This will create an executable named `runge`.

### Running the Application

To run the application, execute the following command from project root directory:

```sh
./build/runge
```

### Running the Tests

To run the tests, execute the following command in build directory:

```sh
ctest
```


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details