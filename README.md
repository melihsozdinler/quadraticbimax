# Quadratic Bimax

A C++17 implementation of the Quadratic Bimax algorithm for finding maximum biclusters in binary matrices.

## Description

This project implements the Quadratic Bimax algorithm, which is used to find maximum biclusters in binary matrices. A bicluster is a submatrix where all elements are 1s, and the algorithm aims to find the largest such submatrix.

## Requirements

- C++17 compatible compiler
- CMake 3.10 or higher
- Git (for cloning the repository)

## Building the Project

1. Clone the repository:
```bash
git clone https://github.com/yourusername/quadraticbimax.git
cd quadraticbimax
```

2. Create a build directory and build the project:
```bash
mkdir build
cd build
cmake ..
cmake --build .
```

The executable will be created in the `build/bin` directory.

## Usage

The program takes the following command-line arguments:

```bash
./QuadraticBimax <input_file> [output_file] [status_file]
```

- `input_file`: Required. The input binary matrix file
- `output_file`: Optional. File to store the results (default: BIMAXResult.txt)
- `status_file`: Optional. File to store status information (default: statusBimax.txt)

### Input File Format

The input file should be a text file containing:
1. Number of rows
2. Number of columns
3. Minimum number of rows for bicluster
4. Minimum number of columns for bicluster
5. The binary matrix (0s and 1s)

Example:
```
4 4 2 2
1 0 1 0
0 1 1 0
1 1 0 1
0 0 1 1
```

### Output Format

The program generates two output files:

1. Result file (default: BIMAXResult.txt):
   - Contains the dimensions and elements of the found biclusters
   - Format: `<rows> <columns>` followed by row and column indices

2. Status file (default: statusBimax.txt):
   - Contains summary statistics
   - Format: `<number_of_biclusters> <max_rows> <max_columns>`

## Implementation Details

The implementation uses:
- Bit vectors for efficient storage and operations
- Recursive divide-and-conquer approach
- C++17 features for modern C++ programming

## License

[Add your license information here]

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Author

[Your Name]

## Acknowledgments

- Original Bimax algorithm authors
- Contributors and maintainers 