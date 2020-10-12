# Sudoku
Sudoku solver

## Build

```
g++ -std=c++17 -O3 sudoku.cpp -o sudoku
```

## Usage

### Normal Sudoku

Pass the grid as a space-separated list of lines. Replace unknowns by any non-numeric symbol.

**NOTE:** Only 9 by 9 grids can be input. Other sizes are currently not implemented. You can hardcode a grid and change the grid size constant for that.

```
./sudoku 86xxxxx2x xxxx6xxxx xx42x16xx 4xx5xx26x xx61x94xx x13xx6xx9 xx54x39xx xxxx1xxxx x8xxxxx74
```

### Killer Sudoku

Add the `-k` switch, assign a symbol to each region (starting at 0, use lowercase letters afterwards) and pass the regions' values as a comma-separated list of numbers.

```
./sudoku -k 001112345 qqrr22345 qqss2gff5 oppshgfe5 ojjihgee6 njkihcdd6 nkkiacc77 nlmaabb77 nlma99988 3,15,22,4,16,15,12,14,17,13,10,15,20,6,17,20,8,17,20,13,6,8,16,27,6,14,25,17,9
```
