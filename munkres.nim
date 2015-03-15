#
type
  Matrix[W, H: static[int]] =
    array[0..W-1, array[0..H-1, int]]
type
  Vector[L: static[int]] =
    array[L, int]

proc step_one[W, H](matrix: var Matrix[W, H]): int =
  var min_in_row = 0
  for r in 0..high(matrix):
    min_in_row = matrix[r][0]
    for c in 0..high(matrix[0]):
      if min_in_row > matrix[r][c]:
        min_in_row = matrix[r][c]
    for c in 0..high(matrix[0]):
      matrix[r][c] -= min_in_row
  result = 2

proc step_two[W, H](matrix: Matrix[W, H], rowCover: var Vector[W], colCover: var Vector[H], m: var Matrix[W, H]): int =
  for r in 0..high(matrix):
    for c in 0..high(matrix[0]):
      if matrix[r][c] == 0 and rowCover[r] == 0 and  colCover[c] == 0:
        m[r][c] = 1
        rowCover[r] = 1
        colCover[c] = 1
  for r in 0..high(matrix):
    rowCover[r] = 0
  for c in 0..high(matrix[0]):
    colCover[c] = 0
  result = 3

proc step_three[W, H](m: Matrix[W, H], colCover: var Vector[H]): int =
  var colcount = 0;
  for r in 0..high(m):
    for c in 0..high(m[0]):
      if m[r][c] == 1:
        colCover[c] = 1
  for c in 0..high(m[0]):
    if colCover[c] == 1:
      colcount += 1
  result = 4
  if colcount >= high(m)+1 or colcount >= high(m[0])+1:
    result = 7

proc find_a_zero[W, H](matrix: Matrix[W, H], rowCover: Vector[W], colCover: Vector[H]): array[2, int] =
  var done = false
  var c = 0
  var r = 0
  var row = -1
  var col = -1
  while done != true:
    c = 0
    while true:
      if matrix[r][c] == 0 and rowCover[r] == 0 and colCover[c] == 0:
        row = r
        col = c
        done = true
      c += 1
      if c >= high(matrix[0])+1 or done:
        break
    r += 1
    if r >= high(matrix)+1:
      done = true
  result = [row, col]

proc star_in_row[W, H](m: Matrix[W, H], row: int): bool =
  result = false
  for c in 0..high(m[0]):
    if m[row][c] == 1:
      result = true

proc find_star_in_row[W, H](m: Matrix[W, H], row: int): int =
  result = -1
  for c in 0..high(m[0]):
    if m[row][c] == 1:
      result = c

proc step_four[W, H](matrix: Matrix[W, H], rowCover: var Vector[W], colCover: var Vector[H], m: var Matrix[W, H]): array[3, int] =
  var prow = -1
  var pcol = -1
  var step = -1
  var row = -1
  var col = -1
  var done = false
  
  while done != true:
    let r = find_a_zero(matrix, rowCover, colCover)
    row = r[0]
    col = r[1]
    if row == -1:
      done = true
      step = 6
    else:
      m[row][col] = 2
      if star_in_row(m, row):
        col = find_star_in_row(m, row)
        rowCover[row] = 1
        colCover[col] = 0
      else:
        done = true
        step = 5
        prow = row
        pcol = col
  result = [prow, pcol, step]

proc find_star_in_col[W, H](m: Matrix[W, H], c: int): int =
  result = -1
  for i in 0..high(m):
    if m[i][c] == 1:
      result = i

proc find_prime_in_row[W, H](m: Matrix[W, H], r: int): int =
  result = -1
  for j in 0..high(m[0]):
    if m[r][j] == 2:
      result = j

proc augment_path[W, H, P, Q](m: var Matrix[W, H], path: var Matrix[P, Q]): int =
  for p in 0..high(path):
    if m[path[p][0]][path[p][1]] == 1:
      m[path[p][0]][path[p][1]] = 0
    else:
      m[path[p][0]][path[p][1]] = 1
  result = 0

proc clear_covers[W, H](rowCover: var Vector[W], colCover: var Vector[H]): int =
  for r in 0..high(rowCover):
    rowCover[r] = 0
  for c in 0..high(colCover):
    colCover[c] = 0
  result = 0

proc erase_primes[W, H](m: var Matrix[W, H]): int = 
  for r in 0..high(m):
    for c in 0..high(m[0]):
      if m[r][c] == 2:
        m[r][c] = 0
  result = 0

proc step_five[W, H, P, Q](rowCover: var Vector[W], colCover: var Vector[H], m: var Matrix[W, H], prow: int, pcol: int, path: var Matrix[P, Q]): int =
  var done = false
  var r = -1
  var c = -1

  var path_count = 1
  path[path_count - 1][0] = prow
  path[path_count - 1][1] = pcol
  while done != true:
    r = find_star_in_col(m, path[path_count - 1][1])
    if r > -1:
      path_count += 1
      path[path_count - 1][0] = r
      path[path_count - 1][1] = path[path_count - 2][1]
    else:
      done = true
    if done != true:
      c = find_prime_in_row(m, path[path_count - 1][0])
      path_count += 1
      path[path_count - 1][0] = path[path_count - 2][0]
      path[path_count - 1][1] = c;
  var v = augment_path(m, path)
  v = clear_covers(rowCover, colCover)
  v = erase_primes(m)
  result = 3

proc find_smallest[W, H](matrix: Matrix[W, H], rowCover: Vector[W], colCover: Vector[H]): int =
  result = 10000
  for r in 0..high(matrix):
    for c in 0..high(matrix[0]):
      if rowCover[r] == 0 and colCover[c] == 0:
        if result > matrix[r][c]:
           result = matrix[r][c] 

proc step_six[W, H](matrix: var Matrix[W, H], rowCover: Vector[W], colCover: Vector[H]): int =
  var minval = find_smallest(matrix, rowCover, colCover)
  for r in 0..high(matrix):
    for c in 0..high(matrix[0]):
      if rowCover[r] == 1:
        matrix[r][c] += minval
      if colCover[c] == 0:
        matrix[r][c] -= minval
  result = 4
  
proc runMunkres[W, H](matrix: var Matrix[W, H], max: bool): Matrix[W, H] =
  if max:
    var maxValue = matrix[0][0]
    for i in 0..high(matrix):
      for j in 0..high(matrix[0]):
        if matrix[i][j] > maxValue:
          maxValue = matrix[i][j]
    for i in 0..high(matrix):
      for j in 0..high(matrix[0]):
        matrix[i][j] = maxValue - matrix[i][j]
  ##
  var done = false
  var step = 1
  let row_count : int = high(matrix) + 1
  let col_count : int = high(matrix[0]) + 1
  var rowCover : Vector[high(matrix) + 1]
  var colCover : Vector[high(matrix[0]) + 1]
  var m : Matrix[high(matrix) + 1, high(matrix[0]) + 1]
  var prow = 0
  var pcol = 0
  ##
  var path : Matrix[high(matrix) + high(matrix[0]) + 2, high(matrix[0]) + 1]
  var iteration = 0
  while done != true:
    iteration += 1
    case step:
      of 1:
        step = step_one(matrix)
      of 2:
        step = step_two(matrix, rowCover, colCover, m)
      of 3:
        step = step_three(m, colCover)
      of 4:
        var r = step_four(matrix, rowCover, colCover, m)
        prow = r[0]
        pcol = r[1]
        step = r[2]
      of 5:
        step = step_five(rowCover, colCover, m, prow, pcol, path)
      of 6:
        step = step_six(matrix, rowCover, colCover)
      else:
        done = true
  result = matrix

##
#var c : Matrix[100, 100] # = [[1, 2, 3], [2, 4, 6], [3, 6, 9]]
#for i in 0..high(c):
#  for j in 0..high(c[0]):
#    c[i][j] = (i + 1) * (j + 1)
var c : Matrix[3, 3] = [[1, 2, 3], [2, 4, 6], [3, 6, 9]]
echo "cost"
echo(c[0][0], " ", c[0][1], " ", c[0][2])
echo(c[1][0], " ", c[1][1], " ", c[1][2])
echo(c[2][0], " ", c[2][1], " ", c[2][2])
let o = runMunkres(c, false)
##
echo "assignment"
echo(o[0][0], " ", o[0][1], " ", o[0][2])
echo(o[1][0], " ", o[1][1], " ", o[1][2])
echo(o[2][0], " ", o[2][1], " ", o[2][2])
##


