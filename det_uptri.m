## Copyright (C) 2019 Heherson Domael <h.domael@outlook.com>
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {[@var{Adet}, @var{At}] =} det_uptri (@var{A})
## Compute the determinant of A using upper triangular matrix method.
##
## @var{Adet} returns the determinant of A.
##
## @var{At} returns a table of column-by-column elementary row operations
## performed to transform A into an upper triangular matrix.
## @seealso{det}
## @end deftypefn

## Author: Heherson Domael <Heherson Domael@LAPTOP-MQALPEGV>
## Created: 2019-01-13

function [Adet, At] = det_uptri (A)
  if rows(A) != columns(A)
    sprintf("Input matrix must be square! Determinant does not exist.")
    Adet = "N/A";
    At = "N/A";
  else
    Adet = 1;
    for col=1:columns(A)
      for row=1:rows(A)
        if row<=col
          continue
        else
          A(row,:) -= A(row,col)*A(col,:)/A(col,col);
        endif
      endfor
      Adet *= A(col,col);
      At(col*rows(A)-rows(A)+1:col*rows(A),:) = A;
    endfor
  endif
endfunction
