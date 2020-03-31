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
## @deftypefn {} {[@var{Adet}, @var{At}] =} det_chioa11 (@var{A})
## Compute the determinant of A of size n x n using Chio's condensation method
## with reference at A(1,1).
##
## @var{Adet} returns the determinant of A.
##
## @var{At} returns the matrix obtained using reduction method using A(1,1) as
## its reference to create n-1 x n-1 matrix.
##
## @seealso{det}
## @end deftypefn

## Author: Heherson Domael <Heherson Domael@LAPTOP-MQALPEGV>
## Created: 2019-03-07

function [Adet, At] = det_chioa11(A)
  if rows(A) != columns(A)
    sprintf("Input matrix must be square! Determinant does not exist.")
    Adet = "N/A";
    At = "N/A";
  else
    col_b = 1;
    for col_a=2:columns(A)
      for row=1:rows(A)-1
        At(row,col_b) = det([A(1,1) A(1,col_a); A(row+1,1) A(row+1,col_a)]);
      endfor
      if col_b > length(A)
        break
      endif
      col_b++;
      row = 1;
    endfor
    Adet = det(At)*A(1,1)^(2-length(A));
  endif
endfunction
