## Copyright (C) 2020 Heherson Domael <h.domael@outlook.com>
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
## @deftypefn {} {[@var{Ainv}, @var{At}] =} inv_mat (@var{A})
## Compute the inverse of matrix A.
##
## @var{Ainv} returns the inverse of A.
##
## @var{At} returns a table of column-by-column elementary row operations
## performed to transform identity matrix into A^-1.
## @seealso{inv}
## @end deftypefn

## Author: Heherson Domael <Heherson Domael@LAPTOP-MQALPEGV>
## Created: 2019-03-06

function [Ainv, At] = inv_mat (A)
  if det(A)==0
    sprintf("Input matrix is singular! Inverse does not exist.")
    Adet = "N/A";
    At = "N/A";
  elseif rows(A) != columns(A)
    sprintf("Input matrix must be square! Inverse does not exist.")
    Adet = "N/A";
    At = "N/A";
  else
    A = [A eye(rows(A))];
    for col=1:(columns(A)/2)
      for row=1:rows(A)
        if row==col
          # A(row,:) = A(row,:)/A(row,col)
          continue
        else
          A(row,:) -= A(row,col)*A(col,:)/A(col,col);
        endif
        if col==rows(A)
          for i=1:rows(A)
            A(i,:) = (1/A(i,i))*A(i,:);
          endfor
        endif
        At((col*rows(A))-(rows(A)-1):(col*rows(A)),:) = A;
      endfor
    endfor
    Ainv=At(rows(At)-(rows(A)-1):rows(At),(columns(At)-(rows(A)-1)):columns(At));
  endif
endfunction
