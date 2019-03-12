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
## @deftypefn {} {[@var{Adet}, @var{At}] =} det_pivotala11(@var{A})
## Compute the determinant of A of size n x n using Pivotal Method
## selecting A(1,1) as pivot element, converting it to 1.
## 
## @var{Adet} returns the determinant of A.
##
## @var{At} returns the matrix obtained using reduction method using A(1,1) as
## its pivot element to create n-1 x n-1 matrix.
##
## @seealso{det}
## @end deftypefn

## Author: Heherson Domael <Heherson Domael@LAPTOP-MQALPEGV>
## Created: 2019-03-12

function [Adet, At] = det_pivotala11 (A)
    if rows(A) != columns(A)
      sprintf("Input matrix must be square! Determinant does not exist.")
      Adet = "N/A";
      At = "N/A";
    else
      rowref  = 1;
      colref  = 1;
      elemref = A(rowref,colref);
      for ctr=1:length(A)
        rownew = 1;
        for row=1:length(A)
          colnew = 1;
          if row==rowref
            continue
          endif
        for col=1:length(A)
          if col==colref
            continue
          endif
          At(rownew,colnew)=A(row,col)-(A(rowref,col)/elemref)*A(row,colref);
          colnew++;
        endfor
        rownew++;
      endfor
    endfor
    endif
    Adet = elemref*(-1)^(rowref+colref)*det(At);
endfunction
