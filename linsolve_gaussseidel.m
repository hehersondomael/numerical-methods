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
## @deftypefn {} {[@var{At}, @var{Aaug}, @var{Adom}] =} gaussseidel (@var{Acoeff}, @var{Aconst}, @var{Ainit})
## 
## Solve using Gauss-Seidel method the systems of linear equations with equal number of equations and
## unknowns. The input matrix must be diagonally dominant, thus, non-diagonally dominant matrix A
## is automatically arranged.
## 
## @var{At} returns the solutions to linear equations with equal number of
## equations and unkowns.
## 
## @var{Aaug} returns the augmented table from using the x values obtained
## from the previous iteration for the present one.
## 
## @var{Adom} returns the arranged diagonally dominant matrix including the
## constants of the arranged systems of linear equations.
## @end deftypefn

## Author: Heherson Domael <Heherson Domael@LAPTOP-MQALPEGV>
## Created: 2020-02-01

function [At, Aaug, Adom] = linsolve_gaussseidel (Acoeff, Aconst, Ainit)
  Acoeff = [Acoeff Aconst];
  ref = ptr = 1;
  diagDomCount  = 0;

  while 1
    sumNonPivot = 0;
    for i=1:rows(Acoeff)
      if i==ref
        continue
      else
        sumNonPivot+=abs(Acoeff(ref,i));
      endif
    endfor
    if abs(Acoeff(ref,ref))>=sumNonPivot
      ptr = 1;
      ref++;
      diagDomCount++;
      if ref>rows(Acoeff)
        break
      endif
      continue
    else
      if ref+ptr>rows(Acoeff)
        break
      endif
      temp         = Acoeff(ref,:);
      Acoeff(ref,:)     = Acoeff(ref+ptr,:);
      Acoeff(ref+ptr,:) = temp;
      ptr++;
      if ref+ptr>rows(Acoeff)+1
        break
      endif
    endif
  endwhile

  if diagDomCount!=rows(Acoeff) || rows(Acoeff)!=columns(Acoeff)-1
    Aconst   = "Can't solve";
    Aaug = "Can't solve";
    Adom = "N/A";
  else
    Adom = Acoeff;
    curr_ctr  = 1;
    curr = Ainit';
    do
      Aaug(curr_ctr,:) = [curr_ctr-1 curr];
      for row=1:rows(Acoeff)
        value = Acoeff(row,columns(Acoeff));
        for col=1:rows(Acoeff)
          if row==col
            continue
          else
            value-=Acoeff(row,col)*curr(1,col);
          endif
        endfor
        Aaug(curr_ctr+1,1) = curr_ctr;
        curr(1,row) = Aaug(curr_ctr+1,row+1) = value/Acoeff(row,row);
      endfor
      curr_ctr++;
    until abs(Aaug(curr_ctr,2:columns(Acoeff))-Aaug(curr_ctr-1,2:columns(Acoeff)))<=0.00005 == ones(1,rows(Acoeff))
    At = Aaug(rows(Aaug),2:columns(Aaug))';
  endif
endfunction