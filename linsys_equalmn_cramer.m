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
## @deftypefn {} {@var{At} =} linsys_equalmn_cramer (@var{Acoeff}, @var{Aconst})
## Solve using Cramer's rule the given system of linear equations with
## equal number of equations and unknowns.
## 
## @seealso{linsolve}
## @end deftypefn

## Author: Heherson Domael <Heherson Domael@LAPTOP-MQALPEGV>
## Created: 2020-01-21

function At = linsys_equalmn_cramer (Acoeff, Aconst)
  if det(Acoeff) == 0
    sprintf("The system has either no nontrivial solutions or an infinite number of solutions.")
    At = "N/A";
  elseif rows(Acoeff)!=columns(Acoeff)
    sprintf("Input coefficient matrix must be square!")
    At = "N/A";
  elseif rows(Aconst) != rows(Acoeff)
    sprintf("Inputs constant matrix and coefficient matrix do not match!")
    At = "N/A";
  else
    Acoeff_t = Acoeff;
    for col=1:columns(Acoeff)
      Acoeff_t(:,col) = Aconst;
      At(col,:)  = det(Acoeff_t)/det(Acoeff);
      Acoeff_t = Acoeff;
    endfor
  endif
endfunction