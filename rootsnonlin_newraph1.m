## Copyright (C) 2020 Heherson Domael
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
## @deftypefn {} {@var{Aroot}, @var{At} =} rootsnonlin_newraph1 (@var{A}, @var{Alow}, @var{Aup})
## Solve for the nearest root of nonlinear equation A on the [Alow, Aup] interval
## using Newton-Raphson 1st method (method of tangents).
## 
## @var{Aroot} returns the nearest root of A on [Alow, Aup].
## 
## @var{At} returns the table generated by computing the nearest root of A on [Alow, Aup].
## 
## @end deftypefn

## Author: Heherson Domael <Heherson Domael@LAPTOP-MQALPEGV>
## Created: 2019-04-29

function [Aroot, At] = rootsnonlin_newraph1 (A, Alow)
  pkg load optim;
  i=1;
  while 1
    Aup = Alow - (A(Alow)/deriv(A, Alow));
    flow=A(Alow);
    fd=deriv(A, Alow);
    fup=A(Aup);
    At(i,:)=[i Alow Aup flow fd fup];
    if abs(flow)<=0.00005
      Aroot=Alow;
      break
    endif
    if abs(fup)<=0.00005
      Aroot=Aup;
      break
    endif
    Alow=Aup;
    i++;
  endwhile
endfunction