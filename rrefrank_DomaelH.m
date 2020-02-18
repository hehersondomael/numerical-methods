## Copyright (C) 2018 Heherson Domael
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
## @deftypefn {} {@var{retval} =} rrefrank_DomaelH (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Heherson Domael <Heherson Domael@LAPTOP-MQALPEGV>
## Created: 2018-11-30

function [rankA, rrefA, At] = rrefrank_DomaelH (A)
  Atempo=A;
  rankcount=bott=size(A)(1,1);
  rowaug=refrow=refer=ctr=k=1;

  for col=1:(size(A)(1,2))
    refer=ctr;

    if refer > size(A)(1,1)
      if k < size(A)(1,1)
        for m=1:size(A)(1,1)-1
          Atempo(m,:) -= (Atempo(m,size(A)(1,2))/Atempo(size(A)(1,1),size(A)(1,2)))*Atempo(size(A)(1,1),:);
        endfor
        Atempo(size(A)(1,1),:) = Atempo(size(A)(1,1),:)/Atempo(size(A)(1,1),size(A)(1,2));
      endif
      for row=1:bott # zero-row placed at the bottom
        if bott<row
          break
        else
          if isequal((abs(Atempo(row,:)) < 1e-10),ones(1,size(A)(1,2)))
            Atempo(row:bott-1,:)=Atempo(row+1:bott,:);
            Atempo(bott,:)=zeros(1,size(A)(1,2));
            bott--;
          endif
        endif
      endfor
      Aug(rowaug:(rowaug+size(A)(1,1)-1),:) = Atempo;
      rowaug += size(A)(1,1);
      if row >= size(A)(1,1) || col >= size(A)(1,2)
        rrefFinal=Atempo;
      endif
      if refrow<size(A)(1,1)
        refrow++;
      endif
      ctr++;
      break
    endif

    for row = (ctr+1):(size(A)(1,1)) # row with highest absolute valued element swapped to reference row
     if abs(Atempo(refer,col)) < abs(Atempo(row,col))
        refer = row;
      endif
    endfor
    if Atempo(refer,col) == 0 # check if the reference element is 0
      for row=1:bott  # zero-row placed at the bottom
        if isequal(Atempo(row,:),zeros(1,size(A)(1,2)))
          Atempo(row:bott-1,:)=Atempo(row+1:bott,:);
          Atempo(bott,:)=zeros(1,size(A)(1,2));
          bott--;
        endif
      endfor
      Aug(rowaug:(rowaug+size(A)(1,1)-1),:) = Atempo; # apply skipping; copy previously augmented matrix
      rowaug += size(A)(1,1);
      if row >= size(A)(1,1) && col>=size(A)(1,2)
        rrefFinal=Atempo;
      endif
      if refrow<size(A)(1,1)
        refrow++;
      endif
      ctr++;
      continue
    else
      temp = Atempo(refrow,:); # row with highest absolute valued element swapped to reference row
      Atempo(refrow,:) = Atempo(refer,:);
      Atempo(refer,:) = temp;
      for row=1:size(A)(1,1)
        if row == col
          continue
        else
          Atempo(row,:) -= (Atempo(row,col)/Atempo(col,col))*Atempo(col,:);
        endif
      endfor
      Atempo(col,:) = Atempo(col,:)/Atempo(col,col); # reference element normalized to 1 immediately
      for row=1:bott  # zero-row placed at the bottom
        if bott<row
          break
        else
        if isequal(Atempo(row,:),zeros(1,size(A)(1,2)))
          Atempo(row:bott-1,:)=Atempo(row+1:bott,:);
          Atempo(bott,:)=zeros(1,size(A)(1,2));
          bott--;
        endif
        endif
      endfor
      Aug(rowaug:(rowaug+size(A)(1,1)-1),:) = Atempo;
      rowaug += size(A)(1,1);
      if row >= size(A)(1,1) || col>=size(A)(1,2)
          rrefFinal=Atempo;
      endif
      if refrow<size(A)(1,1)
        refrow++;
      endif
      ctr++;
      k++;
    endif
  endfor

  for row=1:size(A)(1,1) # counting of rank
    if isequal(rrefFinal(row,:), zeros(1,size(A)(1,2)))
      rankcount--;
    endif
  endfor

  rankA=rankcount;
  At=Aug;
  rrefA=rrefFinal;

endfunction