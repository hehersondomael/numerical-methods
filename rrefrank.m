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
## @deftypefn {} {@var{rankA}, @var{rrefA}, @var{At} =} rrefrank_DomaelH (@var{A})
## Solve for the reduced row echelon form of A.
## 
## @var{rankA} returns the maximum number of linearly independent column (or
## row) vectors in A called rank.
## 
## @var{rrefA} returns the reduced row echelon form of A.
##
## @var{At} returns a table of column-by-column elementary row operations to
## set A into its reduced row echelon form.
## @seealso{rref}
## @end deftypefn

## Author: Heherson Domael <Heherson Domael@LAPTOP-MQALPEGV>
## Created: 2018-03-31

function [rankA, rrefA, At] = rrefrank (A)
  Atemp=A;
  rankA=bottom=rows(A);
  rowAug=ptr=refRow=ctr=1;

  for col=1:columns(A)
    ptr=ctr;

    if ptr>rows(A)
      if ctr<rows(A)
        for row=1:rows(A)-1
          Atemp(row,:)-=(Atemp(row,columns(A))/Atemp(rows(A),columns(A)))*Atemp(rows(A),:);
        endfor
        Atemp(rows(A),:) = Atemp(rows(A),:)/Atemp(rows(A),columns(A));
      endif

      for row=1:bottom # place zero rows at the bottom
        if bottom<row
          break
        else
          if isequal((abs(Atemp(row,:))<1e-10),ones(1,columns(A)))
            Atemp(row:bottom-1,:)=Atemp(row+1:bottom,:);
            Atemp(bottom,:)=zeros(1,columns(A));
            bottom--;
          endif
        endif
      endfor

      At(rowAug:(rowAug+rows(A)-1),:) = Atemp;
      rowAug+=rows(A);
      if row>=rows(A) || col>=columns(A)
        rrefA=Atemp;
      endif
      if refRow<rows(A)
        refRow++;
      endif
      ctr++;
      break
    endif

    for row=(ctr+1):(rows(A))
     if abs(Atemp(ptr,col))<abs(Atemp(row,col))
        ptr=row;
      endif
    endfor

    if Atemp(ptr,col)==0 # check if the reference element is 0
      for row=1:bottom  # place zero rows at the bottom
        if isequal(Atemp(row,:),zeros(1,columns(A)))
          Atemp(row:bottom-1,:)=Atemp(row+1:bottom,:);
          Atemp(bottom,:)=zeros(1,columns(A));
          bottom--;
        endif
      endfor
      At(rowAug:(rowAug+rows(A)-1),:) = Atemp; # skip and copy previously augmented matrix
      rowAug+=rows(A);
      if row>=rows(A) && col>=columns(A)
        rrefA=Atemp;
      endif
      if refRow<rows(A)
        refRow++;
      endif
      ctr++;
      continue
    else
      temp = Atemp(refRow,:);
      Atemp(refRow,:) = Atemp(ptr,:);
      Atemp(ptr,:) = temp;
      for row=1:rows(A)
        if row == col
          continue
        else
          Atemp(row,:) -= (Atemp(row,col)/Atemp(col,col))*Atemp(col,:);
        endif
      endfor
      Atemp(col,:) = Atemp(col,:)/Atemp(col,col); # normalize reference element to 1 immediately
      for row=1:bottom  # place zero rows at the bottom
        if bottom<row
          break
        else
        if isequal(Atemp(row,:),zeros(1,columns(A)))
          Atemp(row:bottom-1,:)=Atemp(row+1:bottom,:);
          Atemp(bottom,:)=zeros(1,columns(A));
          bottom--;
        endif
        endif
      endfor
      At(rowAug:(rowAug+rows(A)-1),:)=Atemp;
      rowAug+=rows(A);
      if row>=rows(A) || col>=columns(A)
          rrefA=Atemp;
      endif
      if refRow<rows(A)
        refRow++;
      endif
      ctr++;
    endif
  endfor

  for row=1:rows(A)
    if isequal(rrefA(row,:),zeros(1,columns(A)))
      rankA--;
    endif
  endfor
endfunction