% {\rtf1\ansi\ansicpg1252\deff0\deflang1033{\fonttbl{\f0\froman\fcharset0 Times New Roman;}{\f1\fswiss\fcharset0 Courier New;}{\f2\fswiss\fcharset0 Arial;}}
% {\*\generator Msftedit 5.41.15.1507;}\viewkind4\uc1\pard\sb100\sa100\tx2520\f0\fs24\par
% \pard\tx2520\f1 
function  c = fchcode2(b, conn, dir) %\par
%FCHCODE Computes the Freeman chain code of a boundary. \par
%   C = FCHCODE(B) computes the 8-connected Freeman chain code of a \par
%   set of 2-D coordinate pairs contained in B, an np-by-2 array. C \par
%   is a structure with the following fields:  \par
% \par
%     c.fcc    = Freeman chain code (1-by-np) \par
%     c.diff   = First difference of code c.fcc (1-by-np) \par
%     c.mm     = Integer of minimum magnitude from c.fcc (1-by-np) \par
%     c.diffmm = First difference of code c.mm (1-by-np) \par
%     c.x0y0   = Coordinates where the code starts (1-by-2)  \par
% \par
%   C = FCHCODE(B, CONN) produces the same outputs as above, but \par
%   with the code connectivity specified in CONN. CONN can be 8 for \par
%   an 8-connected chain code, or CONN can be 4 for a 4-connected \par
%   chain code. Specifying CONN=4 is valid only if the input \par
%   sequence, B, contains transitions with values 0, 2, 4, and 6, \par
%   exclusively. \par
%        \par
%   C = FHCODE(B, CONN, DIR) produces the same outputs as above, but, \par
%   in addition, the desired code direction is specified. Values for \par
%   DIR can be:  \par
% \par
%     'same'      Same as the order of the sequence of points in b. \par
%                 This is the default. \par
% \par
%     'reverse'   Outputs the code in the direction opposite to the  \par
%                 direction of the points in B.  The starting point  \par
%                 for each DIR is the same. \par
% \par
%   The elements of B are assumed to correspond to a 1-pixel-thick, \par
%   fully-connected, closed boundary. B cannot contain duplicate \par
%   coordinate pairs, except in the first and last positions, which \par
%   is a common feature of boundary tracing programs.  \par
% \par
%   FREEMAN CHAIN CODE REPRESENTATION \par
%   The table on the left shows the 8-connected Freeman chain codes  \par
%   corresponding to allowed deltax, deltay pairs. An 8-chain is \par
%   converted to a 4-chain if (1) if conn = 4; and (2) only \par
%   transitions 0, 2, 4, and 6 occur in the 8-code.  Note that \par
%   dividing 0, 2, 4, and 6 by 2 produce the 4-code.  \par
% \par
%       -----------------------  ---------------- \par
%       deltax | deltay | 8-code  corresp 4-code \par
%       -----------------------  ---------------- \par
%         0        1       0            0 \par
%        -1        1       1 \par
%        -1        0       2            1 \par
%        -1       -1       3 \par
%         0       -1       4            2 \par
%         1       -1       5 \par
%         1        0       6            3 \par
%         1        1       7 \par
%       -----------------------  ---------------- \par
% \par
%   The formula z = 4*(deltax + 2) + (deltay + 2) gives the following \par
%   sequence corresponding to rows 1-8 in the preceding table: z = \par
%   11,7,6,5,9,13,14,15. These values can be used as indices into the \par
%   table, improving the speed of computing the chain code. The \par
%   preceding formula is not unique, but it is based on the smallest \par
%   integers (4 and 2) that are powers of 2.  \par
%  \par
%   Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins \par
%   Digital Image Processing Using MATLAB, Prentice-Hall, 2004 \par
%   $Revision: 1.6 $  $Date: 2003/11/21 14:34:49 $ \par
%  \par
% Preliminaries. \par

if nargin == 1  % \par
   dir = 'same';  %\par
   conn = 8; %\par
elseif nargin == 2 % \par
   dir = 'same'; %\par
elseif nargin == 3   % \par
   % Nothing to do here. \par
else % \par
   error('Incorrect number of inputs.') %\par
end %\par
[np, nc] = size(b);% \par
if np < nc % \par
   error('B must be of size np-by-2.');  %\par
end %\par
 %\par
% Some boundary tracing programs, such as boundaries.m, output a \par
% sequence in which the coordinates of the first and last points are \par
% the same. If this is the case, eliminate the last point.  \par
if isequal(b(1, :), b(np, :)) %\par
   np = np -1; 
   %\par
   b = b(1:np, :); 
   %\par
end %\par

 %\par
% Build the code table using the single indices from the formula  \par
% for z given above: \par
C(11)=0; C(7)=1; C(6)=2; C(5)=3; C(9)=4;% \par
C(13)=5; C(14)=6; C(15)=7; %\par
 %\par
% End of Preliminaries. \par
 %\par
% Begin processing. \par
x0 = b(1, 1);
%\par
y0 = b(1, 2); 
%\par
c.x0y0 = [x0, y0];
%\par
% \par
% Make sure the coordinates are organized sequentially: \par
% Get the deltax and deltay between successive points in b. The  \par
% last row of a is the first row of b. \par
a = circshift(b, [-1, 0]); 
 
% DEL = a - b is an nr-by-2 matrix in which the rows contain the \par
% deltax and deltay between successive points in b. The two  \par
% components in the kth row of matrix DEL are deltax and deltay  \par
% between point (xk, yk) and (xk+1, yk+1).  The last row of DEL  \par
% contains the deltax and deltay between (xnr, ynr) and (x1, y1), \par
% (i.e., between the last and first points in b). \par
DEL = a - b; 
 
% If the abs value of either (or both) components of a pair  \par
% (deltax, deltay) is greater than 1, then by definition the curve  \par
% is broken (or the points are out of order), and the program  \par
% terminates. \par
if any(abs(DEL(:, 1)) > 1) | any(abs(DEL(:, 2)) > 1); 
   error('The input curve is broken or points are out of order.') 
end 
 
% Create a single index vector using the formula described above. \par
z = 4*(DEL(:, 1) + 2) + (DEL(:, 2) + 2); 
 
% Use the index to map into the table. The following are \par
% the Freeman 8-chain codes, organized in a 1-by-np array. \par
fcc = C(z); 
 
% Check if direction of code sequence needs to be reversed. \par
if strcmp(dir, 'reverse') 
   fcc = coderev(fcc); % See below for function coderev. \par
end 
 
% If 4-connectivity is specified, check that all components \par
% of fcc are 0, 2, 4, or 6. \par
if conn == 4 
   val = find(fcc == 1 | fcc == 3 | fcc == 5 | fcc ==7 ); 
   if isempty(val) 
      fcc = fcc./2; 
   else 
      warning('The specified 4-connected code cannot be satisfied.') 
   end 
end 
 
% Freeman chain code for structure output. \par
c.fcc = fcc; 
 
% Obtain the first difference of fcc. \par
c.diff = codediff(fcc,conn); % See below for function codediff. \par
 
% Obtain code of the integer of minimum magnitude. \par
c.mm = minmag(fcc);
% See below for function minmag. \par
 
% Obtain the first difference of fcc \par
c.diffmm = codediff(c.mm, conn); 
 
%-------------------------------------------------------------------% \par
function cr = coderev(fcc) 
%   Traverses the sequence of 8-connected Freeman chain code fcc in \par
%   the opposite direction, changing the values of each code \par
%   segment. The starting point is not changed. fcc is a 1-by-np \par
%   array. \par
 
% Flip the array left to right.  This redefines the starting point  \par
% as the last point and reverses the order of "travel" through the  \par
% code. \par
cr = fliplr(fcc); 
 
% Next, obtain the new code values by traversing the code in the  \par
% opposite direction. (0 becomes 4, 1 becomes 5, ... , 5 becomes 1,  \par
% 6 becomes 2, and 7 becomes 3). \par
ind1 = find(0 <= cr & cr <= 3); 
ind2 = find(4 <= cr & cr <= 7); 
cr(ind1) = cr(ind1) + 4; 
cr(ind2) = cr(ind2) - 4; 
 
%-------------------------------------------------------------------% \par
function z = minmag(c) 
%MINMAG Finds the integer of minimum magnitude in a chain code. \par
%   Z = MINMAG(C) finds the integer of minimum magnitude in a given \par
%   4- or 8-connected Freeman chain code, C. The code is assumed to \par
%   be a 1-by-np array. \par
 
% The integer of minimum magnitude starts with min(c), but there  \par
% may be more than one such value. Find them all, \par

I = find(c == min(c)); 
% and shift each one left so that it starts with min(c). \par
J = 0; 
A = zeros(length(I), length(c)); 
for k = I; 
   J = J + 1; 
   A(J, :) = circshift(c,[0 -(k-1)]); 
end 
 
% Matrix A contains all the possible candidates for the integer of \par
% minimum magnitude. Starting with the 2nd column, succesively find \par
% the minima in each column of A. The number of candidates decreases \par
% as the seach moves to the right on A.  This is reflected in the \par
% elements of J.  When length(J)=1, one candidate remains.  This is \par
% the integer of minimum magnitude.   \par
[M, N] = size(A); z=11;
J = (1:M)'; 
for k = 2:N 
   D(1:M, 1) = Inf; 
   D(J, 1) = A(J, k); 
   amin = min(A(J, k)); 
   J = find(D(:, 1) == amin); 
   if length(J)==1 
      z = A(J, :); 
      return 
   end 
end 
     
%-------------------------------------------------------------------% \par
function d = codediff(fcc, conn) 
%CODEDIFF Computes the first difference of a chain code. \par
%   D = CODEDIFF(FCC) computes the first difference of code, FCC. The \par
%   code FCC is treated as a circular sequence, so the last element \par
%   of D is the difference between the last and first elements of \par
%   FCC.  The input code is a 1-by-np vector.  \par
% \par
%   The first difference is found by counting the number of direction \par
%   changes (in a counter-clockwise direction) that separate two \par
%   adjacent elements of the code.  \par
 
sr = circshift(fcc, [0, -1]); % Shift input left by 1 location. \par
delta = sr - fcc; 
d = delta; 
I = find(delta < 0); 
  
type = conn; 
switch type 
case 4 % Code is 4-connected \par
   d(I) = d(I) + 4; 
case 8 % Code is 8-connected \par
   d(I) = d(I) + 8; 
end 
 
%\pard\sb100\sa100\tx2520\f0\line\par
%\pard\tx2520\f2\par
%}
%?