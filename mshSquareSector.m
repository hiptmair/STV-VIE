function mesh = mshSquareSector(N, type)
%+========================================================================+
%|                                                                        |
%|                 OPENMSH - LIBRARY FOR MESH MANAGEMENT                  |
%|           openMsh is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : mshSquare.m                                   |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Build uniform mesh for a square               |
%|  `---'  |                                                              |
%+========================================================================+

C = 1;
if type == 1
X = C* [1 0 0;
        1 1 0;
        0 1 0;
       -1 1 0;
       -1 0 0;
       -1 -1 0;
        0 -1 0;
        1 -1 0;
        0 0 0];

Points = X;
Elts = [1 2 3;
        1 3 9;
        3 4 9;
        4 5 9;
        5 6 7;
        5 7 9;
        7 8 9];

elseif type == 2
X = C* [1 0 0;
        1 1 0;
        0 1 0;
       -1 1 0;
       -1 0 0;
       -1 -1 0;
        0 0 0;
        0 -1 0];

Points = X;
Elts = [1 2 7;
        2 3 7;
        3 4 5;
        3 5 7;
        5 6 7;
        6 8 7];
    
elseif type == 3

X = C* [1 -1 0;
        1 1 0;
       -1 1 0;
       -1 -1 0];

Points = X;
Elts = [4 1 2;
        2 3 4];
end
    
    
    
mesh = msh(Points,Elts);
% 
if N > 0
for j = 1:N
    
    mesh = mesh.refine;
    
end
end


end
