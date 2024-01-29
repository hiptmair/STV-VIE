function mesh = mshU(N)
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

C = 1/4;
X = C* [-2 2 0;
        -2 1 0;
        -2 0 0;
        -2 -1 0;
        -2 -2 0;
        -1 -2 0;
         0 -2 0;
         1 -2 0;
         2 -2 0;
         2 -1 0;
         1 -1 0;
         0 -1 0;
        -1 -1 0;
        -1 0 0;
        -1 1 0;
         0 1 0;
         1 1 0;
         2 1 0;
         2 2 0;
         1 2 0;
         0 2 0;
        -1 2 0];


Points = X;
Elts = [1 2 15;
        2 3 15;
        3 14 15;
        3 13 14;
        3 4 13;
        4 5 13;
        5 6 13;
        6 7 13;
        7 12 13;
        7 11 12;
        7 8 11;
        8 9 11;
        9 10 11;
        17 18 19;
        17 19 20;
        17 20 21;
        16 17 21;
        15 16 21;
        15 21 22;
        1 15 22];

mesh = msh(Points,Elts);
% 
% for j = 1:floor(log(N)/log(4))
%     
%     mesh = mesh.refine;
%     
% end


end
