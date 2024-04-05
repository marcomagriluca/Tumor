
#######
##      NUMBERING
##                      1
##      
##      i-ny-1           i-1             i+ny-1
##      o---------------o---------------o
##      |               |               |
##      |               |               |
##      |               |               |
##      |               |               |
##      | i-ny          | i             | i+ny      
##  3   o---------------o---------------o         4
##      |               |               |
##      |               |               |
##      |               |               |
##      |               |               |
##      | i-ny+1        | i+1           | i+ny+1
##      o---------------o---------------o
##
##                      2
#######


function ret = boundary_nodes (nx, ny, idx)

  N  = nx*ny;
  
  [ix, iy, kk] = ...
  indices (nx, ny);

  ret = [];
  for ii = 1 : numel (idx)
    switch idx(ii)
      case 1
        ret = [ret; kk(iy == 1)];
      case 2
        ret = [ret; kk(iy == ny)];
      case 3
        ret = [ret; kk(ix == 1)];
      case 4
        ret = [ret; kk(ix == nx)];
    endswitch
  endfor
  
endfunction



