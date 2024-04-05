function [ix, iy, ii, hastop, hasbot, haslft, hasrgt] = ...
         indices (nx, ny)

  ## #####
  ##      NUMBERING
  ##      
  ##      i-ny-1           i-1             i+ny-1
  ##      o---------------o---------------o
  ##      |               |               |
  ##      |               |               |
  ##      |               |               |
  ##      |               |               |
  ##      | i-ny          | i             | i+ny
  ##      o---------------o---------------o
  ##      |               |               |
  ##      |               |               |
  ##      |               |               |
  ##      |               |               |
  ##      | i-ny+1        | i+1           | i+ny+1
  ##      o---------------o---------------o
  ##
  ## #####

  
  [ix, iy] = meshgrid (1:nx, 1:ny);
  ii = sub2ind (size(ix), iy, ix);
  ix = ix(:); iy = iy(:); ii = ii(:);

  if (nargout > 3)
    hastop = iy > 1;
    hasbot = iy < ny;

    hasrgt = ix < nx;
    haslft = ix > 1;
  endif
  
endfunction 
