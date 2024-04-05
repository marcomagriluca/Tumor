
#######
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
#######


function A = structured_laplacian (msh, mu, vth)

  N  = msh.nx*msh.ny;
  D  = mu.*vth;
  
  [ix, iy, ii, hastop, hasbot, haslft, hasrgt] = ...
  indices (msh.nx, msh.ny);

  fltop = flbot = flrgt = fllft = zeros (size (ii));

  hx_over_hy = (msh.hx(:) ./ msh.hy(:).').'(:);
  fltop(hastop & hasrgt) += D .* hx_over_hy / 2;
  fltop(hastop & haslft) += D .* hx_over_hy / 2;

  flbot(hasbot & hasrgt) += D .* hx_over_hy / 2;
  flbot(hasbot & haslft) += D .* hx_over_hy / 2;

  flrgt(hasrgt & hastop) += D ./ hx_over_hy / 2;
  flrgt(hasrgt & hasbot) += D ./ hx_over_hy / 2;

  fllft(haslft & hastop) += D ./ hx_over_hy / 2;
  fllft(haslft & hasbot) += D ./ hx_over_hy / 2;

  A = flux_assembly (N, ii,
                     msh.nx, msh.ny,
                     hastop, hasbot,
                     haslft, hasrgt,
                     fltop, flbot,
                     fllft, flrgt);
  
endfunction

%!demo 
%! msh.nx = 201;
%! msh.ny = 201;
%! x  = linspace (0, 1, msh.nx);
%! y  = linspace (0, 1, msh.ny);
%! [X, Y] = meshgrid (x, y);
%! msh.hx = diff (x);
%! msh.hy = diff (y);
%! tic, A = structured_laplacian (msh, 1, 1); assembly_time = toc
%! f = structured_rhs (msh,  1, 2*(X .* (1 - X) + Y .* (1 - Y))(:));
%! dnodes = boundary_nodes (msh.nx, msh.ny, 1:4);
%! inodes = setdiff (1:(msh.nx*msh.ny), dnodes);
%! u = zeros (msh.nx*msh.ny, 1);
%! tic,u(inodes) = A(inodes,inodes) \ f(inodes); solve_time = toc
%!
%! subplot (1, 2, 1)
%! contourf (X, Y,  (X .* (1 - X) .* Y .* (1 - Y)))
%! set (gca, 'fontsize', 14)
%! colorbar ('fontsize', 14)
%! xlabel x
%! ylabel y
%! axis equal
%! axis tight
%! title exact
%!
%! subplot (1, 2, 2)
%! contourf (X, Y, reshape (u, size (X)))
%! set (gca, 'fontsize', 14)
%! colorbar ('fontsize', 14)
%! xlabel x
%! ylabel y
%! axis equal
%! axis tight
%! title numerical
