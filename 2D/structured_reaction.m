
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

function A = structured_reaction (msh, elrhs, ndrhs)

  rhs = structured_rhs (msh, elrhs, ndrhs);
  A = spdiags (rhs (:), 0, numel (rhs), numel (rhs));
  
endfunction

%!demo 
%! msh.nx = 401;
%! msh.ny = 201;
%! x  = linspace (0, 1, msh.nx);
%! y  = linspace (0, 1, msh.ny);
%! [X, Y] = meshgrid (x, y);
%! msh.hx = diff (x);
%! msh.hy = diff (y);
%! tic, A = structured_laplacian (msh, 1e-5, 1); assembly_time = toc
%! tic, A += structured_reaction (msh, 1, 1); assembly_time = toc
%! issparse (A)
%! f = structured_rhs (msh, 1, 1);
%! dnodes = boundary_nodes (msh.nx, msh.ny, 1:4);
%! inodes = setdiff (1:(msh.nx*msh.ny), dnodes);
%! u = zeros (msh.nx*msh.ny, 1);
%! tic,u(inodes) = A(inodes,inodes) \ f(inodes); solve_time = toc
%! contourf (X, Y, reshape (u, size (X)))
%! figure
%! surf (X, Y, reshape (u, size (X)), 'edgecolor', 'none')
%! set (gca, 'fontsize', 14)
%! colorbar ('fontsize', 14)
%! xlabel x
%! ylabel y
%! axis equal
%! axis tight
