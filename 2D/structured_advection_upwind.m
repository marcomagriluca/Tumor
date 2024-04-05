
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

function A = structured_advection_upwind (msh, psi, D)

  N  = msh.nx*msh.ny;
  
  [ix, iy, ii, hastop, hasbot, haslft, hasrgt] = ...
  indices (msh.nx, msh.ny);
  
  fltop = flbot = flrgt = fllft = zeros (size (ii));

  hx_over_hy = (msh.hx(:) ./ msh.hy(:).').'(:);
  B = upwind_difference (psi(ii(hastop & hasrgt)-1) - psi(ii(hastop & hasrgt)));
  fltop(hastop & hasrgt) += D .* B .* hx_over_hy / 2;
  B = upwind_difference (psi(ii(hastop & haslft)-1) - psi(ii(hastop & haslft)));
  fltop(hastop & haslft) += D .* B .* hx_over_hy / 2;

  B = upwind_difference (psi(ii(hasbot & hasrgt)+1) - psi(ii(hasbot & hasrgt)));
  flbot(hasbot & hasrgt) += D .* B .* hx_over_hy / 2;
  B = upwind_difference (psi(ii(hasbot & haslft)+1) - psi(ii(hasbot & haslft)));
  flbot(hasbot & haslft) += D .* B .* hx_over_hy / 2;

  B = upwind_difference (psi(ii(hasrgt & hastop)+msh.ny) - psi(ii(hasrgt & hastop)));
  flrgt(hasrgt & hastop) += D .* B ./ hx_over_hy / 2;
  B = upwind_difference (psi(ii(hasrgt & hasbot)+msh.ny) - psi(ii(hasrgt & hasbot)));
  flrgt(hasrgt & hasbot) += D .* B ./ hx_over_hy / 2;

  B = upwind_difference (psi(ii(haslft & hastop)-msh.ny) - psi(ii(haslft & hastop)));
  fllft(haslft & hastop) += D .* B ./ hx_over_hy / 2;
  B = upwind_difference (psi(ii(haslft & hasbot)-msh.ny) - psi(ii(haslft & hasbot)));
  fllft(haslft & hasbot) += D .* B ./ hx_over_hy / 2;

  A = flux_assembly (N, ii,
                     msh.nx, msh.ny,
                     hastop, hasbot,
                     haslft, hasrgt,
                     fltop, flbot,
                     fllft, flrgt);
  
endfunction

function y = upwind_difference (x)
  y = ifelse (x < 0, -x, 0);
endfunction
  
%!demo 
%! msh.nx = 300;
%! msh.ny = 300;
%! x  = linspace (0, 1, msh.nx);
%! y  = linspace (0, 1, msh.ny);
%! [X, Y] = meshgrid (x, y);
%! msh.hx = diff (x);
%! msh.hy = diff (y);
%!
%! tic, A = structured_advection_diffusion (msh, (X+Y), 1, 1e-1); assembly_time = toc
%! f = structured_rhs (msh, ones((msh.nx-1)*(msh.ny-1),1), ones(msh.nx*msh.ny,1));
%! dnodes = boundary_nodes (msh.nx, msh.ny, 1:4);
%! inodes = setdiff (1:(msh.nx*msh.ny), dnodes);
%! u = zeros (msh.nx*msh.ny, 1);
%! u((X == 0 & Y < .3) | (Y == 0)) = 1;
%! tic,u(inodes) = A(inodes,inodes) \ (f(inodes) - A(inodes, dnodes) * u(dnodes)); solve_time = toc
%!
%! tic, A = structured_advection_upwind (msh, (X+Y), 1); assembly_time = toc
%! tic, A += structured_laplacian (msh, 1, 1e-1); assembly_time = toc
%! f = structured_rhs (msh, ones((msh.nx-1)*(msh.ny-1),1), ones(msh.nx*msh.ny,1));
%! dnodes = boundary_nodes (msh.nx, msh.ny, 1:4);
%! inodes = setdiff (1:(msh.nx*msh.ny), dnodes);
%! uu = zeros (msh.nx*msh.ny, 1);
%! uu((X == 0 & Y < .3) | (Y == 0)) = 1;
%! tic,uu(inodes) = A(inodes,inodes) \ (f(inodes) - A(inodes, dnodes) * uu(dnodes)); solve_time = toc
%!
%! subplot (1, 2, 1)
%! contourf (x, y, reshape (u, msh.nx, msh.ny))
%! set (gca, 'fontsize', 14)
%! xlabel x
%! ylabel y
%! colorbar ('fontsize', 14)
%! axis equal
%! axis tight
%! title "FVSG"
%!
%! subplot (1, 2, 2)
%! contourf (x, y, reshape (uu, msh.nx, msh.ny))
%! set (gca, 'fontsize', 14)
%! xlabel x
%! ylabel y
%! colorbar ('fontsize', 14)
%! axis equal
%! axis tight
%! title "UW"
%!
%! figure
%! contourf (x, y, reshape (abs (uu-u), msh.nx, msh.ny))
%! set (gca, 'fontsize', 14)
%! xlabel x
%! ylabel y
%! colorbar ('fontsize', 14)
%! axis equal
%! axis tight
%! title "difference"

