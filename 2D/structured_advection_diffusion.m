
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

function A = structured_advection_diffusion (msh, psi, mu, vth)

  N  = msh.nx*msh.ny;
  D  = mu.*vth;
  
  [ix, iy, ii, hastop, hasbot, haslft, hasrgt] = ...
  indices (msh.nx, msh.ny);
  
  fltop = flbot = flrgt = fllft = zeros (size (ii));

  hx_over_hy = (msh.hx(:) ./ msh.hy(:).').'(:);
  fltop(hastop & hasrgt) += ...
  (D .* hx_over_hy / 2) .* ...
  bimu_bernoulli ((psi(ii(hastop & hasrgt)-1) - psi(ii(hastop & hasrgt)))/vth);
  fltop(hastop & haslft) += ...
  (D .* hx_over_hy / 2) .* ...
  bimu_bernoulli ((psi(ii(hastop & haslft)-1) - psi(ii(hastop & haslft)))/vth);

  flbot(hasbot & hasrgt) += ...
  (D .* hx_over_hy / 2) .* ...
  bimu_bernoulli ((psi(ii(hasbot & hasrgt)+1) - psi(ii(hasbot & hasrgt)))/vth);
  flbot(hasbot & haslft) += ...
  (D .* hx_over_hy / 2) .* ...
  bimu_bernoulli ((psi(ii(hasbot & haslft)+1) - psi(ii(hasbot & haslft)))/vth);

  flrgt(hasrgt & hastop) += ...
  (D ./ hx_over_hy / 2) .* ...
  bimu_bernoulli ((psi(ii(hasrgt & hastop)+msh.ny) - psi(ii(hasrgt & hastop)))/vth);
  flrgt(hasrgt & hasbot) += ...
  (D ./ hx_over_hy / 2) .* ...
  bimu_bernoulli ((psi(ii(hasrgt & hasbot)+msh.ny) - psi(ii(hasrgt & hasbot)))/vth);

  fllft(haslft & hastop) += ...
  (D ./ hx_over_hy / 2) .* ...
  bimu_bernoulli ((psi(ii(haslft & hastop)-msh.ny) - psi(ii(haslft & hastop)))/vth);
  fllft(haslft & hasbot) += ...
  (D ./ hx_over_hy / 2) .* ...
  bimu_bernoulli ((psi(ii(haslft & hasbot)-msh.ny) - psi(ii(haslft & hasbot)))/vth);

  
  A = flux_assembly (N, ii, msh.nx, msh.ny, hastop, hasbot,
                     haslft, hasrgt, fltop,
                     flbot, fllft, flrgt);
  
endfunction

%!demo 
%! msh.nx = 501;
%! msh.ny = 401;
%! x  = linspace (0, 1, msh.nx);
%! y  = linspace (0, 1, msh.ny);
%! [X, Y] = meshgrid (x, y);
%! msh.hx = diff (x);
%! msh.hy = diff (y);
%! tic, A = structured_advection_diffusion (msh, (X+Y), 1, 1e-6); assembly_time = toc
%! f = structured_rhs (msh, ones((msh.nx-1)*(msh.ny-1),1), ones(msh.nx*msh.ny,1));
%! dnodes = boundary_nodes (msh.nx, msh.ny, 1:4);
%! inodes = setdiff (1:(msh.nx*msh.ny), dnodes);
%! u = zeros (msh.nx*msh.ny, 1);
%! u((X == 0 & Y < .3) | (Y == 0)) = 1;
%! tic,u(inodes) = A(inodes,inodes) \ (f(inodes) - A(inodes, dnodes) * u(dnodes)); solve_time = toc
%! contourf (x, y, reshape (u, msh.ny, msh.nx))
%! set (gca, 'fontsize', 14)
%! xlabel x
%! ylabel y
%! colorbar ('fontsize', 14)
%! axis equal
%! axis tight

