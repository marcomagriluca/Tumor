
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

function [Jx, Jy] = structured_advection_diffusion_flux (msh, n, psi, mu, vth)

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

