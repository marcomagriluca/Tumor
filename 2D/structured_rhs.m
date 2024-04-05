
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

function rhs = structured_rhs (msh, elrhs, ndrhs)
  
  [ix, iy, ii, hastop, hasbot, haslft, hasrgt] = ...
  indices (msh.nx, msh.ny);
  
  rhs = zeros (size (ii));
  hxhy = (msh.hx(:) .* msh.hy(:).').'(:);
  rhs(hastop & hasrgt) += 1*(hxhy / 4) .* elrhs;
  rhs(hastop & haslft) += 1*(hxhy / 4) .* elrhs;
  rhs(hasbot & hasrgt) += 1*(hxhy / 4) .* elrhs;
  rhs(hasbot & haslft) += 1*(hxhy / 4) .* elrhs;

  rhs .*= ndrhs;
  
endfunction


%!demo 
%! msh.nx = 201;
%! msh.ny = 101;
%! x  = linspace (0, 1, msh.nx);
%! y  = linspace (0, 1, msh.ny);
%! [X, Y] = meshgrid (x, y);
%! msh.hx = diff (x);
%! msh.hy = diff (y);
%! tic, A = structured_reaction (msh, 1, 1); assembly_time = toc
%! f = structured_rhs (msh, 1, (X-Y)(:));
%! tic,u = A \ f; solve_time = toc
%! contourf (X, Y, reshape (u, size (X)))
%! figure
%! surf (X, Y, reshape (u, size (X)), 'edgecolor', 'none')
%! set (gca, 'fontsize', 14)
%! colorbar ('fontsize', 14)
%! xlabel x
%! ylabel y
%! axis equal
%! axis tight
