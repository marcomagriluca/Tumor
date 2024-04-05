

## Copyright (C) 2017 Carlo de Falco

clear all
close all
pkg load bim

Lx    = -4;
Rx    = 4;
Ly    = -4;
Ry    = 4;

mu    = 0.5;
nu    = 1;
g     = 30;
PM    = 25;

msh.nx     = 128;
msh.ny     = 128;
N      = msh.nx*msh.ny;
x      = linspace (Lx, Rx, msh.nx);
y      = linspace (Ly, Ry, msh.ny);
[X, Y] = meshgrid (x, y);
msh.hx     = diff (x);
msh.hy     = diff (y);
Xm = .25 * (X(1:end-1, 2:end) + X(1:end-1, 1:end-1) + X(2:end, 2:end) + X(2:end, 1:end-1));
Ym = .25 * (Y(2:end, 1:end-1) + Y(1:end-1, 1:end-1) + Y(2:end, 2:end) + Y(1:end-1, 2:end));


%function um = hmg (u, ii, hastop, hasbot, haslft, hasrgt, Xm)
%  usum = 1./u(ii(hastop & hasrgt)) + 1./u(ii(hastop & haslft)) + ...
%         1./u(ii(hasbot & hasrgt)) + 1./u(ii(hasbot & haslft));
%  um = 4 ./ usum(:);
%endfunction


function um = hmg (u, p, ii, hastop, hasbot, haslft, hasrgt, Xm)
  if p(ii(hasrgt)) >= p(ii(haslft));
    aa=hasrgt;
  else
    aa=haslft;
  endif
  if p(ii(hastop)) >= p(ii(hasbot));
    bb=hastop;
  else
    bb=hasbot;
  endif
  um=u(ii(aa & bb));
endfunction


[ix, iy, ii, hastop, hasbot, haslft, hasrgt] = ...
indices (msh.nx, msh.ny);

hm = @(u, p) hmg (u, p, ii, hastop, hasbot, haslft, hasrgt, Xm);


%%hm0 = @(u) 2 ./  ((1 ./ u(2:end)) + (1 ./ u(1:end-1)));
%%hm = @(u, p) ifelse (p(2:end) >= p(1:end-1), u(2:end), u(1:end-1));

G0  = @(p) PM - p;
G   = @(p) min (max (G0(p), 0), PM);

pressure = @(m, n) ((g + 1)/g) .* (n + m) .^ g;
dpdm     = @(m, n) (g + 1) .* (n + m) .^ (g - 1);
dpdn     = @(m, n) (g + 1) .* (n + m) .^ (g - 1);
diffm    = @(m,n) hm (mu * g * ((g + 1)/g) .* m .* (n + m) .^ (g-1),pressure(m,n));
diffn    = @(m,n) hm (nu * g * ((g + 1)/g) .* n .* (n + m) .^ (g-1),pressure(m,n));
%%diffm0    = @(m,n) hm0 (mu * g * ((g + 1)/g) .* m .* (n + m) .^ (g-1));
%%diffn0    = @(m,n) hm0 (nu * g * ((g + 1)/g) .* n .* (n + m) .^ (g-1));

dG0dm    = @(m, n, p) - dpdm (m, n);
dGdm     = @(m, n, p) ifelse ((G0(p) >= 0) & (G0(p) < PM), dG0dm (m, n, p), 0*p);

dG0dn    = @(m, n, p) - dpdn (m, n);
dGdn     = @(m, n, p) ifelse ((G0(p) >= 0) & (G0(p) < PM), dG0dn (m, n, p), 0*p);



dt    = 1e-3;
T     = 1E-1;

%m  = .1  * exp (-5e-2  * (X.^2 + Y.^2));
%n  = .99 * exp (-5e-3 * (X.^2 + Y.^2));



alpha=pi/4;
x0=0
r=x0+1;
speed =  PM*nu*sqrt(mu)/(r*sqrt(mu)+nu)

n=[];
for i=1:length(x)
  for j=1:length(y)
    x1=x(i)*cos(alpha)+y(j)*sin(alpha);
    if x1 <= x0;
      n(j,i) = 0;
    elseif x1 >= r;
      n(j,i) = 0;
    else
      n(j,i) = ((speed/nu*(r-(x1-x0)))/(g+1)*g)^(1/g);
    endif
  endfor
endfor


m=[];
for i=1:length(x)
  for j=1:length(y)
    x1=x(i)*cos(alpha)+y(j)*sin(alpha);
    if x1<=x0;
      m(j,i) = ((PM + (speed*r /nu - PM)*exp((x1-x0)/sqrt(mu)))/(g+1)*g)^(1/g);
    else
      m(j,i) = 0;
    endif
  endfor
endfor



figure (4)
subplot (2, 1, 1)
contourf (X, Y, reshape (m, size (X))), colorbar
axis equal
subplot (2, 1, 2)
contourf (X, Y, reshape (n, size (X))), colorbar
axis equal
drawnow

mass   = structured_reaction (msh, 1, 1);

stiffm = @(m, n) structured_laplacian (msh, diffm (m,n), 1);
stiffn = @(m, n) structured_laplacian (msh, diffn (m,n), 1);

function [res, jac] = fun (u, mold, nold, nx, ny, N, ...
                           stiffm, stiffn, mass, ...
                           pressure, G, dGdm, dGdn, dt)


  mest = u(1:N);
  nest = u(N+1:end);

  Sm  = stiffm (mest, nest);
  Sn  = stiffn (mest, nest);
  pp  = pressure (mest, nest);
  GG  = G (pp);
  Amm = (mass + dt * (Sm - mass .* sparse (diag (GG))));
  Amn = dt * Sm;
  bm  = mass * mold(:);

  Ann = (mass + dt * Sn);
  Anm = dt * Sn;
  bn  = mass * nold(:);

  A = [Amm, Amn;
       Anm, Ann];

  b = [bm; bn];

  res = A*u(:)-b;

  cfm = mest .* dGdm (mest, nest, pp);
  cfn = mest .* dGdn (mest, nest, pp);
  GGm = mass .* sparse (diag (cfm));
  GGn = mass .* sparse (diag (cfn));
  Amm += (-dt) * GGm;
  Amn += (-dt) * GGn;

  jac = [Amm, Amn;
         Anm, Ann];

endfunction

nt       = 100;
tsave    = linspace (0, T, nt);
ii       = 2;
t        = 0;
tstore   = t;

msave = mvold = mold = m = m(:);
nsave = nvold = nold = n = n(:);

for its = 2 : nt

  while (t < tsave(its))

    if (t + dt > tsave(its))
      dt = tsave(its) - t;
    endif

    printf ("t = %g, dt = %g\n", t, dt)

    mvold = mold;
    nvold = nold;

    mold = m;
    nold = n;

    uguess = [mold;nold];
    if (ii > 3)
      mguess = mvold * (t+dt - t)/(tstore(ii-2) - t) + ...
               mold * (t+dt - tstore(ii-2))/(t - tstore(ii-2));

      nguess = nvold * (t+dt - t)/(tstore(ii-2) - t) + ...
               nold * (t+dt - tstore(ii-2))/(t - tstore(ii-2));

      if (! all ([mguess;nguess] >= 0))
        printf ("negative guess\n")
      endif

      uguess = max (0, [mguess;nguess]);

    endif

    while (true)

      uold = [mold;nold];
      u = fsolve (@(u) ...
                   fun (u(:), mold(:), nold(:), msh.nx, msh.ny, N,
                        stiffm, stiffn, mass, pressure,
                        G, dGdm, dGdn, dt),
                  uguess, optimset ("Jacobian", "on"));

      all_nzero = all (u >= 0);
      incr = norm ((uguess - u)./(u + 1), inf);
      small_err = incr < 1e-3;

      if (all_nzero && small_err)
        break
      else
        dt *= .5;
        printf (["reducing time step : t = %g, ", ...
                 "dt = %g, all_nzero = %d, small_err = %d\n"],
                t, dt, all_nzero, small_err)
        uguess = [mold;nold];
        if (ii > 3)
          mguess = mvold * (t+dt - t)/(tstore(ii-2) - t) + ...
                   mold * (t+dt - tstore(ii-2))/(t - tstore(ii-2));

          nguess = nvold * (t+dt - t)/(tstore(ii-2) - t) + ...
                   nold * (t+dt - tstore(ii-2))/(t - tstore(ii-2));

          if (! all ([mguess;nguess] >= 0))
            printf ("negative guess\n")
          endif
          uguess = max (0, [mguess;nguess]);

        endif

      endif

    endwhile


    t += dt;
    tstore(ii) = t;
    if (incr == 0)
      dt *= 2;
    else
      dtfact = min (sqrt(.38) * sqrt (1e-3 / incr), 2);
      dt = max (min (dtfact * dt, min (diff (tsave)) / 2), 1e-6);
    endif

    m = u(1:numel(m));
    n = u(numel(m)+1:end);

    assert (all (n >= 0))
    assert (all (m >= 0))


    ii += 1;

  endwhile
  msave (:, end+1) = m;
  nsave (:, end+1) = n;

   figure (1)
   contourf (X, Y, reshape (m, size (X))), colorbar
   title (sprintf ("m @ t = %g", t))
   drawnow

   figure(2)
   contourf (X, Y, reshape (n, size (X))), colorbar
   title (sprintf ("n @ t = %g", t))
   drawnow

   figure (3)
   contourf (X, Y, reshape (pressure (m, n)/PM, size (X))), colorbar
   title (sprintf ("p/PM @ t = %g", t))
   drawnow

endfor


for iplot = 1 : columns (msave)

  figure (1)
  title (sprintf ("m @t = %g", tsave(:, iplot)))
  surf (X, Y, reshape (msave(:, iplot), size (X))), colorbar
  axis ([-2 10 -2 10 0 1.4])
  view (56, 26)
  light
  lighting gouraud
  print ("-dpng", sprintf ("frame_m_%4.4d.png", iplot))
  drawnow

  hold on

  surf (X, Y, reshape (nsave(:, iplot), size (X))), colorbar
  title (sprintf ("n @t = %g", tsave(:, iplot)))
  axis ([-2 10 -2 10 0 1.4])
  view (56, 26)
  light
  lighting gouraud
  print ("-dpng", sprintf ("frame_n_%4.4d.png", iplot))
  drawnow

   figure (3)
   surf (X, Y, reshape ( pressure (msave(:, iplot),nsave(:, iplot))/PM, size (X))), colorbar
   title (sprintf ("p/PM t = %g", tsave(:, iplot)))
   axis ([0 45 0 45 0 1])
   view (56, 26)
   light
   lighting gouraud
   print ("-dpng", sprintf ("frame_p_%4.4d.png", iplot))
   drawnow


endfor

