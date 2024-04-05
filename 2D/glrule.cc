/*

Copyright (C) 2018 Carlo de Falco

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see
<https://www.gnu.org/licenses/>.

Author: Carlo de Falco <carlo@guglielmo>
Created: 2018-09-26

*/

#include <octave/oct.h>
#include <sandia_rules.hpp>

using namespace webbur;

DEFUN_DLD(glrule, args, nargout,
          "-*- texinfo -*-\n\
@deftypefn {} {@var{retval} =} glrule (@var{input1}, @var{input2})\n\
@seealso{}\n\
@end deftypefn")
{
  octave_value_list retval;
  int nargin = args.length ();
  if (nargout != 2 || nargin != 1)
    error ("wrong number of inputs or outputs");
  octave_idx_type n = args(0).idx_type_value ();
  
  ColumnVector w(n, 0.), x(n, 0.);
  double *wp = w.fortran_vec (), *xp = x.fortran_vec ();

  laguerre_compute ((int)n, xp, wp);

  retval(1) = w;
  retval(0) = x;  
  
  return retval;
}
