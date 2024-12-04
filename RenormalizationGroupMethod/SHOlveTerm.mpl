# SHOlveTerm --- solve a simple harmonic oscillator, harmonically forced
#
# Copyright 2024 Robert M. Corless 
#
# dsolve can be inefficient at solving y'' + omega^2*y = sum of simple terms
# (probably getting bogged down in type checking or simplification
#  or its own heuristics, or perhaps it's just its generality).
# This routine attempts to do this simpler task more efficiently,
# by giving up on generality.
#
# Input: expr  : right-hand side: polynomial in sine/cos/exp only
#                with coefficients polynomials in var
#        var   : variable (usually t)
#        omega : natural frequency
#
# Processing: map over a sum of terms
#             match the terms to poly(var)*exp( I*(freq*var + phase) )
#                            or  poly(var)*sin( freq*var + phase )
#                            or  poly(var)*cos( freq*var + phase )
#
#          solve the problem y'' + omega^2*y = sum poly(var)* fn( freq*var+phase )
#
# Output: a particular solution to y'' + omega^2 = ...
#         where the rhs is a sum of polynomials in var times sin/cos/exp.
#
# Comments: dsolve would solve this, but for large expressions might take
# a seemingly unreasonable length of time to do so, and may return 
# unevaluated integrals (because it thinks the answers are too complicated).
# However, this simple routine will not succeed for unexpected inputs.
# Take your pick.

# Copyright (c) 2024 Robert M. Corless
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# ``Software''), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED ``AS IS'', WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

SHOlveTerm := proc( expr, var, omega )
	local A, c, fn, freq, la, m, phase, sol, term, zerp;
	description "Solve y''+ omega^2*y=expr";
    if has( expr, sin ) or has( expr, cos ) then
      # convert,exp will make cos(A) --> (exp(I*A)+exp(-I*A))/2, etc
      # expand will expand (A*exp(B) + C*exp(D))*(E*exp(F)+G*exp(H))
      # to A*E*exp(B)*exp(F) + A*G*exp(B)*exp(H) + ... 
      # combine,exp will make exp(B)*exp(F)-->exp(B+F)
      term := combine( expand( convert(expr,exp ) ), exp );
    else
      term := expr;
    end if; 
	if type( term, `+` ) then
      # Apply this procedure to each term in a sum  
	  return map( procname, term, var, omega );
	else 
      # Look for A*exp(I*(freq*var+phase))
      term := combine( term, exp );
      A := eval(term, exp=1);
      fn := normal( term/A );
      # Each term *should* be A*exp(B).  Should be.
      if normal( expr-A*fn ) <> 0 then
        error "Unexpected term type", term;
      end if;
      # Identify frequency and phase
      if match( fn=exp(I*(freq*var+phase)), var, 'la' ) then
        if eval(freq,la)^2 = omega^2 then
          # resonant case
          m := degree(A, var) + 1; # A must be polynomial
        else
          m := degree(A, var);
        end if;
        P := add( c[k]*var^k, k=0..m );
        # y = P(var)*exp(I*(freq*var+phase))
        # for P(var) a polynomial to be found.
        zerp := eval( omega^2*P + diff(P,var,var) + 2*I*freq*diff(P,var) - freq^2*P - A, la);
        sol := solve( {seq( coeff(zerp,var,k), k=0..m)}, {seq(c[k],k=0..m)} );
        # triangular system, guaranteed to have a solution
        # Throw away the redundant terms c[0]*exp(I*(freq*var+phase))
        return eval( eval( eval(P,sol), c[0]=0)*exp(I*(freq*var+phase)), la );
      else
        error "Unexpected term type, after match", term, A, fn;
      end if;
    end if;
end proc:
